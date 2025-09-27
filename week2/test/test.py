# File: week2/test/test.py
# Standalone test script with a more robust aligner test.

from python import sys, os, subprocess, re
from typing import Dict, List, Tuple

# ==============================================================================
# Part 1: Copied Classes and Functions from main_type.py
# ==============================================================================

class Decomposer:
    """Decomposes DNA sequences into simple overlapping k-mers."""
    k: int
    def __init__(self, k: int = 4):
        self.k = k

    def decompose(self, sequence: str) -> List[str]:
        motifs: List[str] = []
        n: int = len(sequence)
        if n < self.k:
            return [sequence] if sequence else []
        for i in range(n - self.k + 1):
            motifs.append(sequence[i:i+self.k])
        return motifs

PRINTABLE_CHARS: str = "!\"#$%&'()*+,-./0123456789:;<=?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"

def build_encoding_maps(all_decomposed_trs: Dict[str, List[str]]) -> Tuple[Dict[str, int], Dict[str, str], Dict[str, str]]:
    """Counts motifs, sorts by frequency, and creates encoding maps."""
    counter: Dict[str, int] = {}
    for name in all_decomposed_trs:
        for motif in all_decomposed_trs[name]:
            counter[motif] = counter.get(motif, 0) + 1
    
    items: List[Tuple[str, int]] = []
    for key in counter:
        items.append((key, counter[key]))

    n = len(items)
    for i in range(n):
        for j in range(0, n-i-1):
            if items[j][1] < items[j+1][1]:
                items[j], items[j+1] = items[j+1], items[j]
    
    sorted_motifs: List[str] = [item[0] for item in items]
    symbol_to_motif: Dict[str, str] = {}
    motif_to_symbol: Dict[str, str] = {}
    
    for i in range(len(sorted_motifs)):
        if i >= len(PRINTABLE_CHARS): break
        motif = sorted_motifs[i]
        char = PRINTABLE_CHARS[i]
        symbol_to_motif[char] = motif
        motif_to_symbol[motif] = char
    return counter, symbol_to_motif, motif_to_symbol

def encode_motifs(decomposed_motifs: List[str], motif_to_symbol: Dict[str, str]) -> str:
    """Encodes a list of motifs into a single string."""
    encoded_str = ""
    for motif in decomposed_motifs:
        encoded_str += motif_to_symbol.get(motif, '?')
    return encoded_str

class MotifAligner:
    """Aligns character-encoded sequences using an external MSA tool."""
    mafft_path: str
    def __init__(self, mafft_path: str = "mafft"):
        self.mafft_path = mafft_path

    def align_encoded_seqs(self, encoded_seqs: Dict[str, str]) -> Dict[str, str]:
        if len(encoded_seqs) == 0:
            return {}
        tmp_input_file: str = "tmp_encoded_sequences.fasta"
        aligned_output: Dict[str, str] = {}
        try:
            with open(tmp_input_file, "w") as f:
                for name, seq in encoded_seqs.items():
                    f.write(f">{name}\n{seq}\n")
            mafft_command = [self.mafft_path, "--quiet", "--text", "--auto", "--inputorder", tmp_input_file]
            result = subprocess.run(mafft_command, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"MAFFT failed. Stderr: {result.stderr}")
                return {}
            
            current_name = ""
            current_seq = ""
            for line in result.stdout.splitlines():
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    if current_name and current_seq:
                        aligned_output[current_name] = current_seq
                    current_name = line[1:]
                    current_seq = ""
                else:
                    current_seq += line
            if current_name and current_seq:
                aligned_output[current_name] = current_seq
            return aligned_output
        finally:
            if os.path.exists(tmp_input_file):
                os.remove(tmp_input_file)

# ==============================================================================
# Part 2: Test Functions
# ==============================================================================

def test_decomposer_logic():
    """Tests the Decomposer class logic."""
    print("-> Running test: Decomposer...")
    decomposer = Decomposer(k=4)
    sequence = "AGATAGAT"
    expected_motifs = ["AGAT", "GATA", "ATAG", "TAGA", "AGAT"]
    actual_motifs = decomposer.decompose(sequence)
    assert actual_motifs == expected_motifs
    print("   Decomposer test PASSED! ✅")

def test_encoding_logic():
    """Tests the motif encoding logic."""
    print("-> Running test: Encoder...")
    sample_decomposed = {
        'seq1': ['ATAT', 'GCGC', 'ATAT', 'CGCG'], # Freq: ATAT=3, GCGC=2, CGCG=1
        'seq2': ['ATAT', 'GCGC']
    }
    _, _, motif_to_symbol = build_encoding_maps(sample_decomposed)
    expected_motif_map = {'ATAT': '!', 'GCGC': '"', 'CGCG': '#'}
    assert motif_to_symbol == expected_motif_map
    motifs_to_encode = ['GCGC', 'CGCG', 'ATAT', 'UNKNOWN']
    expected_encoded_string = '"#!?'
    actual_encoded_string = encode_motifs(motifs_to_encode, motif_to_symbol)
    assert actual_encoded_string == expected_encoded_string
    print("   Encoding logic test PASSED! ✅")

def test_aligner_logic():
    """
    A more robust integration test for MotifAligner.
    Requires mafft to be installed.
    """
    print("-> Running test: MotifAligner...")
    aligner = MotifAligner(mafft_path="mafft")
    
    encoded_sequences = {
        'sample_A': '!#%',  # Shorter sequence
        'sample_B': '!"#$%' # Longer sequence
    }
    
    actual_alignment = aligner.align_encoded_seqs(encoded_sequences)
    
    # --- Start of the new, robust checks ---
    
    # 1. Check if the output contains the same sample names as the input
    assert sorted(list(actual_alignment.keys())) == sorted(list(encoded_sequences.keys())), "Aligner FAILED: Output keys do not match input keys."

    # 2. Check if the alignment produced any output sequences
    aligned_seqs = list(actual_alignment.values())
    assert len(aligned_seqs) > 0, "Aligner FAILED: No sequences were returned."
    
    # 3. THE MOST IMPORTANT CHECK: All output sequences must have the same length
    first_seq_len = len(aligned_seqs[0])
    for i in range(1, len(aligned_seqs)):
        assert len(aligned_seqs[i]) == first_seq_len, f"Aligner FAILED: Sequences have unequal lengths after alignment. Got {first_seq_len} and {len(aligned_seqs[i])}"

    # --- End of the new, robust checks ---
    
    # Optional: Print the actual alignment for manual inspection
    print(f"   - MAFFT actual output: {actual_alignment}")
    
    print("   MotifAligner test PASSED! ✅")


# ==============================================================================
# Part 3: Main execution block
# ==============================================================================
if __name__ == "__main__":
    print("\nStarting Codon tests...")
    test_decomposer_logic()
    test_encoding_logic()
    test_aligner_logic()
    print("\nAll tests finished successfully! 🎉")