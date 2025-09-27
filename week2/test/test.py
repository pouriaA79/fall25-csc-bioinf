
from python import sys, os, subprocess, re
from typing import Dict, List, Tuple
from python import matplotlib.pyplot as plt
from python import matplotlib.colors as mcolors

# ==============================================================================
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

class Visualizer:
    """Visualizes the aligned, encoded tandem repeat data."""
    output_dir: str
    symbol_to_color: Dict[str, Tuple[float, float, float]]
    def __init__(self, output_dir: str = "."):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        self.symbol_to_color = {}

    def assign_colors(self, symbol_to_motif: Dict[str, str]):
        """Assigns a distinct color to each motif symbol."""
        cmap = plt.get_cmap('tab20')
        symbols = sorted(list(symbol_to_motif.keys()))
        num_colors = len(symbols)
        if num_colors == 0: return
        for i, symbol in enumerate(symbols):
            self.symbol_to_color[symbol] = cmap(float(i) / num_colors)[:3]
        self.symbol_to_color['-'] = (1.0, 1.0, 1.0)
        self.symbol_to_color['?'] = (0.5, 0.5, 0.5)

    def trplot(self, aligned_seqs: Dict[str, str], symbol_to_motif: Dict[str, str], filename: str):
        if len(aligned_seqs) == 0: return
        sample_ids = list(aligned_seqs.keys())
        sequences = [s.split() for s in aligned_seqs.values()] # Assuming space-separated symbols if using numbers
        num_samples = len(sample_ids)
        if num_samples == 0: return
        max_len = max(len(s) for s in sequences) if sequences else 0
        color_matrix: List[List[Tuple[float, float, float]]] = []
        for seq in sequences:
            row: List[Tuple[float, float, float]] = []
            for j in range(max_len):
                symbol = seq[j] if j < len(seq) else '-'
                row.append(self.symbol_to_color.get(symbol, (0.0, 0.0, 0.0)))
            color_matrix.append(row)
        fig, ax = plt.subplots(figsize=(max_len * 0.5, num_samples * 0.5))
        ax.imshow(color_matrix, aspect='auto')
        ax.set_yticks([float(i) for i in range(num_samples)])
        ax.set_yticklabels(sample_ids)
        ax.grid(False)
        plt.savefig(os.path.join(self.output_dir, filename), bbox_inches='tight')
        plt.close()

# ==============================================================================
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
    sample_decomposed = {'seq1': ['ATAT', 'GCGC', 'ATAT', 'CGCG'], 'seq2': ['ATAT', 'GCGC']}
    _, symbol_to_motif, motif_to_symbol = build_encoding_maps(sample_decomposed)
    expected_motif_map = {'ATAT': '!', 'GCGC': '"', 'CGCG': '#'}
    assert motif_to_symbol == expected_motif_map
    motifs_to_encode = ['GCGC', 'CGCG', 'ATAT', 'UNKNOWN']
    expected_encoded_string = '"#!?'
    actual_encoded_string = encode_motifs(motifs_to_encode, motif_to_symbol)
    assert actual_encoded_string == expected_encoded_string
    print("   Encoding logic test PASSED! ✅")

def test_aligner_logic():
    """Robust integration test for MotifAligner."""
    print("-> Running test: MotifAligner...")
    aligner = MotifAligner(mafft_path="mafft")
    encoded_sequences = {'sample_A': '!#%', 'sample_B': '!"#$%' }
    actual_alignment = aligner.align_encoded_seqs(encoded_sequences)
    assert sorted(list(actual_alignment.keys())) == sorted(list(encoded_sequences.keys()))
    aligned_seqs = list(actual_alignment.values())
    assert len(aligned_seqs) > 0
    first_seq_len = len(aligned_seqs[0])
    for i in range(1, len(aligned_seqs)):
        assert len(aligned_seqs[i]) == first_seq_len
    print(f"   - MAFFT actual output: {actual_alignment}")
    print("   MotifAligner test PASSED! ✅")

def run_visualization_demo():
    """Runs the full pipeline with test data and saves an image."""
    print("-> Running visualization demo...")
    # 1. Setup sample data
    decomposer = Decomposer(k=4)
    seq_data = {
        'Sample-1': 'AGATAGATGCGC',
        'Sample-2': 'AGATGCGC',
        'Sample-3': 'AGATAGATAGATCGCG'
    }
    decomposed_data = {name: decomposer.decompose(seq) for name, seq in seq_data.items()}
    
    # 2. Encode
    _, symbol_to_motif, motif_to_symbol = build_encoding_maps(decomposed_data)
    encoded_data = {name: encode_motifs(motifs, motif_to_symbol) for name, motifs in decomposed_data.items()}
    
    # 3. Align (Note: This uses single characters, so no space processing is needed)
    aligner = MotifAligner()
    aligned_data = aligner.align_encoded_seqs(encoded_data)
    
    # 4. Visualize
    aligned_data_spaced = {name: ' '.join(list(seq)) for name, seq in aligned_data.items()}

    viz = Visualizer(output_dir=".")
    viz.assign_colors(symbol_to_motif)
    viz.trplot(aligned_data_spaced, symbol_to_motif, filename="test_visualization.png")
    
    print("   Visualization demo complete! Check for 'test_visualization.png' file. 🖼️")

# ==============================================================================
# ==============================================================================
if __name__ == "__main__":
    print("\nStarting Codon tests...")
    test_decomposer_logic()
    test_encoding_logic()
    test_aligner_logic()
    run_visualization_demo() # New visualization step
    print("\nAll tests and demos finished successfully!")



