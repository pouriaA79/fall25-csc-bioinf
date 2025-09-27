from python import os, subprocess, re
from python import matplotlib.pyplot as plt
from python import matplotlib.colors as mcolors # Added for normalization
from typing import Dict, List, Optional, Tuple


class Decomposer:
    """Decomposes DNA sequences into simple overlapping k-mers."""
    k: int

    def __init__(self, k: int = 4):
        self.k = k

    def decompose(self, sequence: str) -> List[str]:
        """Converts a single DNA sequence into overlapping k-mers."""
        motifs: List[str] = []
        n: int = len(sequence)
        if n < self.k:
            return [sequence] if sequence else []
        for i in range(n - self.k + 1):
            motifs.append(sequence[i:i+self.k])
        return motifs

    def decompose_from_fasta(self, filepath: str) -> Dict[str, List[str]]:
        """Reads sequences from a FASTA file and decomposes them into motifs."""
        sequences: Dict[str, List[str]] = {}
        current_name: str = ""
        current_seq: str = ""
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    if current_name and current_seq:
                        sequences[current_name] = self.decompose(current_seq)
                    current_name = line[1:]
                    current_seq = ""
                else:
                    current_seq += line.upper()
            if current_name and current_seq:
                sequences[current_name] = self.decompose(current_seq)
        return sequences


def build_encoding_maps(all_decomposed_trs: Dict[str, List[str]]) -> Tuple[Dict[str, int], Dict[str, str], Dict[str, str]]:
    """
    Counts all unique motifs, sorts them by frequency, and returns encoding maps
    using numeric symbols from 1 to 1000.
    """
    # 1. Count frequencies of all motifs
    counter: Dict[str, int] = {}
    for name in all_decomposed_trs:
        for motif in all_decomposed_trs[name]:
            counter[motif] = counter.get(motif, 0) + 1
            
    # 2. Sort all motifs by frequency
    items: List[Tuple[str, int]] = []
    for motif, count in counter.items():
        items.append((motif, count))

    n = len(items)
    for i in range(n):
        for j in range(0, n-i-1):
            if items[j][1] < items[j+1][1]:
                items[j], items[j+1] = items[j+1], items[j]
    
    sorted_motifs: List[str] = [item[0] for item in items]

    # 3. Assign a numeric symbol to each unique motif up to 1000
    symbol_to_motif: Dict[str, str] = {}
    motif_to_symbol: Dict[str, str] = {}
    
    if len(sorted_motifs) > 1000:
        print("Warning: More than 1000 unique motifs found. Only the top 1000 will be encoded.")
    
    for i in range(len(sorted_motifs)):
        if i >= 1000: break
        motif = sorted_motifs[i]
        symbol = str(i + 1) # Use numbers as symbols
        symbol_to_motif[symbol] = motif
        motif_to_symbol[motif] = symbol
            
    return counter, symbol_to_motif, motif_to_symbol

def encode_motifs(decomposed_motifs: List[str], motif_to_symbol: Dict[str, str]) -> str:
    """Encodes a single list of motifs into a space-separated string of numeric symbols."""
    encoded_parts: List[str] = []
    for motif in decomposed_motifs:
        # If a motif was too infrequent to get a number, it will be marked with '?'
        encoded_parts.append(motif_to_symbol.get(motif, '?'))
    # Join with spaces to handle multi-digit numbers
    return " ".join(encoded_parts)


class MotifAligner:
    """Aligns the character-encoded TR sequences using an external MSA tool."""
    mafft_path: str

    def __init__(self, mafft_path: str = "mafft"):
        self.mafft_path = mafft_path

    def align_encoded_seqs(self, encoded_seqs: Dict[str, str]) -> Dict[str, str]:
        """Runs MAFFT on a dictionary of encoded sequences."""
        if len(encoded_seqs) == 0:
            return {}

        tmp_input_file: str = "tmp_encoded_sequences.fasta"
        tmp_output_file: str = "tmp_aligned_sequences.fasta"
        aligned_output: Dict[str, str] = {}
        try:
            with open(tmp_input_file, "w") as f:
                for name, seq in encoded_seqs.items():
                    f.write(f">{name}\n{seq}\n")

            mafft_command = [
                self.mafft_path, "--quiet", "--text", "--auto", "--inputorder", tmp_input_file
            ]
            
            result = subprocess.run(
                mafft_command, capture_output=True, text=True
            )

            if result.returncode != 0:
                print(f"MAFFT failed with return code: {result.returncode}")
                print(f"MAFFT stderr: {result.stderr}")
                return {}
            
            with open(tmp_output_file, "w") as f_out:
                f_out.write(result.stdout)

            current_name = ""
            current_seq = ""
            with open(tmp_output_file, "r") as f_in:
                for line in f_in:
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

            processed_alignment: Dict[str, str] = {}
            for name, seq in aligned_output.items():
                # This regex ensures that every number (one or more digits) or a hyphen
                # is treated as a distinct "word" by surrounding it with spaces.
                # This makes sure that seq.split() correctly separates them.
                # Example: "12-3" -> " 12 - 3 " -> "12 - 3"
                spaced_seq = re.sub(r'(\d+|-|\?)', r' \1 ', seq) # Also include '?' for safety
                spaced_seq = re.sub(r'\s+', ' ', spaced_seq).strip() # Consolidate multiple spaces
                processed_alignment[name] = spaced_seq
            
            return processed_alignment

        finally:
            if os.path.exists(tmp_input_file):
                os.remove(tmp_input_file)
            if os.path.exists(tmp_output_file):
                os.remove(tmp_output_file)

class Visualizer:
    """Visualizes the aligned, encoded tandem repeat data."""
    output_dir: str
    symbol_to_color: Dict[str, Tuple[float, float, float]]

    def __init__(self, output_dir: str = "results"):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        self.symbol_to_color = {}

    def assign_colors(self, symbol_to_motif: Dict[str, str]):
        """Assigns a distinct color to each motif symbol directly from a colormap."""
        
        # Sort symbols for consistent color assignment
        numeric_symbols = sorted([s for s in symbol_to_motif.keys() if s.isdigit()], key=int)
        other_symbols = sorted([s for s in symbol_to_motif.keys() if not s.isdigit()])
        all_sorted_symbols = numeric_symbols + other_symbols

        num_unique_motifs = len(all_sorted_symbols)
        
        if num_unique_motifs == 0:
            return

        # Choose a colormap based on the number of unique motifs
        cmap_name = 'tab20' if num_unique_motifs <= 20 else 'gist_ncar'
        
        # Create a colormap with the exact number of colors needed
        cmap = plt.get_cmap(cmap_name, num_unique_motifs)
        
        # Directly assign colors from the colormap's color list. This is simpler.
        for i, symbol in enumerate(all_sorted_symbols):
            # For a discrete colormap created with N colors, cmap(i) gives the i-th color.
            # We take the first 3 elements (R, G, B) and ignore Alpha if it exists.
            self.symbol_to_color[symbol] = cmap(i)[:3]

        self.symbol_to_color['-'] = (1.0, 1.0, 1.0) # White for gaps
        if '?' not in self.symbol_to_color:
            self.symbol_to_color['?'] = (0.5, 0.5, 0.5) # Grey for private/unknown motifs


    def trplot(self, aligned_seqs: Dict[str, str], symbol_to_motif: Dict[str, str], title: str = "Tandem Repeat Visualization", filename: Optional[str] = None):
        """The main visualization function, modified to match the target image."""
        if len(aligned_seqs) == 0:
            print("Cannot plot empty alignment.")
            return

        sample_ids = list(aligned_seqs.keys())
        sequences_as_lists: List[List[str]] = [seq.split() for seq in aligned_seqs.values()]

        num_samples = len(sample_ids)
        max_len = max(len(s) for s in sequences_as_lists) if num_samples > 0 else 0

        if max_len == 0:
            print("Cannot plot alignment with zero length.")
            return
            
        consensus_labels: List[str] = []
        for i in range(max_len):
            column_symbols: Dict[str, int] = {}
            for seq_list in sequences_as_lists:
                if i < len(seq_list):
                    symbol = seq_list[i]
                    if symbol != '-':
                        column_symbols[symbol] = column_symbols.get(symbol, 0) + 1
            
            if len(column_symbols) == 0:
                consensus_labels.append('')
                continue

            most_common_symbol = max(column_symbols, key=lambda k: column_symbols[k])
            dna_motif = symbol_to_motif.get(most_common_symbol, '?')
            consensus_labels.append(dna_motif)

        color_matrix: List[List[Tuple[float, float, float]]] = []
        for symbols in sequences_as_lists:
            row: List[Tuple[float, float, float]] = []
            for j in range(max_len):
                symbol = symbols[j] if j < len(symbols) else '-'
                color = self.symbol_to_color.get(symbol, (0.0, 0.0, 0.0)) # Default to black if symbol not found
                row.append(color)
            color_matrix.append(row)

        fig, ax = plt.subplots(figsize=(max_len* 0.2 , num_samples * 0.5)) 
        # ==========================================================

        ax.imshow(color_matrix, aspect='auto', interpolation='nearest')

        ax.imshow(color_matrix, aspect='auto', interpolation='nearest')

        y_ticks: List[float] = [float(i) for i in range(num_samples)]
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(sample_ids)
        
        x_ticks: List[float] = [float(i) for i in range(max_len)]
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(consensus_labels, rotation='vertical', fontsize=6) 
        ax.tick_params(axis='x', length=0)

        ax.set_xlabel("") 
        ax.set_ylabel("")
        ax.set_title(title, pad=20)
        
        ax.grid(False)

        if filename:
            plt.savefig(os.path.join(self.output_dir, filename), bbox_inches='tight', dpi=150)
        plt.close()

    def plot_motif_color_map(self, symbol_to_motif: Dict[str, str], filename: Optional[str] = None):
        """Creates a legend showing which color maps to which motif."""
        if len(symbol_to_motif) == 0:
            return
        
        numeric_symbols = sorted([s for s in symbol_to_motif.keys() if s.isdigit()], key=int)
        other_symbols = sorted([s for s in symbol_to_motif.keys() if not s.isdigit()])
        symbols = numeric_symbols + other_symbols
        
        fig, ax = plt.subplots(figsize=(8, len(symbols) * 0.3))
        
        for i, symbol in enumerate(symbols):
            motif = symbol_to_motif[symbol]
            color = self.symbol_to_color.get(symbol, (0.0, 0.0, 0.0))
            ax.add_patch(plt.Rectangle((0, i), 1, 1, color=color))
            ax.text(1.2, i + 0.5, f"'{symbol}' = {motif}", va='center', fontsize=12)
            
        ax.set_ylim(-0.5, len(symbols))
        ax.set_xlim(0, 10)
        ax.axis('off')
        ax.set_title("Motif Color Legend")

        if filename:
            plt.savefig(os.path.join(self.output_dir, filename), bbox_inches='tight')
        plt.close()
# ==============================================================================
# Main execution logic
# ==============================================================================
def main() -> None:
    fasta_file: str = "SORL1_tr_sequence.fa"
    results_dir: str = "."

    if not os.path.exists(fasta_file):
        print(f"Error: Input file '{fasta_file}' not found.")
        return

    os.makedirs(results_dir, exist_ok=True)

    # 1. Decompose all sequences into motifs
    # Changed k to 3 based on common practice, adjust if your motifs are k=4
    decomposer = Decomposer(k=3) # Changed back to k=4 as it was previously
    decomposed_trs = decomposer.decompose_from_fasta(fasta_file)
    print("Step 1: Decomposition complete.")

    # 2. Build encoding maps from ALL motifs and then encode each sequence
    motif_counter, symbol_to_motif, motif_to_symbol = build_encoding_maps(decomposed_trs)
    
    encoded_seqs: Dict[str, str] = {}
    for name, motifs in decomposed_trs.items():
        encoded_seqs[name] = encode_motifs(motifs, motif_to_symbol)
    
    print("Step 2: Motif encoding complete.")
    print(f" - {len(symbol_to_motif)} motifs encoded.")
    # print(" - Symbol Map:", symbol_to_motif) # This can be very long, so commented out for cleaner output


    # 3. Align the character-encoded sequences
    aligner = MotifAligner(mafft_path="mafft")
    aligned_encoded_seqs = aligner.align_encoded_seqs(encoded_seqs)
    print("Step 3: Alignment of encoded sequences complete.")
    # For debugging, you can uncomment this to see the processed MAFFT output
    # for name, seq in aligned_encoded_seqs.items():
    #      print(f" - {name}: {seq}")

    # 4. Visualize the results
    viz = Visualizer(output_dir=results_dir)
    viz.assign_colors(symbol_to_motif)
    
    viz.trplot(aligned_encoded_seqs,
               symbol_to_motif,
               title="Visualization of Tandem Repeats in SORL1",
               filename="tr_visualization.png")

    viz.plot_motif_color_map(symbol_to_motif,
                             filename="motif_legend.png")
    print(f"\nStep 4: Visualization complete. Results saved in '{results_dir}' directory.")

if __name__ == "__main__":

    main()
