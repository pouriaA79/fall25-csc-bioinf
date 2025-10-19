import sys
import time
from python import os
from typing import Dict

from .. src.global_align_codon import global_align
from .. src.local_codon import local_align
from .. src.semi_global_codon import semi_global_align
from .. src.affine_codon import affine_align


def read_fasta(path: str) -> str:
    seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            cleaned = "".join(ch for ch in line if ch.isalpha() or ch == "-").upper()
            if cleaned:
                seq.append(cleaned)
    return "".join(seq)

from typing import Dict, Optional

def read_multifasta(path: str) -> Dict[str, str]:
    seqs = {}  
    current_name_opt: Optional[str] = None
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                parts = line[1:].split()
                if len(parts) > 0:
                    current_name = parts[0]
                    seqs[current_name] = ""
                    current_name_opt = current_name
            else:
                clean = "".join(ch for ch in line if ch.isalpha() or ch == "-").upper()
                if current_name_opt is not None and current_name_opt in seqs:
                    seqs[current_name_opt] += clean
    return seqs


def test_pair(method_name: str, func, file_a: str, file_b: str, base_path: str):
    seq1 = read_fasta(os.path.join(base_path, file_a))
    seq2 = read_fasta(os.path.join(base_path, file_b))
    start = time.time()
    _, _, score = func(seq1, seq2)
    end = time.time()
    ms = int((end - start) * 1000)
    print(f"{method_name:<20} codon      {ms:>6}ms   score={score}")


def test_pair_q(method_name: str, func, seq1: str, seq2: str):
    start = time.time()
    _, _, score = func(seq1, seq2)
    end = time.time()
    ms = int((end - start) * 1000)
    print(f"{method_name:<20} codon      {ms:>6}ms   score={score}")


def resolve_data_folder() -> str:
    argv = sys.argv
    if len(argv) >= 2:
        path = argv[1]
        if os.path.isdir(path):
            return os.path.abspath(path)

    here = os.path.dirname(argv[0])
    candidate = os.path.abspath(os.path.join(here, "..", "data"))
    if os.path.isdir(candidate):
        return candidate

    print("❌ Cannot find data folder. Please pass it as argument.")
    sys.exit(1)


def main():
    data_path = resolve_data_folder()

    print("✅ Using data folder:", data_path)
    print("Method              Language    Runtime")
    print("--------------------------------------")

    # ---- Main tests ----
    test_pair("global-mt_human", global_align, "MT-human.fa", "MT-orang.fa", data_path)
    test_pair("affine-mt_human", affine_align, "MT-human.fa", "MT-orang.fa", data_path)
    test_pair("semi-global-mt_human", semi_global_align, "MT-human.fa", "MT-orang.fa", data_path)
    test_pair("local-mt_human", local_align, "MT-human.fa", "MT-orang.fa", data_path)

    # ---- q1..q5 vs t1..t5 ----
    qs = read_multifasta(os.path.join(data_path, "q1.fa"))
    ts = read_multifasta(os.path.join(data_path, "t1.fa"))

    for i in range(1, 6):
        qkey, tkey = f"q{i}", f"t{i}"
        if qkey in qs and tkey in ts:
            test_pair_q(f"global-{qkey}", global_align, qs[qkey], ts[tkey])
            test_pair_q(f"local-{qkey}", local_align, qs[qkey], ts[tkey])
            test_pair_q(f"semi-{qkey}", semi_global_align, qs[qkey], ts[tkey])
            test_pair_q(f"affine-{qkey}", affine_align, qs[qkey], ts[tkey])
        else:
            print(f"⚠️ Missing {qkey} or {tkey} in multi-FASTA files.")


if __name__ == "__main__":
    main()


