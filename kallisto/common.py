# common.py — Codon-safe helpers

def read_fastq(path: str) -> list[str]:
    reads: list[str] = []
    with open(path) as f:
        i = 0
        seq = ""
        for line in f:
            i = (i + 1) % 4
            if i == 1:
                # header
                pass
            elif i == 2:
                seq = line.strip().upper().replace('U', 'T')
                reads.append(seq)
            # i==3 plus و i==0 کیفیت—برای MVP اهمیتی ندارد
    return reads

def effective_length(L: int, frag_len: int) -> float:
    e = float(L - frag_len + 1)
    if e < 1.0:
        e = 1.0
    return e

# common.py



def write_abundance(path: str, rows: list[tuple[str,int,float,float,float]]) -> None:
    from python import os
    path = f"{path}".strip()
    with open(path, "w") as w:
        w.write("target_id\tlength\teff_length\test_counts\ttpm\n")
        for t, L, effL, cnt, tpm in rows:
            w.write(f"{t}\t{L}\t{effL:.6f}\t{cnt:.6f}\t{tpm:.6f}\n")

