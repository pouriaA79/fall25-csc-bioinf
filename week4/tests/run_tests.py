# tests/run_tests.py  (فقط حذف type hints از امضاها)
import sys
import time
from pathlib import Path

ROOT = Path(__file__).parent.parent
SRC_DIR = ROOT / "src"
DATA_DIR = ROOT / "data"
sys.path.insert(0, str(SRC_DIR))

from global_align import align as global_align
from local import align as local_align
from semi_global import align as semi_global_align
from affine import affine_align as affine_align

def read_fasta(path):
    seq = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            cleaned = "".join(ch for ch in line if ch.isalpha() or ch == "-").upper()
            if cleaned:
                seq.append(cleaned)
    return "".join(seq)

def read_multifasta(path):
    seqs = {}
    current_name = None
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_name = line[1:].split()[0]
                seqs[current_name] = ""
            else:
                clean = "".join(ch for ch in line if ch.isalpha() or ch == "-").upper()
                if current_name:
                    seqs[current_name] += clean
    return seqs


def test_pair(method_name, func, seq1, seq2):
    start = time.time()
    _, _, score = func(seq1, seq2)
    end = time.time()
    ms = int((end - start) * 1000)
    print(f"{method_name:<20} python      {ms:>6}ms   score={score}")

def test_pair_q1(method_name, func, seq1, seq2):
    """اجرای تست و چاپ زمان اجرا (می‌گیرد سکانس‌ها، نه اسم فایل‌ها)"""
    start = time.time()
    _, _, score = func(seq1, seq2)
    end = time.time()
    ms = int((end - start) * 1000)
    print(f"{method_name:<20} python      {ms:>6}ms   score={score}")

def main():
    print("Method              Language    Runtime")
    print("--------------------------------------")
    human = read_fasta(DATA_DIR/ "MT-human.fa")
    orang = read_fasta(DATA_DIR/ "MT-orang.fa")
    # تست‌های اصلی
    test_pair("global-mt_human", global_align, human, orang)
    test_pair("affine-mt_human", affine_align, human, orang)
    test_pair("semiglobal-mt_human", semi_global_align, human, orang)
    test_pair("local-mt_human", local_align, human, orang)

    # q1..q5 در برابر t1..t5
    qs = read_multifasta(DATA_DIR / "q1.fa")
    # print(qs)
    ts = read_multifasta(DATA_DIR / "t1.fa")

    for i in range(1, 6):
        qkey, tkey = f"q{i}", f"t{i}"
        if qkey in qs and tkey in ts:
            print(qs[qkey])
            test_pair_q1(f"global-{qkey}", global_align, qs[qkey], ts[tkey])
            test_pair_q1(f"local-{qkey}", local_align, qs[qkey], ts[tkey])
            test_pair_q1(f"semi-{qkey}", semi_global_align, qs[qkey], ts[tkey])
            test_pair_q1(f"affine-{qkey}", affine_align, qs[qkey], ts[tkey])
        else:
            print(f"⚠️  Missing {qkey} or {tkey} in multi-FASTA files.")

if __name__ == "__main__":
    main()
