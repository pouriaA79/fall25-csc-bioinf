# from dbg_type import DBG
# from utils_type import read_data
# from python import sys
# from python import os

# sys.setrecursionlimit(100000)


# if __name__ == "__main__":
#     import sys
#     argv = sys.argv
#     data_path = f"../../data/{argv[1]}"
#     # data_path = os.path.abspath(data_path)
#     short1, short2, long1 = read_data(data_path)

#     k = 25
#     # output_dir = f"../data/{argv[1]}"
#     output_file_path = f"{data_path}/contigs.fasta"
#     dbg = DBG(k=k, data_list=[short1, short2, long1])
#     # dbg.show_count_distribution()
#     # try:
#     #     os.mkdir(output_dir)
#     # except OSError:
#     #     pass # اگر دایرکتوری از قبل وجود داشته باشد، مشکلی نیست

#     with open(output_file_path, 'w') as f:
#         for i in range(20):
#             c = dbg.get_longest_contig()
#             if c is None:
#                 break
#             print(i, len(c))
#             f.write(f">contig_{i}\n")
#             f.write(c + "\n")
#     print(f"WROTE: {output_file_path}")

# from dbg_type import DBG
# from utils_type import read_data
# from python import os, sys  # در Codon بهتره os/sys را از bridge بیاریم

# # اختیاری: روی Codon اثر خاصی نداره
# sys.setrecursionlimit(100000)

# def resolve_data_path(arg: str) -> str:
#     """
#     اگر arg مسیر کامل/نسبیِ موجود باشد همان را برمی‌گرداند.
#     اگر فقط نام دیتاست باشد (data1, data2, ...)، آن را زیر cwd/data می‌سازد.
#     """
#     p = os.path.abspath(arg)
#     if os.path.isdir(p):
#         return p
#     base = os.getcwd()  # در CI == week1
#     return os.path.abspath(os.path.join(base, "data", arg))

# if __name__ == "__main__":
#     if len(sys.argv) < 2:
#         print("usage: codon run -plugin seq main_type.py <dataset-path | dataset-name>")
#         sys.exit(1)

#     data_path = resolve_data_path(sys.argv[1])
#     short1, short2, long1 = read_data(data_path)

#     k = 25
#     dbg = DBG(k=k, data_list=[short1, short2, long1])

#     out_fa = os.path.join(data_path, "contigs.fasta")
#     with open(out_fa, "w") as f:
#         for i in range(20):
#             c = dbg.get_longest_contig()
#             if c is None:
#                 break
#             print(i, len(c))
#             f.write(f">contig_{i}\n{c}\n")

#     # برای evaluate.sh
#     print(f"WROTE: {out_fa}")




# week1/codes/codon/main_type.py

from dbg_type import DBG
from python import os, sys  # از پل پایتون فقط os/sys را می‌گیریم؛ __file__ استفاده نمی‌کنیم

def resolve_data_path(arg: str) -> str:
    """
    اگر arg مسیر موجود باشد همان را برمی‌گرداند.
    اگر فقط نام دیتاست باشد (مثلاً data1)، مسیر cwd/data/arg را می‌سازد.
    (در CI، cwd = week1 است، پس week1/data/data1 ساخته می‌شود.)
    """
    p = os.path.abspath(arg)
    if os.path.isdir(p):
        return p
    base = os.getcwd()
    return os.path.abspath(os.path.join(base, "data", arg))

def read_fasta(fp: str):
    """ساده و مطمئن: کل FASTA را به لیست توالی‌ها تبدیل می‌کند."""
    seqs = []
    cur = []
    f = open(fp, "r")
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
                cur = []
        else:
            cur.append(line)
    if cur:
        seqs.append("".join(cur))
    f.close()
    return seqs

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: codon run -plugin seq main_type.py <dataset-path | dataset-name>")
        sys.exit(1)

    data_path = resolve_data_path(sys.argv[1])
    print(f"[CODON] cwd={os.getcwd()}  arg={sys.argv[1]}  data_path={data_path}")

    # فایل‌ها را مستقیم بخوان (بدون utils_type)
    s1 = read_fasta(os.path.join(data_path, "short_1.fasta"))
    s2 = read_fasta(os.path.join(data_path, "short_2.fasta"))
    l1 = read_fasta(os.path.join(data_path, "long.fasta"))

    # پیام‌های سازگار با لاگ Python برای دیباگ و evaluate.sh
    print(f"short_1.fasta {len(s1)} {len(s1[0]) if s1 else 0}")
    print(f"short_2.fasta {len(s2)} {len(s2[0]) if s2 else 0}")
    print(f"long.fasta {len(l1)} {len(l1[0]) if l1 else 0}")

    k = 25
    dbg = DBG(k=k, data_list=[s1, s2, l1])

    out_fa = os.path.join(data_path, "contigs.fasta")
    f = open(out_fa, "w")
    for i in range(20):
        c = dbg.get_longest_contig()
        if c is None:
            break
        print(i, len(c))
        f.write(f">contig_{i}\n{c}\n")
    f.close()

    print(f"WROTE: {out_fa}")











