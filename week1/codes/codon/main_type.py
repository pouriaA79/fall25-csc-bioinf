from dbg_type import DBG
from utils_type import read_data
from python import sys
from python import os

sys.setrecursionlimit(100000)


if __name__ == "__main__":
    import sys
    argv = sys.argv
    data_path = f"data/{argv[1]}"
    print(data_path, 546)
    # data_path = os.path.abspath(data_path)
    short1, short2, long1 = read_data(data_path)

    k = 25
    # output_dir = f"../data/{argv[1]}"
    output_file_path = f"{data_path}/contigs.fasta"
    dbg = DBG(k=k, data_list=[short1, short2, long1])
    # dbg.show_count_distribution()
    # try:
    #     os.mkdir(output_dir)
    # except OSError:
    #     pass # اگر دایرکتوری از قبل وجود داشته باشد، مشکلی نیست

    with open(output_file_path, 'w') as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            print(i, len(c))
            f.write(f">contig_{i}\n")
            f.write(c + "\n")
    print(f"WROTE: {output_file_path}")

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



# from dbg_type import DBG
# from python import os, sys, traceback

# # ---------- path utils ----------
# def resolve_data_path(arg: str) -> str:
#     """
#     اگر arg دایرکتوری باشد همان را برمی‌گرداند.
#     اگر فقط نام دیتاست باشد (مثل data1)، مسیر cwd/data/arg را می‌سازد.
#     (در CI، cwd = week1 است)
#     """
#     p = os.path.abspath(arg)
#     if os.path.isdir(p):
#         return p
#     base = os.getcwd()
#     return os.path.abspath(os.path.join(base, "data", arg))

# # ---------- FASTA I/O (ساده و سازگار با Codon) ----------
# def read_fasta_file(path: str):
#     seqs = []
#     buf = []
#     with open(path, "r") as fh:
#         for raw in fh:
#             if not raw:
#                 continue
#             if raw[0] == ">":
#                 if buf:
#                     seqs.append("".join(buf))
#                     buf = []
#             else:
#                 buf.append(raw.strip())
#         if buf:
#             seqs.append("".join(buf))
#     return seqs

# def read_data_local(dp: str):
#     s1p = os.path.join(dp, "short_1.fasta")
#     s2p = os.path.join(dp, "short_2.fasta")
#     lp  = os.path.join(dp, "long.fasta")

#     # وجود فایل‌ها را چک کنیم تا اگر نبود، پیام واضح بدهیم
#     for p in (s1p, s2p, lp):
#         if not os.path.isfile(p):
#             raise FileNotFoundError(f"missing file: {p}")

#     s1 = read_fasta_file(s1p)
#     s2 = read_fasta_file(s2p)
#     l1 = read_fasta_file(lp)

#     def lensafe(v): return len(v[0]) if v else 0
#     print(f"short_1.fasta {len(s1)} {lensafe(s1)}")
#     print(f"short_2.fasta {len(s2)} {lensafe(s2)}")
#     print(f"long.fasta {len(l1)} {lensafe(l1)}")
#     return s1, s2, l1

# # ---------- main ----------
# if __name__ == "__main__":
#     try:
#         if len(sys.argv) < 2:
#             print("usage: codon run -plugin seq main_type.py <dataset-path | dataset-name>")
#             sys.exit(1)

#         data_path = resolve_data_path(sys.argv[1])
#         print(f"[CODON] data_path = {data_path}")

#         short1, short2, long1 = read_data_local(data_path)

#         k = 25
#         dbg = DBG(k=k, data_list=[short1, short2, long1])

#         out_fa = os.path.join(data_path, "contigs.fasta")
#         with open(out_fa, "w") as f:
#             for i in range(20):
#                 c = dbg.get_longest_contig()
#                 if c is None:
#                     break
#                 print(i, len(c))
#                 f.write(f">contig_{i}\n{c}\n")

#         print(f"WROTE: {out_fa}")

#     except BaseException as e:
#         # traceback کامل تا در CI بفهمیم دقیقاً کجا گیر کرده
#         print("[CODON] ERROR:", repr(e))
#         traceback.print_exc()
#         sys.exit(2)




