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

from dbg_type import DBG
from utils_type import read_data
from python import os, sys  # در Codon بهتره os/sys را از bridge بیاریم

# اختیاری: روی Codon اثر خاصی نداره
# sys.setrecursionlimit(100000)

def resolve_data_path(arg: str) -> str:
    """
    اگر arg مسیر کامل/نسبیِ موجود باشد همان را برمی‌گرداند.
    اگر فقط نام دیتاست باشد (data1, data2, ...)، آن را زیر cwd/data می‌سازد.
    """
    p = os.path.abspath(arg)
    if os.path.isdir(p):
        return p
    base = os.getcwd()  # در CI == week1
    return os.path.abspath(os.path.join(base, "data", arg))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: codon run -plugin seq main_type.py <dataset-path | dataset-name>")
        sys.exit(1)

    data_path = resolve_data_path(sys.argv[1])
    short1, short2, long1 = read_data(data_path)

    k = 25
    dbg = DBG(k=k, data_list=[short1, short2, long1])

    out_fa = os.path.join(data_path, "contigs.fasta")
    with open(out_fa, "w") as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            print(i, len(c))
            f.write(f">contig_{i}\n{c}\n")

    # برای evaluate.sh
    print(f"WROTE: {out_fa}")













