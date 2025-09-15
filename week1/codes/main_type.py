from dbg_type import DBG
from utils_type import read_data
from python import sys
import os

sys.setrecursionlimit(100000)


if __name__ == "__main__":
    import sys
    argv = sys.argv
    data_path = f"../data/{argv[1]}"
    short1, short2, long1 = read_data(data_path)

    k = 25
    output_dir = f"../data/{argv[1]}"
    output_file_path = f"{output_dir}/contig.fasta"
    dbg = DBG(k=k, data_list=[short1, short2, long1])
    # dbg.show_count_distribution()
    try:
        os.mkdir(output_dir)
    except OSError:
        pass # اگر دایرکتوری از قبل وجود داشته باشد، مشکلی نیست

    with open(output_file_path, 'w') as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            print(i, len(c))
            f.write(f">contig_{i}\n")
            f.write(c + "\n")