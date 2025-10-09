import sys, os
base_dir = os.path.dirname(os.path.abspath(__file__))

src_dir = os.path.join(base_dir, "..", "src")

if src_dir not in sys.path:
    sys.path.append(src_dir)


import time
from phylo_python import upgma

def test_distances():
    # D = [
    #     [0.0, 5.0, 9.0, 9.0, 8.0],
    #     [5.0, 0.0, 10.0, 10.0, 9.0],
    #     [9.0, 10.0, 0.0, 8.0, 7.0],
    #     [9.0, 10.0, 8.0, 0.0, 3.0],
    #     [8.0, 9.0, 7.0, 3.0, 0.0],
    # ]
    D = [[abs(i - j) * 1.1 for j in range(100)] for i in range(100)]

    tree = upgma(D)
    print("Distance(0,0) =", tree.get_distance(0, 0))
    print("Distance(0,1) =", tree.get_distance(0, 1))
    print("Distance(3,4) =", tree.get_distance(3, 4))
    print("Topological distance(3,4) =", tree.get_distance(3, 4, topological=True))

if __name__ == "__main__":
    t0 = time.perf_counter()
    for _ in range(20):
        test_distances()
    # test_distances()
    print(f"Runtime: {int((time.perf_counter() - t0)*1000)}ms")


