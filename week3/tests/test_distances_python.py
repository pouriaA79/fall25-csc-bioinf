import sys, os
import time
from phylo_python import upgma
import numpy as np
base_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(base_dir, "..", "src")
if src_dir not in sys.path:
    sys.path.append(src_dir)

def test_distances():
    dist_path = os.path.join(base_dir, "distances.txt")

    D = np.loadtxt(dist_path, dtype=float)
    tree = upgma(D)
    print("Distance(0,0) =", tree.get_distance(0, 0))
    print("Distance(0,1) =", tree.get_distance(0, 1))
    print("Distance(3,4) =", tree.get_distance(3, 4))
    print("Topological distance(3,4) =", tree.get_distance(3, 4, topological=True))

if __name__ == "__main__":
    t0 = time.perf_counter()
    for _ in range(1):
        test_distances()
    print(f"Runtime: {int((time.perf_counter() - t0)*1000)}ms")
