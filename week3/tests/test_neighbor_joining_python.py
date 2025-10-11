import sys, os
import time
from phylo_python import neighbor_joining
import numpy as np
base_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(base_dir, "..", "src")
if src_dir not in sys.path:
    sys.path.append(src_dir)

def test_neighbor_joining():
    dist_path = os.path.join(base_dir, "distances.txt")
    
    D = np.loadtxt(dist_path, dtype=float)
    tree = neighbor_joining(D)
    d34 = tree.get_distance(3,4)
    d01 = tree.get_distance(0,1)
    print("NeighborJoining: dist(3,4)=%.3f, dist(0,1)=%.3f" %
          (tree.get_distance(3, 4), tree.get_distance(0, 1)))

if __name__ == "__main__":
    t0 = time.perf_counter()
    for _ in range(1):
        test_neighbor_joining()
    print(f"Runtime: {int((time.perf_counter() - t0)*1000)}ms")
