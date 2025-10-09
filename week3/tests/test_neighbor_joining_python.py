import sys, os

base_dir = os.path.dirname(os.path.abspath(__file__))

src_dir = os.path.join(base_dir, "..", "src")

if src_dir not in sys.path:
    sys.path.append(src_dir)



import time
from phylo_python import neighbor_joining

def test_neighbor_joining():
    # D = [
    #     [0.0, 5.0, 9.0, 9.0, 8.0],
    #     [5.0, 0.0, 10.0, 10.0, 9.0],
    #     [9.0, 10.0, 0.0, 8.0, 7.0],
    #     [9.0, 10.0, 8.0, 0.0, 3.0],
    #     [8.0, 9.0, 7.0, 3.0, 0.0],
    # ]
    D = [[abs(i - j) * 1.1 for j in range(100)] for i in range(100)]

    tree = neighbor_joining(D)
    d34 = tree.get_distance(3,4)
    d01 = tree.get_distance(0,1)
    print("NeighborJoining: dist(3,4)=%.3f, dist(0,1)=%.3f" %
          (tree.get_distance(3, 4), tree.get_distance(0, 1)))

if __name__ == "__main__":
    t0 = time.perf_counter()
    # test_neighbor_joining()
    for _ in range(20):
        test_neighbor_joining()
    print(f"Runtime: {int((time.perf_counter() - t0)*1000)}ms")

