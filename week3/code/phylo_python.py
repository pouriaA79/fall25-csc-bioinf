from math import inf

# ============================================================
# TreeNode
# ============================================================
class TreeNode:
    def __init__(self, children=None, distances=None, index=-1):
        if children is None:
            children = []
        if distances is None:
            distances = []
        assert len(children) == len(distances)
        self.children = children
        self.distances = distances
        self.index = index

    def is_leaf(self):
        return len(self.children) == 0

    def _path_to(self, target, path):
        if self.is_leaf():
            if self.index == target:
                path.append((self, -1, 0.0))
                return True
            return False
        for ci, child in enumerate(self.children):
            if child._path_to(target, path):
                path.append((self, ci, self.distances[ci]))
                return True
        return False


# ============================================================
# Tree
# ============================================================
class Tree:
    def __init__(self, root):
        self.root = root
        self.leaves = []
        self._collect_leaves(root)
        self.leaves.sort(key=lambda n: n.index)
        # Ensure leaves are indexed 0..n-1
        for i, leaf in enumerate(self.leaves):
            assert leaf.index == i

    def _collect_leaves(self, node):
        if node.is_leaf():
            self.leaves.append(node)
        else:
            for ch in node.children:
                self._collect_leaves(ch)

    def get_distance(self, i1, i2, topological=False):
        if i1 == i2:
            return 0.0

        def trace_path(node, target, path):
            if node.is_leaf():
                if node.index == target:
                    path.append((node, 0.0))
                    return True
                return False
            for ci, child in enumerate(node.children):
                if trace_path(child, target, path):
                    path.append((node, node.distances[ci]))
                    return True
            return False

        path1 = []
        path2 = []

        assert trace_path(self.root, i1, path1)
        assert trace_path(self.root, i2, path2)

        path1.reverse()
        path2.reverse()

        k = 0
        while k < len(path1) and k < len(path2) and path1[k][0] is path2[k][0]:
            k += 1

        def branch_length(path, start):
            total = 0.0
            for j in range(start, len(path)):
                total += (1.0 if topological else path[j][1])
            return total

        return branch_length(path1, k - 1) + branch_length(path2, k - 1)


# ============================================================
# UPGMA
# ============================================================
def upgma(D):
    n = len(D)
    if n == 0:
        raise ValueError("Empty distance matrix")

    clusters = [[i] for i in range(n)]
    heights = [0.0 for _ in range(n)]
    nodes = [TreeNode(index=i) for i in range(n)]
    active = [True for _ in range(n)]

    def cluster_distance(a, b):
        total = 0.0
        count = 0
        for i in clusters[a]:
            for j in clusters[b]:
                total += D[i][j]
                count += 1
        return total / float(count)

    active_count = n
    while active_count > 1:
        best_i, best_j, best_d = -1, -1, 1e9
        for i in range(n):
            if not active[i]:
                continue
            for j in range(i + 1, n):
                if not active[j]:
                    continue
                d = cluster_distance(i, j)
                if d < best_d:
                    best_d = d
                    best_i = i
                    best_j = j

        new_height = best_d / 2.0

        left_len = new_height - heights[best_i]
        right_len = new_height - heights[best_j]
        if left_len < 0.0:
            left_len = 0.0
        if right_len < 0.0:
            right_len = 0.0

        parent = TreeNode([nodes[best_i], nodes[best_j]],
                          [left_len, right_len],
                          -1)

        nodes[best_i] = parent
        heights[best_i] = new_height
        clusters[best_i].extend(clusters[best_j])
        active[best_j] = False
        active_count -= 1

    for i in range(n):
        if active[i]:
            return Tree(nodes[i])

    raise RuntimeError("No active cluster found in UPGMA")


# ============================================================
# Neighbor Joining (Saitou & Nei, 1987)
# ============================================================
def neighbor_joining(dist):
    n0 = len(dist)
    assert n0 > 1 and all(len(r) == n0 for r in dist)

    active_ids = [i for i in range(n0)]
    nodes = [TreeNode(index=i) for i in range(n0)]

    def dij(i, j):
        return dist[active_ids[i]][active_ids[j]]

    while len(active_ids) > 2:
        m = len(active_ids)
        r = [0.0 for _ in range(m)]
        for i in range(m):
            for j in range(m):
                if i != j:
                    r[i] += dij(i, j)

        best_i, best_j, best_q = -1, -1, inf
        for i in range(m):
            for j in range(i + 1, m):
                q = (m - 2.0) * dij(i, j) - r[i] - r[j]
                if q < best_q:
                    best_q, best_i, best_j = q, i, j

        di_j = dij(best_i, best_j)
        li = 0.5 * di_j + (r[best_i] - r[best_j]) / (2.0 * (m - 2.0))
        lj = di_j - li
        li = max(0.0, li)
        lj = max(0.0, lj)

        u = TreeNode([nodes[active_ids[best_i]],
                      nodes[active_ids[best_j]]],
                     [li, lj],
                     -1)
        new_id = len(nodes)
        nodes.append(u)

        new_row = []
        for k in range(m):
            if k == best_i or k == best_j:
                continue
            duk = 0.5 * (dij(best_i, k) + dij(best_j, k) - di_j)
            new_row.append(duk)

        needed = len(nodes)
        if needed > len(dist):
            for row in dist:
                row.extend([0.0] * (needed - len(row)))
            for _ in range(needed - len(dist)):
                dist.append([0.0] * needed)

        survivors = []
        for idx in range(m):
            if idx not in (best_i, best_j):
                survivors.append(active_ids[idx])

        for p, oldpos in enumerate(range(m)):
            if oldpos in (best_i, best_j):
                continue
            survivor_id = active_ids[oldpos]
            dval = new_row[len(survivors) - (m - 2)]
            dist[new_id][survivor_id] = dval
            dist[survivor_id][new_id] = dval

        active_ids = survivors + [new_id]

    a_id, b_id = active_ids[0], active_ids[1]
    dab = dist[a_id][b_id]
    root = TreeNode([nodes[a_id], nodes[b_id]], [dab / 2.0, dab / 2.0], -1)
    return Tree(root)

