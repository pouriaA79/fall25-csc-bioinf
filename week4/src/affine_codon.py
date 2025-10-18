from typing import List, Tuple

def affine_align(seq1: str, seq2: str, match: float = 3.0, mismatch: float = -3.0, gap_open: float = -5.0, gap_extend: float = -1.0) -> Tuple[str, str, float]:
    n, m = len(seq1), len(seq2)
    M: List[List[float]] = [[0.0] * (m + 1) for _ in range(n + 1)]
    X: List[List[float]] = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
    Y: List[List[float]] = [[float('-inf')] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        X[i][0] = gap_open + (i - 1) * gap_extend
        M[i][0] = X[i][0]
    for j in range(1, m + 1):
        Y[0][j] = gap_open + (j - 1) * gap_extend
        M[0][j] = Y[0][j]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i - 1] == seq2[j - 1] else mismatch
            X[i][j] = max(M[i - 1][j] + gap_open, X[i - 1][j] + gap_extend)
            Y[i][j] = max(M[i][j - 1] + gap_open, Y[i][j - 1] + gap_extend)
            M[i][j] = max(M[i - 1][j - 1] + s, X[i][j], Y[i][j])

    i, j = n, m
    a1, a2 = [], []
    while i > 0 and j > 0:
        s = match if seq1[i - 1] == seq2[j - 1] else mismatch
        if M[i][j] == M[i - 1][j - 1] + s:
            a1.append(seq1[i - 1])
            a2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif M[i][j] == X[i][j]:
            a1.append(seq1[i - 1])
            a2.append('-')
            i -= 1
        else:
            a1.append('-')
            a2.append(seq2[j - 1])
            j -= 1

    return ''.join(reversed(a1)), ''.join(reversed(a2)), M[n][m]
