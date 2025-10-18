from typing import List, Tuple

def global_align(seq1: str, seq2: str, match: int = 3, mismatch: int = -3, gap: int = -2) -> Tuple[str, str, int]:
    n, m = len(seq1), len(seq2)
    dp: List[List[int]] = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + gap
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i - 1] == seq2[j - 1] else mismatch
            dp[i][j] = max(dp[i - 1][j - 1] + s, dp[i - 1][j] + gap, dp[i][j - 1] + gap)

    i, j = n, m
    a1, a2 = [], []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch):
            a1.append(seq1[i - 1])
            a2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + gap:
            a1.append(seq1[i - 1])
            a2.append('-')
            i -= 1
        else:
            a1.append('-')
            a2.append(seq2[j - 1])
            j -= 1

    return ''.join(reversed(a1)), ''.join(reversed(a2)), dp[n][m]
