from typing import List, Tuple

def semi_global_align(seq1: str, seq2: str, match: int = 3, mismatch: int = -3, gap_open: int = -5, gap_extend: int = -1) -> Tuple[str, str, int]:
    n, m = len(seq1), len(seq2)
    dp: List[List[int]] = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = 0
    for j in range(1, m + 1):
        dp[0][j] = 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i - 1] == seq2[j - 1] else mismatch
            gap_up = dp[i - 1][j] + (gap_extend if dp[i - 1][j] != 0 else gap_open)
            gap_left = dp[i][j - 1] + (gap_extend if dp[i][j - 1] != 0 else gap_open)
            dp[i][j] = max(dp[i - 1][j - 1] + s, gap_up, gap_left)

    i, j = n, max(range(m + 1), key=lambda x: dp[n][x])
    a1, a2 = [], []
    while i > 0 and j > 0:
        s = match if seq1[i - 1] == seq2[j - 1] else mismatch
        if dp[i][j] == dp[i - 1][j - 1] + s:
            a1.append(seq1[i - 1])
            a2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif dp[i][j] == dp[i - 1][j] + gap_extend or dp[i][j] == dp[i - 1][j] + gap_open:
            a1.append(seq1[i - 1])
            a2.append('-')
            i -= 1
        else:
            a1.append('-')
            a2.append(seq2[j - 1])
            j -= 1

    return ''.join(reversed(a1)), ''.join(reversed(a2)), dp[n][max(range(m + 1), key=lambda x: dp[n][x])]
