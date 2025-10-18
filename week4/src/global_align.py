# # week4/global_align.py  (pure Python)
# def align(x, y, match=3, mismatch=-3, gap=-2):
#     m, n = len(x), len(y)
#     dp = [[0]*(n+1) for _ in range(m+1)]
#     bt = [['']*(n+1) for _ in range(m+1)]

#     for i in range(1, m+1):
#         dp[i][0] = i * gap
#         bt[i][0] = 'U'
#     for j in range(1, n+1):
#         dp[0][j] = j * gap
#         bt[0][j] = 'L'

#     def s(a, b):
#         return match if a == b else mismatch

#     for i in range(1, m+1):
#         xi = x[i-1]
#         for j in range(1, n+1):
#             yj = y[j-1]
#             diag = dp[i-1][j-1] + s(xi, yj)
#             up   = dp[i-1][j] + gap
#             left = dp[i][j-1] + gap
#             best = diag
#             bt[i][j] = 'D'
#             if up > best:
#                 best = up
#                 bt[i][j] = 'U'
#             if left > best:
#                 best = left
#                 bt[i][j] = 'L'
#             dp[i][j] = best

#     i, j = m, n
#     ax, ay = [], []
#     while i > 0 or j > 0:
#         if i > 0 and j > 0 and bt[i][j] == 'D':
#             ax.append(x[i-1]); ay.append(y[j-1])
#             i -= 1; j -= 1
#         elif i > 0 and bt[i][j] == 'U':
#             ax.append(x[i-1]); ay.append('-')
#             i -= 1
#         else:
#             ax.append('-'); ay.append(y[j-1])
#             j -= 1

#     ax.reverse(); ay.reverse()
#     return ''.join(ax), ''.join(ay), dp[m][n]



def global_align(seq1: str, seq2: str, match=1, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    prev = [j * gap for j in range(m + 1)]
    curr = [0] * (m + 1)

    for i in range(1, n + 1):
        curr[0] = i * gap
        for j in range(1, m + 1):
            diag = prev[j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            delete = prev[j] + gap
            insert = curr[j - 1] + gap
            curr[j] = max(diag, delete, insert)
        prev, curr = curr, [0] * (m + 1)

    score = prev[m]
    return "", "", score


