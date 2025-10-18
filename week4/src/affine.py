# # week4/affine.py  (pure Python)
# NEG_INF = -10**15

# def align(x, y, match=3, mismatch=-3, gap_open=-5, gap_extend=-1):
#     m, n = len(x), len(y)
#     M  = [[NEG_INF]*(n+1) for _ in range(m+1)]
#     Ix = [[NEG_INF]*(n+1) for _ in range(m+1)]  # gap in x (advance j)
#     Iy = [[NEG_INF]*(n+1) for _ in range(m+1)]  # gap in y (advance i)

#     TBM  = [[-1]*(n+1) for _ in range(m+1)]
#     TBIx = [[-1]*(n+1) for _ in range(m+1)]
#     TBIy = [[-1]*(n+1) for _ in range(m+1)]

#     M[0][0] = 0
#     for j in range(1, n+1):
#         Ix[0][j] = gap_open + (j-1)*gap_extend
#         TBIx[0][j] = 1 if j > 1 else 0
#     for i in range(1, m+1):
#         Iy[i][0] = gap_open + (i-1)*gap_extend
#         TBIy[i][0] = 2 if i > 1 else 0

#     def s(a, b):
#         return match if a == b else mismatch

#     for i in range(1, m+1):
#         xi = x[i-1]
#         for j in range(1, n+1):
#             yj = y[j-1]

#             prev, src = max(
#                 (M[i-1][j-1], 0),
#                 (Ix[i-1][j-1], 1),
#                 (Iy[i-1][j-1], 2),
#                 key=lambda t: t[0]
#             )
#             M[i][j] = prev + s(xi, yj)
#             TBM[i][j] = src

#             open_x = M[i][j-1] + gap_open
#             ext_x  = Ix[i][j-1] + gap_extend
#             if open_x >= ext_x:
#                 Ix[i][j] = open_x
#                 TBIx[i][j] = 0
#             else:
#                 Ix[i][j] = ext_x
#                 TBIx[i][j] = 1

#             open_y = M[i-1][j] + gap_open
#             ext_y  = Iy[i-1][j] + gap_extend
#             if open_y >= ext_y:
#                 Iy[i][j] = open_y
#                 TBIy[i][j] = 0
#             else:
#                 Iy[i][j] = ext_y
#                 TBIy[i][j] = 2

#     score, state = max((M[m][n], 0), (Ix[m][n], 1), (Iy[m][n], 2), key=lambda t: t[0])

#     i, j = m, n
#     ax, ay = [], []
#     while i > 0 or j > 0:
#         if state == 0:
#             ax.append(x[i-1]); ay.append(y[j-1])
#             src = TBM[i][j]
#             i -= 1; j -= 1
#             state = src
#         elif state == 1:
#             ax.append('-'); ay.append(y[j-1])
#             src = TBIx[i][j]
#             j -= 1
#             state = 0 if src == 0 else 1
#         else:
#             ax.append(x[i-1]); ay.append('-')
#             src = TBIy[i][j]
#             i -= 1
#             state = 0 if src == 0 else 2

#     ax.reverse(); ay.reverse()
#     return ''.join(ax), ''.join(ay), score


# affine_align_fast.py
# نسخه‌ی بهینه‌شده‌ی الگوریتم Affine Gap Penalty Alignment
# سازگار با Python و Codon

def affine_align(seq1, seq2, match=1, mismatch=-1, gap_open=-5, gap_extend=-1):
    n = len(seq1)
    m = len(seq2)

    NEG_INF = -10**9  # جایگزین float('-inf') چون Codon ممکنه مشکل داشته باشه

    # --- تخصیص ماتریس‌ها (با استفاده از list comprehension سریع‌تر از append) ---
    M = [[0] * (m + 1) for _ in range(n + 1)]
    X = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

    # --- مقداردهی اولیه ---
    M[0][0] = 0
    for i in range(1, n + 1):
        M[i][0] = NEG_INF
        X[i][0] = gap_open + (i - 1) * gap_extend
        Y[i][0] = NEG_INF
    for j in range(1, m + 1):
        M[0][j] = NEG_INF
        X[0][j] = NEG_INF
        Y[0][j] = gap_open + (j - 1) * gap_extend

    # --- حلقه‌ی اصلی ---
    for i in range(1, n + 1):
        s1 = seq1[i - 1]
        Mi = M[i]
        Xi = X[i]
        Yi = Y[i]
        Mprev = M[i - 1]
        Xprev = X[i - 1]
        Yprev = Y[i - 1]

        for j in range(1, m + 1):
            s2 = seq2[j - 1]
            score = match if s1 == s2 else mismatch

            # match/mismatch
            best_prev = Mprev[j - 1]
            if Xprev[j - 1] > best_prev:
                best_prev = Xprev[j - 1]
            if Yprev[j - 1] > best_prev:
                best_prev = Yprev[j - 1]
            Mi[j] = best_prev + score

            # gap در seq1 (یعنی حذف در seq2)
            open_gap = Mi[j - 1] + gap_open
            extend_gap = Xi[j - 1] + gap_extend
            if extend_gap > open_gap:
                Xi[j] = extend_gap
            else:
                Xi[j] = open_gap

            # gap در seq2 (یعنی حذف در seq1)
            open_gap2 = Mprev[j] + gap_open
            extend_gap2 = Yprev[j] + gap_extend
            if extend_gap2 > open_gap2:
                Yi[j] = extend_gap2
            else:
                Yi[j] = open_gap2

    # --- محاسبه‌ی نهایی ---
    final_score = M[n][m]
    if X[n][m] > final_score:
        final_score = X[n][m]
    if Y[n][m] > final_score:
        final_score = Y[n][m]

    return "", "", final_score
