NEG_INF = -10**15

def align(x, y, match=3, mismatch=-3, gap_open=-5, gap_extend=-1):
    m, n = len(x), len(y)
    M  = [[NEG_INF]*(n+1) for _ in range(m+1)]
    Ix = [[NEG_INF]*(n+1) for _ in range(m+1)]  # gap in x (advance j)
    Iy = [[NEG_INF]*(n+1) for _ in range(m+1)]  # gap in y (advance i)

    TBM  = [[-1]*(n+1) for _ in range(m+1)]
    TBIx = [[-1]*(n+1) for _ in range(m+1)]
    TBIy = [[-1]*(n+1) for _ in range(m+1)]

    M[0][0] = 0
    for j in range(1, n+1):
        Ix[0][j] = gap_open + (j-1)*gap_extend
        TBIx[0][j] = 1 if j > 1 else 0
    for i in range(1, m+1):
        Iy[i][0] = gap_open + (i-1)*gap_extend
        TBIy[i][0] = 2 if i > 1 else 0

    def s(a, b):
        return match if a == b else mismatch

    for i in range(1, m+1):
        xi = x[i-1]
        for j in range(1, n+1):
            yj = y[j-1]

            prev, src = max(
                (M[i-1][j-1], 0),
                (Ix[i-1][j-1], 1),
                (Iy[i-1][j-1], 2),
                key=lambda t: t[0]
            )
            M[i][j] = prev + s(xi, yj)
            TBM[i][j] = src

            open_x = M[i][j-1] + gap_open
            ext_x  = Ix[i][j-1] + gap_extend
            if open_x >= ext_x:
                Ix[i][j] = open_x
                TBIx[i][j] = 0
            else:
                Ix[i][j] = ext_x
                TBIx[i][j] = 1

            open_y = M[i-1][j] + gap_open
            ext_y  = Iy[i-1][j] + gap_extend
            if open_y >= ext_y:
                Iy[i][j] = open_y
                TBIy[i][j] = 0
            else:
                Iy[i][j] = ext_y
                TBIy[i][j] = 2

    score, state = max((M[m][n], 0), (Ix[m][n], 1), (Iy[m][n], 2), key=lambda t: t[0])

    i, j = m, n
    ax, ay = [], []
    while i > 0 or j > 0:
        if state == 0:
            ax.append(x[i-1]); ay.append(y[j-1])
            src = TBM[i][j]
            i -= 1; j -= 1
            state = src
        elif state == 1:
            ax.append('-'); ay.append(y[j-1])
            src = TBIx[i][j]
            j -= 1
            state = 0 if src == 0 else 1
        else:
            ax.append(x[i-1]); ay.append('-')
            src = TBIy[i][j]
            i -= 1
            state = 0 if src == 0 else 2

    ax.reverse(); ay.reverse()
    return ''.join(ax), ''.join(ay), score
