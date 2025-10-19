NEG_INF = -10**15

def align(x, y, match=3, mismatch=-3, gap_open=-5, gap_extend=-1):
    m, n = len(x), len(y)
    M = [[NEG_INF]*(n+1) for _ in range(m+1)]
    Ix = [[NEG_INF]*(n+1) for _ in range(m+1)]  # gap in x (advance j)
    Iy = [[NEG_INF]*(n+1) for _ in range(m+1)]  # gap in y (advance i)

    TBM  = [[-1]*(n+1) for _ in range(m+1)]
    TBIx = [[-1]*(n+1) for _ in range(m+1)]
    TBIy = [[-1]*(n+1) for _ in range(m+1)]

    M[0][0] = 0
    for j in range(1, n+1):
        M[0][j] = 0  # free to skip prefix of y

    Iy[1][0] = gap_open
    TBIy[1][0] = 2
    M[1][0] = NEG_INF
    Ix[1][0] = NEG_INF
    for i in range(2, m+1):
        Iy[i][0] = Iy[i-1][0] + gap_extend
        TBIy[i][0] = 2

    def s(a, b):
        return match if a == b else mismatch

    for i in range(1, m+1):
        xi = x[i-1]
        for j in range(1, n+1):
            yj = y[j-1]

            prevM, srcM = max(
                (M[i-1][j-1], 0),
                (Ix[i-1][j-1], 1),
                (Iy[i-1][j-1], 2),
                key=lambda t: t[0]
            )
            M[i][j] = prevM + s(xi, yj)
            TBM[i][j] = srcM

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

    end_candidates = []
    for j in range(n+1):
        end_candidates.append((M[m][j], j, 0))
        end_candidates.append((Ix[m][j], j, 1))
        end_candidates.append((Iy[m][j], j, 2))
    best_score, j, state = max(end_candidates, key=lambda t: t[0])

    i = m
    ax, ay = [], []
    while i > 0:
        if state == 0:  # M
            ax.append(x[i-1]); ay.append(y[j-1])
            src = TBM[i][j]
            i -= 1; j -= 1
            state = src
        elif state == 1:  # Ix
            ax.append('-'); ay.append(y[j-1])
            src = TBIx[i][j]
            j -= 1
            state = 0 if src == 0 else 1
        else:            # Iy
            ax.append(x[i-1]); ay.append('-')
            src = TBIy[i][j]
            i -= 1
            state = 0 if src == 0 else 2

    ax.reverse(); ay.reverse()
    return ''.join(ax), ''.join(ay), best_score

