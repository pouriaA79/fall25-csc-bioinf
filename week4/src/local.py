# week4/local.py  (pure Python)
def align(x, y, match=3, mismatch=-3, gap=-2):
    m, n = len(x), len(y)
    dp = [[0]*(n+1) for _ in range(m+1)]
    bt = [['']*(n+1) for _ in range(m+1)]

    def s(a, b):
        return match if a == b else mismatch

    best_i, best_j, best_score = 0, 0, 0
    for i in range(1, m+1):
        xi = x[i-1]
        for j in range(1, n+1):
            yj = y[j-1]
            diag = dp[i-1][j-1] + s(xi, yj)
            up   = dp[i-1][j] + gap
            left = dp[i][j-1] + gap
            score = max(0, diag, up, left)
            dp[i][j] = score
            if score == 0:
                bt[i][j] = 'Z'
            elif score == diag:
                bt[i][j] = 'D'
            elif score == up:
                bt[i][j] = 'U'
            else:
                bt[i][j] = 'L'
            if score > best_score:
                best_score = score
                best_i, best_j = i, j

    i, j = best_i, best_j
    ax, ay = [], []
    while i > 0 and j > 0 and dp[i][j] > 0:
        if bt[i][j] == 'D':
            ax.append(x[i-1]); ay.append(y[j-1])
            i -= 1; j -= 1
        elif bt[i][j] == 'U':
            ax.append(x[i-1]); ay.append('-')
            i -= 1
        elif bt[i][j] == 'L':
            ax.append('-'); ay.append(y[j-1])
            j -= 1
        else:
            break

    ax.reverse(); ay.reverse()
    return ''.join(ax), ''.join(ay), best_score
