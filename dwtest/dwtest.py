import numpy as np
import pasty
import warnings


def _pan(a, m, n=15, c=0):
    nu = np.argmax(a > a[0])
    if a[1] > a[m]:
        h = m
        k = -1
        i = 1
    else:
        h = 1
        k = 1
        i = m
    x = a[0]
    for nu in range(h, i+k, k):
        if a[nu] >= x:
            if nu == h and c <= 0:
                gradsol = 0
                return gradsol

            if k == 1:
                nu = nu - 1
            h = m - nu
            if c == 0:
                y = h - nu
            else:
                y = c * (a[1] - a[m])

            if y >= 0:
                d = 2
                h = nu
                k = -k
                j1 = 0
                j2 = 2
                j3 = 3
                j4 = 1
            else:
                d  = -2
                nu = nu + 1
                j1 = m - 2
                j2 = m - 1
                j3 = m + 1
                j4 = m
            pin = np.pi / 2 / n
            sum = 0.5 * (k + 1)
            sgn = k / n
            n2 = n + n - 1

            # first integrals
            for l1 in range(h-2*int(h/2), 0-1, -1):
                for l2 in range(j2, nu+d, d):
                    sum1 = a[j4]
                    prod = a[l2]
                    u = 0.5 * (sum1 + prod)
                    v = 0.5 * (sum1 - prod)
                    sum1 = 0
                    for i in range(1, n2+2, 2):
                        y = u - v * np.cos(i * pin)
                        num = y - x
                        prod = np.exp(-c / num)
                        for k in range(1, j1+1):
                            prod = prod * num / (y - a[k])
                        for k in range(j3, m+1):
                            prod = prod * num / (y - a[k])
                        sum1 = sum1 + np.sqrt(np.abs(prod))
                    sgn = -sgn
                    sum = sum + sgn * sum1
                    j1 = j1 + d
                    j3 = j3 + d
                    j4 = j4 + d

                if d == 2:
                    j3 = j3 - 1
                else:
                    j1 = j1 + 1
                j2 = 0
                nu = 0

            gradsol = sum
            return gradsol[0]

    if c >= 0:
        gradsol = 1
        return gradsol
        
    return 1


# altervative: ("greater", "two_sided", "less")
def dwtest(formula, data, order_by = None, alternative = "greater", exact = None, tol = 1e-10, iterations=15):    
    # ensure alternative to be legal
    assert alternative in ("greater", "two_sided", "less"), "argument \"alternative\" should be one of \"greater\", \"two_sided\", \"less\""

    result = smf.ols(formula, data=data).fit()
    y, X = patsy.dmatrices(formula, data, return_type='matrix')
    
    n = result.nobs
    if exact is None:
        exact = (n < 100)
    k = result.df_model
    
    res = result.resid
    dw = sum(np.diff(res) ** 2) / sum(res ** 2)
    q, r = np.linalg.qr(X)
    Q1 = np.linalg.inv(r.T @ r)
    
    if n < 3:
        warnings.warn('not enough observations for computing a p value, set to 1', UserWarning)
        p_value = 1
        
    else:
        if exact:
            A = np.diag(np.concatenate([[1],np.repeat(2, n-2),[1]]))
            row = np.arange(n-1)
            idx = np.vstack([row, row+1]).transpose()
            idx = np.vstack([idx, idx[:, [1,0]]]).astype(int)
            for i, j in idx:
                A[i, j] = -1
            MA = np.diag(np.repeat(1, n)) - X @ Q1 @ X.T
            MA = MA @ A
            ev = np.linalg.eigvals(MA)[:int(n-k-1)]
            if np.any(np.array([complex(e).imag for e in ev]) > tol):
                warnings.warn('imaginary parts of eigenvalues discarded', UserWarning)
            ev = np.array([complex(e).real for e in ev])
            ev = ev[ev > tol]
            pdw = _pan(np.concatenate([np.array([[dw]]), np.array(sorted(ev)).reshape(-1, 1)]), ev.shape[0], iterations)
            print(pdw)
            if alternative == 'two_sided':
                pval = 2 * min(pdw, 1 - pdw)
            elif alternative == 'less':
                pval = 1 - pdw
            elif alternative == 'greater':
                pval = pdw

            if np.isnan(pdw) or pdw > 1 or pdw < 0:
                warnings.warn('exact p value cannot be computed (not in [0,1]), approximate p value will be used', UserWarning)
                exact = False

        if not exact:
            if n < max(5, k):
                warnings.warn('not enough observations for computing an approximate p value, set to 1', UserWarning)
                pval = 1
            else:
                AX = (-pd.DataFrame(X).shift(periods=-1) + 2*pd.DataFrame(X) - pd.DataFrame(X).shift(periods=1)).values
                AX[0] = X[0] - X[1]
                AX[int(n-1)] = X[int(n-1)] - X[int(n-2)]
                XAXQ = X.T @ AX @ Q1
                P = 2 * (n-1) - np.sum(np.diag(XAXQ))
                Q = 2 * (3*n - 4) - 2 * np.sum(np.diag(AX.T @ AX @ Q1)) + np.sum(np.diag(XAXQ @ XAXQ))
                dmean = P / (n-k)
                dvar = 2 / ((n-k) * (n-k+2)) * (Q - P*dmean)
                if alternative == 'two_sided':
                    pval = 2 * stats.norm.sf(np.abs(dw-dmean), 0, np.sqrt(dvar))
                elif alternative == 'less':
                    pval = stats.norm.sf(dw, dmean, np.sqrt(dvar))
                elif alternative == 'greater':
                    pval = stats.norm.cdf(dw, dmean, np.sqrt(dvar))
        
    return dw, pval