# Scipy's L-BFGS-B, run on many independent problems in lockstep (ADR-019).
#
# PyDESeq2 fits the dispersion of every gene, and the log fold changes of the
# genes whose IRLS diverges, with scipy's L-BFGS-B, one gene at a time. The
# vectorised per-study fit ('vec_deseq2.py') needs the same answers for all
# genes at once. This module drives scipy's compiled optimiser
# ('scipy.optimize._lbfgsb.setulb') exactly as scipy's own
# 'minimize(method="L-BFGS-B")' drives it, one state per problem, and answers
# the optimiser's requests for f and the gradient for all problems in one
# array evaluation. The optimiser's steps, line searches and stopping decisions
# are scipy's own, so each problem ends where PyDESeq2's per-gene call ends;
# only the objective is batched.
#
# 'setulb' is a private function and its signature has changed between scipy
# releases. 'available()' runs a small problem through the lockstep driver and
# compares it with 'scipy.optimize.minimize'; the engine is offered only when
# they agree. The call below matches scipy 1.15-1.17.

import numpy as np
from scipy.special import gammaln

try:
    from scipy.optimize import _lbfgsb
    _SETULB = _lbfgsb.setulb
except Exception:  # noqa: BLE001 -- a private module: any failure means unavailable
    _SETULB = None

M_CORR = 10
FACTR = 2.2204460492503131e-09 / np.finfo(float).eps   # scipy's default ftol as factr
PGTOL = 1e-5
MAXLS = 20
MAXITER = 15000
MAXFUN = 15000

_AVAILABLE = None


def objective(y, s, X, disp, beta, min_mu, ridge):
    """NB negative log-likelihood and its gradient for the genes in the columns of y.

    PyDESeq2's bounded fallback of 'irls_solver' minimises
    ``nll(y, max(s * exp(X beta), min_mu), disp) + 0.5 * ridge * sum(beta ** 2)``.

    Parameters
    ----------
    y : ndarray (n, G)
    s : ndarray (n,)
    X : ndarray (n, p)
    disp : ndarray (G,)
    beta : ndarray (G, p)

    Returns
    -------
    (ndarray (G,), ndarray (G, p))
    """
    mu = np.maximum(s[:, None] * np.exp(X @ beta.T), min_mu)             # (n, G)
    a = 1.0 / disp
    f = (a * np.log(disp) - (gammaln(y + a) - gammaln(y + 1.0) - gammaln(a))
         + (y + a) * np.log(a + mu) - y * np.log(mu)).sum(axis=0) \
        + 0.5 * ridge * (beta ** 2).sum(axis=1)
    grad = (-(X.T @ y) + X.T @ ((a + y) * mu / (a + mu))).T + ridge * beta
    return f, grad


def lbfgsb_lockstep(evaluate, x0, lo, hi):
    """Scipy's L-BFGS-B (default options, one box ``lo..hi`` on every coordinate) on G problems.

    Parameters
    ----------
    evaluate : callable(idx, x) -> (f, g)
        Answers the optimiser's requests: for the problems ``idx`` at the points
        ``x`` (k, p), the objective (k,) and its gradient (k, p).
    x0 : ndarray (G, p)
        Starting points, clipped into the box.
    lo, hi : float

    Returns
    -------
    (ndarray (G, p), ndarray (G,) of bool)
        The solutions and whether the optimiser reported convergence for each.
    """
    if _SETULB is None:
        raise RuntimeError("scipy.optimize._lbfgsb.setulb is not available")
    n_prob, p = x0.shape
    m = M_CORR
    low = np.full(p, lo, dtype=np.float64)
    upp = np.full(p, hi, dtype=np.float64)
    nbd = np.full(p, 2, dtype=np.int32)
    x = np.ascontiguousarray(np.clip(x0, lo, hi), dtype=np.float64)
    f = np.zeros(n_prob)
    g = np.zeros((n_prob, p))
    wa = [np.zeros(2 * m * p + 5 * p + 11 * m * m + 8 * m) for _ in range(n_prob)]
    iwa = [np.zeros(3 * p, dtype=np.int32) for _ in range(n_prob)]
    task = [np.zeros(2, dtype=np.int32) for _ in range(n_prob)]
    ln_task = [np.zeros(2, dtype=np.int32) for _ in range(n_prob)]
    lsave = [np.zeros(4, dtype=np.int32) for _ in range(n_prob)]
    isave = [np.zeros(44, dtype=np.int32) for _ in range(n_prob)]
    dsave = [np.zeros(29) for _ in range(n_prob)]
    n_it = np.zeros(n_prob, dtype=int)
    nfev = np.zeros(n_prob, dtype=int)
    status = np.full(n_prob, -1)
    active = list(range(n_prob))
    while active:
        need = []
        for k in active:
            while True:
                _SETULB(m, x[k], low, upp, nbd, float(f[k]), g[k], FACTR, PGTOL, wa[k],
                        iwa[k], task[k], lsave[k], isave[k], dsave[k], MAXLS, ln_task[k])
                t = task[k][0]
                if t == 3:                       # the optimiser wants f and g at x[k]
                    need.append(k)
                    break
                if t == 1:                       # a new iterate
                    n_it[k] += 1
                    if n_it[k] >= MAXITER:
                        task[k][0], task[k][1] = 5, 504
                    elif nfev[k] > MAXFUN:
                        task[k][0], task[k][1] = 5, 502
                    continue
                status[k] = 0 if t == 4 else (1 if (nfev[k] > MAXFUN or n_it[k] >= MAXITER) else 2)
                break
        if need:
            idx = np.array(need)
            f[idx], g[idx] = evaluate(idx, x[idx])
            nfev[idx] += 1
        active = need
    return x, status == 0


def lbfgsb_batch(y, s, X, disp, beta_init, min_mu=0.5, lo=-30.0, hi=30.0, ridge=1e-6):
    """The bounded NB-GLM fallback of PyDESeq2's ``irls_solver`` for the genes in y's columns.

    Parameters
    ----------
    y : ndarray (n, G)
    s : ndarray (n,)
    X : ndarray (n, p)
    disp : ndarray (G,)
    beta_init : ndarray (G, p)

    Returns
    -------
    (ndarray (G, p), ndarray (G,) of bool)
    """
    return lbfgsb_lockstep(
        lambda idx, b: objective(y[:, idx], s, X, disp[idx], b, min_mu, ridge),
        beta_init, lo, hi)


def available():
    """True when the lockstep driver reproduces ``scipy.optimize.minimize``.

    Runs once per process on a small bounded problem with a known solution and
    compares the result with scipy's own L-BFGS-B. False when scipy's private
    interface is absent or has changed, in which case the vectorised engine is
    not offered.
    """
    global _AVAILABLE
    if _AVAILABLE is None:
        _AVAILABLE = _self_check()
    return _AVAILABLE


def _self_check():
    if _SETULB is None:
        return False
    try:
        from scipy.optimize import minimize
        rng = np.random.default_rng(0)
        a = rng.uniform(0.5, 3.0, size=(6, 3))
        c = rng.uniform(-2.0, 2.0, size=(6, 3))

        def evaluate(idx, x):
            d = x - c[idx]
            return (a[idx] * d ** 2).sum(axis=1) + 0.1 * (x ** 4).sum(axis=1), \
                2 * a[idx] * d + 0.4 * x ** 3

        x0 = np.zeros((6, 3))
        x, ok = lbfgsb_lockstep(evaluate, x0, -1.0, 1.0)
        for k in range(6):
            ref = minimize(
                lambda z, k=k: tuple(v[0] for v in evaluate(np.array([k]), z[None])),
                x0[k], jac=True, method='L-BFGS-B', bounds=[(-1.0, 1.0)] * 3)
            if not (ok[k] and np.array_equal(x[k], ref.x)):
                return False
        return True
    except Exception:  # noqa: BLE001 -- any failure of a private API means unavailable
        return False
