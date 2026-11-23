"""
Statistics used in the MissMap analyses. Small, dependency-light, and tested,
so the numbers in the thesis can be traced to code. scipy is used for
distribution tails when available.

  wilson_ci            proportion with a 95% interval that behaves near 0 and 1
  mcnemar_exact        paired binary comparison of two strategies/models
  cochrans_q           >2 models on the same items
  holm                 multiple-comparison correction
  fleiss_kappa         agreement of >2 raters (e.g. authorities) on categories
  jaccard / mean_pairwise_jaccard   set overlap, run-to-run stability
  bootstrap_ci         percentile CI for any statistic (medians of skewed counts)
"""

import math
import random
from itertools import combinations


def wilson_ci(successes, n, z=1.959964):
    if n == 0:
        return (float("nan"), float("nan"), float("nan"))
    p = successes / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (p, max(0.0, centre - half), min(1.0, centre + half))


def _binom_two_sided(k, n):
    """Exact two-sided binomial test p-value for p=0.5."""
    if n == 0:
        return 1.0
    k = min(k, n - k)
    tail = sum(math.comb(n, i) for i in range(0, k + 1)) / 2 ** n
    return min(1.0, 2 * tail)


def mcnemar_exact(a, b):
    """
    a, b: equal-length sequences of 0/1 (e.g. synonym recovered by model A / B).
    Returns dict(b01, b10, p) where b01 = A missed & B found.
    """
    if len(a) != len(b):
        raise ValueError("paired sequences must be the same length")
    b01 = sum(1 for x, y in zip(a, b) if not x and y)
    b10 = sum(1 for x, y in zip(a, b) if x and not y)
    return {"b01": b01, "b10": b10, "p": _binom_two_sided(b01, b01 + b10)}


def _chi2_sf(x, df):
    try:
        from scipy.stats import chi2
        return float(chi2.sf(x, df))
    except ImportError:  # Wilson-Hilferty approximation
        if x <= 0:
            return 1.0
        z = ((x / df) ** (1 / 3) - (1 - 2 / (9 * df))) / math.sqrt(2 / (9 * df))
        return 0.5 * math.erfc(z / math.sqrt(2))


def cochrans_q(matrix):
    """
    matrix: list of rows (items), each a list of 0/1 per treatment (model).
    Returns dict(Q, df, p). Rows that are all 0 or all 1 carry no information.
    """
    k = len(matrix[0]) if matrix else 0
    if k < 2:
        raise ValueError("need at least two treatments")
    col = [sum(r[j] for r in matrix) for j in range(k)]
    row = [sum(r) for r in matrix]
    n_total = sum(row)
    num = (k - 1) * (k * sum(c * c for c in col) - n_total ** 2)
    den = k * n_total - sum(r * r for r in row)
    if den == 0:
        return {"Q": 0.0, "df": k - 1, "p": 1.0}
    q = num / den
    return {"Q": q, "df": k - 1, "p": _chi2_sf(q, k - 1)}


def holm(pvalues):
    """Holm-Bonferroni adjusted p-values, same order as input."""
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: pvalues[i])
    adj = [0.0] * m
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (m - rank) * pvalues[i]))
        adj[i] = running
    return adj


def fleiss_kappa(table):
    """
    table: list of items, each a list of counts per category (raters per category).
    Every item must have the same number of raters.
    """
    n_items = len(table)
    if n_items == 0:
        return float("nan")
    n_raters = sum(table[0])
    if any(sum(r) != n_raters for r in table) or n_raters < 2:
        raise ValueError("each item needs the same number (>=2) of ratings")
    k = len(table[0])
    p_j = [sum(r[j] for r in table) / (n_items * n_raters) for j in range(k)]
    p_i = [(sum(c * c for c in r) - n_raters) / (n_raters * (n_raters - 1)) for r in table]
    p_bar = sum(p_i) / n_items
    p_e = sum(p * p for p in p_j)
    if p_e == 1:
        return 1.0
    return (p_bar - p_e) / (1 - p_e)


def jaccard(a, b):
    a, b = set(a), set(b)
    u = a | b
    return 1.0 if not u else len(a & b) / len(u)


def mean_pairwise_jaccard(sets):
    sets = [set(s) for s in sets]
    pairs = list(combinations(sets, 2))
    if not pairs:
        return float("nan")
    return sum(jaccard(x, y) for x, y in pairs) / len(pairs)


def bootstrap_ci(values, stat=None, n_boot=10000, alpha=0.05, seed=1):
    """Percentile bootstrap CI. Default statistic is the median."""
    vals = list(values)
    if not vals:
        return (float("nan"), float("nan"), float("nan"))

    def median(v):
        s = sorted(v)
        m = len(s) // 2
        return s[m] if len(s) % 2 else (s[m - 1] + s[m]) / 2

    stat = stat or median
    rng = random.Random(seed)
    boots = sorted(stat([rng.choice(vals) for _ in vals]) for _ in range(n_boot))
    lo = boots[int((alpha / 2) * n_boot)]
    hi = boots[min(n_boot - 1, int((1 - alpha / 2) * n_boot))]
    return (stat(vals), lo, hi)


def precision_recall(predicted, reference):
    p, r = set(predicted), set(reference)
    tp = len(p & r)
    prec = tp / len(p) if p else float("nan")
    rec = tp / len(r) if r else float("nan")
    if p and r and (prec + rec) > 0:
        f1 = 2 * prec * rec / (prec + rec)
    else:
        f1 = float("nan")
    return {"tp": tp, "n_pred": len(p), "n_ref": len(r), "precision": prec, "recall": rec, "f1": f1}
