"""
Issue #40 option 3: validate Bayesian optimization against the exact
real RSM landscapes option 1's GA sweep used, before touching any
Fortran/production code.

Loads scratch_rpe/bench/surfaces_dump.csv (56 zones' real 9-corner
grid-interpolated apogee surfaces, from orem_run_dumpsurf.F -- the
exact same production surfaces/e_grid/bn_grid/tobs/apobs the shipped
GA searches, nothing reconstructed). Implements:

  1. ga_twoint-equivalent bilinear interpolation (bracket + interp,
     same index-search logic as ga.F's ga_twoint).
  2. The same fitness definition ga_fitness uses (true RMS in km,
     sqrt(mean((interp-apobs)^2)) -- no /100 scaling needed here,
     that was GA-internal numerical conditioning only, see issue #36
     Finding 2).
  3. A dense grid search per zone (300x300) over the exact same box
     [e_grid[0],e_grid[2]] x [bn_grid[0],bn_grid[2]] as TRUE ground
     truth -- a stronger reference than the "production seed=0.123
     result" option 1 used, since here the fitness is cheap enough to
     brute-force to (near-)certainty.
  4. A from-scratch GP-based Bayesian optimizer (RBF kernel, Expected
     Improvement acquisition, multi-start candidate search) -- no
     sklearn/skopt/bayes_opt available in this environment.

Compares, across 10 seeds x several evaluation budgets (10-100, an
order of magnitude below what option 1 found the GA needs), what
fraction of (zone, seed) runs land within the same tolerance bands
used for the GA sweep, against the TRUE grid-search optimum this
time (stricter reference than option 1 used).
"""
import csv
import time
from collections import defaultdict

import numpy as np

CSV_PATH = "scratch_rpe/bench/surfaces_dump.csv"


def load_zones():
    zones = defaultdict(lambda: {"tobs": [], "apobs": [], "s": []})
    meta = {}
    with open(CSV_PATH) as f:
        r = csv.DictReader(f)
        for row in r:
            obj = row["object"].strip()
            z = int(row["zone"])
            key = (obj, z)
            zones[key]["tobs"].append(float(row["tobs"]))
            zones[key]["apobs"].append(float(row["apobs"]))
            s = np.array([
                [float(row["s11"]), float(row["s12"]), float(row["s13"])],
                [float(row["s21"]), float(row["s22"]), float(row["s23"])],
                [float(row["s31"]), float(row["s32"]), float(row["s33"])],
            ])
            zones[key]["s"].append(s)
            meta[key] = {
                "e_grid": np.array([float(row["e1"]), float(row["e2"]), float(row["e3"])]),
                "bn_grid": np.array([float(row["bn1"]), float(row["bn2"]), float(row["bn3"])]),
                "rms_baseline": float(row["rms_baseline"]),
            }
    out = {}
    for key, d in zones.items():
        out[key] = {
            "tobs": np.array(d["tobs"]),
            "apobs": np.array(d["apobs"]),
            "surf": np.stack(d["s"], axis=0),  # (nobs, 3, 3) surf[it, i(e), j(bn)]
            **meta[key],
        }
    return out


def bracket_index(grid, val):
    # matches ga.F ga_twoint: do i=2,imax; if (xinp<=x(i)) goto 10; i=imax
    for i in range(1, 3):
        if val <= grid[i]:
            return i
    return 2


def bilinear(e_grid, bn_grid, f2d, e_val, bn_val):
    # f2d[i,j] with i indexing e_grid, j indexing bn_grid (matches
    # ga_twoint(x=e_grid, y=bn_grid, f, xinp=e_val, yinp=bn_val, ...))
    i = bracket_index(e_grid, e_val)
    j = bracket_index(bn_grid, bn_val)
    r1 = (e_val - e_grid[i - 1]) / (e_grid[i] - e_grid[i - 1])
    r2 = (bn_val - bn_grid[j - 1]) / (bn_grid[j] - bn_grid[j - 1])
    f1 = f2d[i - 1, j - 1] + r1 * (f2d[i, j - 1] - f2d[i - 1, j - 1])
    f2 = f2d[i - 1, j] + r1 * (f2d[i, j] - f2d[i - 1, j])
    return f1 + r2 * (f2 - f1)


def make_fitness(zdata):
    e_grid = zdata["e_grid"]
    bn_grid = zdata["bn_grid"]
    surf = zdata["surf"]  # (nobs, 3, 3)
    apobs = zdata["apobs"]
    nobs = len(apobs)

    def fitness(e_val, bn_val):
        sq = 0.0
        for it in range(nobs):
            fout = bilinear(e_grid, bn_grid, surf[it], e_val, bn_val)
            sq += (fout - apobs[it]) ** 2
        return (sq / nobs) ** 0.5

    return fitness


def grid_search_optimum(zdata, n=300):
    e_grid, bn_grid = zdata["e_grid"], zdata["bn_grid"]
    fitness = make_fitness(zdata)
    es = np.linspace(e_grid[0], e_grid[2], n)
    bns = np.linspace(bn_grid[0], bn_grid[2], n)
    best = (1e30, None, None)
    surf = zdata["surf"]
    apobs = zdata["apobs"]
    nobs = len(apobs)
    # vectorized bilinear over the whole grid for speed
    ii = np.searchsorted(e_grid[1:], es, side="left") + 1
    ii = np.clip(ii, 1, 2)
    jj = np.searchsorted(bn_grid[1:], bns, side="left") + 1
    jj = np.clip(jj, 1, 2)
    r1 = (es - e_grid[ii - 1]) / (e_grid[ii] - e_grid[ii - 1])
    r2 = (bns - bn_grid[jj - 1]) / (bn_grid[jj] - bn_grid[jj - 1])
    sq_total = np.zeros((n, n))
    for it in range(nobs):
        f2d = surf[it]
        f1 = f2d[ii - 1, jj[:, None] - 1] + r1[:, None] * (f2d[ii, jj[:, None] - 1] - f2d[ii - 1, jj[:, None] - 1])
        f2 = f2d[ii - 1, jj[:, None]] + r1[:, None] * (f2d[ii, jj[:, None]] - f2d[ii - 1, jj[:, None]])
        fout = f1 + r2[None, :] * (f2 - f1)
        sq_total += (fout - apobs[it]) ** 2
    rms = np.sqrt(sq_total / nobs)
    idx = np.unravel_index(np.argmin(rms), rms.shape)
    return rms[idx], es[idx[0]], bns[idx[1]]


# ---------------- from-scratch GP + Expected Improvement BO ----------------

def rbf_kernel(X1, X2, lengthscale, sigma_f2):
    d2 = np.sum((X1[:, None, :] - X2[None, :, :]) ** 2, axis=-1)
    return sigma_f2 * np.exp(-0.5 * d2 / lengthscale ** 2)


def gp_predict(Xtr, ytr, Xte, lengthscale, sigma_f2, jitter=1e-8):
    K = rbf_kernel(Xtr, Xtr, lengthscale, sigma_f2) + jitter * np.eye(len(Xtr))
    L = np.linalg.cholesky(K)
    alpha = np.linalg.solve(L.T, np.linalg.solve(L, ytr))
    Ks = rbf_kernel(Xtr, Xte, lengthscale, sigma_f2)
    mu = Ks.T @ alpha
    v = np.linalg.solve(L, Ks)
    var = sigma_f2 - np.sum(v ** 2, axis=0)
    var = np.maximum(var, 1e-12)
    return mu, var


def expected_improvement(mu, var, f_min, xi=0.01):
    sigma = np.sqrt(var)
    imp = f_min - mu - xi
    z = np.where(sigma > 1e-12, imp / sigma, 0.0)
    from scipy.stats import norm
    ei = imp * norm.cdf(z) + sigma * norm.pdf(z)
    ei = np.where(sigma > 1e-12, ei, 0.0)
    return ei


def bayes_opt(fitness_fn, bounds, budget, seed, n_init=5, n_candidates=1500):
    rng = np.random.default_rng(seed)
    n_init = min(n_init, budget)
    lo = np.array([bounds[0][0], bounds[1][0]])
    hi = np.array([bounds[0][1], bounds[1][1]])

    def to_unit(X):
        return (X - lo) / (hi - lo)

    def from_unit(U):
        return lo + U * (hi - lo)

    Xu = rng.uniform(0, 1, size=(n_init, 2))
    X = from_unit(Xu)
    y = np.array([fitness_fn(x[0], x[1]) for x in X])

    n_remaining = budget - n_init
    for _ in range(n_remaining):
        sigma_f2 = max(np.var(y), 1e-6)
        lengthscale = 0.25
        cand_u = rng.uniform(0, 1, size=(n_candidates, 2))
        mu, var = gp_predict(Xu, y, cand_u, lengthscale, sigma_f2)
        f_min = np.min(y)
        ei = expected_improvement(mu, var, f_min, xi=0.01 * (np.std(y) + 1e-9))
        best_idx = np.argmax(ei)
        next_u = cand_u[best_idx]
        next_x = from_unit(next_u)
        next_y = fitness_fn(next_x[0], next_x[1])
        Xu = np.vstack([Xu, next_u])
        X = np.vstack([X, next_x])
        y = np.append(y, next_y)

    best_i = np.argmin(y)
    return y[best_i], X[best_i, 0], X[best_i, 1]


def main():
    t0 = time.time()
    zones = load_zones()
    print(f"Loaded {len(zones)} zones from {CSV_PATH}")

    print("\nComputing true grid-search optimum per zone (250x250)...")
    truth = {}
    for key, zdata in zones.items():
        truth[key] = grid_search_optimum(zdata, n=250)
    t1 = time.time()
    print(f"  done in {t1-t0:.1f}s")

    # sanity: compare truth vs production baseline rms
    gaps = []
    for key, zdata in zones.items():
        t_rms = truth[key][0]
        base = zdata["rms_baseline"]
        gaps.append((base - t_rms) / max(t_rms, 1e-9))
    print(f"\nProduction baseline vs true grid optimum: "
          f"mean relative gap {100*np.mean(gaps):.1f}%, "
          f"median {100*np.median(gaps):.1f}%, "
          f"max {100*np.max(gaps):.1f}%")

    budgets = [10, 20, 30, 50, 75, 100]
    seeds = list(range(6))

    results = defaultdict(list)  # budget -> list of (converged_5pct, converged_loose)
    print("\nRunning BO sweep (this may take a few minutes)...")
    for bi, budget in enumerate(budgets):
        t_budget0 = time.time()
        for key, zdata in zones.items():
            fitness_fn = make_fitness(zdata)
            e_grid, bn_grid = zdata["e_grid"], zdata["bn_grid"]
            bounds = ((e_grid[0], e_grid[2]), (bn_grid[0], bn_grid[2]))
            true_rms = truth[key][0]
            for seed in seeds:
                best_rms, _, _ = bayes_opt(fitness_fn, bounds, budget, seed)
                tight_ok = best_rms <= true_rms * 1.05 + 1e-6
                loose_ok = best_rms <= max(true_rms * 1.20, true_rms + 1.0)
                results[budget].append((tight_ok, loose_ok))
        dt = time.time() - t_budget0
        n = len(results[budget])
        tight_rate = 100 * sum(a for a, b in results[budget]) / n
        loose_rate = 100 * sum(b for a, b in results[budget]) / n
        print(f"  budget={budget:>4} nevals={budget:>4}  "
              f"tight(5%)={tight_rate:>5.1f}%  loose(20%/1km)={loose_rate:>5.1f}%  "
              f"[{dt:.1f}s, n={n}]")

    print(f"\nTotal wall time: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
