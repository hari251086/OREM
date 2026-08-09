"""
Issue #40 option 1: aggregate scratch_rpe/bench/gasweep_results.csv
(56 real zones x 42 pop/gen combos x 10 seeds each) into a per-combo
convergence-reliability table, to find the minimum GA evaluation
budget that reliably reproduces the production pop=20/gen=500 optimum
on real RSM landscapes.
"""
import csv
from collections import defaultdict

rows = []
with open("scratch_rpe/bench/gasweep_results.csv") as f:
    r = csv.DictReader(f)
    for row in r:
        row["object"] = row["object"].strip()
        for k in ("nzone",):
            row[k] = int(row[k])
        for k in ("zone",):
            row[k] = int(row[k])
        for k in ("pop", "gen", "nevals", "nconverged"):
            row[k] = int(row[k])
        for k in ("e_span", "bn_span", "rms_baseline", "mean_rms",
                   "min_rms", "max_rms", "mean_bn", "std_bn"):
            row[k] = float(row[k])
        rows.append(row)

nzones = len(set((r["object"], r["zone"]) for r in rows))
print(f"Loaded {len(rows)} rows, {nzones} zones\n")

# --- Per (pop,gen) aggregate across all zones ---
agg = defaultdict(list)
for r in rows:
    agg[(r["pop"], r["gen"])].append(r)

print(f"{'pop':>4} {'gen':>5} {'nevals':>7} {'conv_rate':>10} "
      f"{'zones_100%':>11} {'zones_0%':>9}")
combo_stats = []
for (pop, gen), rs in sorted(agg.items(), key=lambda kv: kv[0][0]*(kv[0][1]+1)):
    nevals = pop * (gen + 1)
    total_conv = sum(r["nconverged"] for r in rs)
    total_possible = len(rs) * 10
    conv_rate = total_conv / total_possible
    zones_full = sum(1 for r in rs if r["nconverged"] == 10)
    zones_zero = sum(1 for r in rs if r["nconverged"] == 0)
    combo_stats.append((pop, gen, nevals, conv_rate, zones_full, zones_zero, len(rs)))
    print(f"{pop:>4} {gen:>5} {nevals:>7} {conv_rate*100:>9.1f}% "
          f"{zones_full:>6}/{len(rs):<4} {zones_zero:>6}/{len(rs):<4}")

# --- Minimum-eval combo reaching >=90% / >=95% / 100% convergence rate ---
print()
for thresh in (0.90, 0.95, 1.00):
    candidates = [c for c in combo_stats if c[3] >= thresh]
    if candidates:
        best = min(candidates, key=lambda c: c[2])
        print(f">= {thresh*100:.0f}% conv_rate: min-eval combo pop={best[0]} "
              f"gen={best[1]} (nevals={best[2]}, actual conv_rate={best[3]*100:.1f}%)")
    else:
        print(f">= {thresh*100:.0f}% conv_rate: NOT REACHED by any combo tested")

# --- Production baseline evals for reference ---
prod_nevals = 20 * (500 + 1)
print(f"\nProduction baseline (pop=20, gen=500): nevals={prod_nevals}")

# --- pop=4 "range-invariant seed artifact" check: does std_bn stay
#     near-zero (seeds agree) while conv_rate stays near-zero (wrong
#     answer) -- the specific issue #12 signature, as opposed to normal
#     under-convergence (seeds disagree, gradually improving) ---
print("\n--- pop=4 signature check (mean over 56 zones, by gen) ---")
print(f"{'gen':>5} {'mean_std_bn':>12} {'mean_bn_span':>13} "
      f"{'std_bn/span':>12} {'conv_rate':>10}")
for gen in sorted(set(r["gen"] for r in rows)):
    rs = [r for r in rows if r["pop"] == 4 and r["gen"] == gen]
    mean_std_bn = sum(r["std_bn"] for r in rs) / len(rs)
    mean_span = sum(r["bn_span"] for r in rs) / len(rs)
    conv = sum(r["nconverged"] for r in rs) / (len(rs) * 10)
    print(f"{gen:>5} {mean_std_bn:>12.4f} {mean_span:>13.4f} "
          f"{mean_std_bn/mean_span:>12.4f} {conv*100:>9.1f}%")

print("\n--- same check at pop=20 (production pop) for contrast ---")
for gen in sorted(set(r["gen"] for r in rows)):
    rs = [r for r in rows if r["pop"] == 20 and r["gen"] == gen]
    mean_std_bn = sum(r["std_bn"] for r in rs) / len(rs)
    mean_span = sum(r["bn_span"] for r in rs) / len(rs)
    conv = sum(r["nconverged"] for r in rs) / (len(rs) * 10)
    print(f"{gen:>5} {mean_std_bn:>12.4f} {mean_span:>13.4f} "
          f"{mean_std_bn/mean_span:>12.4f} {conv*100:>9.1f}%")

# --- Per-zone worst case at a few candidate budgets ---
print("\n--- Per-zone conv_rate at candidate combos (worst zones shown) ---")
for pop, gen in [(8, 40), (10, 75), (14, 75), (20, 40), (20, 75), (20, 125)]:
    rs = [r for r in rows if r["pop"] == pop and r["gen"] == gen]
    rs_sorted = sorted(rs, key=lambda r: r["nconverged"])
    nevals = pop * (gen + 1)
    worst = rs_sorted[:5]
    print(f"\npop={pop} gen={gen} (nevals={nevals}):")
    for r in worst:
        print(f"  {r['object']:>8} Z{r['zone']} nconv={r['nconverged']:>2}/10 "
              f"bn_span={r['bn_span']:>8.2f} rms_base={r['rms_baseline']:>8.3f}")
