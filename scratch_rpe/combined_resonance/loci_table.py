"""Where can each driver be stationary in the HEO regime? (a, e, i) loci + a sanity check of the
critical inclinations (the s = j = 0 members must reproduce 46.38/56.06/63.43/69.01/73.15 deg)."""
import math, sys
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from drivers import *   # noqa

A_MIN, A_MAX = 12000.0, 50000.0     # km, HEO / GTO / Molniya / high-altitude debris range of interest


def crit_inclinations():
    out = {}
    for (lab, k, m, s, j, pump) in family():
        if s == 0 and j == 0:
            hits = []
            for i10 in range(0, 1800):
                i = i10 / 10.0
                g0, g1 = G(i, k, m), G(i + 0.1, k, m)
                if g0 == 0 or g0 * g1 < 0:
                    # refine
                    lo, hi = i, i + 0.1
                    for _ in range(40):
                        mid = 0.5 * (lo + hi)
                        if G(lo, k, m) * G(mid, k, m) <= 0:
                            hi = mid
                        else:
                            lo = mid
                    hits.append(round(0.5 * (lo + hi), 2))
            out[lab] = sorted(set(hits))
    return out


def main():
    print("Critical inclinations from the s=j=0 members (deg):")
    for lab, v in crit_inclinations().items():
        print(f"  {lab:>8}: {v}")

    print("\nReachable resonances (a in %.0f-%.0f km) at e=0.70, selected inclinations:" % (A_MIN, A_MAX))
    for i in (50.0, 56.06, 60.0, 63.4, 66.0, 69.0, 73.15, 80.0, 90.0):
        row = []
        for (lab, k, m, s, j, pump) in family():
            if s == 0 and j == 0:
                continue
            a = a_resonant(0.70, i, k, m, s, j)
            if a and A_MIN <= a <= A_MAX:
                row.append(f"{lab}@{a:,.0f}{'*' if pump else ''}")
        print(f"  i={i:6.2f}: " + ("  ".join(row) if row else "(none)"))
    print("  (* = directly pumps e at quadrupole order, k != 0)")

    print("\nMolniya-like state a=26562 km, e=0.70, i=63.4 deg: drivers sorted by |Psi_dot| (deg/day)")
    st = []
    for (lab, k, m, s, j, pump) in family():
        st.append((abs(psi_dot_j2(26562.0, 0.70, 63.4, k, m, s, j)), lab, psi_dot_j2(26562.0, 0.70, 63.4, k, m, s, j), pump))
    for av, lab, v, pump in sorted(st)[:8]:
        print(f"  {lab:>28}  {v:+.4f}  {'(pumps e)' if pump else ''}")
    print(f"  lam(26562, 0.70) = {lam_deg_per_day(26562.0, 0.70):.4f} deg/day;  |Omega_L_dot| = {abs(OMEGA_L_DOT):.4f} deg/day")


if __name__ == '__main__':
    main()
