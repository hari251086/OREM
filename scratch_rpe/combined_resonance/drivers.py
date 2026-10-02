"""Combined resonant-azimuth framework (issue #56 follow-up): slow-angle family and J2 rates.

Quadrupole lunisolar model (Celletti, Gales, Pucacco & Rosengren 2016, full text): after averaging
over the satellite mean anomaly (a is then an exact constant) the slow angles are
        Psi = k*omega + m*Omega + s*Omega_L        k in {2, 0, -2}  (= 2-2p),  m, s integers
with Omega_L the (regressing) lunar node. Because the disturbing function is even in Psi (cos Psi),
k=-2 is the same as k=+2 with (m, s) negated, so the family is k = 2 (m in -2..2, s in -2..2) and k = 0
(m in 1..2, s in -2..2). Semi-secular Sun members add -j*n_S (j integer); Wang & Gurfil's solar apsidal
U-turn is the member 2*omega + 2*Omega - 2*lambda_S, i.e. (k, m, s, j) = (2, 2, 0, 2).

J2 secular rates (standard):
        omega_dot = lam * (5 c^2 - 1) / 2 ,   Omega_dot = -lam * c ,   c = cos i
        lam(a, e) = 1.5 * J2 * n * (R/p)^2 ,  p = a (1 - e^2)
so   Psi_dot = lam * G(i) + s * Omega_L_dot - j * n_S ,   G = k (5c^2-1)/2 - m c.
Resonance (stationary Psi): lam * G = j*n_S - s*Omega_L_dot.  As drag shrinks a, lam ~ a^(-7/2) grows
and Psi_dot sweeps through zero -> U-turn of the angle -> monotonic e (quadrupole: de/dt ~ -dR/dpsi * k).
"""
import math

MU = 398600.4415
RE = 6378.1363
J2 = 1.08263e-3
OMEGA_L_DOT = -1934.136261 / 36525.0       # deg/day, regression of the lunar node (ecliptic)
OMEGA_L0 = 125.04452                       # deg at J2000
N_SUN = 360.0 / 365.2422                   # deg/day, Sun's mean motion (semi-secular term)
RAD = math.pi / 180.0


def lam_deg_per_day(a, e):
    """J2 drift scale lam = 1.5 J2 n (R/p)^2 in deg/day."""
    n = math.sqrt(MU / a ** 3) * 86400.0 / RAD            # deg/day
    p = a * (1.0 - e * e)
    return 1.5 * J2 * n * (RE / p) ** 2


def G(i_deg, k, m):
    c = math.cos(i_deg * RAD)
    return k * (5.0 * c * c - 1.0) / 2.0 - m * c


def psi_dot_j2(a, e, i_deg, k, m, s, j=0):
    return lam_deg_per_day(a, e) * G(i_deg, k, m) + s * OMEGA_L_DOT - j * N_SUN


def a_resonant(e, i_deg, k, m, s, j=0):
    """Semi-major axis (km) at which Psi is stationary for given e, i, or None."""
    g = G(i_deg, k, m)
    target = j * N_SUN - s * OMEGA_L_DOT                    # lam * g = target
    if abs(g) < 1e-9:
        return None
    lam_res = target / g
    if lam_res <= 0:
        return None
    # lam = 1.5 J2 sqrt(mu) R^2 a^(-7/2) (1-e^2)^(-2) in rad/s -> deg/day
    k0 = 1.5 * J2 * math.sqrt(MU) * RE ** 2 / (1.0 - e * e) ** 2 * 86400.0 / RAD
    return (k0 / lam_res) ** (2.0 / 7.0)


def family():
    """List of (label, k, m, s, j, direct_e_pumping) members of the quadrupole slow-angle family."""
    out = []
    for m in range(-2, 3):
        for s in range(-2, 3):
            if m == 0 and s == 0:
                out.append(("2w", 2, 0, 0, 0, True))
                continue
            out.append((f"2w{m:+d}O{s:+d}L", 2, m, s, 0, True))
    for m in (1, 2):
        for s in range(-2, 3):
            if s == 0:
                out.append((f"{m}O", 0, m, 0, 0, False))
            else:
                out.append((f"{m}O{s:+d}L", 0, m, s, 0, False))
    out.append(("2w+2O-2lamS (solar U-turn)", 2, 2, 0, 2, True))
    return out


def theta_L(jd):
    """Lunar ascending-node longitude (deg, ecliptic) at JD."""
    return OMEGA_L0 + OMEGA_L_DOT * (jd - 2451545.0)
