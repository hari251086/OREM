"""Doubly-averaged quadrupole lunisolar drive on the perigee radius, evaluated from TLE mean elements.

For a circular perturber P (Sun, Moon) with orbit normal h_P, the doubly-averaged quadrupole disturbing
function is  R = C [ (2+3e^2)(3cos^2 psi - 1) + 15 e^2 sin^2 psi cos 2w' ],  C = mu_P a^2 / (16 a_P^3),
with psi the mutual inclination and w' the argument of perigee measured from the mutual node (Lidov 1962,
Kozai 1962; same model as Celletti et al. 2016 / Alessi et al. 2021).  Lagrange's equation gives
        de/dt = (15/8) (mu_P / (a_P^3 n)) e sqrt(1-e^2) sin^2 psi sin 2w'
and, since a is constant in the averaged motion, dq/dt = -a de/dt for the perigee radius q = a(1-e).
Nothing here assumes which angle combination is resonant: summing the exact drive and averaging it over a
window leaves only the slowly varying (resonant) part -- the generalization of the solar-azimuth U-turn.
"""
import math
import numpy as np

MU = 398600.4415
MU_SUN = 1.32712440018e11
MU_MOON = 4902.800066
A_SUN = 149597870.7
A_MOON = 384400.0
EPS = math.radians(23.4392911)      # obliquity
I_MOON = math.radians(5.145)        # lunar orbit inclination to the ecliptic
OMEGA_L_DOT = -1934.136261 / 36525.0
OMEGA_L0 = 125.04452


def _unit_normals(jd):
    """Orbit-normal unit vectors (equatorial frame) of the Sun and Moon at times jd (array)."""
    jd = np.asarray(jd, float)
    n = len(jd)
    # Sun: ecliptic normal (0,0,1) -> equatorial
    hs = np.tile(np.array([0.0, -math.sin(EPS), math.cos(EPS)]), (n, 1))
    # Moon: normal in ecliptic frame, node Omega_L(t), then rotate ecliptic -> equatorial about x by EPS
    OL = np.radians(OMEGA_L0 + OMEGA_L_DOT * (jd - 2451545.0))
    hx = math.sin(I_MOON) * np.sin(OL)
    hy = -math.sin(I_MOON) * np.cos(OL)
    hz = np.full(n, math.cos(I_MOON))
    hm = np.stack([hx, hy * math.cos(EPS) - hz * math.sin(EPS), hy * math.sin(EPS) + hz * math.cos(EPS)], axis=1)
    return hs, hm


def drive_dqdt(jd, a, e, i_deg, raan_deg, argp_deg):
    """Return dq/dt in km/day for (Sun, Moon) as two arrays, from mean elements at each epoch."""
    jd = np.asarray(jd, float)
    a = np.asarray(a, float); e = np.asarray(e, float)
    i = np.radians(i_deg); O = np.radians(raan_deg); w = np.radians(argp_deg)
    h = np.stack([np.sin(i) * np.sin(O), -np.sin(i) * np.cos(O), np.cos(i)], axis=1)
    ev = np.stack([np.cos(O) * np.cos(w) - np.sin(O) * np.sin(w) * np.cos(i),
                   np.sin(O) * np.cos(w) + np.cos(O) * np.sin(w) * np.cos(i),
                   np.sin(w) * np.sin(i)], axis=1)
    n = np.sqrt(MU / a ** 3)                                   # rad/s
    out = []
    for hp, mu_p, a_p in zip(_unit_normals(jd), (MU_SUN, MU_MOON), (A_SUN, A_MOON)):
        cospsi = np.sum(h * hp, axis=1)
        sin2psi = 1.0 - cospsi ** 2
        node = np.cross(hp, h)
        nn = np.linalg.norm(node, axis=1)
        nn = np.where(nn < 1e-12, np.nan, nn)
        node = node / nn[:, None]
        cw = np.sum(node * ev, axis=1)
        sw = np.sum(np.cross(node, ev) * h, axis=1)
        sin2w = 2.0 * sw * cw
        dedt = 1.875 * (mu_p / a_p ** 3) / n * e * np.sqrt(1 - e ** 2) * sin2psi * sin2w   # 1/s
        out.append(-a * dedt * 86400.0)                          # km/day
    return out[0], out[1]


def angles(jd, a, e, i_deg, raan_deg, argp_deg):
    """Per-perturber azimuth omega'_P (perigee measured from the mutual node, deg, wrapped), sin^2(psi_P),
    and the dq/dt amplitude A_P (km/day) so that dq/dt = -A_P-signed * sin(2 omega'_P):
        dq/dt = -a * (15/8)(mu_P/(a_P^3 n)) e sqrt(1-e^2) sin^2(psi) sin(2 omega')."""
    jd = np.asarray(jd, float); a = np.asarray(a, float); e = np.asarray(e, float)
    i = np.radians(i_deg); O = np.radians(raan_deg); w = np.radians(argp_deg)
    h = np.stack([np.sin(i) * np.sin(O), -np.sin(i) * np.cos(O), np.cos(i)], axis=1)
    ev = np.stack([np.cos(O) * np.cos(w) - np.sin(O) * np.sin(w) * np.cos(i),
                   np.sin(O) * np.cos(w) + np.cos(O) * np.sin(w) * np.cos(i),
                   np.sin(w) * np.sin(i)], axis=1)
    n = np.sqrt(MU / a ** 3)
    res = []
    for hp, mu_p, a_p in zip(_unit_normals(jd), (MU_SUN, MU_MOON), (A_SUN, A_MOON)):
        cospsi = np.sum(h * hp, axis=1)
        sin2psi = 1.0 - cospsi ** 2
        node = np.cross(hp, h)
        nn = np.linalg.norm(node, axis=1)
        nn = np.where(nn < 1e-12, np.nan, nn)
        node = node / nn[:, None]
        cw = np.sum(node * ev, axis=1)
        sw = np.sum(np.cross(node, ev) * h, axis=1)
        wp = np.degrees(np.arctan2(sw, cw))
        amp = a * 1.875 * (mu_p / a_p ** 3) / n * e * np.sqrt(1 - e ** 2) * sin2psi * 86400.0
        res.append((wp, sin2psi, amp))
    return res
