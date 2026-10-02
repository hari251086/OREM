"""Doubly-averaged secular propagator: J2 + Sun + Moon quadrupole, a = const (no drag).

Lagrange's planetary equations (R positive convention) with
    R = R_J2 + sum_P C_P [ (2+3e^2)(3cos^2 psi_P - 1) + 15 e^2 sin^2 psi_P cos 2w'_P ],
    R_J2 = mu J2 Re^2 / (4 a^3 (1-e^2)^1.5) * (3cos^2 i - 1),   C_P = mu_P a^2 / (16 a_P^3),
partial derivatives by central differences.  The Moon's orbit normal precesses with the 18.6-yr nodal regression.
No free parameters; state = (e, i, Omega, omega) from a TLE mean-element set.
"""
import math
import numpy as np
from quadrupole_drive import (MU, MU_SUN, MU_MOON, A_SUN, A_MOON, EPS, I_MOON, OMEGA_L_DOT, OMEGA_L0)

J2 = 1.08263e-3
RE = 6378.1363
_hs = np.array([0.0, -math.sin(EPS), math.cos(EPS)])


def _h_moon(jd):
    OL = math.radians(OMEGA_L0 + OMEGA_L_DOT * (jd - 2451545.0))
    hx = math.sin(I_MOON) * math.sin(OL); hy = -math.sin(I_MOON) * math.cos(OL); hz = math.cos(I_MOON)
    return np.array([hx, hy * math.cos(EPS) - hz * math.sin(EPS), hy * math.sin(EPS) + hz * math.cos(EPS)])


def R_fun(e, i, O, w, jd, a, use=(True, True, True)):
    r = 0.0
    if use[0]:
        r += MU * J2 * RE ** 2 / (4 * a ** 3 * (1 - e * e) ** 1.5) * (3 * math.cos(i) ** 2 - 1)
    h = np.array([math.sin(i) * math.sin(O), -math.sin(i) * math.cos(O), math.cos(i)])
    ev = np.array([math.cos(O) * math.cos(w) - math.sin(O) * math.sin(w) * math.cos(i),
                   math.sin(O) * math.cos(w) + math.cos(O) * math.sin(w) * math.cos(i),
                   math.sin(w) * math.sin(i)])
    for flag, hp, mu_p, a_p in ((use[1], _hs, MU_SUN, A_SUN), (use[2], _h_moon(jd), MU_MOON, A_MOON)):
        if not flag:
            continue
        cp = float(np.dot(h, hp)); s2 = 1 - cp * cp
        node = np.cross(hp, h); nn = float(np.linalg.norm(node))
        if nn < 1e-12:
            continue
        node /= nn
        cw = float(np.dot(node, ev)); sw = float(np.dot(np.cross(node, ev), h))
        c2w = cw * cw - sw * sw
        C = mu_p * a * a / (16 * a_p ** 3)
        r += C * ((2 + 3 * e * e) * (3 * cp * cp - 1) + 15 * e * e * s2 * c2w)
    return r


def rates(x, jd, a, use):
    e, i, O, w = x
    n = math.sqrt(MU / a ** 3)
    d = [1e-6, 1e-6, 1e-6, 1e-6]
    def dR(k):
        xp = list(x); xm = list(x); xp[k] += d[k]; xm[k] -= d[k]
        return (R_fun(*xp, jd, a, use) - R_fun(*xm, jd, a, use)) / (2 * d[k])
    Re_, Ri, RO, Rw = dR(0), dR(1), dR(2), dR(3)
    s = math.sqrt(1 - e * e); si, ci = math.sin(i), math.cos(i)
    de = -s / (n * a * a * e) * Rw
    di = (ci * Rw - RO) / (n * a * a * s * si)
    dO = Ri / (n * a * a * s * si)
    dw = s / (n * a * a * e) * Re_ - ci / (n * a * a * s * si) * Ri
    return np.array([de, di, dO, dw])          # per second


def propagate(e0, i0_deg, O0_deg, w0_deg, a, jd0, years, dt_days=10.0, use=(True, True, True), hp_stop=None):
    x = np.array([e0, math.radians(i0_deg), math.radians(O0_deg), math.radians(w0_deg)])
    n_steps = int(years * 365.25 / dt_days)
    dt = dt_days * 86400.0
    out = []
    jd = jd0
    for k in range(n_steps + 1):
        out.append((jd, x[0], x[1], x[2], x[3]))
        if hp_stop is not None and a * (1 - x[0]) - RE < hp_stop:
            break
        k1 = rates(x, jd, a, use)
        k2 = rates(x + 0.5 * dt * k1, jd + 0.5 * dt_days, a, use)
        k3 = rates(x + 0.5 * dt * k2, jd + 0.5 * dt_days, a, use)
        k4 = rates(x + dt * k3, jd + dt_days, a, use)
        x = x + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        jd += dt_days
    arr = np.array(out)
    return arr   # columns: jd, e, i(rad), Omega(rad), omega(rad)
