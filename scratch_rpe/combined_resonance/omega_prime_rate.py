"""Closed form for the rate of the perturber azimuth w'_P and its U-turn (stationarity) locus.

Geometry (derived here, checked numerically below).  Satellite normal h(i, Omega), perturber normal h_P(I, Omega_P),
Delta = Omega - Omega_P, psi = mutual inclination.  The mutual node, measured from the satellite's ascending node
along the orbit, is
        w_N = atan2( sin I sin D ,  sin i cos I - cos i sin I cos D ),     w' = w - w_N,
        sin^2 psi = X^2 + Y^2  with X = sin i cos I - cos i sin I cos D,  Y = sin I sin D,
and the partials are
        dw_N/dD = sin I (sin i cos I cos D - cos i sin I) / sin^2 psi
        dw_N/di = - sin I sin D cos psi / sin^2 psi
        dw_N/dI =   sin i sin D / sin^2 psi.
With the J2 secular rates (w_dot = lam (5c^2-1)/2, Omega_dot = -lam c, i_dot = 0):
        w'_dot = lam * f(i, D, I) + g ,
        f = (5c^2-1)/2 + c * dw_N/dD ,
        g = dw_N/dD * Omega_P_dot - dw_N/dI * I_dot          (g = 0 for the Sun: fixed plane).
U-turn (stationary azimuth):   lam(a,e) = -g / f      (Moon),      f(i, D) = 0   (Sun; independent of a, e).
"""
import math
import numpy as np
from drivers import lam_deg_per_day, OMEGA_L_DOT, theta_L
from quadrupole_drive import EPS, I_MOON, _unit_normals

RAD = math.pi / 180.0


def plane_params(h):
    """Inclination I and node Omega_P (rad) of an orbit normal h in the equatorial frame."""
    I = math.acos(max(-1, min(1, h[2])))
    OmP = math.atan2(h[0], -h[1])
    return I, OmP


def moon_plane_rates(jd):
    """dI/dt and dOmega_P/dt (rad/day) of the lunar plane (finite difference of the exact geometry)."""
    d = 1.0
    h1 = _unit_normals(np.array([jd - d]))[1][0]; h2 = _unit_normals(np.array([jd + d]))[1][0]
    I1, O1 = plane_params(h1); I2, O2 = plane_params(h2)
    dO = (O2 - O1 + math.pi) % (2 * math.pi) - math.pi
    return (I2 - I1) / (2 * d), dO / (2 * d)


def partials(i, D, I):
    ci, si, cI, sI = math.cos(i), math.sin(i), math.cos(I), math.sin(I)
    X = si * cI - ci * sI * math.cos(D); Y = sI * math.sin(D)
    s2 = X * X + Y * Y
    cpsi = ci * cI + si * sI * math.cos(D)
    dND = sI * (si * cI * math.cos(D) - ci * sI) / s2
    dNi = -sI * math.sin(D) * cpsi / s2
    dNI = si * math.sin(D) / s2
    return dND, dNi, dNI, s2


def wprime_rate(a, e, i_deg, O_deg, jd, who):
    """Return (w'_dot, lam*f, g) in deg/day for who in {'sun','moon'} from J2 + perturber-plane precession."""
    i = i_deg * RAD; O = O_deg * RAD
    lam = lam_deg_per_day(a, e)
    hs, hm = _unit_normals(np.array([jd]))
    if who == 'sun':
        I, OmP = plane_params(hs[0]); Idot = 0.0; OmPdot = 0.0
    else:
        I, OmP = plane_params(hm[0]); Idot, OmPdot = moon_plane_rates(jd)
    D = O - OmP
    dND, dNi, dNI, s2 = partials(i, D, I)
    c = math.cos(i)
    f = (5 * c * c - 1) / 2 + c * dND
    g = (dND * OmPdot - dNI * Idot) / RAD
    return lam * f + g, lam * f, g, f, lam


def w_prime_deg(i_deg, O_deg, w_deg, jd, who):
    hs, hm = _unit_normals(np.array([jd]))
    I, OmP = plane_params((hs if who == 'sun' else hm)[0])
    i = i_deg * RAD; D = O_deg * RAD - OmP
    wN = math.atan2(math.sin(I) * math.sin(D), math.sin(i) * math.cos(I) - math.cos(i) * math.sin(I) * math.cos(D))
    return (w_deg - wN / RAD)


if __name__ == '__main__':
    # --- verify the closed form against the vector-geometry azimuth of quadrupole_drive.angles ---
    from quadrupole_drive import angles
    rng = np.random.default_rng(0)
    err = []
    for _ in range(200):
        i = rng.uniform(46, 90); O = rng.uniform(0, 360); w = rng.uniform(0, 360); jd = rng.uniform(2440000, 2460000)
        (wps, _, _), (wpm, _, _) = angles(np.array([jd]), np.array([26500.0]), np.array([0.7]), np.array([i]), np.array([O]), np.array([w]))
        for who, ref in (('sun', wps[0]), ('moon', wpm[0])):
            mine = w_prime_deg(i, O, w, jd, who)
            err.append(abs((mine - ref + 180) % 360 - 180))
    print(f"closed-form w' vs vector geometry: max |diff| = {max(err):.2e} deg over {len(err)} cases")
    # --- verify the rate against finite differences of w' evolved with the J2 rates + plane precession ---
    a, e, i, O, w, jd = 26560.0, 0.71, 64.5, 288.0, 276.0, 2444000.0
    lam = lam_deg_per_day(a, e); c = math.cos(i * RAD)
    dt = 0.5
    for who in ('sun', 'moon'):
        w1 = w_prime_deg(i, O, w, jd, who)
        w2 = w_prime_deg(i, O + (-lam * c) * dt, w + (lam * (5 * c * c - 1) / 2) * dt, jd + dt, who)
        w0 = w_prime_deg(i, O - (-lam * c) * dt, w - (lam * (5 * c * c - 1) / 2) * dt, jd - dt, who)
        fd = ((w2 - w0 + 180) % 360 - 180) / (2 * dt)
        cf = wprime_rate(a, e, i, O, jd, who)[0]
        print(f"{who:>5}: finite-difference w'_dot = {fd:+.5f} deg/day   closed form = {cf:+.5f} deg/day")


def wprime_rate_full(a, e, i_deg, O_deg, w_deg, jd, who='moon', use=(True, True, True)):
    """Complete (J2 + Sun + Moon quadrupole) rate of the azimuth w'_P in deg/day.
        w'_dot = w_dot - dw_N/dD (Omega_dot - Omega_P_dot) - dw_N/di i_dot - dw_N/dI I_dot,
    with w_dot, Omega_dot, i_dot from the full secular Lagrange equations (secular_propagator.rates)."""
    from secular_propagator import rates
    x = np.array([e, i_deg * RAD, O_deg * RAD, w_deg * RAD])
    de, di, dO, dw = rates(x, jd, a, use) * 86400.0            # rad/day
    hs, hm = _unit_normals(np.array([jd]))
    if who == 'sun':
        I, OmP = plane_params(hs[0]); Idot = 0.0; OmPdot = 0.0
    else:
        I, OmP = plane_params(hm[0]); Idot, OmPdot = moon_plane_rates(jd)
    D = O_deg * RAD - OmP
    dND, dNi, dNI, s2 = partials(x[1], D, I)
    return (dw - dND * (dO - OmPdot) - dNi * di - dNI * Idot) / RAD
