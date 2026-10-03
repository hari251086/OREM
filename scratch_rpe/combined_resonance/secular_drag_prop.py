"""Secular propagator WITH drag: state (a, e, i, Omega, omega); J2 + Sun + Moon quadrupole (secular_propagator.rates) plus orbit-averaged drag
(secular_drag.drag_rates).  Adaptive step by perigee altitude; stops at perigee altitude < 80 km (OREM's re-entry threshold) = decay."""
import math
import numpy as np
from secular_propagator import rates, RE
from secular_drag import drag_rates, density_kg_m3

HP_DECAY = 80.0          # perigee altitude (km): OREM's re-entry threshold
HA_DECAY = 150.0         # apogee altitude (km): an orbit lying wholly below this circularizes inside dense atmosphere and decays within days


def _step_days(hp):
    if hp > 800.0: return 5.0
    if hp > 400.0: return 2.0
    if hp > 250.0: return 1.0
    if hp > 160.0: return 0.25
    return 0.05


def _f(x, jd, B):
    a, e, i, O, w = x
    if a <= 0.0 or e <= 0.0 or e >= 1.0 or a * (1.0 - e) - RE < 0.0:      # invalid trial state (only after an overshoot)
        return np.zeros(5)
    r4 = rates(np.array([e, i, O, w]), jd, a, (True, True, True))                    # per second (de, di, dO, dw)
    hp = a * (1.0 - e) - RE
    da = de = 0.0
    if hp < 1500.0:
        da, de = drag_rates(a, e, B)
    return np.array([da, r4[0] + de, r4[1], r4[2], r4[3]])


def propagate_drag(a0, e0, i0_deg, O0_deg, w0_deg, jd0, B, years=25.0):
    """Returns (decay_jd or nan, trajectory array of rows jd, a, e, hp)."""
    x = np.array([a0, e0, math.radians(i0_deg), math.radians(O0_deg), math.radians(w0_deg)])
    jd = jd0
    out = []
    t_end = jd0 + years * 365.25
    while jd < t_end:
        a, e = x[0], x[1]
        hp = a * (1.0 - e) - RE
        ha = a * (1.0 + e) - RE
        out.append((jd, a, e, hp))
        if hp < HP_DECAY or ha < HA_DECAY or e <= 0.0:          # decay (apogee collapse ends with e -> 0 at low altitude)
            return jd, np.array(out)
        if e >= 0.995:
            return float('nan'), np.array(out)
        dt_d = _step_days(hp)
        k1 = _f(x, jd, B)
        # limit the step so the perigee altitude changes by at most 10 km (drag rates explode in the last hours)
        dhp_dt = abs(k1[0] * (1.0 - e) - a * k1[1]) * 86400.0                   # km/day
        if dhp_dt > 0.0:
            dt_d = max(1e-4, min(dt_d, 10.0 / dhp_dt))
        dt = dt_d * 86400.0
        k2 = _f(x + 0.5 * dt * k1, jd + 0.5 * dt_d, B)
        k3 = _f(x + 0.5 * dt * k2, jd + 0.5 * dt_d, B)
        k4 = _f(x + dt * k3, jd + dt_d, B)
        x = x + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        jd += dt_d
    return float('nan'), np.array(out)
