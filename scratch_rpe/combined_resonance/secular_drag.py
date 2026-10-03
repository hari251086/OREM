"""Orbit-averaged atmospheric drag for the secular propagator (issue #56 follow-up: make the event decay itself).

Drag acceleration a_D = -(1/2) B rho v^2 v_hat  (B = Cd A / m, atmosphere at rest -- co-rotation neglected).
Instantaneous  da/dt = (2 a^2 / mu) v . a_D = -(a^2/mu) B rho v^3,   dh/dt = -(1/2) B rho v h   (h = angular momentum),
and with h^2 = mu a (1 - e^2):  de/dt = [ (1 - e^2) da/dt - 2 h dh/dt / mu ] / (2 a e).
Both are averaged over the orbit by quadrature in eccentric anomaly (dt = (1 - e cos E)/n dE), on a grid concentrated at perigee.
Density: OREM's own input/ATM.DAT (60-500 km: scale height and density; a static, low-activity atmosphere), extended above 500 km with a
linearly growing scale height.  Units: km, s, kg; B in km^2/kg, rho in kg/km^3.
"""
import math
from pathlib import Path
import numpy as np

MU = 398600.4415
RE = 6378.1363
_ATM = Path('E:/GitHub/OREM/input/ATM.DAT')


def _load_atm():
    lines = _ATM.read_text().splitlines()

    def parse(ls):
        out = []
        for l in ls:
            for k in range(0, len(l), 10):
                t = l[k:k + 10].strip()
                if t:
                    out.append(float(t))
        return out
    sch = parse(lines[0:8])[:61] + parse(lines[8:37])
    den = parse(lines[37:45])[:61] + parse(lines[45:74])
    alt = [60.0 + i for i in range(141)] + [200.0 + 2.0 * (i + 1) for i in range(150)]
    return np.array(alt), np.array(den), np.array(sch)


_ALT, _DEN, _SCH = _load_atm()
_LOGDEN = np.log(_DEN)
_H500, _RHO500, _HS500 = _ALT[-1], _DEN[-1], _SCH[-1]
_K = 0.25                                  # d(scale height)/d(altitude) above 500 km
_U = np.linspace(-1.0, 1.0, 401)
_E_GRID = math.pi * np.sinh(4.0 * _U) / math.sinh(4.0)       # eccentric anomaly, dense near perigee (E = 0)


def density_kg_m3(h_km):
    """Atmospheric density (kg/m^3) at altitude h (km), vectorised."""
    h = np.asarray(h_km, float)
    below = np.exp(np.interp(h, _ALT, _LOGDEN))
    d = np.maximum(h - _H500, 0.0)
    above = _RHO500 * (1.0 + _K * d / _HS500) ** (-1.0 / _K)
    below_range = np.where(h < _ALT[0], _DEN[0] * np.exp((_ALT[0] - h) / _SCH[0]), below)
    return np.where(h > _H500, above, below_range)


def drag_rates(a, e, B):
    """Orbit-averaged (da/dt [km/s], de/dt [1/s]) for ballistic coefficient B = Cd A / m in km^2/kg."""
    E = _E_GRID
    r = a * (1.0 - e * np.cos(E))
    h = r - RE
    rho = density_kg_m3(h) * 1.0e9                         # kg/km^3
    v2 = MU * (2.0 / r - 1.0 / a)
    v = np.sqrt(np.maximum(v2, 0.0))
    w = (1.0 - e * np.cos(E)) / (2.0 * math.pi)             # dt/P weight in dE
    dadt_i = -(a * a / MU) * B * rho * v ** 3
    hang = math.sqrt(MU * a * (1.0 - e * e))
    dhdt_i = -0.5 * B * rho * v * hang
    dadt = np.trapezoid(dadt_i * w, E)
    dhdt = np.trapezoid(dhdt_i * w, E)
    dedt = ((1.0 - e * e) * dadt - 2.0 * hang * dhdt / MU) / (2.0 * a * e)
    return float(dadt), float(dedt)
