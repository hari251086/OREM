"""Orbit-averaged atmospheric drag for the secular propagator (issue #56 follow-up: make the event decay itself).

Drag acceleration a_D = -(1/2) B rho v^2 v_hat  (B = Cd A / m, atmosphere at rest -- co-rotation neglected).
Instantaneous  da/dt = (2 a^2 / mu) v . a_D = -(a^2/mu) B rho v^3,   dh/dt = -(1/2) B rho v h   (h = angular momentum),
and with h^2 = mu a (1 - e^2):  de/dt = [ (1 - e^2) da/dt - 2 h dh/dt / mu ] / (2 a e).
Both are averaged over the orbit by quadrature in eccentric anomaly (dt = (1 - e cos E)/n dE), on a grid concentrated at perigee.
Density: OREM's own input/ATM.DAT (60-500 km: scale height and density; a static, low-activity atmosphere), extended above 500 km with a
linearly growing scale height.  Units: km, s, kg; B in km^2/kg, rho in kg/km^3.
ATM.DAT matches the J71 table input/ATM2D.DAT at an exospheric temperature of ~620 K at 300-500 km, i.e. deep solar
minimum, night side.  set_atmosphere('hp') switches to Harris-Priester for mean solar activity instead (Montenbruck &
Gill, Satellite Orbits, Table 3.8, after Long et al. 1989), orbit-averaged over the bulge as (rho_min + rho_max) / 2
(the mean of cos^2(psi/2)); ln(rho) interpolated linearly in altitude between the tabulated heights (THEORY section 22).
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


# Harris-Priester, mean solar activity: height (km), rho_min, rho_max (g/km^3 = 1e-12 kg/m^3)
_HP = np.array([
    (100, 497400.0, 497400.0), (120, 24900.0, 24900.0), (130, 8377.0, 8710.0), (140, 3899.0, 4059.0),
    (150, 2122.0, 2215.0), (160, 1263.0, 1344.0), (170, 800.8, 875.8), (180, 528.3, 601.0), (190, 361.7, 429.7),
    (200, 255.7, 316.2), (210, 183.9, 239.6), (220, 134.1, 185.3), (230, 99.49, 145.5), (240, 74.88, 115.7),
    (250, 57.09, 93.08), (260, 44.03, 75.55), (270, 34.30, 61.82), (280, 26.97, 50.95), (290, 21.39, 42.26),
    (300, 17.08, 35.26), (320, 10.99, 25.11), (340, 7.214, 18.19), (360, 4.824, 13.37), (380, 3.274, 9.955),
    (400, 2.249, 7.492), (420, 1.558, 5.684), (440, 1.091, 4.355), (460, 0.7701, 3.362), (480, 0.5474, 2.612),
    (500, 0.3916, 2.042), (520, 0.2819, 1.605), (540, 0.2042, 1.267), (560, 0.1488, 1.005), (580, 0.1092, 0.7997),
    (600, 0.08070, 0.6390), (620, 0.06012, 0.5123), (640, 0.04519, 0.4121), (660, 0.03430, 0.3325),
    (680, 0.02632, 0.2691), (700, 0.02043, 0.2185), (720, 0.01607, 0.1779), (740, 0.01281, 0.1452),
    (760, 0.01036, 0.1190), (780, 0.008496, 0.09776), (800, 0.007069, 0.08059), (840, 0.004680, 0.05741),
    (880, 0.003200, 0.04210), (920, 0.002210, 0.03130), (960, 0.001560, 0.02360), (1000, 0.001150, 0.01810)])
_HP_H, _HP_LNMIN, _HP_LNMAX = _HP[:, 0], np.log(_HP[:, 1] * 1e-12), np.log(_HP[:, 2] * 1e-12)
_ATMOSPHERE = 'atm_dat'


def set_atmosphere(name):
    """'atm_dat' (default: OREM's static input/ATM.DAT) or 'hp' (Harris-Priester, mean solar activity)."""
    global _ATMOSPHERE
    if name not in ('atm_dat', 'hp'):
        raise ValueError(name)
    _ATMOSPHERE = name


def _hp_density(h):
    def ln_rho(lnr):
        lo = (lnr[1] - lnr[0]) / (_HP_H[1] - _HP_H[0])          # end slopes extend the table both ways
        hi = (lnr[-1] - lnr[-2]) / (_HP_H[-1] - _HP_H[-2])
        return np.where(h < _HP_H[0], lnr[0] + lo * (h - _HP_H[0]),
                        np.where(h > _HP_H[-1], lnr[-1] + hi * (h - _HP_H[-1]), np.interp(h, _HP_H, lnr)))
    return 0.5 * (np.exp(ln_rho(_HP_LNMIN)) + np.exp(ln_rho(_HP_LNMAX)))


def density_kg_m3(h_km):
    """Atmospheric density (kg/m^3) at altitude h (km), vectorised."""
    h = np.asarray(h_km, float)
    if _ATMOSPHERE == 'hp':
        return _hp_density(h)
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
