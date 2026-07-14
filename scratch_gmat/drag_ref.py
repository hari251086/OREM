"""Independent drag reference: two-body + exponential-atmosphere drag, RK4.

Matches propagate_ks's own drag assumptions exactly:
  rho(r) = rho_p * exp(-(r - rp0)/H), H = const = table scale height at perigee
  a_drag = -rho * V * vec_v / (2 BN)   (with atmosphere rotation OFF, F factor off)
Orbit: N13 test case (SMA=9636.790207, ECC=0.321216, BN=80), 7 days.
Compare 7-day apogee-radius drop against propagate_ks (-16.24 km w/ new table)
and the GMAT JR reference (-62.26 km).
"""
import math

mu = 398600.4415  # km^3/s^2 (value irrelevant at this precision)
a0, e0 = 9636.790207, 0.321216
BN = 80.0          # kg/m^2
RE = 6378.1363

# new-table values at the perigee altitude
rp0 = a0 * (1 - e0)          # km
hp0 = rp0 - RE               # ~159.6 km
# from the generated table: interpolate DEN/SCH at hp0
lines = open(r"C:\Users\hari2\OneDrive\Documents\GitHub\OREM\input\ATM.DAT").read().splitlines()
def fields(ls, n, w=10):
    vals = []
    for li, l in enumerate(ls):
        for i in range(0, len(l), w):
            s = l[i:i+w].strip()
            if s: vals.append(float(s))
        if len(vals) >= n: return vals[:n], ls[li+1:]
rest = lines
sch1, rest = fields(rest, 61); sch2, rest = fields(rest, 230)
den1, rest = fields(rest, 61); den2, rest = fields(rest, 230)
den = den1 + den2; sch = sch1 + sch2
alt = [60.0 + i for i in range(141)] + [202.0 + 2*i for i in range(150)]
def interp(h, ys):
    for i in range(len(alt)-1):
        if alt[i] <= h <= alt[i+1]:
            f = (h - alt[i])/(alt[i+1] - alt[i])
            return ys[i] + f*(ys[i+1] - ys[i])
rho_p = interp(hp0, den)     # kg/m^3
H     = interp(hp0, sch)     # km
print(f"hp0={hp0:.2f} km  rho_p={rho_p:.3e} kg/m3  H={H:.2f} km")

def acc(state):
    x, y, vx, vy = state
    r = math.hypot(x, y)
    ax, ay = -mu*x/r**3, -mu*y/r**3
    rho = rho_p * math.exp(-(r - rp0)/H)      # kg/m^3
    v = math.hypot(vx, vy)                    # km/s
    # a = rho * V^2/(2BN): SI -> km/s^2 gives factor 1e3
    f = -rho * v * 1e3 / (2*BN)               # 1/s * (km/s)/... => km/s^2 per (km/s)
    return ax + f*vx, ay + f*vy

def rk4(state, dt):
    def deriv(s):
        ax, ay = acc(s)
        return (s[2], s[3], ax, ay)
    k1 = deriv(state)
    k2 = deriv([state[i] + 0.5*dt*k1[i] for i in range(4)])
    k3 = deriv([state[i] + 0.5*dt*k2[i] for i in range(4)])
    k4 = deriv([state[i] + dt*k3[i] for i in range(4)])
    return [state[i] + dt/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(4)]

# start at perigee, planar orbit
rp = a0*(1-e0); vp = math.sqrt(mu*(2/rp - 1/a0))
s = [rp, 0.0, 0.0, vp]
T = 2*math.pi*math.sqrt(a0**3/mu)
dt = 5.0  # s
t, tend = 0.0, 7*86400.0
def elems(s):
    x, y, vx, vy = s
    r = math.hypot(x, y); v2 = vx*vx + vy*vy
    a = 1/(2/r - v2/mu)
    hz = x*vy - y*vx
    e = math.sqrt(max(0.0, 1 - hz*hz/(mu*a)))
    return a, e, a*(1+e)
a_, e_, ra0 = elems(s)
print(f"initial apogee radius = {ra0:.3f} km")
while t < tend:
    s = rk4(s, dt); t += dt
a_, e_, ra7 = elems(s)
print(f"final   apogee radius = {ra7:.3f} km")
print(f"7-day apogee drop = {ra7 - ra0:+.2f} km")
print(f"propagate_ks (new table): -16.24 km   GMAT JR ref: -62.26 km")
