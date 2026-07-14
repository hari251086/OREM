"""Independent drag reference: two-body + exponential-atmosphere drag, RK4.

Matches propagate_ks's own drag assumptions exactly:
  rho(r) = rho_p * exp(-(r - rp0)/H), H = const = table scale height at perigee
  a_drag = -F * rho * V * vec_v / (2 BN)
with, optionally matched to propagate_ks's conventions:
  - F = (1 - rp*wE*cos(i)/Vp)^2  atmosphere co-rotation factor
  - rho_p evaluated at the OBLATE perigee altitude hp = rp - RE*(1 - eps*sin^2(i)*sin^2(aop))
Orbit: N13/N14 cross-validation case (SMA=9636.790207, ECC=0.321216, i=19.1713,
AOP=273.5646), duration = 35 orbital periods (matching test_gmat's nrev_test=35 --
NOTE: that is ~3.8 days, NOT the 7 days GMAT's reference numbers span; the
historical N13/N14 "2x discrepancy" was exactly this duration mismatch, see #25).

History: the first version of this script (uncorrected duration, no F, spherical
perigee) gave -37.02 km over 7 days and seeded issue #25's "2x drag deficit"
hypothesis. With duration/F/oblateness matched, propagate_ks agrees with this
reference to ~1% (BN=80: -16.05 vs -16.24; BN=160: -8.06 vs -8.12), exonerating
the drag implementation at the revolution level.
"""
import math

mu = 398600.4415  # km^3/s^2
a0, e0 = 9636.790207, 0.321216
inc = math.radians(19.1713)
aop = math.radians(273.564599)
RE = 6378.1363
wE = 7.2921150e-5  # rad/s
eps_f = 3.35281066e-3
NREV = 35          # match test_gmat nrev_test

# --- table lookup at the OBLATE perigee altitude (as propagate_ks does)
rp0 = a0 * (1 - e0)
RREQ = RE * (1.0 - eps_f * math.sin(inc)**2 * math.sin(aop)**2)
hp_obl = rp0 - RREQ
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
rho_p = interp(hp_obl, den)   # kg/m^3
H     = interp(hp_obl, sch)   # km

# --- co-rotation factor F at perigee
vp = math.sqrt(mu*(2/rp0 - 1/a0))
F = (1.0 - rp0*wE*math.cos(inc)/vp)**2
print(f"hp(oblate)={hp_obl:.2f} km  rho_p={rho_p:.3e} kg/m3  H={H:.2f} km  F={F:.4f}")

def run(BN):
    def acc(state):
        x, y, vx, vy = state
        r = math.hypot(x, y)
        ax, ay = -mu*x/r**3, -mu*y/r**3
        rho = rho_p * math.exp(-(r - rp0)/H)
        v = math.hypot(vx, vy)
        f = -F * rho * v * 1e3 / (2*BN)
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
    def elems(s):
        x, y, vx, vy = s
        r = math.hypot(x, y); v2 = vx*vx + vy*vy
        a = 1/(2/r - v2/mu)
        hz = x*vy - y*vx
        e = math.sqrt(max(0.0, 1 - hz*hz/(mu*a)))
        return a*(1+e)
    s = [rp0, 0.0, 0.0, vp]
    T = 2*math.pi*math.sqrt(a0**3/mu)
    ra0 = elems(s)
    t, tend, dt = 0.0, NREV*T, 5.0
    while t < tend:
        s = rk4(s, dt); t += dt
    ra1 = elems(s)
    return ra0, ra1

for BN in (80.0, 160.0):
    ra0, ra1 = run(BN)
    print(f"BN={BN:6.1f}: apogee {ra0:.3f} -> {ra1:.3f} km, "
          f"drop over {NREV} revs = {ra1-ra0:+.2f} km")
