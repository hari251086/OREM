"""Sanity checks for secular_drag (run: python test_secular_drag.py)."""
import math, sys
sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
from secular_drag import drag_rates, density_kg_m3, MU, RE

# 1) density table anchors (ATM.DAT) and monotone decrease
assert abs(density_kg_m3(300.0) - 4.99e-12) / 4.99e-12 < 1e-6
assert abs(density_kg_m3(150.0) - 1.53e-9) / 1.53e-9 < 1e-6
hs = [100, 200, 300, 450, 500, 600, 800, 1200]
d = [float(density_kg_m3(h)) for h in hs]
assert all(d[i] > d[i + 1] for i in range(len(d) - 1)), d

# 2) circular-orbit limit: da/dt = -B rho sqrt(mu a)
a = RE + 300.0; B = 0.002 * 1e-3 * 2.2 / 1.0   # arbitrary B in km^2/kg
da, de = drag_rates(a, 1e-3, B)
exact = -B * float(density_kg_m3(300.0)) * 1e9 * math.sqrt(MU * a)
assert abs(da - exact) / abs(exact) < 0.02, (da, exact)

# 3) eccentric orbit: drag lowers a and e; effect is concentrated at perigee (higher perigee -> weaker)
da1, de1 = drag_rates(26550.0, 0.74, 1.2e-3 * 1e-3)    # perigee ~ 300 km
da2, de2 = drag_rates(26550.0 , 0.70, 1.2e-3 * 1e-3)   # perigee ~ 1,200 km? (higher)
assert da1 < 0 and de1 < 0
assert abs(da2) < abs(da1)
print('secular_drag sanity checks passed')
