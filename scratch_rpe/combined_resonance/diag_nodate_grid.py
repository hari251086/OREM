"""Model-minus-real decay error vs B multiplier for selected no-date cases (issue #56, THEORY section 16)."""
import sys

import numpy as np

sys.path.insert(0, 'E:/GitHub/OREM/scratch_rpe/combined_resonance')
import diag_nodate as dn  # noqa: E402

CASES = [(int(a), float(b)) for a, b in (arg.split(':') for arg in sys.argv[1:])] or [(27906, 0.5), (21935, 1.0)]
for n, L in CASES:
    o = dn.objects[n]
    rec = dn.load_full_record(f"{dn.BASE}/{o['tle_file']}")
    jd, a, e, inc, O, w = np.array(dn.clean_rows(rec)).T
    decay = dn.cal2jd(*o['decay_ymd'])
    p = dn.get_object_params(n, dn.cfg); B0 = p.cd * p.area_m2 / p.mass_kg * 1e-6
    st, jd0 = dn.start_state(jd, a, e, inc, O, w, decay - L * 365.25)
    print(f"{n} lead {L}: start hp={st[0]*(1-st[1])-dn.RE:.0f} km ha={st[0]*(1+st[1])-dn.RE:.0f} km, real decay {decay-jd0:.0f} d after start")
    for m in [1, 8, 64, 512, 4096]:
        dj, traj = dn.decay_with(st, jd0, B0 * m, L + 4.0)
        hpmin = traj[:, 3].min()
        t_hpmin = traj[np.argmin(traj[:, 3]), 0] - jd0
        msg = f"decay {dj-decay:+.0f} d vs real" if dj == dj else "no decay in horizon"
        print(f"   B x{m:5d} (BN {p.mass_kg/(p.cd*p.area_m2)/m:7.2f}): {msg}; model perigee min {hpmin:.0f} km at +{t_hpmin:.0f} d")
