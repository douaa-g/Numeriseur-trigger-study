"""Shared models: the four-channel CFD (bit-exact fixed point) and the cycle model of rtl/trigger4.v.
Used by system.py (study) and viewer_qt.py (GUI)."""
from cfd_core import cfd_fixed_stream

TS = 4.0                                 # ns per sample (250 MS/s)
AMAX = 1500                              # LSB, full-scale amplitude used by the study
D, K, THR, HOLD = 3, 5, 90, 16           # CFD: f = 5/8, delay 3 samples (12 ns); arming at 6 sigma (SNR 100)

def cfd_hits(xi):
    """Per channel: list of (n, y_prev, y) from the bit-exact fixed-point CFD model (same as rtl/cfd.v)."""
    return [cfd_fixed_stream(xi[c], D, K, THR, HOLD) for c in range(len(xi))]

def t_interp(h):
    """Sub-sample hit time in ns from (n, y[n-1], y[n]) by linear interpolation (done in software, as in the RTL split)."""
    n, yp, yc = h
    return ((n - 1) + yp / (yp - yc)) * TS

def trigger_model(hits, tt, win, n_edges):
    """Cycle model of rtl/trigger4.v.  hits: per-channel lists from cfd_hits.  Edge t consumes the hit mask
    that was visible after edge t-1.  Returns trigger indices (visible after edge t) and scalers [A, B, C, D, trig]."""
    H = {}
    for c in range(4):
        for n, _, _ in hits[c]:
            H[n] = H.get(n, 0) | (1 << c)
    win = max(win, 1)
    opn, gc, sticky = False, 0, 0
    cnt, ctrig, trig_idx = [0] * 4, 0, []
    for t in range(n_edges):
        hit = H.get(t - 1, 0)
        acc = sticky | hit
        if not opn:
            if hit:
                opn, sticky, gc = True, hit, win
        else:
            sticky = acc
            if gc == 1:
                opn, sticky = False, 0
                if (tt >> acc) & 1:
                    trig_idx.append(t); ctrig += 1
            else:
                gc -= 1
        for c in range(4):
            cnt[c] += (hit >> c) & 1
    return trig_idx, cnt + [ctrig]
