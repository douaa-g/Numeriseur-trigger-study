"""Core models: pulse generator, leading-edge and constant-fraction timing (float),
and a bit-exact fixed-point CFD model that mirrors rtl/cfd.v."""
import numpy as np

TAU_R, TAU_D = 2.0, 20.0          # ns, illustrative detector pulse (not the real Lohengrin detector)
T_PRE, T_WIN = 40.0, 200.0        # ns of baseline before the pulse, total window
A_MIN = 0.3                       # pulse amplitudes uniform in [A_MIN, 1] x Amax

_tp = TAU_R * TAU_D / (TAU_D - TAU_R) * np.log(TAU_D / TAU_R)
_NORM = np.exp(-_tp / TAU_D) - np.exp(-_tp / TAU_R)

def pulse(t, A):
    """Bi-exponential pulse, peak = A. t, A broadcastable. Zero for t<0."""
    tt = np.clip(t, 0, None)
    return np.where(t >= 0, A * (np.exp(-tt / TAU_D) - np.exp(-tt / TAU_R)) / _NORM, 0.0)

def make_events(n, fs, snr, rng, a_min=A_MIN):
    """n pulses sampled at fs (GS/s), amplitude spread, white noise sigma = 1/snr (Amax = 1).
    Random sub-sample arrival phase. Returns x (n,L), t0 (n,) true arrival time [ns], A (n,)."""
    ts = 1.0 / fs
    L = int(T_WIN * fs)
    tk = np.arange(L) * ts
    t0 = T_PRE + rng.uniform(0, ts, n)
    A = a_min + (1 - a_min) * rng.uniform(0, 1, n)
    x = pulse(tk[None, :] - t0[:, None], A[:, None]) + rng.normal(0, 1.0 / snr, (n, L))
    return x, t0, A

def _first_true(mask):
    idx = np.argmax(mask, axis=1)
    ok = mask[np.arange(mask.shape[0]), idx]
    return idx, ok

def leading_edge(x, thr, ts):
    """Fixed-threshold timing with linear interpolation. Returns t_est [ns], valid mask."""
    idx, ok = _first_true(x >= thr)
    ok &= idx > 0
    i = np.where(ok, idx, 1)
    r = np.arange(x.shape[0])
    x0, x1 = x[r, i - 1], x[r, i]
    frac = (thr - x0) / (x1 - x0)
    return (i - 1 + frac) * ts, ok

def cfd(x, f, d, arm, ts):
    """Digital CFD: y[n] = f*x[n] - x[n-d]; zero crossing (+ to -) after arming, linear interpolation."""
    n_ev, L = x.shape
    xd = np.zeros_like(x)
    xd[:, d:] = x[:, :L - d]
    y = f * x - xd
    a, ok_a = _first_true(x >= arm)
    n = np.arange(L)[None, :]
    cond = np.zeros_like(x, dtype=bool)
    cond[:, 1:] = (y[:, :-1] > 0) & (y[:, 1:] <= 0)
    cond &= n > a[:, None]
    idx, ok = _first_true(cond)
    ok &= ok_a
    i = np.where(ok, idx, 1)
    r = np.arange(n_ev)
    y0, y1 = y[r, i - 1], y[r, i]
    frac = y0 / np.where(y1 == y0, 1, (y0 - y1))
    return (i - 1 + frac) * ts, ok

def stats(t_est, ok, t0):
    """RMS timing error (constant offset removed, i.e. calibrated) in ps, and efficiency."""
    if ok.sum() < 10:
        return np.nan, ok.mean()
    e = (t_est - t0)[ok]
    return 1000.0 * e.std(), ok.mean()

# ---------------- fixed-point model (mirrors rtl/cfd.v, cycle for cycle) -------------
def cfd_fixed_stream(x, D, K, THR, HOLD=16):
    """x: int sequence (12-bit signed samples). y = ((K*x)>>>3) - x[n-D].
    arm when x>=THR (not while locked-out); fire when armed & y_prev>0 & y<=0 (n>arm sample);
    after firing the unit stays locked out until x<THR for HOLD consecutive samples.  Returns list of (n, y_prev, y)."""
    out = []
    dl = [0] * D
    armed = False
    lock = False
    low = 0
    y_prev = 0
    for n, xv in enumerate(x):
        xd = dl[D - 1]
        y = ((K * int(xv)) >> 3) - xd
        fire = armed and (y_prev > 0) and (y <= 0)
        if fire:
            out.append((n, y_prev, y))
        # state update (same priority order as the RTL)
        if lock:
            if xv < THR:
                if low == HOLD - 1:
                    lock = False
                    low = 0
                else:
                    low += 1
            else:
                low = 0
        elif armed:
            if fire:
                armed = False
                lock = True
        else:
            if xv >= THR:
                armed = True
        y_prev = y
        dl = [int(xv)] + dl[:-1]
    return out
