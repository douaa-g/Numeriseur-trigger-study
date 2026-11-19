"""250 MS/s front-end study: what timing resolution survives when a pulse is sampled every 4 ns?

Pulses (same illustrative bi-exponential as cfd_core) are generated on a 16 GS/s grid, passed through a
2nd-order Butterworth low-pass (the board's anti-alias / shaping stage, cut-off fc), sampled at 250 MS/s
with a random sub-sample arrival phase, and ADC noise (white, sigma = 1/SNR, peak = 1) is added.
Timing = digital CFD, zero crossing found by linear or cubic (4-point Lagrange) interpolation.

Writes figures/fig_frontend.png and results_frontend250.json.  Illustrative model, not the LPSC detectors/board."""
import json, itertools, numpy as np
from scipy.signal import butter, sosfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter
from cfd_core import pulse, T_PRE, T_WIN, A_MIN

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 130})
OUT = "../figures/"
REF_1GSPS_PS = 46.7        # github.com/douaa-g/CFD-timing-study, results.json (held-out events)
FS = 0.25                  # GS/s  -> 4 ns per sample
TS = 1 / FS
FS_H = 16.0                # GS/s  fine grid for the analogue part
DEC = int(FS_H / FS)       # 64 fine samples per ADC sample
SNR = 100
N = 3000
FR = [k / 8 for k in range(1, 8)]
D_SAMPLES = [1, 2, 3, 4, 5]
ARM = max(6 / SNR, 0.05)

def frontend_scale(fc):
    """Peak of a unit pulse after the filter, so that Amax = 1 at the ADC input for every fc."""
    if fc is None:
        return None, 1.0
    sos = butter(2, fc * 1e-3 / (FS_H / 2), output="sos")
    tk = np.arange(int(T_WIN * FS_H)) / FS_H
    return sos, 1.0 / sosfilt(sos, pulse(tk - T_PRE, 1.0)).max()

def make_events_250(n, fc, rng, snr=SNR):
    sos, scale = frontend_scale(fc)
    Lh = int(T_WIN * FS_H)
    tk = np.arange(Lh) / FS_H
    t0 = T_PRE + rng.uniform(0, TS, n)
    A = A_MIN + (1 - A_MIN) * rng.uniform(0, 1, n)
    xh = pulse(tk[None, :] - t0[:, None], A[:, None])
    if sos is not None:
        xh = sosfilt(sos, xh, axis=1) * scale
    x = xh[:, ::DEC] + rng.normal(0, 1.0 / snr, (n, Lh // DEC))
    return x, t0, A

def _first(mask):
    idx = np.argmax(mask, axis=1)
    return idx, mask[np.arange(mask.shape[0]), idx]

def cfd_crossing(x, f, d, arm=ARM):
    """CFD signal y = f*x[n] - x[n-d]; first +->- crossing after the arming sample. Returns y, idx, ok."""
    n_ev, L = x.shape
    xd = np.zeros_like(x); xd[:, d:] = x[:, :L - d]
    y = f * x - xd
    a, ok_a = _first(x >= arm)
    cond = np.zeros_like(x, dtype=bool)
    cond[:, 1:] = (y[:, :-1] > 0) & (y[:, 1:] <= 0)
    cond &= np.arange(L)[None, :] > a[:, None]
    idx, ok = _first(cond)
    return y, np.where(ok & ok_a, idx, 2), ok & ok_a

def interp_linear(y, idx):
    r = np.arange(y.shape[0]); y0, y1 = y[r, idx - 1], y[r, idx]
    return (idx - 1 + y0 / np.where(y1 == y0, 1, y0 - y1)) * TS

def interp_cubic(y, idx):
    """Root of the 4-point Lagrange cubic through y[i-2..i+1], searched in [i-1, i] by bisection."""
    r = np.arange(y.shape[0]); L = y.shape[1]
    ym, y0, y1, y2 = y[r, idx - 2], y[r, idx - 1], y[r, idx], y[r, np.minimum(idx + 1, L - 1)]
    def p(u):
        return (ym * (-u * (u - 1) * (u - 2) / 6) + y0 * ((u + 1) * (u - 1) * (u - 2) / 2)
                + y1 * (-(u + 1) * u * (u - 2) / 2) + y2 * ((u + 1) * u * (u - 1) / 6))
    lo, hi = np.zeros(len(idx)), np.ones(len(idx))
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        pos = p(mid) > 0
        lo, hi = np.where(pos, mid, lo), np.where(pos, hi, mid)
    return (idx - 1 + 0.5 * (lo + hi)) * TS

def timing(x, f, d, interp):
    y, idx, ok = cfd_crossing(x, f, d)
    t = interp_cubic(y, idx) if interp == "cubic" else interp_linear(y, idx)
    return t, ok

def sigma_ps(t, ok, t0):
    return 1000 * (t - t0)[ok].std() if ok.sum() > 50 else np.nan

def best_cfg(x, t0, interp):
    res = {}
    for f, d in itertools.product(FR, D_SAMPLES):
        t, ok = timing(x, f, d, interp)
        if ok.mean() >= 0.99:
            res[(f, d)] = sigma_ps(t, ok, t0)
    k = min(res, key=res.get)
    return k, res[k]

def main():
    fcs = [None, 150, 100, 60, 40, 25]
    out = {"fs_MSps": 250, "snr": SNR, "events_per_set": N, "interp": {}, "fc_MHz": [fc if fc else "none" for fc in fcs]}
    for interp in ("linear", "cubic"):
        sig, cfg = [], []
        for j, fc in enumerate(fcs):
            xc, t0c, _ = make_events_250(N, fc, np.random.default_rng(10 + j))     # calibration events
            k, _ = best_cfg(xc, t0c, interp)
            xh, t0h, _ = make_events_250(N, fc, np.random.default_rng(500 + j))   # held-out events
            t, ok = timing(xh, k[0], k[1], interp)
            sig.append(float(sigma_ps(t, ok, t0h))); cfg.append([k[0], k[1]])
            print(f"{interp:6s} fc={fc}: f={k[0]:.3f} d={k[1]} -> {sig[-1]:.1f} ps (held-out, eff {ok.mean():.3f})")
        out["interp"][interp] = {"sigma_ps": sig, "best_f_d_samples": cfg}
    out["ref_1GSps_cfd_ps"] = REF_1GSPS_PS      # from the CFD-timing-study repository (1 GS/s, SNR 100, no filter)
    json.dump(out, open("../results_frontend250.json", "w"), indent=1)

    fig, ax = plt.subplots(figsize=(4.8, 3.1))
    xs = np.arange(len(fcs))
    ax.plot(xs, out["interp"]["linear"]["sigma_ps"], "o-", label="CFD, linear interpolation")
    ax.plot(xs, out["interp"]["cubic"]["sigma_ps"], "s-", label="CFD, cubic interpolation")
    ax.axhline(out["ref_1GSps_cfd_ps"], color="k", ls="--", lw=0.9, label="1 GS/s reference (47 ps)")
    ax.set_xticks(xs); ax.set_xticklabels(["none" if fc is None else str(fc) for fc in fcs])
    ax.set_xlabel("input low-pass cut-off fc [MHz] (2nd order)"); ax.set_ylabel("single-channel timing rms [ps]")
    ax.set_title("250 MS/s, SNR = 100: front-end bandwidth", fontsize=9)
    ax.set_yscale("log"); ax.set_yticks([30, 50, 100, 200, 400]); ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter()); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(OUT + "fig_frontend.png"); plt.close(fig)

if __name__ == "__main__":
    main()
