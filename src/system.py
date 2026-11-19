"""System-level study: 4 channels at 250 MS/s -> 4 digital CFDs -> configurable coincidence trigger.

 1. builds a 160 us, 4-channel stream with known event topologies (which channels fire, true arrival times, amplitudes)
 2. compiles text equations to truth tables (trigger_eq.py) and runs the Verilog (rtl/numeriseur_top.v) in Icarus
 3. checks the RTL bit-for-bit against a cycle-accurate Python model (CFD hits, gate decisions, scalers)
 4. checks the number of triggers against the ground truth of the injected events
 5. software side of the offer's GUI list: amplitude histograms, channel-to-channel time-difference histogram, oscillogram
Writes results_system.json, figures/fig_viewer.png, figures/fig_trigger.png.  Illustrative pulses, not the LPSC detectors."""
import json, os, subprocess, numpy as np
from scipy.signal import sosfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cfd_core import pulse
from models import TS, AMAX, D, K, THR, HOLD, cfd_hits, trigger_model, t_interp
from frontend250 import frontend_scale, FS_H, DEC
from trigger_eq import compile_equation, describe

plt.rcParams.update({"font.size": 8, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 130})
FIG = "../figures/"
SIM = "../sim/"
NS = 40000                     # samples per channel (160 us)
SNR, FC_MHZ = 100, 60
WIN = 5                                  # coincidence gate: hits up to 5 cycles (20 ns) after the first one
EXTRA = 16
CH = "ABCD"

# ----------------------------------------------------------------------------- stream generation
MASKS = [0b0111, 0b0011, 0b0101, 0b0110, 0b1011, 0b0001, 0b0010, 0b0100, 0b1000, 0b1100, 0b1111]
PROBS = np.array([0.20, 0.16, 0.06, 0.06, 0.10, 0.09, 0.08, 0.08, 0.05, 0.07, 0.05]); PROBS /= PROBS.sum()

def gen_events(rng):
    ev, t = [], 1000.0
    while t < NS * TS - 1500:
        mask = int(rng.choice(MASKS, p=PROBS))
        e = {"T": t, "mask": mask, "t0": {}, "amp": {}}
        for c in range(4):
            if (mask >> c) & 1:
                e["t0"][c] = t + rng.uniform(-6, 6)                       # channels arrive within +-6 ns
                a = {0: rng.choice([0.45, 0.75]) + rng.normal(0, 0.03),   # A: two "lines" (energy peaks)
                     1: rng.uniform(0.3, 1.0),
                     2: rng.normal(0.6, 0.05),
                     3: rng.uniform(0.3, 1.0)}[c]
                e["amp"][c] = float(np.clip(a, 0.25, 1.0))
        ev.append(e)
        t += 600 + rng.exponential(400)
    return ev

def gen_stream(events, rng):
    sos, scale = frontend_scale(FC_MHZ)
    Lh = NS * DEC
    x = np.zeros((4, NS))
    for c in range(4):
        hi = np.zeros(Lh)
        for e in events:
            if c in e["t0"]:
                k0 = int(e["t0"][c] * FS_H)
                tk = (np.arange(k0, k0 + int(200 * FS_H))) / FS_H
                hi[k0:k0 + len(tk)] += pulse(tk - e["t0"][c], e["amp"][c])
        x[c] = (sosfilt(sos, hi) * scale)[::DEC] * AMAX + rng.normal(0, AMAX / SNR, NS)
    return np.clip(np.round(x), -2048, 2047).astype(int)

# ----------------------------------------------------------------------------- RTL run
def build_rtl():
    cmd = (f"iverilog -g2005 -P tb_numeriseur.NS={NS} -P tb_numeriseur.D={D} -P tb_numeriseur.K={K} "
           f"-P tb_numeriseur.THR={THR} -P tb_numeriseur.HOLD={HOLD} -o tb_num "
           f"../tb/tb_numeriseur.v ../rtl/numeriseur_top.v ../rtl/cfd.v ../rtl/trigger4.v")
    r = subprocess.run(cmd, shell=True, cwd=SIM, capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(r.stderr)

def write_input(xi):
    with open(SIM + "input.hex", "w") as f:
        for n in range(NS):
            f.write("%03x%03x%03x%03x\n" % tuple(int(xi[c, n]) & 0xFFF for c in (3, 2, 1, 0)))

def run_rtl(tt, win):
    r = subprocess.run(f"vvp tb_num +tt={tt:04x} +win={win}", shell=True, cwd=SIM, capture_output=True, text=True)
    hits = [[] for _ in range(4)]; trig = []; hitmask = {}; scal = None
    for line in open(SIM + "rtl_log.txt"):
        p = line.split()
        if p[0] == "C":   hits[int(p[1])].append((int(p[2]), int(p[3]), int(p[4])))
        elif p[0] == "T": trig.append(int(p[1]))
        elif p[0] == "H": hitmask[int(p[1])] = int(p[2])
        elif p[0] == "S": scal = list(map(int, p[1:]))
    return hits, trig, scal

# ----------------------------------------------------------------------------- analysis helpers
def associate(hits, events):
    """hits[c] -> {event index: hit}; also counts extra / missing hits."""
    Ts = np.array([e["T"] for e in events])
    got = [dict() for _ in range(4)]; extra = 0
    for c in range(4):
        for h in hits[c]:
            j = int(np.searchsorted(Ts, t_interp(h))) - 1
            if j in got[c] or j < 0 or c not in events[j]["t0"]:
                extra += 1
            else:
                got[c][j] = h
    missing = sum(1 for j, e in enumerate(events) for c in e["t0"] if j not in got[c])
    return got, extra, missing

def main():
    rng = np.random.default_rng(20261119)
    events = gen_events(rng)
    xi = gen_stream(events, rng)
    os.makedirs(SIM, exist_ok=True)
    write_input(xi); build_rtl()
    os.makedirs('../data', exist_ok=True)
    np.savez_compressed('../data/capture.npz', xi=xi.astype(np.int16), fs_MSps=250, amax=AMAX)
    n_pulses = sum(len(e["t0"]) for e in events)
    print(f"{len(events)} events, {n_pulses} pulses on 4 channels, {NS} samples/channel")

    # python CFD model (same fixed-point model as the CFD repo) once; the RTL is compared against it below
    py_hits = cfd_hits(xi)
    equations = ["A & B", "A & B & C", "(A & B) & !D", "MULT(2, A, B, C, D)", "A XOR B", "(A & B) | (C & D)", "B & !A"]
    rows, all_exact = [], True
    for eq in equations:
        tt = compile_equation(eq)
        r_hits, r_trig, r_scal = run_rtl(tt, WIN)
        m_trig, m_scal = trigger_model(py_hits, tt, WIN, NS + EXTRA)
        exact = (r_hits == py_hits) and (r_trig == m_trig) and (r_scal == m_scal)
        all_exact &= exact
        expected = sum((tt >> e["mask"]) & 1 for e in events)
        rows.append({"equation": eq, "truth_table": f"0x{tt:04X}", "fires_on": describe(tt),
                     "expected_from_event_labels": expected, "rtl_triggers": len(r_trig),
                     "bit_exact_vs_python_model": bool(exact)})
        print(f"{eq:22s} tt=0x{tt:04X}  expected {expected:3d}  RTL {len(r_trig):3d}  bit-exact: {exact}")

    # ---------------- hit quality (CFD hits from the RTL run of the last equation == py_hits)
    got, extra, missing = associate(py_hits, events)
    print(f"hits: {sum(len(h) for h in py_hits)} for {n_pulses} pulses; extra {extra}; missing {missing}")

    # ---------------- time-difference B - A, software interpolation on the RTL outputs
    dt_est, dt_true = [], []
    for j, e in enumerate(events):
        if j in got[0] and j in got[1]:
            dt_est.append(t_interp(got[1][j]) - t_interp(got[0][j])); dt_true.append(e["t0"][1] - e["t0"][0])
    res = np.array(dt_est) - np.array(dt_true); res -= res.mean()
    sig_dt = 1000 * res.std()
    print(f"time difference B-A: {len(res)} coincident events, residual rms {sig_dt:.0f} ps")

    # ---------------- amplitude spectra (software MCA on the captured samples around each hit)
    amp = {c: [] for c in range(4)}
    for c in range(4):
        for j, h in got[c].items():
            n = h[0]
            amp[c].append(xi[c, max(0, n - 4):n + 4].max() / AMAX)

    # ---------------- window scan (python model, bit-exact with the RTL as verified above)
    tt3 = compile_equation("A & B & C")
    exp3 = sum((tt3 >> e["mask"]) & 1 for e in events)
    scan = {w: len(trigger_model(py_hits, tt3, w, NS + EXTRA)[0]) for w in range(1, 8)}

    # ---------------- figures
    ev = next(e for e in events if e["mask"] == 0b0111)
    tt_v = compile_equation("A & B & C")
    r_hits, r_trig, _ = run_rtl(tt_v, WIN)
    t_lo, t_hi = ev["T"] - 60, ev["T"] + 190
    fig = plt.figure(figsize=(7.4, 5.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.5, 1])
    ax = fig.add_subplot(gs[0, :])
    s0, s1 = int(t_lo / TS), int(t_hi / TS)
    tn = np.arange(s0, s1) * TS
    for c in range(4):
        off = 1.25 * (3 - c)
        ax.step(tn, xi[c, s0:s1] / AMAX + off, where="post", lw=0.9, label=f"channel {CH[c]}")
        ax.text(t_lo + 2, off + 0.55, CH[c], fontsize=8, fontweight="bold")
        for h in r_hits[c]:
            th = t_interp(h)
            if t_lo <= th <= t_hi:
                ax.plot([th, th], [off - 0.15, off + 0.1], color="k", lw=1.2)
    first = min(h[0] for c in range(4) for h in r_hits[c] if t_lo <= h[0] * TS <= t_hi)
    ttrig = [i for i in r_trig if t_lo <= i * TS <= t_hi][0] * TS
    ax.axvline(ttrig, color="tab:red", ls="--", lw=1.0)
    ax.text(ttrig + 3, 4.95, f"trigger (A & B & C), {ttrig - first * TS:.0f} ns after the first hit", color="tab:red", fontsize=7, va="top")
    ax.set_xlim(t_lo, t_hi); ax.set_ylim(-0.4, 5.0)
    ax.set_xlabel("t [ns]"); ax.set_yticks([])
    ax.set_title("Oscillogram of the four digitised channels (250 MS/s) with CFD hit marks (black) and trigger", fontsize=8)
    a1 = fig.add_subplot(gs[1, 0])
    a1.hist(amp[0], bins=np.linspace(0.25, 1.0, 40), color="tab:blue", alpha=0.8, label="A")
    a1.hist(amp[2], bins=np.linspace(0.25, 1.0, 40), color="tab:green", alpha=0.7, label="C")
    a1.set_xlabel("amplitude [fraction of Amax]"); a1.set_ylabel("counts"); a1.legend(fontsize=7)
    a1.set_title("Amplitude histograms", fontsize=8)
    a2 = fig.add_subplot(gs[1, 1])
    a2.hist(1000 * res, bins=30, color="tab:purple", alpha=0.85)
    a2.set_xlabel("measured - injected (B - A) [ps]"); a2.set_ylabel("counts")
    a2.set_title(f"Time difference B-A: rms {sig_dt:.0f} ps", fontsize=8)
    a3 = fig.add_subplot(gs[1, 2])
    dti, dtm = np.array(dt_true), np.array(dt_est) - np.mean(np.array(dt_est) - np.array(dt_true))
    a3.scatter(dti, dtm, s=7, color="tab:purple", alpha=0.8)
    a3.plot([-12, 12], [-12, 12], color="0.4", lw=0.8)
    a3.set_xlabel("injected t_B - t_A [ns]"); a3.set_ylabel("measured [ns]"); a3.set_title("Measured vs injected", fontsize=8)
    fig.tight_layout(); fig.savefig(FIG + "fig_viewer.png"); plt.close(fig)
    from PIL import Image                      # histogram row only, for the note (the Qt screenshot replaces the oscillogram)
    im = Image.open(FIG + "fig_viewer.png"); w, h = im.size
    im.crop((0, int(0.585 * h), w, h)).save(FIG + "fig_viewer_hist.png")

    fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.9))
    xs = np.arange(len(rows)); wbar = 0.38
    ax[0].bar(xs - wbar / 2, [r["expected_from_event_labels"] for r in rows], wbar, label="expected (event labels)", color="0.65")
    ax[0].bar(xs + wbar / 2, [r["rtl_triggers"] for r in rows], wbar, label="RTL triggers", color="tab:blue")
    ax[0].set_xticks(xs); ax[0].set_xticklabels([r["equation"].replace("MULT(2, A, B, C, D)", "MULT(2,ABCD)") for r in rows],
                                                rotation=35, ha="right", fontsize=6.5)
    ax[0].set_ylabel("triggers in 160 us"); ax[0].legend(fontsize=7); ax[0].set_title("Same bitstream, 7 equations", fontsize=8)
    ax[1].plot(list(scan), [scan[w] for w in scan], "o-")
    ax[1].axhline(exp3, color="0.5", ls="--", lw=0.9); ax[1].text(5.0, exp3 - 3.0, f"expected = {exp3}", fontsize=7, color="0.4")
    ax[1].set_xlabel("gate length WIN [cycles of 4 ns]"); ax[1].set_ylabel("A & B & C triggers")
    ax[1].set_title("Coincidence gate vs arrival spread (+-6 ns)", fontsize=8)
    fig.tight_layout(); fig.savefig(FIG + "fig_trigger.png"); plt.close(fig)

    lat_cycles = (int(round(ttrig / TS)) - first)
    out = {"stream": {"events": len(events), "pulses": n_pulses, "samples_per_channel": NS, "fs_MSps": 250,
                      "front_end_lowpass_MHz": FC_MHZ, "snr": SNR},
           "cfd": {"f": K / 8, "delay_samples": D, "arming_LSB": THR, "hold_samples": HOLD},
           "gate_win_cycles": WIN, "trigger_latency_from_first_hit_cycles": lat_cycles,
           "hits": {"total": int(sum(len(h) for h in py_hits)), "extra": extra, "missing": missing},
           "rtl_bit_exact_all_equations": bool(all_exact), "equations": rows,
           "time_difference_B_minus_A": {"coincident_events": len(res), "residual_rms_ps": float(sig_dt)},
           "window_scan_ABC": {"expected": exp3, "triggers_vs_win": {str(w): v for w, v in scan.items()}}}
    json.dump(out, open("../results_system.json", "w"), indent=1)
    print("all equations bit-exact:", all_exact, "| latency from first hit:", lat_cycles, "cycles")
    print("window scan A&B&C:", scan, "expected", exp3)

if __name__ == "__main__":
    main()
