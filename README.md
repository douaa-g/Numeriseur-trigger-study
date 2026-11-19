# 4-channel digitiser + configurable trigger at 250 MS/s: a preliminary study

Author: Douaa Guermassi. Simulation (Python) and behavioural RTL (Verilog, verified with Icarus Verilog).

A study of the *principles* behind a 4-channel, 250 MS/s digitiser-trigger board: what timing resolution survives a 4 ns
sampling step, and how a trigger defined by a text equation can be implemented so that changing the equation does not
change the firmware. It extends my digital CFD study (github.com/douaa-g/CFD-timing-study) and was written to understand the
design space of configurable trigger/digitiser boards of the UCTM / UCTM2 family (O. Bourrion et al., arXiv:1603.06816).
Pulses are simulated (bi-exponential, 2 ns rise / 20 ns decay, illustrative), so the numbers are not those of any real detector or board.

## 1. Timing resolution at 250 MS/s (`src/frontend250.py`)

Pulses are generated at 16 GS/s, low-pass filtered (2nd-order Butterworth, cut-off fc), sampled at 250 MS/s with a random
sub-sample phase, white ADC noise added (SNR = 100). Timing by digital CFD, zero crossing by linear or cubic (4-point)
interpolation. Best (f, d) picked on one event set, reported on an independent one (3000 events each).

| fc [MHz] | none | 150 | 100 | 60 | 40 | 25 |
|---|---|---|---|---|---|---|
| CFD, linear [ps rms] | 370 | 273 | 186 | 134 | 135 | 164 |
| CFD, cubic [ps rms] | 358 | 239 | 158 | 118 | 129 | 169 |

Reference: 47 ps at 1 GS/s without filter (CFD repository). At 250 MS/s the input bandwidth matters more than the
interpolation: an optimum near 60 MHz (about 0.24 x fs); cubic interpolation gains about 12 % there.

## 2. Four channels and a configurable trigger (`rtl/`, `src/system.py`)

```
x0..x3 (12 bit, 250 MS/s) -> 4 x cfd.v -> hit[3:0] -> trigger4.v -> trig
                                                       gate + truth table + scalers
```
* `rtl/cfd.v`: fixed-point CFD (from the CFD repository), one per channel.
* `rtl/trigger4.v`: the first hit opens a gate; hits arriving up to `WIN` cycles later are accumulated; when the gate closes, the
  pattern {D,C,B,A} addresses a 16-bit truth table. Deciding once per gate makes vetoes (`A & B & !D`) and XOR correct, at the price
  of a fixed latency (`WIN`+1 cycles). 32-bit scalers count hits and triggers.
* `rtl/numeriseur_top.v`: the four CFDs, the trigger block and a tiny register interface (`tt_we`, `win_we`); the table is
  loaded at run start, so a new equation needs no new synthesis (same idea as UCTM/UCTM2).
* `src/trigger_eq.py`: equation compiler. Operators AND OR XOR NAND NOR XNOR NOT (or `& | ^ !`), parentheses, `MULT(k, A, B, ...)`.
  Example: `(A & B) & !D` -> `0x0088`; `MULT(2, A, B, C, D)` -> `0xFEE8`.

### Verification (`python3 system.py`)
A 160 us stream (40 000 samples per channel) with 162 events of known topology (which channels, arrival times within +-6 ns,
amplitudes) = 349 pulses, front-end low-pass 60 MHz, CFD f = 5/8, d = 3 samples, gate `WIN` = 5 cycles.

| Check | Result |
|---|---|
| CFD hits vs injected pulses | 349 / 349, no extra, none missing |
| RTL vs cycle-accurate Python model (hits, triggers, scalers), 7 equations | identical, bit for bit |
| Triggers vs event labels (e.g. `A & B & C`: 40 expected) | equal for all 7 equations |
| Negative test (wrong table / wrong window in the model) | comparison fails, as it should |
| Gate length for +-6 ns arrival spread | 3 cycles needed (1 cycle: 21 of 40 coincidences) |
| Time difference B - A (85 coincidences) | 193 ps rms, consistent with sqrt(2) x 134 ps |
| Trigger latency | 24 ns after the first hit (`WIN` = 5) |

Static version of the viewer outputs (oscillogram, amplitude histograms, time-difference histogram and a measured-vs-injected check): `figures/fig_viewer.png`.

### Qt viewer (`src/viewer_qt.py`, PySide6 + pyqtgraph)
```
python3 src/viewer_qt.py                      # opens data/capture.npz (the 160 us capture used by the study)
```
* equation entry with live compilation (truth table and the input patterns it fires on, error message if the syntax is wrong), gate length, scalers
* oscillogram of the four channels with CFD hit marks and trigger marks; list of triggers and previous/next navigation
* amplitude histograms and a time-difference histogram between two chosen channels, for all hits or only for gates that fired
* it uses the same models as the study (`src/models.py`, checked bit for bit against the Verilog); a scripted test reproduced the trigger counts of the 7 equations above, 21 vs 40 coincidences for gate 1 vs 5
Screenshots: `figures/qt_osc.png`, `qt_amp.png`, `qt_dt.png` (made with `QT_QPA_PLATFORM=offscreen python3 src/viewer_qt.py --screenshots figures`).
Equations tested: `A & B`, `A & B & C`, `(A & B) & !D`, `MULT(2, A, B, C, D)`, `A XOR B`, `(A & B) | (C & D)`, `B & !A`.

## Layout
```
src/frontend250.py   front-end bandwidth x interpolation study     -> figures/fig_frontend.png, results_frontend250.json
src/trigger_eq.py    text equation -> truth table (self-test: python3 trigger_eq.py)
src/system.py        stream generation, Icarus runs, bit-exact checks, viewer figures -> results_system.json, fig_viewer/fig_trigger.png
src/models.py        4-channel CFD hits and the cycle model of trigger4.v (used by system.py and the viewer)
src/viewer_qt.py     Qt viewer
src/cfd_core.py      pulse model and fixed-point CFD model (shared with the CFD repository)
data/capture.npz     the 4 x 40 000 sample capture used by the study and the viewer
rtl/                 cfd.v, trigger4.v, numeriseur_top.v          tb/tb_numeriseur.v   file-based testbench
note_numeriseur.pdf  2-page note (French)
```

## Reproduce
```
pip install numpy scipy matplotlib PySide6 pyqtgraph   # and: apt-get install iverilog
./run_all.sh
```

## Limitations
Illustrative pulses and white noise; identical shaping and CFD settings on all channels; no ADC defects (INL, clock jitter); the
sub-sample interpolation is done in software (the RTL outputs the sample index and y[n-1], y[n]); the RTL is behavioural and was only
simulated (no synthesis, no timing analysis); the viewer works on a recorded capture (no live acquisition); the gate adds a fixed trigger latency.
