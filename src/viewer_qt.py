"""Qt viewer (PySide6 + pyqtgraph) for a 4-channel 250 MS/s capture.

Mirrors the software functions listed for a configurable digitiser board:
  * trigger equation entry (live compile to a truth table), gate length, scalers
  * oscillogram of the four channels with CFD hit marks and trigger marks, jump from trigger to trigger
  * amplitude histograms (optionally only for hits belonging to a gate that fired the trigger)
  * histogram of time differences between two chosen channels
The CFD and trigger are the Python models that were checked bit for bit against the Verilog (models.py).

Usage:   python3 viewer_qt.py [../data/capture.npz]
Headless screenshots (used for the README):   QT_QPA_PLATFORM=offscreen python3 viewer_qt.py --screenshots ../figures"""
import argparse, os, sys, numpy as np
from PySide6 import QtCore, QtWidgets
import pyqtgraph as pg
from models import TS, AMAX, cfd_hits, t_interp, trigger_model
from trigger_eq import compile_equation, describe

pg.setConfigOption("background", "w")
pg.setConfigOption("foreground", "k")
CH = "ABCD"
COL = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
PAIR_WINDOW_NS = 30.0
EXTRA = 16


class Viewer(QtWidgets.QMainWindow):
    def __init__(self, xi, equation="A & B & C", win=5):
        super().__init__()
        self.setWindowTitle("Numériseur 4 voies, 250 MS/s: viewer")
        self.xi = xi
        self.ns = xi.shape[1]
        self.t = np.arange(self.ns) * TS
        self.hits = cfd_hits(xi)                                    # CFD model, computed once
        self.hit_t = [np.array([t_interp(h) for h in hs]) for hs in self.hits]
        self.hit_n = [np.array([h[0] for h in hs]) for hs in self.hits]
        self.trig_idx, self.tt, self.win, self.cur = [], 0, win, -1
        self.marks = []

        # ---------------------------------------------------------------- controls
        ctrl = QtWidgets.QWidget()
        form = QtWidgets.QVBoxLayout(ctrl)
        form.addWidget(QtWidgets.QLabel("<b>Trigger equation</b> (inputs A B C D)"))
        self.eq = QtWidgets.QLineEdit(equation)
        self.eq.setPlaceholderText("e.g. (A & B) & !D   or   MULT(2, A, B, C, D)")
        form.addWidget(self.eq)
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        form.addWidget(self.status)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Gate length (cycles of 4 ns)"))
        self.win_spin = QtWidgets.QSpinBox(); self.win_spin.setRange(1, 64); self.win_spin.setValue(win)
        row.addWidget(self.win_spin)
        form.addLayout(row)
        self.apply = QtWidgets.QPushButton("Apply")
        form.addWidget(self.apply)
        form.addWidget(QtWidgets.QLabel("<b>Scalers</b>"))
        self.scal = QtWidgets.QTableWidget(5, 1)
        self.scal.setVerticalHeaderLabels(["Hits A", "Hits B", "Hits C", "Hits D", "Triggers"])
        self.scal.horizontalHeader().setVisible(False)
        self.scal.horizontalHeader().setStretchLastSection(True)
        self.scal.setMaximumHeight(160)
        form.addWidget(self.scal)
        self.only_trig = QtWidgets.QCheckBox("Histograms: only gates that fired")
        self.only_trig.setChecked(True)
        form.addWidget(self.only_trig)
        form.addWidget(QtWidgets.QLabel("<b>Triggers</b> (click to display)"))
        self.tlist = QtWidgets.QListWidget()
        form.addWidget(self.tlist, 1)
        nav = QtWidgets.QHBoxLayout()
        self.b_prev, self.b_next, self.b_full = (QtWidgets.QPushButton(s) for s in ("< trigger", "trigger >", "Full capture"))
        for b in (self.b_prev, self.b_next, self.b_full):
            nav.addWidget(b)
        form.addLayout(nav)

        # ---------------------------------------------------------------- oscillogram
        self.glw = pg.GraphicsLayoutWidget()
        self.plots = []
        for c in range(4):
            p = self.glw.addPlot(row=c, col=0)
            p.showGrid(x=True, y=True, alpha=0.25)
            p.setLabel("left", f"{CH[c]} [ADC]")
            p.setDownsampling(auto=True, mode="peak"); p.setClipToView(True)
            p.plot(self.t, xi[c], pen=pg.mkPen(COL[c], width=1))
            sc = pg.ScatterPlotItem(x=self.hit_t[c], y=[xi[c, int(n)] for n in self.hit_n[c]], symbol="t1", size=9,
                                    brush=pg.mkBrush("k"), pen=None)
            p.addItem(sc)
            if c < 3:
                p.hideAxis("bottom")
            else:
                p.setLabel("bottom", "t [ns]")
            if c:
                p.setXLink(self.plots[0])
            self.plots.append(p)

        # ---------------------------------------------------------------- amplitude histograms
        self.amp_plot = pg.PlotWidget(title="Amplitude histograms (peak sample near the CFD hit)")
        self.amp_plot.addLegend(); self.amp_plot.showGrid(x=True, y=True, alpha=0.25)
        self.amp_plot.setLabel("bottom", "amplitude [fraction of full scale]"); self.amp_plot.setLabel("left", "counts")
        self.amp_bins = QtWidgets.QSpinBox(); self.amp_bins.setRange(10, 200); self.amp_bins.setValue(50)
        amp_tab = QtWidgets.QWidget(); al = QtWidgets.QVBoxLayout(amp_tab)
        r = QtWidgets.QHBoxLayout(); r.addWidget(QtWidgets.QLabel("Bins")); r.addWidget(self.amp_bins); r.addStretch(1)
        al.addLayout(r); al.addWidget(self.amp_plot)

        # ---------------------------------------------------------------- time-difference histogram
        self.dt_plot = pg.PlotWidget(title="Time difference between two channels")
        self.dt_plot.showGrid(x=True, y=True, alpha=0.25); self.dt_plot.setLabel("left", "counts")
        self.cx, self.cy = QtWidgets.QComboBox(), QtWidgets.QComboBox()
        self.cx.addItems(list(CH)); self.cy.addItems(list(CH)); self.cx.setCurrentIndex(0); self.cy.setCurrentIndex(1)
        self.dt_bins = QtWidgets.QSpinBox(); self.dt_bins.setRange(10, 200); self.dt_bins.setValue(60)
        self.dt_info = QtWidgets.QLabel()
        dt_tab = QtWidgets.QWidget(); dl = QtWidgets.QVBoxLayout(dt_tab)
        r = QtWidgets.QHBoxLayout()
        for w in (QtWidgets.QLabel("t(Y) - t(X):  X"), self.cx, QtWidgets.QLabel("Y"), self.cy, QtWidgets.QLabel("Bins"), self.dt_bins):
            r.addWidget(w)
        r.addStretch(1); dl.addLayout(r); dl.addWidget(self.dt_plot); dl.addWidget(self.dt_info)

        self.tabs = QtWidgets.QTabWidget()
        self.tabs.addTab(self.glw, "Oscillogram"); self.tabs.addTab(amp_tab, "Amplitudes"); self.tabs.addTab(dt_tab, "Time differences")
        split = QtWidgets.QSplitter()
        split.addWidget(ctrl); split.addWidget(self.tabs); split.setStretchFactor(1, 1); split.setSizes([300, 900])
        self.setCentralWidget(split)
        self.statusBar().showMessage(f"{self.ns} samples per channel ({self.ns * TS / 1000:.0f} us), "
                                     f"hits per channel: " + ", ".join(f"{CH[c]} {len(self.hits[c])}" for c in range(4)))

        # ---------------------------------------------------------------- wiring
        self.eq.textChanged.connect(self.check_equation)
        self.eq.returnPressed.connect(self.run)
        self.apply.clicked.connect(self.run)
        self.win_spin.valueChanged.connect(self.run)
        self.only_trig.toggled.connect(self.refresh_histograms)
        for w in (self.amp_bins, self.dt_bins, self.cx, self.cy):
            (w.valueChanged if hasattr(w, "valueChanged") else w.currentIndexChanged).connect(self.refresh_histograms)
        self.tlist.currentRowChanged.connect(self.goto)
        self.b_prev.clicked.connect(lambda: self.goto(self.cur - 1))
        self.b_next.clicked.connect(lambda: self.goto(self.cur + 1))
        self.b_full.clicked.connect(lambda: self.plots[0].setXRange(0, self.ns * TS, padding=0))
        self.check_equation()
        self.run()

    # -------------------------------------------------------------------------- logic
    def check_equation(self):
        try:
            tt = compile_equation(self.eq.text())
        except ValueError as e:
            self.status.setText(f"<span style='color:#b00020'>Error: {e}</span>")
            return None
        self.status.setText(f"<span style='color:#1b5e20'>truth table 0x{tt:04X}</span><br>fires on: " + ", ".join(describe(tt)))
        return tt

    def run(self, *_):
        tt = self.check_equation()
        if tt is None:
            return
        self.tt, self.win = tt, self.win_spin.value()
        self.trig_idx, scal = trigger_model(self.hits, tt, self.win, self.ns + EXTRA)
        for i, v in enumerate(scal):
            self.scal.setItem(i, 0, QtWidgets.QTableWidgetItem(str(v)))
        self.tlist.blockSignals(True)
        self.tlist.clear()
        self.tlist.addItems([f"t = {i * TS:,.0f} ns" for i in self.trig_idx])
        self.tlist.blockSignals(False)
        for line in self.marks:
            line[0].removeItem(line[1])
        self.marks = []
        for i in self.trig_idx:
            for p in self.plots:
                ln = pg.InfiniteLine(pos=i * TS, angle=90, pen=pg.mkPen("r", width=1, style=QtCore.Qt.DashLine))
                p.addItem(ln); self.marks.append((p, ln))
        self.cur = -1
        if self.trig_idx:
            self.goto(0)
        self.refresh_histograms()

    def goto(self, k):
        if not self.trig_idx or not 0 <= k < len(self.trig_idx):
            return
        self.cur = k
        tc = self.trig_idx[k] * TS
        self.plots[0].setXRange(tc - 120, tc + 120, padding=0)
        self.tlist.blockSignals(True); self.tlist.setCurrentRow(k); self.tlist.blockSignals(False)

    def fired_mask(self, c):
        """Boolean mask over the hits of channel c: True if the hit belongs to a gate that produced a trigger.
        A trigger visible at index t closes a gate that collected the hits with sample index in [t-1-win, t-1]."""
        n = self.hit_n[c]
        ok = np.zeros(len(n), bool)
        for t in self.trig_idx:
            ok |= (n >= t - 1 - self.win) & (n <= t - 1)
        return ok

    def refresh_histograms(self, *_):
        sel = (lambda c: self.fired_mask(c)) if self.only_trig.isChecked() else (lambda c: np.ones(len(self.hit_n[c]), bool))
        # amplitudes
        self.amp_plot.clear(); self.amp_plot.addLegend()
        edges = np.linspace(0.0, 1.0, self.amp_bins.value() + 1)
        for c in range(4):
            m = sel(c)
            amps = np.array([self.xi[c, max(0, int(n) - 4):int(n) + 4].max() / AMAX for n in self.hit_n[c][m]])
            if len(amps):
                cnt, _ = np.histogram(amps, bins=edges)
                self.amp_plot.plot(edges, cnt, stepMode="center", pen=pg.mkPen(COL[c], width=2), name=f"{CH[c]} (n={len(amps)})")
        # time differences
        self.dt_plot.clear()
        x, y = self.cx.currentIndex(), self.cy.currentIndex()
        if x == y:
            self.dt_info.setText("Choose two different channels."); return
        tx, ty = self.hit_t[x][sel(x)], self.hit_t[y][sel(y)]
        d = []
        for t0 in tx:
            j = np.argmin(np.abs(ty - t0)) if len(ty) else None
            if j is not None and abs(ty[j] - t0) <= PAIR_WINDOW_NS:
                d.append(ty[j] - t0)
        d = np.array(d)
        self.dt_plot.setLabel("bottom", f"t({CH[y]}) - t({CH[x]}) [ns]")
        if len(d) < 2:
            self.dt_info.setText("Not enough coincident pairs."); return
        cnt, edges = np.histogram(d, bins=self.dt_bins.value())
        self.dt_plot.plot(edges, cnt, stepMode="center", fillLevel=0, brush=(148, 103, 189, 150), pen=pg.mkPen("#6a3d9a", width=2))
        self.dt_info.setText(f"{len(d)} pairs within +-{PAIR_WINDOW_NS:.0f} ns,  mean {d.mean():+.2f} ns,  rms {d.std():.2f} ns")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture", nargs="?", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "capture.npz"))
    ap.add_argument("--equation", default="A & B & C")
    ap.add_argument("--screenshots", metavar="DIR", help="save one PNG per tab into DIR and exit (use QT_QPA_PLATFORM=offscreen)")
    a = ap.parse_args()
    xi = np.load(a.capture)["xi"].astype(int)
    app = QtWidgets.QApplication(sys.argv[:1])
    v = Viewer(xi, a.equation)
    v.resize(1250, 760)
    v.show()
    if a.screenshots:
        os.makedirs(a.screenshots, exist_ok=True)
        for k, name in enumerate(["osc", "amp", "dt"]):
            v.tabs.setCurrentIndex(k)
            app.processEvents()
            v.grab().save(os.path.join(a.screenshots, f"qt_{name}.png"))
        return
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
