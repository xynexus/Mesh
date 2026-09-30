"""Static engineering plots (PNG) for the TX diplexer model."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from . import circuit as ckt  # noqa: E402
from . import efr32  # noqa: E402
from . import txdiplexer as tx  # noqa: E402

# Validated categorical slots 1-2 (dataviz reference palette, light mode)
# plus text/grid tokens. Series keep the same colour in every chart.
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"
BAND = "#efeeea"
C_VHF = "#2a78d6"
C_UHF = "#eb6834"


def _style(ax, title, xlabel, ylabel):
    ax.set_facecolor(SURFACE)
    ax.set_title(title, loc="left", color=TEXT, fontsize=11)
    ax.set_xlabel(xlabel, color=TEXT_2)
    ax.set_ylabel(ylabel, color=TEXT_2)
    ax.tick_params(colors=TEXT_2, labelsize=9)
    ax.grid(True, color=GRID, linewidth=0.8)
    for s in ax.spines.values():
        s.set_color(GRID)


def _label_end(ax, x, y, text, color):
    ax.annotate(text, (x[-1], y[-1]), xytext=(4, 0), textcoords="offset points",
                color=TEXT, fontsize=9, va="center")
    ax.plot(x[-1], y[-1], "o", color=color, ms=4)


def harmonic_leakage(vhf, uhf, parts, sw, path, limit_dbm):
    """Conducted power at each SMA while VHF transmits (calibrated source,
    1/n harmonic extrapolation), across 100 MHz-1.4 GHz."""
    cal = tx.calibrate(169e6, parts)
    f0 = np.array([151.64e6, 173.65e6])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True, facecolor=SURFACE)
    for ax, fc in zip(axes, f0):
        hs = np.arange(1, 9)
        fr = ckt.frequency(fc * hs)
        n = tx.build(fr, vhf, uhf, parts, sw, vhf_selected=True)
        pu = ckt.delivered_power(n, 1.0, efr32.pa_source_z(fr.f))
        v2 = np.array([cal.vn(int(h)) ** 2 for h in hs])
        f_mhz = fr.f / 1e6
        for port, color, name in ((1, C_VHF, "VHF SMA"), (2, C_UHF, "UHF SMA")):
            y = ckt.dbm(pu[:, port] * v2)
            ax.plot(f_mhz, y, "-o", color=color, lw=2, ms=6, label=name)
            _label_end(ax, f_mhz, y, name, color)
        ax.axvspan(915, 928, color=BAND, zorder=0)
        ax.text(921.5, 10, "UHF\nband", ha="center", color=TEXT_2, fontsize=8)
        ax.axhline(limit_dbm, color=TEXT_2, lw=1, ls="--")
        ax.text(f_mhz[0], limit_dbm + 2, f"limit {limit_dbm:.0f} dBm", color=TEXT_2, fontsize=8)
        top = ax.secondary_xaxis("top")
        top.set_xticks(f_mhz, [f"{h}f" for h in hs])
        top.tick_params(colors=TEXT_2, labelsize=8, length=0)
        top.spines["top"].set_color(GRID)
        _style(ax, f"VHF TX at {fc / 1e6:.2f} MHz", "Frequency (MHz)", "Conducted power (dBm)")
        ax.set_ylim(-135, 30)
        ax.set_xlim(0, f_mhz[-1] * 1.12)
    axes[0].legend(frameon=False, loc="upper right", fontsize=9, labelcolor=TEXT)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def pa_load(vhf, uhf, parts, sw, path):
    """|Gamma| of the die-side impedance against the optimum PA load."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True, facecolor=SURFACE)
    spans = ((140e6, 190e6, True, tx.ZT_VHF, tx.VHF_F, C_VHF, "VHF selected"),
             (890e6, 950e6, False, tx.ZT_UHF, tx.UHF_F, C_UHF, "UHF selected"))
    for ax, (lo, hi, vsel, zt, band, color, name) in zip(axes, spans):
        fr = ckt.frequency(np.linspace(lo, hi, 201))
        n = tx.build(fr, vhf, uhf, parts, sw, vhf_selected=vsel)
        g = np.abs(ckt.gamma_to_target(ckt.zin(n), zt))
        x = fr.f / 1e6
        ax.plot(x, g, color=color, lw=2, label=name)
        _label_end(ax, x, g, name, color)
        for a, b in ((tx.VHF_F[0], tx.VHF_F[2]), (tx.VHF_F[3], tx.VHF_F[5])) if vsel else ((band[0], band[-1]),):
            ax.axvspan(a / 1e6, b / 1e6, color=BAND, zorder=0)
        _style(ax, f"PA load vs optimum ({zt.real:.1f}{zt.imag:+.1f}j Ω), {name}",
               "Frequency (MHz)", "|Γ| to optimum load")
        ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
