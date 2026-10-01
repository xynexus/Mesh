#!/usr/bin/env python3
"""Shared-TX-pin VHF/UHF diplexer model for the EFR32FG28 (SUBG_O1).

  python tx_diplexer.py validate            # reproduce AN923.2 reference matches
  python tx_diplexer.py evaluate            # Rev A values: report, plots, .s3p
  python tx_diplexer.py evaluate --seed     # un-optimized seed values
  python tx_diplexer.py optimize            # re-tune from the seed (~5 min)

Antenna loads: --vhf-ant / --uhf-ant take a one-port Touchstone (e.g. from
openEMS/NEC) instead of the default 50 ohm.
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import asdict

import numpy as np
import skrf as rf

from fe import circuit as ckt
from fe import efr32
from fe import txdiplexer as tx
from fe.parts import PartsModel

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

# Rev A candidate: tx.optimize() from the seed with PartsModel() and Switch()
# defaults (spdt+shunt), snapped to E24.
REV_A_VHF = tx.VhfArm(L1=10e-9, C6=0.91e-12, C1=75e-12, L2=30e-9, C2=33e-12,
                      L3=47e-9, C3=7.5e-12, L3t=33e-9)
REV_A_UHF = tx.UhfArm(Lb=27e-9, Cb=39e-12, L1=0.75e-9, C6=3.6e-12, L2=0.43e-9,
                      C2=16e-12, L3=8.2e-9, C3=3.0e-12)


def _fmt_z(z):
    return f"{z.real:6.2f}{z.imag:+7.2f}j"


def _fmt_arm(a):
    parts = []
    for k, v in asdict(a).items():
        parts.append(f"{k}={v * 1e9:.3g} nH" if k.startswith("L") else f"{k}={v * 1e12:.3g} pF")
    return f"{type(a).__name__}: " + ", ".join(parts)


def cmd_validate(args):
    print("AN923.2 reference TX matches (+20 dBm HPA), die-side impedance")
    print("model = inferred Fig. 4.1 ladder + 1.2 nH bond + 0.75 nH trace\n")
    print(f"{'freq':>8}  {'target':>16}  {'lossless':>16}  {'lossy parts':>16}   |Z| at 3f (target)")
    for f0, zt in efr32.ZIN_TARGET.items():
        fr = ckt.frequency([f0, 3 * f0])
        bom = efr32.REF_BOM_20DBM[f0]
        z_ll = ckt.zin(efr32.reference_tx(fr, bom))
        z_lo = ckt.zin(efr32.reference_tx(fr, bom, PartsModel()))
        z3t = efr32.ZIN_3RD_TARGET.get(f0)
        z3 = f"{abs(z_ll[1]):.0f} ({abs(z3t):.0f})" if z3t else f"{abs(z_ll[1]):.0f} (n/a)"
        print(f"{f0 / 1e6:6.0f} MHz  {_fmt_z(zt):>16}  {_fmt_z(z_ll[0]):>16}  {_fmt_z(z_lo[0]):>16}   {z3}")
    print("\nSource calibration (fitted to measured conducted power of the reference match):")
    for src in tx.CAL_SOURCES:
        for f0 in (169e6, 915e6):
            cal = tx.calibrate(f0, PartsModel(), src)
            v_ideal = efr32.class_d_v1(efr32.VDD_HPA_868_915 if f0 > 800e6 else efr32.VDD_HPA)
            hs = ", ".join(f"{h}f {20 * np.log10(r):.1f}" for h, r in sorted(cal.ratio.items()))
            print(f"  {src:8s} {f0 / 1e6:3.0f} MHz: V1 = {cal.v1:.2f} V pk "
                  f"({20 * np.log10(cal.v1 / v_ideal):+.1f} dB vs ideal square wave); dBc: {hs}")


def _load(path):
    return rf.Network(path) if path else None


def report(vhf, uhf, parts, sw, limit, vhf_ant=None, uhf_ant=None, cal="an923") -> list[str]:
    r = tx.evaluate(vhf, uhf, parts, sw, vhf_ant, uhf_ant, cal)
    lines = [_fmt_arm(vhf), _fmt_arm(uhf),
             f"Switch: {sw.kind}, Ron {sw.ron} Ω, Coff {sw.coff * 1e12:.2f} pF; "
             f"parts: {'lossless' if parts.lossless else f'Q_L {parts.q_thinfilm:g}/{parts.q_wirewound:g}, Q_C {parts.q_cap:g}'}; "
             f"source calibration: {cal}",
             ""]
    for band in ("VHF", "UHF"):
        lines.append(f"== {band} selected, transmitting ==")
        lines.append(f"{'f MHz':>9} {'Z die':>16} {'|Γopt|':>6} {'|Z| 3f':>6} "
                     f"{'P VHF':>6} {'P UHF':>6}  (dBm, calibrated to the reference match)")
        worst = {}
        for row in r[f"{band} TX"]:
            lines.append(f"{row['f'] / 1e6:9.3f} {_fmt_z(row['z_die']):>16} {row['gamma_opt']:6.2f} "
                         f"{abs(row['z_die_3f']):6.0f} {row['p_vhf_dbm']:6.1f} {row['p_uhf_dbm']:6.1f}")
            for h, (_, pv, pu) in row["harmonics"].items():
                w = worst.setdefault(h, [-999, -999])
                w[0], w[1] = max(w[0], pv), max(w[1], pu)
        lines.append("  worst harmonic, conducted dBm at [VHF SMA / UHF SMA]:")
        for h, (pv, pu) in worst.items():
            flag = "  <-- over limit" if max(pv, pu) > limit else ""
            lines.append(f"    {h}f: {pv:7.1f} / {pu:7.1f}{flag}")
        lines.append("")
    lines.append("== Receive: loss at each node from the idle TX path (PA off) ==")
    for node, rows in r["RX"].items():
        lines.append(f"  {node}: " + ", ".join(f"{f / 1e6:.1f} MHz {loss:.2f} dB" for f, _, loss in rows))
    lines.append("")
    lines.append(f"Harmonics without reference data are extrapolated (1/n from the fitted 2nd/3rd); "
                 f"limit = {limit} dBm.")
    return lines


def export(vhf, uhf, parts, sw, tag):
    os.makedirs(OUT, exist_ok=True)
    fr = ckt.frequency(np.linspace(50e6, 3e9, 2951))
    for sel, name in ((True, "vhf_selected"), (False, "uhf_selected")):
        n = tx.build(fr, vhf, uhf, parts, sw, vhf_selected=sel)
        n.name = f"txdiplexer_{tag}_{name}"
        n.write_touchstone(os.path.join(OUT, n.name))


def cmd_evaluate(args):
    vhf, uhf = (tx.VhfArm(), tx.UhfArm()) if args.seed else (REV_A_VHF, REV_A_UHF)
    parts = PartsModel(lossless=args.lossless)
    sw = tx.Switch(kind=args.switch, coff=args.coff * 1e-12, ron=args.ron)
    lines = report(vhf, uhf, parts, sw, args.limit, _load(args.vhf_ant), _load(args.uhf_ant), args.cal)
    print("\n".join(lines))
    if not args.no_files:
        from fe import plots
        tag = ("seed" if args.seed else "revA") + ("" if args.cal == "an923" else f"_{args.cal}")
        os.makedirs(OUT, exist_ok=True)
        with open(os.path.join(OUT, f"report_{tag}.txt"), "w") as fh:
            fh.write("\n".join(lines) + "\n")
        plots.harmonic_leakage(vhf, uhf, parts, sw, os.path.join(OUT, f"harmonics_{tag}.png"), args.limit,
                               args.cal)
        plots.pa_load(vhf, uhf, parts, sw, os.path.join(OUT, f"pa_load_{tag}.png"))
        export(vhf, uhf, parts, sw, tag)
        print(f"\nwrote report, plots and .s3p files to {OUT}")


def cmd_optimize(args):
    parts = PartsModel(lossless=args.lossless)
    sw = tx.Switch(kind=args.switch, coff=args.coff * 1e-12, ron=args.ron)
    v, u, sol = tx.optimize(tx.VhfArm(), tx.UhfArm(), parts, sw)
    print(f"least_squares: status {sol.status}, cost {sol.cost:.3f}, {sol.nfev} evaluations\n")
    print("continuous:\n  " + _fmt_arm(v) + "\n  " + _fmt_arm(u))
    v, u = tx.snapped(v), tx.snapped(u)
    print("E24:\n  " + _fmt_arm(v) + "\n  " + _fmt_arm(u) + "\n")
    print("\n".join(report(v, u, parts, sw, args.limit, cal=args.cal)))
    print("\nPaste into REV_A_VHF / REV_A_UHF in tx_diplexer.py to adopt:")
    print(f"  {v!r}\n  {u!r}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("validate", help="reproduce AN923.2 reference impedances")
    for name in ("evaluate", "optimize"):
        s = sub.add_parser(name)
        s.add_argument("--switch", default="spdt+shunt", choices=["spdt+shunt", "spdt", "shunt"])
        s.add_argument("--coff", type=float, default=0.25, help="switch off capacitance, pF")
        s.add_argument("--ron", type=float, default=1.5, help="switch on resistance, ohm")
        s.add_argument("--limit", type=float, default=-36.0, help="spurious limit, dBm")
        s.add_argument("--lossless", action="store_true", help="ideal components")
        s.add_argument("--cal", default="an923", choices=tx.CAL_SOURCES,
                       help="harmonic calibration source: AN923.2 xG23 boards (conservative) "
                            "or FG28 datasheet typicals")
        if name == "evaluate":
            s.add_argument("--seed", action="store_true", help="evaluate seed values, not Rev A")
            s.add_argument("--vhf-ant", help="one-port Touchstone for the VHF SMA load")
            s.add_argument("--uhf-ant", help="one-port Touchstone for the UHF SMA load")
            s.add_argument("--no-files", action="store_true", help="print only")
    args = p.parse_args(argv)
    {"validate": cmd_validate, "evaluate": cmd_evaluate, "optimize": cmd_optimize}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
