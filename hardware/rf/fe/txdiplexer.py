"""Shared-TX-pin diplexer: one EFR32FG28 Class-D output (SUBG_O1) feeding a
VHF arm (150.79-174 MHz) and a UHF arm (915-928 MHz), each ending on its own
SMA node.

VHF harmonics land in and around the UHF passband (5th-7th: 754-1218 MHz),
so while VHF is selected the UHF arm must be kept off the UHF SMA. Two
options are modelled (``Switch.kind``):

* ``"shunt"``: SPST from the UHF node to ground, closed for VHF.
* ``"spdt"``: UHF arm -> SPDT common; throw A to the UHF SMA, throw B to a
  50 ohm termination. VHF selected = throw B.
* ``"spdt+shunt"`` (default): the SPDT plus an SPST from the UHF SMA node to
  ground, closed for VHF. The model needs this second stage to keep VHF
  harmonics at the UHF SMA below about -40 dBm.

The assembled network is a 3-port: 0 = PA die (before bond wire), 1 = VHF
node, 2 = UHF node (the SMA side of the switch). Antenna loads default to
50 ohm; pass impedances or one-port Networks via ``loads``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields, replace

import numpy as np
from scipy.optimize import least_squares

from . import circuit as ckt
from . import efr32
from .parts import Capacitor, Inductor, Parallel, PartsModel, Resistor, Series, snap_e24

MHZ = 1e6
VHF_F = np.array([150.7875, 151.640625, 152.49375, 173.29375, 173.646875, 174.0]) * MHZ
UHF_F = np.array([915.0, 921.5, 928.0]) * MHZ
VHF_HARMONICS = range(2, 9)  # the 8th lands on GPS L2 / Galileo E5b
UHF_HARMONICS = range(2, 4)
ZT_VHF = efr32.ZIN_TARGET[169e6]  # nearest published +20 dBm optimum
ZT_UHF = efr32.ZIN_TARGET[915e6]


@dataclass(frozen=True)
class VhfArm:
    """Series-L-first ladder, seeded from the 169 MHz +20 dBm reference
    scaled to 162 MHz, plus a 2nd-harmonic trap (L3t in series with C3) and
    a small C6 across L1 to raise its impedance near 3 x 915 MHz."""

    L1: float = 12.5e-9
    C6: float = 0.3e-12
    C1: float = 64.7e-12
    L2: float = 41.7e-9
    C2: float = 28.2e-12
    L3: float = 53.2e-9
    C5: float = 220e-12  # DC block (fixed)
    C3: float = 12.6e-12
    L3t: float = 19.0e-9  # C3-L3t series resonance ~ 325 MHz (2nd harmonic)


@dataclass(frozen=True)
class UhfArm:
    """VHF-blocking parallel tank (Lb || Cb, ~160 MHz) followed by the 915 MHz
    +20 dBm reference ladder (L1 || C6 3rd-harmonic trap, L2, C2, L3, C3)."""

    Lb: float = 33e-9
    Cb: float = 30e-12
    L1: float = 1.5e-9
    C6: float = 2.0e-12
    L2: float = 1.3e-9
    C2: float = 7.6e-12
    L3: float = 13e-9
    C5: float = 220e-12  # DC block (fixed)
    C3: float = 1.3e-12


@dataclass(frozen=True)
class Switch:
    """UHF-node switch and RX-tap off-state loading. Generic small SOI/pHEMT
    values; replace with the chosen part's data (AN923.2 §3.9 lists
    candidates: UPG2214TB-A, CG2179M2-C4, AS179-92LF, AS213-92LF,
    SKY13323-378LF)."""

    kind: str = "spdt+shunt"  # "spdt+shunt", "spdt" or "shunt"
    ron: float = 1.5
    l_via: float = 0.5e-9
    coff: float = 0.25e-12
    rx_tap_coff: float = 0.25e-12


FIXED = {"C5"}


def vhf_arm(freq, a: VhfArm, parts: PartsModel):
    ind, cap = parts.inductor, parts.capacitor
    l1 = Parallel((ind(a.L1), cap(a.C6))) if a.C6 else ind(a.L1)
    return ckt.cascade(
        ckt.series(freq, l1),
        ckt.shunt(freq, cap(a.C1)),
        ckt.series(freq, ind(a.L2)),
        ckt.shunt(freq, cap(a.C2)),
        ckt.series(freq, ind(a.L3)),
        ckt.series(freq, cap(a.C5)),
        ckt.shunt(freq, Series((cap(a.C3), ind(a.L3t)))),
    )


def uhf_arm(freq, a: UhfArm, parts: PartsModel):
    ind, cap = parts.inductor, parts.capacitor
    return ckt.cascade(
        ckt.series(freq, Parallel((ind(a.Lb), cap(a.Cb)))),
        ckt.series(freq, Parallel((ind(a.L1), cap(a.C6)))),
        ckt.series(freq, ind(a.L2)),
        ckt.shunt(freq, cap(a.C2)),
        ckt.series(freq, ind(a.L3)),
        ckt.series(freq, cap(a.C5)),
        ckt.shunt(freq, cap(a.C3)),
    )


def uhf_switch(freq, sw: Switch, vhf_selected: bool):
    """Two-port from the UHF arm output to the UHF SMA node."""
    on = Series((Resistor(sw.ron), Inductor(sw.l_via)))
    off = Capacitor(sw.coff)
    if sw.kind == "shunt":
        return ckt.shunt(freq, on if vhf_selected else off)
    if sw.kind in ("spdt", "spdt+shunt"):
        term_on = Series((on, Resistor(ckt.Z0)))
        term_off = Series((off, Resistor(ckt.Z0)))
        if vhf_selected:  # arm -> termination, SMA isolated
            n = ckt.cascade(ckt.shunt(freq, term_on), ckt.series(freq, off))
        else:
            n = ckt.cascade(ckt.shunt(freq, term_off), ckt.series(freq, Resistor(sw.ron)))
        if sw.kind == "spdt+shunt":  # extra SPST to ground on the SMA side
            n = ckt.cascade(n, ckt.shunt(freq, on if vhf_selected else off))
        return n
    raise ValueError(f"unknown switch kind {sw.kind!r}")


def build(freq, vhf: VhfArm, uhf: UhfArm, parts: PartsModel, sw: Switch,
          vhf_selected: bool):
    """3-port [die, VHF node, UHF node] for the given band-select state."""
    rx_tap = ckt.shunt(freq, Capacitor(sw.rx_tap_coff))
    v = ckt.cascade(vhf_arm(freq, vhf, parts), rx_tap)
    u = ckt.cascade(uhf_arm(freq, uhf, parts), uhf_switch(freq, sw, vhf_selected), rx_tap)
    pin = ckt.series(freq, Inductor(efr32.TX_BOND_L + efr32.TX_TRACE_L))
    return ckt.junction(freq, pin, v, u)


# --- Calibration against AN923.2 measured conducted power -----------------

@dataclass(frozen=True)
class Calibration:
    """Effective Class-D source harmonics, fitted so the single-band reference
    match reproduces measured conducted power at each harmonic for which the
    chosen source (efr32.REF_MEASURED_DBM) has data. Missing harmonics are
    extrapolated as 1/n from the fitted 2nd (even) or 3rd (odd) -- an
    assumption, not measured."""

    v1: float
    ratio: dict  # n -> |Vn| / |V1|

    def vn(self, n: int) -> float:
        if n == 1:
            return self.v1
        if n in self.ratio:
            return self.v1 * self.ratio[n]
        base = 2 if n % 2 == 0 else 3
        return self.v1 * self.ratio[base] * base / n


CAL_SOURCES = tuple(efr32.REF_MEASURED_DBM)


def calibrate(f0: float, parts: PartsModel, source: str = "an923") -> Calibration:
    meas = efr32.REF_MEASURED_DBM[source][f0]
    orders = sorted(meas)
    fr = ckt.frequency([h * f0 for h in orders])
    n = efr32.reference_tx(fr, efr32.REF_BOM_20DBM[f0], parts)
    p_unit = ckt.delivered_power(n, 1.0, efr32.pa_source_z(fr.f))[:, 1]
    v = np.sqrt(1e-3 * 10 ** (np.array([meas[h] for h in orders]) / 10) / p_unit)
    return Calibration(v1=float(v[0]), ratio={h: float(v[i] / v[0]) for i, h in enumerate(orders) if h > 1})


# --- Evaluation -----------------------------------------------------------

def eval_frequency():
    fs = [VHF_F, UHF_F]
    fs += [VHF_F * n for n in VHF_HARMONICS] + [UHF_F * n for n in UHF_HARMONICS]
    return ckt.frequency(np.concatenate(fs))


def _idx(freq, f):
    return np.array([np.argmin(np.abs(freq.f - x)) for x in np.atleast_1d(f)])


def _loads(vhf_ant, uhf_ant):
    loads = {}
    if vhf_ant is not None:
        loads[1] = vhf_ant
    if uhf_ant is not None:
        loads[2] = uhf_ant
    return loads


def evaluate(vhf: VhfArm, uhf: UhfArm, parts: PartsModel, sw: Switch = Switch(),
             vhf_ant=None, uhf_ant=None, cal: str = "an923") -> dict:
    """Per-band TX metrics (die impedance, calibrated port powers, harmonics
    at both SMAs) plus RX-mode loading of each node by the idle TX path."""
    cals = {"VHF": calibrate(169e6, parts, cal), "UHF": calibrate(915e6, parts, cal)}
    fr = eval_frequency()
    zs = efr32.pa_source_z(fr.f)
    loads = _loads(vhf_ant, uhf_ant)
    out = {}
    for band, fs, harms, zt, vhf_sel in (("VHF", VHF_F, VHF_HARMONICS, ZT_VHF, True),
                                          ("UHF", UHF_F, UHF_HARMONICS, ZT_UHF, False)):
        cal = cals[band]
        n = build(fr, vhf, uhf, parts, sw, vhf_selected=vhf_sel)
        z = ckt.zin(n, 0, loads)
        p_unit = ckt.delivered_power(n, 1.0, zs, 0, loads)  # W per V^2 (peak)
        rows = []
        for f in fs:
            i1, i3 = _idx(fr, f)[0], _idx(fr, 3 * f)[0]
            p1 = p_unit[i1] * cal.v1 ** 2
            rows.append(dict(
                f=f, z_die=z[i1], z_die_3f=z[i3],
                gamma_opt=abs(ckt.gamma_to_target(z[i1], zt)),
                p_in_dbm=ckt.dbm(p1[0]), p_vhf_dbm=ckt.dbm(p1[1]), p_uhf_dbm=ckt.dbm(p1[2]),
                harmonics={h: (h * f,
                               ckt.dbm(p_unit[_idx(fr, h * f)[0], 1] * cal.vn(h) ** 2),
                               ckt.dbm(p_unit[_idx(fr, h * f)[0], 2] * cal.vn(h) ** 2))
                           for h in harms},
            ))
        out[f"{band} TX"] = rows
    out["RX"] = rx_loading(vhf, uhf, parts, sw)
    out["calibration"] = cals
    return out


def rx_loading(vhf, uhf, parts, sw) -> dict:
    """Receive loss at each SMA node caused by the idle TX path (PA off),
    with 50 ohm antennas and a 50 ohm RX tap. The band-select state matches
    the band being received."""
    fr = ckt.frequency(np.concatenate([VHF_F, UHF_F]))
    pa_off = efr32.pa_source_z(fr.f, on=False)
    res = {}
    for name, port, band, vhf_sel in (("VHF node", 1, VHF_F, True), ("UHF node", 2, UHF_F, False)):
        n = build(fr, vhf, uhf, parts, sw, vhf_selected=vhf_sel)
        z_tx = ckt.zin(n, port, {0: pa_off})[_idx(fr, band)]
        res[name] = list(zip(band, z_tx, _rx_loss_db(z_tx)))
    return res


def _rx_loss_db(z_tx):
    zp = ckt.Z0 * z_tx / (ckt.Z0 + z_tx)
    return -20 * np.log10(np.abs(2 * zp / (ckt.Z0 + zp)))


# --- Optimization ---------------------------------------------------------

def _vars(vhf: VhfArm, uhf: UhfArm):
    keys = [("vhf", f.name) for f in fields(vhf) if f.name not in FIXED and getattr(vhf, f.name)]
    keys += [("uhf", f.name) for f in fields(uhf) if f.name not in FIXED]
    return keys


def _apply(vhf, uhf, keys, x):
    dv = {k: float(v) for (arm, k), v in zip(keys, x) if arm == "vhf"}
    du = {k: float(v) for (arm, k), v in zip(keys, x) if arm == "uhf"}
    return replace(vhf, **dv), replace(uhf, **du)


def optimize(vhf: VhfArm, uhf: UhfArm, parts: PartsModel, sw: Switch = Switch(),
             span: float = 3.0):
    """Tune all non-fixed values on a log scale within seed/span .. seed*span.

    Residuals: reflection against the optimum PA load across both VHF
    sub-bands and the UHF band; arm loss (power into the network that does
    not reach the selected SMA); VHF 2nd-harmonic transfer at least 10 dB
    below the 169 MHz reference match; |Z_die| >= 250 ohm at the UHF 3rd
    harmonic (AN923.2 §3.10.1); and RX loss from the idle TX path.
    """
    keys = _vars(vhf, uhf)
    x0 = np.log([getattr(vhf if a == "vhf" else uhf, k) for a, k in keys])
    fr = eval_frequency()
    zs = efr32.pa_source_z(fr.f)
    iv, iu = _idx(fr, VHF_F), _idx(fr, UHF_F)
    iv2, iu3 = _idx(fr, 2 * VHF_F), _idx(fr, 3 * UHF_F)
    fr_rx = ckt.frequency(np.concatenate([VHF_F, UHF_F]))
    pa_off = efr32.pa_source_z(fr_rx.f, on=False)

    ref_fr = ckt.frequency([169e6, 338e6])
    ref = efr32.reference_tx(ref_fr, efr32.REF_BOM_20DBM[169e6], parts)
    ref_pu = ckt.delivered_power(ref, 1.0, efr32.pa_source_z(ref_fr.f))[:, 1]
    h2_goal_db = 10 * np.log10(ref_pu[1] / ref_pu[0]) - 10

    def resid(x):
        v, u = _apply(vhf, uhf, keys, np.exp(x))
        r = []
        for vhf_sel, idx, zt, port in ((True, iv, ZT_VHF, 1), (False, iu, ZT_UHF, 2)):
            n = build(fr, v, u, parts, sw, vhf_selected=vhf_sel)
            z = ckt.zin(n)
            g = ckt.gamma_to_target(z[idx], zt)
            r += [g.real, g.imag]
            pu = ckt.delivered_power(n, 1.0, zs)
            r.append(2 * (1 - pu[idx, port] / pu[idx, 0]))
            if vhf_sel:
                h2 = 10 * np.log10(pu[iv2, 1] / pu[iv, 1])
                r.append(np.maximum(0, h2 - h2_goal_db) / 20)
            else:
                r.append(np.maximum(0, 1 - np.abs(z[iu3]) / efr32.Z3_MIN_3V3))
            n_rx = build(fr_rx, v, u, parts, sw, vhf_selected=vhf_sel)
            z_tx = ckt.zin(n_rx, port, {0: pa_off})[_idx(fr_rx, VHF_F if vhf_sel else UHF_F)]
            r.append(_rx_loss_db(z_tx) / 2)
        return np.concatenate([np.ravel(a) for a in r])

    sol = least_squares(resid, x0, bounds=(x0 - np.log(span), x0 + np.log(span)))
    v, u = _apply(vhf, uhf, keys, np.exp(sol.x))
    return v, u, sol


def snapped(a):
    """Copy of an arm with every value moved to the nearest E24 value."""
    return replace(a, **{k: snap_e24(val) for k, val in asdict(a).items() if val})
