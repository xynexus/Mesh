"""EFR32xG23/xG28 sub-GHz front-end models and reference data from
Silicon Labs AN923.2 (Rev. 0.8), "EFR32 Series 2 sub-GHz Matching Guide".

Section/table numbers refer to that document. Values measured on EFR32xG23
boards; AN923.2 notes xG28 values "may slightly differ" (package, layout).
"""

from __future__ import annotations

import numpy as np

from . import circuit as ckt
from .parts import LOSSLESS, Capacitor, Inductor, Parallel, PartsModel, Resistor, Series

# --- Chip-side parasitics -------------------------------------------------
TX_BOND_L = 1.2e-9  # §3.5: TX bond wire estimate
TX_TRACE_L = 0.75e-9  # §3.10.1 step 2: external PCB trace estimate
RX_BOND_L = 1.5e-9  # §3.7
RX_DIRECT_TIE_TRACE_L = 4e-9  # §3.10.2: long inductive RX trace on radio boards

# --- PA (Table 3.1) and LNA (Table 3.4) as parallel RC ---------------------
HPA_ON = (1.2, 4e-12)
HPA_OFF = (17e3, 2e-12)
LNA_RX = (8e3, 1.2e-12)  # 10 k -> 6 k across sub-GHz; midpoint used
LNA_TX_R = 10.0  # internal switch to GND (series 10 ohm) during TX

# Regulated HPA supply (§3.4): 3.0 V, except 3.15 V in the 868/915 bands.
VDD_HPA = 3.0
VDD_HPA_868_915 = 3.15


def pa_source_z(f, on: bool = True) -> np.ndarray:
    r, c = HPA_ON if on else HPA_OFF
    return Parallel((Resistor(r), Capacitor(c))).z(f)


def class_d_v1(vdd: float) -> float:
    """Fundamental (peak) of an ideal 0..VDD square wave."""
    return 2 * vdd / np.pi


# --- Optimum PA load at the die, before the bond wire (Table 3.3) -----------
ZIN_TARGET = {
    169e6: 6.9 + 1j,  # +20 dBm HPA, 3.3 V PAVDD
    868e6: 5.7 - 2.5j,  # +20 dBm HPA
    915e6: 4.8 + 1j,  # +20 dBm HPA
}
ZIN_3RD_TARGET = {169e6: 116j, 915e6: 338j}
# §3.10.1: target |Zin| at the 3rd harmonic, 3.3 V supply.
Z3_MIN_3V3 = 250.0

# --- Reference BOMs, EFR32xG28 +20 dBm (Table 4.2) -------------------------
# Designators per Fig. 4.1. The ladder order (L1-C1-L2-C2-L3-C5-C3 with C6
# across L1, L4/C4 = RX match) was inferred by fitting Table 3.3; see
# tests/test_reference.py. None = not mounted.
REF_BOM_20DBM = {
    169e6: dict(L1=12e-9, L2=40e-9, L3=51e-9, L4=220e-9, C1=62e-12, C2=27e-12,
                C3=16e-12, C4=2.6e-12, C5=220e-12, C6=None),
    915e6: dict(L1=1.5e-9, L2=1.3e-9, L3=13e-9, L4=18e-9, C1=None, C2=7.6e-12,
                C3=1.3e-12, C4=None, C5=220e-12, C6=2.0e-12),
}
REF_BOM_20DBM[868e6] = REF_BOM_20DBM[915e6]

# Conducted measurements, +20 dBm (Table 4.3): (P_fund, P_2nd, P_3rd) dBm.
REF_MEASURED_DBM = {
    169e6: (20.4, -36.1, -52.1),
    915e6: (20.3, -39.6, -50.6),
}


def reference_tx(freq, bom: dict, parts: PartsModel = LOSSLESS):
    """AN923.2 single-band TX-RX direct-tie match in TX mode.

    Two-port: port 1 = PA die (before bond wire), port 2 = antenna (50 ohm).
    The RX match hangs off the antenna node with the LNA switched to 10 ohm.
    """
    ind, cap = parts.inductor, parts.capacitor
    l1 = ind(bom["L1"])
    if bom.get("C6"):
        l1 = Parallel((l1, cap(bom["C6"])))
    chain = [ckt.series(freq, Inductor(TX_BOND_L + TX_TRACE_L)), ckt.series(freq, l1)]
    for name, kind in (("C1", "sh"), ("L2", "se"), ("C2", "sh"), ("L3", "se"),
                       ("C5", "se"), ("C3", "sh")):
        v = bom.get(name)
        if v is None:
            continue
        part = ind(v) if name.startswith("L") else cap(v)
        chain.append(ckt.series(freq, part) if kind == "se" else ckt.shunt(freq, part))
    rx = rx_branch_tx_mode(bom["L4"], bom.get("C4"), parts)
    chain.append(ckt.shunt(freq, rx))
    return ckt.cascade(*chain)


def rx_branch_tx_mode(l4: float, c4: float | None, parts: PartsModel = LOSSLESS):
    """Direct-tie RX match as seen from the antenna node while transmitting."""
    lna = Resistor(LNA_TX_R)
    if c4:
        lna = Parallel((lna, parts.capacitor(c4)))
    return Series((Inductor(RX_DIRECT_TIE_TRACE_L + RX_BOND_L), parts.inductor(l4), lna))
