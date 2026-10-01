import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fe import circuit as ckt  # noqa: E402
from fe import efr32  # noqa: E402
from fe import txdiplexer as tx  # noqa: E402
from fe.parts import Capacitor, Inductor, PartsModel, Resistor, snap_e24  # noqa: E402


def test_series_shunt_l_network_matches_hand_calculation():
    f = 100e6
    fr = ckt.frequency([f])
    n = ckt.cascade(ckt.series(fr, Inductor(50e-9)), ckt.shunt(fr, Capacitor(20e-12)))
    w = 2 * np.pi * f
    zc = 1 / (1j * w * 20e-12)
    expected = 1j * w * 50e-9 + zc * 50 / (zc + 50)
    assert ckt.zin(n)[0] == pytest.approx(expected, rel=1e-9)


def test_junction_keeps_branch_order():
    fr = ckt.frequency([100e6])
    n = ckt.junction(fr, ckt.thru(fr), ckt.series(fr, Resistor(1000.0)),
                     ckt.thru(fr), ckt.series(fr, Resistor(200.0)))
    s = np.abs(n.s[0, 0])
    assert s[1] < s[3] < s[2]  # 1000 ohm branch, 200 ohm branch, thru


def test_power_is_conserved_in_lossless_network():
    fr = ckt.frequency([150e6, 915e6])
    n = tx.build(fr, tx.VhfArm(), tx.UhfArm(), PartsModel(lossless=True),
                 tx.Switch(kind="shunt", ron=0.0, l_via=0.0), vhf_selected=False)
    p = ckt.delivered_power(n, 1.0, 5.0)
    assert p[:, 1] + p[:, 2] == pytest.approx(p[:, 0], rel=1e-9)


def test_loaded_zin_matches_series_load():
    fr = ckt.frequency([100e6])
    n = ckt.series(fr, Resistor(10.0))
    assert ckt.zin(n, 0, {1: 30 + 5j})[0] == pytest.approx(40 + 5j)


@pytest.mark.parametrize("f0", [169e6, 868e6, 915e6])
def test_reference_match_reproduces_an923_pa_load(f0):
    """The inferred Fig. 4.1 ladder reproduces AN923.2 Table 3.3 within ~1.5 ohm."""
    fr = ckt.frequency([f0])
    z = ckt.zin(efr32.reference_tx(fr, efr32.REF_BOM_20DBM[f0]))[0]
    assert abs(z - efr32.ZIN_TARGET[f0]) < 1.5


@pytest.mark.parametrize("source", tx.CAL_SOURCES)
@pytest.mark.parametrize("f0", [169e6, 915e6])
def test_calibration_reproduces_every_measured_harmonic(source, f0):
    parts = PartsModel()
    cal = tx.calibrate(f0, parts, source)
    meas = efr32.REF_MEASURED_DBM[source][f0]
    fr = ckt.frequency([h * f0 for h in sorted(meas)])
    n = efr32.reference_tx(fr, efr32.REF_BOM_20DBM[f0], parts)
    pu = ckt.delivered_power(n, 1.0, efr32.pa_source_z(fr.f))[:, 1]
    p = [ckt.dbm(pu[i] * cal.vn(h) ** 2) for i, h in enumerate(sorted(meas))]
    assert p == pytest.approx([meas[h] for h in sorted(meas)], abs=1e-6)


def test_uhf_arm_blocks_vhf_and_vhf_arm_blocks_uhf():
    fr = ckt.frequency([151.6e6, 173.6e6, 915e6])
    parts = PartsModel(lossless=True)
    z_u = ckt.zin(tx.uhf_arm(fr, tx.UhfArm(), parts))
    z_v = ckt.zin(tx.vhf_arm(fr, tx.VhfArm(), parts))
    assert np.all(np.abs(z_u[:2]) > 100)  # vs ~7 ohm PA load
    assert abs(z_v[2]) > 40  # vs ~5 ohm PA load


def test_snap_e24():
    assert snap_e24(4.75e-9) == pytest.approx(4.7e-9)
    assert snap_e24(9.8e-12) == pytest.approx(10e-12)
