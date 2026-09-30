"""Small helpers on top of scikit-rf for building lumped ladder networks.

Everything is an ``skrf.Network`` referenced to Z0 = 50 ohm, so vendor
S-parameter files and EM-simulated antenna loads can be dropped in directly.
"""

from __future__ import annotations

import numpy as np
import skrf as rf
from skrf.media import DefinedGammaZ0
from skrf.network import a2s, connect

Z0 = 50.0


def frequency(f_hz) -> rf.Frequency:
    return rf.Frequency.from_f(np.sort(np.unique(np.atleast_1d(f_hz))), unit="Hz")


def series(freq: rf.Frequency, part) -> rf.Network:
    """Two-port with ``part`` in series between port 1 and port 2."""
    z = part.z(freq.f)
    a = np.zeros((len(freq), 2, 2), dtype=complex)
    a[:, 0, 0] = 1
    a[:, 0, 1] = z
    a[:, 1, 1] = 1
    return rf.Network(frequency=freq, s=a2s(a, Z0), z0=Z0)


def shunt(freq: rf.Frequency, part) -> rf.Network:
    """Two-port with ``part`` from the through line to ground."""
    y = 1 / part.z(freq.f)
    a = np.zeros((len(freq), 2, 2), dtype=complex)
    a[:, 0, 0] = 1
    a[:, 1, 0] = y
    a[:, 1, 1] = 1
    return rf.Network(frequency=freq, s=a2s(a, Z0), z0=Z0)


def thru(freq: rf.Frequency) -> rf.Network:
    return DefinedGammaZ0(frequency=freq, z0=Z0).thru()


def cascade(*ntwks: rf.Network) -> rf.Network:
    out = ntwks[0]
    for n in ntwks[1:]:
        out = out ** n
    return out


def junction(freq: rf.Frequency, left: rf.Network, *branches: rf.Network) -> rf.Network:
    """Ideal junction: ``left`` port 2 meets port 1 of each branch.

    Returns an (1 + len(branches))-port: [left port 1, branch_i port 2 ...].
    """
    media = DefinedGammaZ0(frequency=freq, z0=Z0)
    node = media.splitter(1 + len(branches))
    out = connect(left, 1, node, 0)  # ports: left1, node1..nodeN
    # connect() puts the branch's far port where the consumed port was, so
    # node port i+1 becomes branch i's far port.
    for i, br in enumerate(branches):
        out = connect(out, 1 + i, br, 0)
    return out


def load_z(freq: rf.Frequency, load) -> np.ndarray:
    """Impedance array from a scalar/array impedance or a one-port Network
    (e.g. an EM-simulated antenna .s1p, interpolated onto ``freq``)."""
    if isinstance(load, rf.Network):
        s = load.interpolate(freq).s[:, 0, 0]
        return Z0 * (1 + s) / (1 - s)
    return np.broadcast_to(np.asarray(load, dtype=complex), (len(freq),))


def _gammas(ntwk: rf.Network, loads: dict | None) -> np.ndarray:
    g = np.zeros((len(ntwk), ntwk.nports), dtype=complex)
    for port, load in (loads or {}).items():
        z = load_z(ntwk.frequency, load)
        g[:, port] = (z - Z0) / (z + Z0)
    return g


def solve(ntwk: rf.Network, src: int, vs_peak, zs, loads: dict | None = None):
    """Port waves (a, b) for a Thevenin source on ``src`` and one-port loads
    on the other ports (missing ports = Z0). Peak-amplitude normalization:
    P = 0.5 * (|a|^2 - |b|^2) into the network."""
    n = len(ntwk)
    zs = np.broadcast_to(np.asarray(zs, dtype=complex), (n,))
    vs = np.broadcast_to(np.asarray(vs_peak, dtype=complex), (n,))
    g = _gammas(ntwk, loads)
    g[:, src] = (zs - Z0) / (zs + Z0)
    c = np.zeros((n, ntwk.nports), dtype=complex)
    c[:, src] = vs * np.sqrt(Z0) / (zs + Z0)  # wave launched by the source
    # a = G b + c,  b = S a   =>   (I - S G) b = S c
    s = ntwk.s
    lhs = np.eye(ntwk.nports)[None] - s * g[:, None, :]
    b = np.linalg.solve(lhs, np.einsum("fij,fj->fi", s, c)[..., None])[..., 0]
    a = g * b + c
    return a, b


def zin(ntwk: rf.Network, port: int = 0, loads: dict | None = None) -> np.ndarray:
    """Input impedance at ``port`` with the other ports loaded (default Z0)."""
    a, b = solve(ntwk, port, 1.0, Z0, {k: v for k, v in (loads or {}).items() if k != port})
    gin = b[:, port] / a[:, port]
    return Z0 * (1 + gin) / (1 - gin)


def delivered_power(ntwk: rf.Network, vs_peak, zs, src: int = 0,
                    loads: dict | None = None) -> np.ndarray:
    """Power (W) at every port for a Thevenin source on ``src``.

    Returns shape (nfreq, nports): the source column is the power delivered
    into the network, the others the power absorbed by each port's load.
    """
    a, b = solve(ntwk, src, vs_peak, zs, loads)
    p = 0.5 * (np.abs(b) ** 2 - np.abs(a) ** 2)
    p[:, src] *= -1
    return p


def gamma_to_target(z, z_target) -> np.ndarray:
    """Power-wave reflection of ``z`` against the optimum load ``z_target``."""
    return (z - z_target) / (z + np.conj(z_target))


def dbm(p_w) -> np.ndarray:
    return 10 * np.log10(np.maximum(np.asarray(p_w, dtype=float), 1e-30) / 1e-3)
