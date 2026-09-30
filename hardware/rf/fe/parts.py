"""Lumped component models.

Every part exposes ``z(f)``: its complex impedance (ohms) at an array of
frequencies (Hz). The loss models are deliberately simple placeholders
(constant Q, optional self-resonance / ESL). Swap in vendor data with
``TouchstonePart`` once real part numbers are chosen.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import skrf as rf


def _w(f):
    return 2 * np.pi * np.asarray(f, dtype=float)


@dataclass(frozen=True)
class PartsModel:
    """Loss assumptions applied to every inductor and capacitor.

    Defaults are rough typical values, not vendor data:
      * 0201 thin-film inductors (e.g. the LQP03 family AN923.2 uses): Q ~ 20
      * 0402 wire-wound inductors for the large VHF values: Q ~ 45
      * C0G capacitors: Q ~ 300
    Inductors at or above ``wirewound_from`` use the wire-wound Q.
    """

    q_thinfilm: float = 20.0
    q_wirewound: float = 45.0
    wirewound_from: float = 22e-9
    q_cap: float = 300.0
    cap_esl: float = 0.0
    lossless: bool = False

    def inductor(self, L: float, srf: float | None = None) -> "Inductor":
        if self.lossless:
            return Inductor(L, q=None, srf=srf)
        q = self.q_wirewound if L >= self.wirewound_from else self.q_thinfilm
        return Inductor(L, q=q, srf=srf)

    def capacitor(self, C: float) -> "Capacitor":
        if self.lossless:
            return Capacitor(C, q=None, esl=0.0)
        return Capacitor(C, q=self.q_cap, esl=self.cap_esl)


LOSSLESS = PartsModel(lossless=True)


@dataclass(frozen=True)
class Inductor:
    L: float
    q: float | None = None  # constant-Q loss; None = lossless
    srf: float | None = None  # self-resonant frequency (Hz); None = ignore

    def z(self, f):
        w = _w(f)
        z = 1j * w * self.L
        if self.q:
            z = z + w * self.L / self.q
        if self.srf:
            cp = 1 / ((2 * np.pi * self.srf) ** 2 * self.L)
            z = parallel(z, 1 / (1j * w * cp))
        return z


@dataclass(frozen=True)
class Capacitor:
    C: float
    q: float | None = None
    esl: float = 0.0

    def z(self, f):
        w = _w(f)
        z = 1 / (1j * w * self.C)
        if self.q:
            z = z + 1 / (w * self.C * self.q)
        if self.esl:
            z = z + 1j * w * self.esl
        return z


@dataclass(frozen=True)
class Resistor:
    R: float

    def z(self, f):
        return np.full(np.shape(f), self.R, dtype=complex)


@dataclass(frozen=True)
class Series:
    """Parts in series, e.g. a series-resonant trap branch."""

    parts: tuple

    def z(self, f):
        return sum(p.z(f) for p in self.parts)


@dataclass(frozen=True)
class Parallel:
    """Parts in parallel, e.g. a parallel-resonant tank."""

    parts: tuple

    def z(self, f):
        y = sum(1 / p.z(f) for p in self.parts)
        return 1 / y


@dataclass(frozen=True)
class TouchstonePart:
    """Two-terminal part taken from a vendor .s2p measured in series-thru.

    Z is extracted as Z = 2*Z0*(1 - S21)/S21 and interpolated onto the
    requested frequencies.
    """

    path: str

    def z(self, f):
        ntwk = rf.Network(self.path)
        freq = rf.Frequency.from_f(np.atleast_1d(f), unit="Hz")
        s21 = ntwk.interpolate(freq).s[:, 1, 0]
        z0 = ntwk.z0[0, 0].real
        return 2 * z0 * (1 - s21) / s21


def parallel(z1, z2):
    return z1 * z2 / (z1 + z2)


def snap_e24(x: float) -> float:
    """Nearest E24 value (component values are optimized continuously)."""
    e24 = np.array([1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0,
                    3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1])
    decade = 10 ** np.floor(np.log10(x))
    cands = np.concatenate([e24 * decade, [10 * decade]])
    return float(cands[np.argmin(np.abs(np.log(cands / x)))])
