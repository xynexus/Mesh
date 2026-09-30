# RF front-end models

Circuit models for the EFR32FG28B322F1024IM68 sub-GHz front end, built with
[scikit-rf](https://scikit-rf.org). On this dual-band OPN the sub-GHz radio has
**one TX pin (SUBG_O1) and one RX pin**, per AN923.2 §3.3. VHF (150.79–174 MHz)
and UHF (915–928 MHz) therefore share a single Class-D PA output.

The first model is the **shared-TX-pin diplexer** (`tx_diplexer.py`). The RX
side (SP3T selector → wideband LNA → broadband pin match) is not modelled yet.

```
                     ┌─ VHF arm: series-L-first LPF + 2nd-harmonic trap ───────────────── VHF SMA
SUBG_O1 ─ bond/trace ┤
                     └─ UHF arm: VHF-blocking tank ─ 915 ladder ─ SPDT ─┬─ shunt SPST ─ UHF SMA
                                                                        └─ 50 Ω term.
```

Whenever VHF is selected, the SPDT routes the UHF arm into the 50 Ω termination
and the shunt SPST shorts the UHF SMA node. The VHF 5th–7th harmonics
(754–1218 MHz) sit near the UHF passband and pass straight through the UHF arm.
The model says a single switch stage gives only about 20–30 dB of isolation,
which is not enough (see Results).

## Usage

```sh
pip install -r requirements.txt
python tx_diplexer.py validate            # reproduce the AN923.2 reference matches
python tx_diplexer.py evaluate            # Rev A values → out/: report, plots, .s3p
python tx_diplexer.py evaluate --seed     # un-optimized seed values
python tx_diplexer.py evaluate --switch spdt --coff 0.1    # compare isolation options
python tx_diplexer.py evaluate --vhf-ant ant_vhf.s1p       # EM-simulated antenna load
python tx_diplexer.py optimize            # re-tune from the seed (~5–10 min)
python -m pytest tests
```

The `.s3p` exports have ports [PA die, VHF SMA, UHF SMA], one file per
band-select state, all referenced to 50 Ω. To use them in another simulator,
renormalize the die port to the PA model.

## What is modelled, and where the numbers come from

All figures are from Silicon Labs AN923.2 Rev. 0.8. `fe/efr32.py` gives the
section references.

- **PA:** HPA on = 1.2 Ω ∥ 4 pF, off = 17 kΩ ∥ 2 pF (Table 3.1). TX bond wire
  1.2 nH, PCB trace 0.75 nH (§3.5, §3.10.1).
- **Optimum PA load at the die:** 6.9+1j Ω at 169 MHz and 4.8+1j Ω at 915 MHz
  for +20 dBm (Table 3.3). The VHF arm targets the 169 MHz value across both
  VHF sub-bands.
- **Seeds:**
  - VHF arm: the FG28 169 MHz +20 dBm BOM (Table 4.2), scaled to 162 MHz.
  - UHF arm: the 915 MHz +20 dBm BOM, including its C6 3rd-harmonic trap.
- **Reference topology:** the ladder order in Fig. 4.1 isn't legible in the text,
  so it was *inferred*. The order L1(∥C6)–C1–L2–C2–L3–C5–C3 reproduces the
  Table 3.3 impedances to within about 1.3 Ω (`validate`, `tests/`).
- **Source calibration:** the Class-D source harmonic content is fitted so the
  single-band reference match reproduces the measured conducted fundamental,
  2nd and 3rd harmonic in Table 4.3 (169 MHz: 20.4 / −36.1 / −52.1 dBm; 915 MHz:
  20.3 / −39.6 / −50.6 dBm). Powers are therefore relative to the Silicon Labs
  reference board.
- **Receive loading:** with the PA off, how much each SMA node is loaded by the
  idle TX path. This is the direct-tie condition from §3.8, with a 50 Ω RX tap
  assumed.

### Known limitations: read before trusting a number

- **Harmonics above the 3rd are extrapolated** as 1/n from the fitted 2nd
  (even) and 3rd (odd). They are not measured. The 5th–7th decide the switch
  requirement, so measure them on Rev A.
- **Component models are placeholders:** constant Q (thin-film Q = 20, wire-wound
  Q = 45 from 22 nH up, C0G Q = 300), no self-resonance, no pad or via
  parasitics. At harmonics these dominate. Swap in vendor S-parameters
  (`fe.parts.TouchstonePart`) once part numbers are chosen, or EM-simulate the
  layout.
- **Some published values aren't reproduced.** The model gives |Z| ≈ 39 Ω at the
  3rd harmonic of the 169 MHz reference, where AN923.2 quotes 116j Ω. Board
  parasitics probably account for the difference, so VHF 3rd-harmonic numbers
  are the least trustworthy.
- **The switch model is generic:** Ron = 1.5 Ω, Coff = 0.25 pF, 0.5 nH via.
  Replace it with the chosen part's data. AN923.2 §3.9 lists SKY13323-378LF,
  AS179-92LF, UPG2214TB-A and others.
- Conducted only. Radiated harmonics also depend on each antenna's efficiency
  at the harmonic, and on trace and PAVDD coupling (§3.10.1).

## Results (Rev A values, default parts and switch)

From `python tx_diplexer.py evaluate`. These are E24 values, optimized from the
seed with the default spdt+shunt switch.

| Band selected | PA load \|Γ\| vs optimum | Conducted fundamental | RX loss from idle TX path |
|---|---|---|---|
| VHF, 150.79–152.49 MHz | 0.03 | 20.6–20.7 dBm | 0.28–0.31 dB |
| VHF, 173.29–174.0 MHz | 0.10–0.11 | 20.2–20.3 dBm | 0.34–0.35 dB |
| UHF, 915–928 MHz | 0.07–0.11 | 20.6 dBm | 0.53–0.57 dB |

Fundamental figures are relative to the AN923.2 reference board, which measures
20.4 dBm at 169 MHz and 20.3 dBm at 915 MHz. So to first order, sharing the pin
costs roughly nothing at the fundamental.

Worst conducted harmonics (dBm), VHF transmitting:

| | 2f | 3f | 4f | 5f | 6f | 7f | 8f |
|---|---|---|---|---|---|---|---|
| VHF SMA | −38.7 | −41.9 | −69.7 | −65.4 | −88.9 | −88.4 | −122 |
| UHF SMA (spdt+shunt) | −73 | −54 | −63 | −43 | −58 | −47 | −67 |

What the model says so far:

1. **A single isolation switch is not enough.** With a plain SPDT (Coff
   0.25 pF), the VHF 5th harmonic reaches about −19 dBm at the UHF SMA; with only
   a shunt SPST, the 7th reaches about −15 dBm. Adding the shunt SPST after the
   SPDT gives about −43 dBm worst case, and −51 dBm if the switch Coff is
   0.1 pF. Compare options with `--switch` / `--coff`.
2. **The VHF 2nd and 3rd harmonics on the VHF SMA are the next margin to
   watch:** −38.7 and −41.9 dBm. The narrowband 169 MHz reference measures −36.1
   and −52.1 dBm. Covering 14 % bandwidth plus the UHF arm hanging on the pin
   costs 3rd-harmonic suppression. Expect to need a stronger 2nd/3rd-harmonic
   trap once real part data is in.
3. **The UHF 3rd-harmonic die impedance** is 185–225 Ω, against the 250 Ω target
   in AN923.2 §3.10.1. The single-band reference model reaches 289 Ω.

Next steps: put vendor S-parameters on the switch and the large VHF inductors,
feed antenna `.s1p` files from EM simulation, and measure the 4th–7th
harmonics on Rev A hardware to replace the extrapolation.
