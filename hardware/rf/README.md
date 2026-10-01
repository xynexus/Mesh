# RF front-end models

Circuit models for the EFR32FG28B322F1024IM68 sub-GHz front end, built with
[scikit-rf](https://scikit-rf.org). On this dual-band OPN the sub-GHz radio has
**one TX pin and one RX pin**: SUBG_I1 is pin 21 and SUBG_O1 is pin 23, while
pins 22/24 are RF2G4_I0/O0 (FG28 datasheet §6, AN923.2 §3.3). VHF
(150.79–174 MHz) and UHF (915–928 MHz) therefore share a single Class-D PA
output.

> **Alternative part.** The **EFR32FG28B320F1024IM68** is the same +20 dBm
> QFN68 die without the 2.4 GHz radio. On that part pins 22/24 become a second
> sub-GHz pair (SUBG_I0/SUBG_O0), and AN923.2 §3.1/§3.3 says the internal
> switch selects between the two pairs. VHF and UHF could then each have their
> own pin and match, and this whole diplexer and its harmonic-isolation switch
> would no longer be needed. The cost is losing 2.4 GHz. **Decision: keep the
> B322.** 2.4 GHz has to work, but only over about 3 m, so its antenna can be
> very poor. The 2.4 GHz port has its own pins and uses the datasheet's
> reference match (Table 5.8, with a BGS12WN6 SPDT tying the TX and RX pins
> together); it doesn't interact with the sub-GHz diplexer.
>
> All four functions are **time-shared on one transceiver**: VHF, UHF, GPS L2
> on RF1 (SUBG_I1/O1), and BLE on RF0 (RF2G4_I0/O0), per the FG28 datasheet
> §3. That has three consequences:
> - Only one path is ever active, so a single band-select state sets every
>   switch, LNA enable and the bias-tee.
> - No path receives while another band transmits. The only cross-band
>   concerns are TX harmonics leaving through idle connectors (modelled here)
>   and idle inputs surviving coupled TX power.
> - BLE can share the UHF antenna through a simple 915 MHz / 2.44 GHz
>   diplexer, because it is never active while the UHF shunt switch is
>   closed.

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
Whether one switch stage is enough depends on how clean the PA really is; the
two sources for that disagree (see Results).

## Usage

```sh
pip install -r requirements.txt
python tx_diplexer.py validate            # reproduce the AN923.2 reference matches
python tx_diplexer.py evaluate            # Rev A values → out/: report, plots, .s3p
python tx_diplexer.py evaluate --seed     # un-optimized seed values
python tx_diplexer.py evaluate --switch spdt --coff 0.1    # compare isolation options
python tx_diplexer.py evaluate --cal fg28-ds               # FG28 datasheet harmonic levels
python tx_diplexer.py evaluate --vhf-ant ant_vhf.s1p       # EM-simulated antenna load
python tx_diplexer.py optimize            # re-tune from the seed (~5–10 min)
python -m pytest tests
```

The `.s3p` exports have ports [PA die, VHF SMA, UHF SMA], one file per
band-select state, all referenced to 50 Ω. To use them in another simulator,
renormalize the die port to the PA model.

## What is modelled, and where the numbers come from

Figures are from Silicon Labs AN923.2 Rev. 0.8 and the EFR32FG28 datasheet.
`fe/efr32.py` gives the section references.

- **PA:** HPA on = 1.2 Ω ∥ 4 pF, off = 17 kΩ ∥ 2 pF (Table 3.1). TX bond wire
  1.2 nH, PCB trace 0.75 nH (§3.5, §3.10.1).
- **Optimum PA load at the die:** 6.9+1j Ω at 169 MHz and 4.8+1j Ω at 915 MHz
  for +20 dBm (Table 3.3). The VHF arm targets the 169 MHz value across both
  VHF sub-bands.
- **Seeds:**
  - VHF arm: the FG28 169 MHz +20 dBm BOM (AN923.2 Table 4.2, identical to
    datasheet Table 5.7), scaled to 162 MHz.
  - UHF arm: the 915 MHz +20 dBm BOM, including its C6 3rd-harmonic trap.
- **Reference topology:** the ladder order in Fig. 4.1 isn't legible in the text,
  so it was *inferred*. The order L1(∥C6)–C1–L2–C2–L3–C5–C3 reproduces the
  Table 3.3 impedances to within about 1.3 Ω (`validate`, `tests/`).
- **Source calibration (`--cal`):** the Class-D source harmonic content is
  fitted so the single-band reference match reproduces measured conducted
  power, so all powers are relative to the Silicon Labs reference board. Two
  sources are available, and they disagree strongly at 169 MHz even though the
  BOM is identical:

  | +20 dBm reference match | Fundamental / 2f / 3f, dBm |
  |---|---|
  | `an923` (default): AN923.2 Table 4.3, measured on xG23 boards | 169 MHz: 20.4 / −36.1 / −52.1; 915 MHz: 20.3 / −39.6 / −50.6 |
  | `fg28-ds`: FG28 datasheet Tables 4.22 and 4.15, typical | 169 MHz: 21.1 / ≤ −52.9 / ≤ −89.3; 915 MHz: 20.5 / ≤ −32.0 / −50.4 |

  The datasheet only gives the worst harmonic in each regulatory group. So only
  2f and 3f are fitted; fitting the other group members to the same bound
  gives impossible sources. As a check, `fg28-ds` predicts 169 MHz 6f–8f well
  below the datasheet's −78.3 dBm "above 1 GHz" figure. Treat `an923` as
  conservative and `fg28-ds` as optimistic.
- **Receive loading:** with the PA off, how much each SMA node is loaded by the
  idle TX path. This is the direct-tie condition from §3.8, with a 50 Ω RX tap
  assumed.

### Known limitations: read before trusting a number

- **Harmonics above the 3rd are extrapolated** as 1/n from the fitted 2nd
  (even) and 3rd (odd). They are not measured. The 5th–7th decide the switch
  requirement, so measure them on Rev A.
- **The FG28 datasheet gives no sensitivity or noise figure at 1227.6 MHz**
  (GPS L2), so the L2 receive budget rests on your own measurements.
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

## Results (Rev A values, default parts and switch, `an923` calibration)

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

Worst conducted harmonic (dBm) under each calibration and switch option:

| Calibration | Switch | VHF TX → VHF SMA | VHF TX → UHF SMA | UHF TX → UHF SMA |
|---|---|---|---|---|
| `an923` | spdt+shunt | −38.7 (2f) | −43.0 (5f) | −43.1 (2f) |
| `an923` | spdt | −38.7 (2f) | −18.8 (5f) | −42.6 (2f) |
| `an923` | shunt | −41.1 (2f) | −15.9 (5f) | −42.4 (2f) |
| `fg28-ds` | spdt+shunt | −55.5 (2f) | −74.9 (6f) | −35.5 (2f) |
| `fg28-ds` | spdt | −55.5 (2f) | −51.1 (6f) | −35.0 (2f) |
| `fg28-ds` | shunt | −57.9 (2f) | −41.0 (2f) | −34.8 (2f) |

What the model says so far:

1. **Whether one isolation switch is enough depends on the PA's real
   harmonics.**
   - With the AN923.2 (xG23) levels, a single stage is not enough. A plain
     SPDT lets the VHF 5th harmonic reach about −19 dBm at the UHF SMA. SPDT
     plus shunt gives about −43 dBm, or −51 dBm with a 0.1 pF Coff switch.
   - With the FG28 datasheet levels, a plain SPDT already gives −51 dBm.
   - Lay out SPDT plus shunt, and leave the shunt unpopulated if Rev A
     measurements allow.
2. **VHF 2nd/3rd harmonics on the VHF SMA:** −38.7 / −41.9 dBm with `an923`, but
   about −55 dBm with `fg28-ds`. With the conservative levels, covering 14 %
   bandwidth with the UHF arm on the pin costs harmonic suppression, and
   a stronger trap would be needed.
   - The UHF 2nd harmonic in the `fg28-ds` rows (about −35 dBm) mostly reflects
     the datasheet's −52.5 dBc worst-case figure. That is a non-restricted-band
     limit under FCC §15.247 (−20 dBc), so check it against your own region's
     limits.
3. **The UHF 3rd-harmonic die impedance** is 185–225 Ω, against the 250 Ω target
   in AN923.2 §3.10.1. The single-band reference model reaches 289 Ω.

Next steps: put vendor S-parameters on the switch and the large VHF inductors,
feed antenna `.s1p` files from EM simulation, and measure the 4th–7th
harmonics on Rev A hardware to replace the extrapolation.
