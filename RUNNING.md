---
status: draft
tags: [skill, slicer, nins, device]
---

How to run `nins_stent_constraints.py`, which turns a few seed clicks into the NINS stent design constraint table.

**Version 0.2.0.** That is exactly what this round is for: run it on real scans with real seeds and tell me where it breaks. Do not put these numbers in a deck yet.

## What question it answers

Does a transnasal stent surface get within **5 mm** of the sphenopalatine foramen, and what are the stent's dimensions? Nothing else. It is separate from the repeatability study, shares no data with it, and writes nothing into it.

## What you place

Point fiducials, named **exactly** as below. Put each one **anywhere inside** the named structure. These are seeds, not landmarks. None of them has to sit on an edge, and the tool refines from there.

| Node name | Where | Needed for |
|---|---|---|
| `SPF_R`, `SPF_L` | anywhere inside the sphenopalatine foramen | the headline number. Required |
| `piriform_R`, `piriform_L` | inside the nasal airway at the piriform aperture, near the midline | corridor profile. **Optional**, see below |
| `choana_R`, `choana_L` | inside the airway at the posterior choana | corridor profile. **Optional**, see below |
| `vidian_R`, `vidian_L` | inside the vidian canal | safety clearance. Optional |
| `FR_R`, `FR_L` | inside foramen rotundum | safety clearance. Optional |

One side alone is fine. Missing optional seeds are reported as absent, never guessed.

**One seed is enough.** With only `SPF_R` present the tool derives the corridor endpoints itself, from the posterior edge of the septum and the anterior edge of the maxilla: per Dr. Iloreta both sit in the midline almost always, so both are findable on the symmetry plane. On the reference scan that gave a 56.0 mm corridor and a smooth 53-station profile from a single click. It is flagged `corridor_derived_*` every time, because it is validated on one scan by eye only.

**If you do place them:** keep `piriform` and `choana` in the nasal airway proper, within roughly a centimetre of the septum. Placed too far laterally they land in the maxillary sinus, which is continuous with the nasal cavity through the ostium and will drag the whole profile with it. Seeds always override the derivation.

## Running it

In the Slicer Python Interactor:

```python
exec(open("/path/to/nins_stent_constraints.py").read())
report = run()
print_report(report)
```

Options: `run(volume_node_name="...")` to pick a specific volume, `write_json_dir="/path"` to save the full record, `skip_symmetry=True` to skip the midline search (faster, but then no corridor profile: see below).

Takes about 20 s on a cropped sinus CT.

## Reading the output

```
  midline: x = -14.57 mm, yaw -1.15 deg, Dice 0.556 (sharp)

  --- side R ---
    SPF aperture      axial 3.63 mm   coronal 4.47 mm
    CONDUCTION GAP    3.49 mm   (budget 5.0 mm -> WITHIN)
    at SPF station    stent dia 5.70 mm   (airway component 94 mm2)
    corridor length   56.00 mm   (53 stations, derived)
    stent diameter    min 2.24  max 13.17 mm
    binding point     2.24 mm at 10.0 mm from the piriform seed
    compressible      3.76 mm (bony minus mucosal)
```

- **CONDUCTION GAP** is the headline: mucosal surface to the SPF centre, against the 5 mm budget.
- **binding point** is the narrowest station and where it sits along the corridor. This is the stent's limiting dimension.
- **compressible** is bony minus mucosal, the margin a stent can recruit by compressing mucosa.

## Checks you must not ignore

- `midline ... (FLAT, do not trust)` — the symmetry search found no sharp optimum. Everything side-dependent is unreliable on that scan. Usually means the volume is cropped too tightly or the head is very asymmetric. Try running on the uncropped series.
- `corridor_no_midline_*` — **the corridor profile is refused outright** when there is no trustworthy midline. Without it the two nasal cavities merge across the septum and the tool reports roughly double the true area. Refusing beats reporting wrong.
- `[JUMPY: n discontinuities]` — n stations changed area by more than 2x against their neighbour, so lumen selection switched structures. Inspect the profile before using the dimensions.
- `spf_centre_moved_*` — the refined SPF centre is more than 3 mm from your seed. Usually the seed was in an adjacent air space rather than the foramen.

## What is derived and what you still place

**Derived, no input, verified on a real scan:**
- the midsagittal plane, by reflective symmetry of the bone mask. Reproducible to **SD 0.13 mm** across 4 bone thresholds and 2 downsample factors, with a sharp optimum.
- the SPF centre and both apertures, direction-free. On the reference scan **3.63 mm** against a published SPF width of 3.79 +/- 0.35 mm, where a hand-drawn chord read 6.35 mm.
- the lumen cross section and the tracked corridor centreline.

**Not derived, and this is settled by testing rather than argument.** Seven approaches to automated anatomical identification have now failed. The most recent was deriving midline bony landmarks by argmax on the symmetry plane: it put ANS at the chin and PNS on the sphenoid floor, and read the hard palate at 66.5 mm against a published 45 to 55 mm. Geometry automates. Recognition does not. Hence the seeds.

## Known limits in v0.1.0

- Validated on **one side of one scan**, with synthetic seeds. No clinician has placed a seed for it yet.
- The corridor endpoint derivation is **n=1, checked by eye**. An unbounded version of it failed outright, putting the anterior point at the chin and the posterior point behind the sphenoid, so the working version is bounded to a height band around the SPF and to 30 mm posterior / 45 mm anterior of it. Those bounds are priors, and unusual anatomy will break them.
- Runs in about 20 s on a cropped volume, about 85 s on a full head, most of it the symmetry search.
- The bony contour is a **star-shaped approximation** about the centreline, so it under-reports a re-entrant recess. It returns nothing rather than a guess when too many rays are unbounded, which is why bony diameters are often blank in the anterior corridor.
- `NASAL_HALF_WIDTH_MM = 16.0` is a **prior, not a measurement**: it is what keeps the centreline out of the maxillary sinus. If a subject's anatomy is unusually wide this will clip the corridor, so it is worth revisiting once there are real cases.
- A handful of stations can drop out with `no_lumen_near_axis`. The profile reports how many of the requested stations survived.
- Mucosal disease filling the airway has not been tested and is the failure mode that would most directly perturb the headline number.

## What to send back

For each scan you run: the printed block, the JSON if you saved it, and a one-line note on anything that looked wrong on screen. Most useful of all are the cases where it **fails**, especially deviated septa, concha bullosa, dental amalgam, and opacified sinuses.

## Related

- `nins_stent_constraints.py`
- `06 System/(C) Autoseg Script Overview.md`, the separate QC arm for the repeatability study
- `00 Protocol & Setup/(C) Appendix A Revision Proposal.md`
- `03 Data & QC/(C) Auto Arm QC Log.md`
- [[SPG Anatomy Study MOC]]
