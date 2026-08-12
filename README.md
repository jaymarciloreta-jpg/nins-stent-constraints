# NINS stent design constraints

Derives transnasal stent design constraints from CT sinus scans, in 3D Slicer, from **one seed click per side**.

> **v0.2.0. Not validated. Do not use these numbers in a deck, a filing, or a design freeze.**
>
> It has never been run against a deviated septum, concha bullosa, dental amalgam, or mucosal disease filling the airway. That last one is the failure mode that would most directly corrupt the headline number.
>
> This round exists to find out where it breaks. Failures are the useful output.

## The question it answers

Does a transnasal stent surface get within **5 mm** of the sphenopalatine foramen, the distance over which microcurrent is expected to conduct through tissue to the neural structures exiting the SPF? And what are the stent's dimensions?

Nothing else. It is separate from the SPG anatomy repeatability study, shares no data with it, and writes nothing into it.

## Quick start

Open a CT sinus scan in 3D Slicer, place one point fiducial named `SPF_R` anywhere inside the right sphenopalatine foramen, then in the Python Interactor:

```python
exec(open("/path/to/nins_stent_constraints.py").read())
report = run()
print_report(report)
```

Full seed list, options and QC flags are in [RUNNING.md](RUNNING.md).

## What it outputs

```
  midline: x = -12.85 mm, yaw -1.15 deg, Dice 0.728 (sharp)
  --- side R ---
    SPF aperture      axial 3.63 mm   coronal 4.47 mm
    CONDUCTION GAP    3.49 mm   (budget 5.0 mm -> WITHIN)
    corridor length   56.00 mm   (53 stations, derived)
    stent diameter    min 2.24  max 13.17 mm
    binding point     2.24 mm at 10.0 mm
    compressible      3.76 mm (bony minus mucosal)
```

| Constraint | Sets |
|---|---|
| Conduction gap, mucosal surface to SPF | whether the 5 mm budget is met. The headline |
| Lumen profile, diameter vs position | stent diameter, taper, binding point |
| Corridor length | stent length |
| Compressible margin, bony minus mucosal | how much a stent can recruit by compressing mucosa |
| Vidian and V2 clearance | current-spread safety surrogates |

## What is automatic and what is not

**Derived with no operator input, verified on a real scan:**

- **Midsagittal plane**, by reflective symmetry of the bone mask. Reproducible to **SD 0.13 mm** across 4 bone thresholds and 2 downsample factors, with a sharp optimum, and confirmed correct by rendering it.
- **SPF centre and both apertures**, by direction-free minimum feret: one seed, no direction, so the operator's angle cannot influence the result. On the reference scan **3.63 mm** against a published SPF width of 3.79 +/- 0.35 mm, where a hand-drawn chord read **6.35 mm**.
- **Corridor endpoints**, from the posterior edge of the septum and the anterior edge of the maxilla, which sit in the midline almost always and so are findable on the symmetry plane.

**Not derived: anatomical identification.** Seven approaches have been tried and all failed. The most recent was deriving midline bony landmarks by unbounded argmax on the symmetry plane, which put ANS at the chin and PNS on the sphenoid floor and read the hard palate at 66.5 mm against a published 45 to 55 mm.

The working rule: **geometry automates, recognition does not, refinement does.** Hence one seed.

## Design decisions worth knowing

- **The corridor is sampled continuously**, every 1 mm, rather than at named anatomical stations. A stent is a continuous device, so taper and binding point are what matter. This also sidesteps a live ambiguity: two CT studies disagree by 4.3 mm on where the basal lamella landmark sits.
- **The stent diameter is an axis-centred circle on a tracked centreline**, not the largest circle in the airway. The maxillary sinus is genuinely continuous with the nasal cavity through the ostium, so a global maximum reported a 25 mm stent diameter on the reference scan.
- **`NASAL_HALF_WIDTH_MM = 16.0` is a prior, not a measurement.** It is what keeps the centreline out of the maxillary sinus, and unusual anatomy will break it.
- **The bony contour is a star-shaped approximation** about the centreline, so it under-reports a re-entrant recess. It returns nothing rather than a guess when too many rays are unbounded.
- **The corridor is refused outright when there is no trustworthy midline**, because without one the two nasal cavities merge across the septum and the tool would report roughly double the true area. Refusing beats reporting wrong.

## Requirements

3D Slicer 5.6+ with its bundled Python (numpy, scipy, vtk). No external packages, no network, no deep learning.

Thresholds self-calibrate per scan from its own histogram, and scans outside the method's competence are refused rather than measured. Run `test_calibration.py` to verify both.

Runs in about 20 s on a cropped sinus CT, about 85 s on a full head, most of it the symmetry search.

## Data handling

The script reads whatever volume is loaded in Slicer and writes only what you ask it to via `write_json_dir`. **It contains no patient data and none should be committed here.** `.gitignore` excludes DICOM, NRRD, NIfTI, MRB and JSON for that reason. Keep scan data and per-subject output under the same access controls as the source DICOM.

## Feedback

For each scan: the printed block, the JSON if saved, and a one-line note on anything that looked wrong on screen. Cases where it fails are worth more than cases where it works.

## Provenance

Generated with Claude Code. Source of truth is the NINS SPG Anatomy Study vault, `03 Projects/NINS/SPG Anatomy Study/`, where this file is `06 System/(C) NINS_Stent_Constraints_0811.py` and RUNNING.md is `07 Skills/(C) Run Stent Constraints Tool.md`. The build record, including the four defects found and fixed during development, is in `03 Data & QC/(C) Auto Arm QC Log.md`.
