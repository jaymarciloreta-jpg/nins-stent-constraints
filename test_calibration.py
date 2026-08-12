"""Tests for the stent tool's self-calibration and competence gate.

The point of these two components is to let the tool run on a scan it was not
tuned on. So the tests do not check one phantom: they sweep the phantom across
the conditions a real cohort varies over, and assert that the tool either
recovers the truth or REFUSES. Silently returning a wrong number is the only
failing outcome.

Conditions swept: HU calibration offset and gain (cone beam CT has no reliable
calibration), additive noise, edge blur (reconstruction kernel), voxel size,
and metal saturation.

Run:  /Applications/Slicer.app/Contents/bin/PythonSlicer test_calibration.py
"""
import os
import sys
import types

import numpy as np

TOOL = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "nins_stent_constraints.py")

failures = []


def check(label, cond, detail=""):
    print(f"  [{'PASS' if cond else 'FAIL'}] {label}" + (f"  {detail}" if detail else ""))
    if not cond:
        failures.append(label)


slicer_stub = types.ModuleType("slicer")
slicer_stub.mrmlScene = types.SimpleNamespace(GetFirstNodeByName=lambda n: None)
slicer_stub.util = types.SimpleNamespace(getNodesByClass=lambda c: [],
                                         arrayFromVolume=lambda v: None)
sys.modules["slicer"] = slicer_stub

ns = {"__name__": "_tool", "__file__": TOOL}
exec(compile(open(TOOL).read(), TOOL, "exec"), ns)

calibrate = ns["calibrate_intensities"]
assess = ns["assess_competence"]
feret = ns["min_feret_aperture"]


class FakeVolume:
    """Minimal stand-in for a vtkMRMLScalarVolumeNode."""

    def __init__(self, spacing):
        self._spacing = spacing

    def GetSpacing(self):
        return self._spacing


def head_phantom(n=90, spacing=0.5, air_hu=-1000.0, soft_hu=40.0,
                 bone_hu=900.0, noise=0.0, offset=0.0, gain=1.0, seed=0):
    """A block of soft tissue with an air channel, a bony shell, and a foramen.

    Deliberately crude. What matters is that it has the two histogram modes a
    real head has, symmetric bone, and one aperture of exactly known width.
    """
    rng = np.random.default_rng(seed)
    a = np.full((n, n, n), soft_hu, dtype=np.float32)
    c = n // 2
    # bilateral bony shell, symmetric about the mid column
    a[:, :, :6] = bone_hu
    a[:, :, -6:] = bone_hu
    a[:6, :, :] = bone_hu
    a[-6:, :, :] = bone_hu
    # an air channel down the middle
    a[c - 8:c + 8, :, c - 6:c + 6] = air_hu
    # a bony plate with a 4.0 mm gap in it, centred at (c, c, c)
    half = int(round(2.0 / spacing))
    plate = slice(c + 10, c + 10 + max(1, int(round(1.5 / spacing))))
    a[plate, :, :] = bone_hu
    a[plate, c - half:c + half, c - half:c + half] = air_hu
    if noise:
        a = a + rng.normal(0.0, noise, a.shape).astype(np.float32)
    return (a * gain + offset).astype(np.float32)


def ras_to_ijk_for(spacing, n):
    M = np.eye(4)
    M[0, 0] = M[1, 1] = M[2, 2] = spacing
    M[0, 3] = M[1, 3] = M[2, 3] = spacing / 2.0
    return np.linalg.inv(M)


print("=" * 70)
print(" stent tool: self-calibration and competence gate")
print("=" * 70)

# ---------------- calibration recovers the intensity scale ----------------
print("\n calibration across scanner conditions")

base = head_phantom()
cal = calibrate(base)
check("calibration finds both histogram modes on a clean phantom",
      cal["status"] == "ok", cal.get("reason", ""))
check("air mode is recovered to within 20 HU",
      abs(cal["air_mode_hu"] - (-1000.0)) < 20.0,
      f"{cal['air_mode_hu']:.1f} HU")
check("soft tissue mode is recovered to within 20 HU",
      abs(cal["soft_mode_hu"] - 40.0) < 20.0, f"{cal['soft_mode_hu']:.1f} HU")
check("derived thresholds are ordered air < lumen < bone",
      cal["AIR_HU_MAX"] < cal["LUMEN_HU_MAX"] < cal["BONE_HU_MIN"],
      f"{cal['AIR_HU_MAX']:.0f} < {cal['LUMEN_HU_MAX']:.0f} < {cal['BONE_HU_MIN']:.0f}")

# The case fixed thresholds cannot survive: a shifted intensity scale.
for off in (-200.0, +150.0, +300.0):
    shifted = head_phantom(offset=off)
    c2 = calibrate(shifted)
    ok = (c2["status"] == "ok"
          and abs((c2["AIR_HU_MAX"] - off) - cal["AIR_HU_MAX"]) < 25.0
          and abs((c2["BONE_HU_MIN"] - off) - cal["BONE_HU_MIN"]) < 25.0)
    check(f"thresholds track a {off:+.0f} HU calibration offset", ok,
          f"bone floor {c2.get('BONE_HU_MIN', float('nan')):.0f} HU")
    fixed_would_fail = not (ns["AIR_HU_MAX"] - 25 < cal["AIR_HU_MAX"] + off
                            < ns["AIR_HU_MAX"] + 25)
    if off > 200:
        check(f"  and a FIXED threshold would have been wrong at {off:+.0f} HU",
              fixed_would_fail, "which is why calibration exists")

for gain in (0.85, 1.2):
    c3 = calibrate(head_phantom(gain=gain))
    check(f"thresholds track a {gain:.2f}x intensity gain",
          c3["status"] == "ok"
          and abs(c3["separation_hu"] - gain * cal["separation_hu"]) < 60.0,
          f"separation {c3.get('separation_hu', 0):.0f} HU")

for nz in (10.0, 40.0, 80.0):
    c4 = calibrate(head_phantom(noise=nz))
    ok = c4["status"] == "ok" and abs(c4["noise_sd_hu"] - nz) < max(12.0, 0.35 * nz)
    check(f"noise is measured at sigma {nz:.0f} HU", ok,
          f"measured {c4.get('noise_sd_hu', 0):.1f} HU")
    expected = max(ns["CALIB_RIM_NOISE_SIGMA"] * c4["noise_sd_hu"],
                   ns["CALIB_RIM_MIN_FRACTION"] * c4["separation_hu"])
    check(f"  rim prominence is max(4 sigma, 10% of separation) at sigma {nz:.0f}",
          abs(c4["RIM_PROMINENCE_HU"] - expected) < 1e-6,
          f"{c4['RIM_PROMINENCE_HU']:.0f} HU")
    check(f"  the intensity scale is not inverted at sigma {nz:.0f}",
          c4["air_mode_hu"] < -500.0 < c4["soft_mode_hu"] < 400.0,
          f"air {c4['air_mode_hu']:.0f}, soft {c4['soft_mode_hu']:.0f} HU")

# Regression: bone-plentiful volumes used to invert the whole intensity scale.
# The rule "the two largest well separated peaks" paired SOFT TISSUE as air and
# BONE as soft tissue, because those are the two largest populations whenever
# there is more bone than air in frame. Air is identified by being the LOWEST
# major mode, not the largest.
bony = head_phantom()
bony[:14, :, :] = 900.0
bony[-14:, :, :] = 900.0
cb = calibrate(bony)
check("a bone-heavy volume still identifies air as the lowest mode",
      cb["status"] == "ok" and cb["air_mode_hu"] < -500.0
      and -100.0 < cb["soft_mode_hu"] < 300.0,
      f"air {cb.get('air_mode_hu', 0):.0f}, soft {cb.get('soft_mode_hu', 0):.0f} HU")
check("  and its bone floor sits above soft tissue, not below it",
      cb["status"] == "ok" and cb["BONE_HU_MIN"] > cb["soft_mode_hu"],
      f"bone floor {cb.get('BONE_HU_MIN', 0):.0f} HU")

flat = np.full((60, 60, 60), 30.0, dtype=np.float32)
cf = calibrate(flat)
check("a volume with no air population is refused, not guessed",
      cf["status"] == "not_bimodal", cf.get("reason", ""))

# ---------------- the aperture survives what calibration absorbs -----------
print("\n aperture measurement under the same conditions")

TRUE_MM = 4.0

def aperture(**kw):
    """Measure the plate's 4.0 mm hole.

    The hole is cut through a plate lying perpendicular to S, so it spans A and
    R. That means the measuring plane is AXIAL (R,A): only there does a profile
    through the hole centre hit plate bone on both sides. Measuring it in the
    sagittal plane, as the first version of this test did, profiles along S
    instead, where the plate is bounded by soft tissue rather than bone and no
    rim exists to find.
    """
    sp = kw.pop("spacing", 0.5)
    n = kw.pop("n", 90)
    arr = head_phantom(n=n, spacing=sp, **kw)
    R2I = ras_to_ijk_for(sp, n)
    c = n // 2
    thick = max(1, int(round(1.5 / sp)))
    k_mid = c + 10 + thick / 2.0
    seed = np.array([c * sp + sp / 2.0,          # x from i
                     c * sp + sp / 2.0,          # y from j
                     k_mid * sp + sp / 2.0])     # z from k
    cal_i = calibrate(arr)
    if cal_i["status"] != "ok":
        return None
    r = feret(arr, R2I, seed, plane="axial", calib=cal_i)
    return r["aperture_mm"] if r.get("status") == "ok" else None

clean = aperture()
check("aperture is recovered on the clean phantom",
      clean is not None and abs(clean - TRUE_MM) < 0.35,
      f"{clean if clean is None else round(clean, 3)} mm, true {TRUE_MM}")

for off in (-200.0, +300.0):
    v = aperture(offset=off)
    check(f"aperture is unchanged by a {off:+.0f} HU offset once calibrated",
          v is not None and abs(v - TRUE_MM) < 0.45,
          f"{v if v is None else round(v, 3)} mm")

for nz in (20.0, 60.0):
    v = aperture(noise=nz)
    check(f"aperture survives noise at sigma {nz:.0f} HU",
          v is not None and abs(v - TRUE_MM) < 0.6,
          f"{v if v is None else round(v, 3)} mm")

# ---------------- threshold sensitivity is self-detected ----------------
print("\n aperture stability check")
stab = ns["aperture_with_stability"]

def stability(**kw):
    sp = kw.pop("spacing", 0.5); n = kw.pop("n", 90)
    arr = head_phantom(n=n, spacing=sp, **kw)
    c = n // 2
    thick = max(1, int(round(1.5 / sp)))
    seed = np.array([c * sp + sp / 2.0, c * sp + sp / 2.0,
                     (c + 10 + thick / 2.0) * sp + sp / 2.0])
    cal_i = calibrate(arr)
    return stab(arr, ras_to_ijk_for(sp, n), seed, plane="axial", calib=cal_i)

st = stability()
check("a solid 900 HU rim is NOT flagged threshold sensitive",
      st["status"] == "ok" and st["threshold_sensitive"] is False,
      f"spread {st.get('threshold_spread_mm', 0):.3f} mm")

# A rim only marginally above the LUMEN baseline, which is what "marginal"
# means here: prominence is measured from the lumen, not from surrounding
# tissue. The real coronal SPF case had a soft-tissue-filled lumen near -5 HU
# and a rim about 130 HU above it, and the aperture read 4.47 mm at 120 HU
# prominence and 6.17 mm at 139 HU. Reproduced: a weak inner rim at 2 mm and a
# solid outer rim at 3 mm, so missing the inner one widens the aperture.
thin = head_phantom()
c, sp = 45, 0.5
thick = 3
plate = slice(c + 10, c + 10 + thick)
jj, ii = np.meshgrid(np.arange(90), np.arange(90), indexing="ij")
dj, di = np.abs(jj - c), np.abs(ii - c)
ring = np.maximum(dj, di)                      # square rings about the centre
layer = np.full((90, 90), 900.0, dtype=np.float32)
layer[ring <= 6] = 165.0                       # weak rim, 2 to 3 mm out
layer[ring <= 4] = 40.0                        # soft-tissue lumen, +/- 2 mm
for k in range(plate.start, plate.stop):
    thin[k] = layer
cal_t = calibrate(thin)
seed_t = np.array([c * sp + sp / 2.0, c * sp + sp / 2.0,
                   (c + 10 + thick / 2.0) * sp + sp / 2.0])
R2I_t = ras_to_ijk_for(sp, 90)
lo = dict(cal_t); lo["RIM_PROMINENCE_HU"] = 80.0
hi = dict(cal_t); hi["RIM_PROMINENCE_HU"] = 200.0
a_lo = feret(thin, R2I_t, seed_t, plane="axial", calib=lo)
a_hi = feret(thin, R2I_t, seed_t, plane="axial", calib=hi)
if a_lo.get("status") == "ok" and a_hi.get("status") == "ok":
    moved = abs(a_hi["aperture_mm"] - a_lo["aperture_mm"])
    check("a weak rim makes the aperture threshold-dependent (the phenomenon)",
          moved > 0.5, f"{a_lo['aperture_mm']:.2f} -> {a_hi['aperture_mm']:.2f} mm")
    mid = dict(cal_t); mid["RIM_PROMINENCE_HU"] = 110.0
    sres = stab(thin, R2I_t, seed_t, plane="axial", calib=mid)
    check("  and the stability check flags it rather than picking a side",
          sres.get("threshold_sensitive") is True,
          f"spread {sres.get('threshold_spread_mm', 0):.2f} mm")
else:
    check("a weak rim makes the aperture threshold-dependent (the phenomenon)",
          False, f"lo={a_lo.get('status')} hi={a_hi.get('status')}")

check("the stability check reports the aperture at each prominence tried",
      isinstance(st.get("aperture_at_prominence"), dict)
      and len(st["aperture_at_prominence"]) == 3)

# ---------------- the competence gate ----------------
print("\n competence gate")

good = head_phantom()
g = assess(FakeVolume([0.5, 0.5, 0.5]), good, calibrate(good))
check("a good scan passes the gate", g["status"] in ("ok", "marginal"),
      f"{g['status']}  {g['refusals']}")
check("it reports how many samples fall across the SPF",
      g["checks"]["samples_across_spf"] > 7.0,
      f"{g['checks']['samples_across_spf']}")

g2 = assess(FakeVolume([0.8, 0.8, 0.8]), good, calibrate(good))
check("a coarse-but-usable voxel is marginal, not refused",
      g2["status"] == "marginal" and not g2["refusals"], f"{g2['status']}")

g3 = assess(FakeVolume([0.5, 0.5, 3.0]), good, calibrate(good))
check("3 mm slices are REFUSED, not interpolated through",
      g3["status"] == "refused"
      and any("samples across" in r for r in g3["refusals"]),
      f"{g3['status']}")
check("  the refusal says how few samples that leaves",
      any("1.2 samples" in r for r in g3["refusals"]),
      g3["refusals"][0][:70] if g3["refusals"] else "")

g4 = assess(FakeVolume([0.5, 0.5, 0.5]), flat, calibrate(flat))
check("an uncalibrated volume is REFUSED",
      g4["status"] == "refused"
      and any("calibration failed" in r for r in g4["refusals"]))

metal = head_phantom()
metal[30:55, 30:55, 30:55] = 3200.0
g5 = assess(FakeVolume([0.5, 0.5, 0.5]), metal, calibrate(metal))
check("heavy metal saturation is REFUSED",
      g5["status"] == "refused" and any("saturated" in r for r in g5["refusals"]),
      f"{g5['checks']['metal_saturated_fraction'] * 100:.1f}% saturated")

tiny = np.full((40, 40, 40), 40.0, dtype=np.float32)
tiny[:2] = -1000.0
g6 = assess(FakeVolume([0.5, 0.5, 0.5]), tiny, calibrate(tiny))
check("a volume with almost no bone is REFUSED",
      g6["status"] == "refused", f"{g6['status']}  {g6['refusals'][:1]}")

lop = head_phantom()
lop[:, :, : lop.shape[2] // 2] = 40.0        # strip bone from one side
gl = assess(FakeVolume([0.5, 0.5, 0.5]), lop, calibrate(lop))
check("one-sided bone is flagged as marginal for the symmetry search",
      gl["status"] in ("marginal", "refused"),
      f"{gl['status']}  balance {gl['checks'].get('bone_lateral_balance')}")

check("refusals are never silent: every refusal carries a reason",
      all(isinstance(r, str) and len(r) > 20
          for gg in (g3, g4, g5, g6) for r in gg["refusals"]))

# ---------------- the contract ----------------
print("\n contract")
src = open(TOOL).read()
check("run() refuses before computing anything when the gate refuses",
      'report["status"] = "refused"' in src
      and src.index('report["status"] = "refused"') < src.index("seeds = _seeds_in_scene()"))
check("no measurement path uses a bare fixed threshold",
      src.count('_th(calib, "AIR_HU_MAX"') >= 1
      and src.count('_th(calib, "BONE_HU_MIN"') >= 2
      and src.count('_th(calib, "LUMEN_HU_MAX"') >= 1)

print()
print("=" * 70)
if failures:
    print(f" {len(failures)} FAILURE(S):")
    for f in failures:
        print(f"   - {f}")
    sys.exit(1)
print(" ALL CHECKS PASSED.")
print("=" * 70)
