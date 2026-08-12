"""NINS stent design constraints from a minimal seed set.

Run inside 3D Slicer (Python Interactor or the slicer MCP).

WHAT THIS IS
------------
A focused tool for one question: does a transnasal stent surface get within
5 mm of the sphenopalatine foramen, and what are the stent's dimensions?

It is NOT the repeatability study. It shares no data with the manual arm and
writes nothing into it. It exists to turn a handful of seed clicks into the
device design constraint table.

WHAT THE OPERATOR PLACES
------------------------
Point fiducials, named exactly. Anywhere INSIDE the named structure: these are
seeds, not landmarks, and none of them has to sit on an edge.

  REQUIRED for the headline number
    SPF_R, SPF_L          one point inside the sphenopalatine foramen

  REQUIRED for length and the lumen profile
    piriform_R, piriform_L    inside the nasal airway at the piriform aperture
    choana_R,   choana_L      inside the airway at the posterior choana

  OPTIONAL, safety clearances
    vidian_R, vidian_L        inside the vidian canal
    FR_R,     FR_L            inside foramen rotundum

Only one side need be present. Missing optional seeds are reported as absent,
never guessed.

WHAT IS DERIVED, AND WHAT IS NOT
--------------------------------
Derived with no operator input, and verified on a real scan:
  - the midsagittal plane, by reflective symmetry of the bone mask.
    Reproducible to SD 0.13 mm across 4 bone thresholds and 2 downsample
    factors, with a sharp optimum (Dice falls 0.037 at 1 mm off-plane).
  - the SPF centre and both apertures, by direction-free minimum feret.
  - the lumen cross section at any plane.

NOT derived, established by testing rather than argument:
  - anatomical identification. Seven approaches have now failed, most recently
    deriving midline bony landmarks by argmax on the symmetry plane: it put ANS
    at the chin and PNS on the sphenoid floor, and read the hard palate at
    66.5 mm against a published 45 to 55 mm. "Most anterior midline bony point
    OF THE MAXILLA" needs to know what the maxilla is. That is recognition.
    Geometry automates. Recognition does not. Hence the seeds above.

HOW TO RUN
----------
    exec(open("(C) NINS_Stent_Constraints_0811.py").read())
    report = run()                       # or run(volume_node_name="...")
    print_report(report)

Add write_json_dir="/path" to persist the full record.
"""

import json
import os
import time
import datetime
import hashlib

import numpy as np
import vtk
import slicer
from scipy import ndimage

SCRIPT_VERSION = "0.1.0"

# ---- tunables. Every uppercase global is captured into the provenance record.
AIR_HU_MAX = -350.0          # airway lumen, mucosal surface
BONE_HU_MIN = 300.0          # cortical bone
RIM_PROMINENCE_HU = 120.0    # a rim must rise this far above the lumen baseline
RIM_SEARCH_MM = 8.0          # how far outward to look for a bounding rim
LUMEN_HU_MAX = 200.0         # above this a sample point is not lumen
FERET_N_ANGLES = 90          # 2 degree steps over a half turn
PROFILE_STEP_MM = 1.0        # lumen sampling interval along the corridor
BONY_N_RAYS = 180            # rays for the bony corridor contour
SYMMETRY_DOWNSAMPLE = 4
SYMMETRY_BONE_HU = 300.0
CONDUCTION_BUDGET_MM = 5.0   # the design budget this study tests
# A PRIOR, not a measurement: the nasal airway proper lies within this distance
# of the midline, while the maxillary sinus lies lateral to it. Needed because
# the sinus is genuinely continuous with the cavity through the maxillary
# ostium, and it has a far larger inscribed circle, so an unconstrained
# centreline tracker migrates into it and reports the sinus as the corridor.
NASAL_HALF_WIDTH_MM = 16.0
AXIS_STEP_SEARCH_MM = 3.0    # how far the centreline may move per 1 mm station
# Bounds for deriving the corridor endpoints, relative to the SPF centre.
CORRIDOR_POST_SEARCH_MM = 30.0   # posterior septal edge lies within this of the SPF
CORRIDOR_ANT_SEARCH_MM = 45.0    # anterior maxillary edge likewise

SEED_ROLES = ("SPF", "piriform", "choana", "vidian", "FR")
REQUIRED_SEEDS = ("SPF",)
CORRIDOR_SEEDS = ("piriform", "choana")


# --------------------------------------------------------------- provenance

def collect_params():
    return {k: v for k, v in globals().items()
            if k.isupper() and not k.startswith("_")
            and isinstance(v, (int, float, str, bool, tuple, list))}


def config_hash(params):
    blob = json.dumps(params, sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


# ------------------------------------------------------------ volume access

def _volume_arrays(volume_node):
    array = slicer.util.arrayFromVolume(volume_node).astype(np.float32)
    mat = vtk.vtkMatrix4x4()
    volume_node.GetIJKToRASMatrix(mat)
    M = np.array([[mat.GetElement(r, c) for c in range(4)] for r in range(4)])
    return array, M, np.linalg.inv(M)


def _sample(array, ras_to_ijk, points):
    """Trilinear sample at RAS points (N,3)."""
    pts = np.atleast_2d(np.asarray(points, float))
    ijk = (np.c_[pts, np.ones(len(pts))] @ ras_to_ijk.T)[:, :3]
    k, j, i = ijk[:, 2], ijk[:, 1], ijk[:, 0]
    k0 = np.floor(k).astype(int); j0 = np.floor(j).astype(int)
    i0 = np.floor(i).astype(int)
    out = np.zeros(len(pts))
    K, J, I = array.shape
    for dk in (0, 1):
        for dj in (0, 1):
            for di in (0, 1):
                w = (np.abs(1 - dk - (k - k0)) * np.abs(1 - dj - (j - j0))
                     * np.abs(1 - di - (i - i0)))
                out += w * array[np.clip(k0 + dk, 0, K - 1),
                                 np.clip(j0 + dj, 0, J - 1),
                                 np.clip(i0 + di, 0, I - 1)]
    return out


def _profile(array, ras_to_ijk, origin, direction, length_mm, step_mm=0.1):
    direction = np.asarray(direction, float)
    direction = direction / np.linalg.norm(direction)
    t = np.arange(0.0, length_mm + 1e-9, step_mm)
    pts = np.asarray(origin, float)[None, :] + t[:, None] * direction[None, :]
    return t, _sample(array, ras_to_ijk, pts)


# ------------------------------------------------- midsagittal by symmetry

def find_midsagittal(array, ijk_to_ras, downsample=SYMMETRY_DOWNSAMPLE,
                     bone_hu=SYMMETRY_BONE_HU):
    """Midsagittal plane as the reflection that best maps the bone mask onto itself.

    Bone is used rather than the whole image because soft tissue, the table and
    the head holder are not symmetric. Scored by Dice on a soft occupancy mask,
    which is bounded and gives the search a gradient to follow.

    Verified on a real scan: SD 0.13 mm across 4 thresholds x 2 downsample
    factors, and a sharp optimum. Returns None if the optimum is flat, which
    means the scan is too asymmetric or too cropped to trust.
    """
    d = (array > bone_hu).astype(np.float32)[::downsample, ::downsample,
                                             ::downsample]
    if d.max() <= 0:
        return None
    d = d / d.max()
    K, J, I = d.shape
    jc, kc = J / 2.0, K / 2.0
    kk, jj, ii = np.meshgrid(np.arange(K), np.arange(J), np.arange(I),
                             indexing="ij")

    def score(x0, yaw, roll):
        xp = x0 + yaw * (jj - jc) + roll * (kk - kc)
        mi = 2.0 * xp - ii
        valid = (mi >= 0) & (mi <= I - 1)
        if not valid.any():
            return -1.0
        mic = np.clip(mi, 0, I - 1)
        f0 = np.floor(mic).astype(int)
        f1 = np.minimum(f0 + 1, I - 1)
        frac = mic - f0
        mirror = d[kk, jj, f0] * (1 - frac) + d[kk, jj, f1] * frac
        a, b = d[valid], mirror[valid]
        if a.sum() < 1 or b.sum() < 1:
            return -1.0
        return float(2 * (a * b).sum() / (a.sum() + b.sum()))

    best = None
    for x0 in np.arange(I * 0.35, I * 0.65, 1.0):
        s = score(x0, 0.0, 0.0)
        if best is None or s > best[0]:
            best = (s, x0, 0.0, 0.0)
    for _ in range(3):
        s0, x0, yaw, roll = best
        for dx in (-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0):
            for dy in (-0.06, -0.02, 0.0, 0.02, 0.06):
                for dr in (-0.06, -0.02, 0.0, 0.02, 0.06):
                    s = score(x0 + dx, yaw + dy, roll + dr)
                    if s > best[0]:
                        best = (s, x0 + dx, yaw + dy, roll + dr)
    s, x0, yaw, roll = best

    # Sharpness. A flat optimum means the estimate is arbitrary, so refuse it.
    spacing_ds = abs(ijk_to_ras[0, 0]) * downsample or 1.0
    drop_1mm = s - max(score(x0 - 1.0 / spacing_ds, yaw, roll),
                       score(x0 + 1.0 / spacing_ds, yaw, roll))
    point = ijk_to_ras @ np.array([x0 * downsample, J * downsample / 2.0,
                                   K * downsample / 2.0, 1.0])
    return {"dice": s,
            "x_ras_mm": float(point[0]),
            "yaw_deg": float(np.degrees(np.arctan(yaw))),
            "roll_deg": float(np.degrees(np.arctan(roll))),
            "dice_drop_at_1mm": float(drop_1mm),
            "sharp_enough": bool(drop_1mm > 0.01),
            "downsample": int(downsample), "bone_hu": float(bone_hu)}


# ------------------------------------------------- direction-free aperture

def _lumen_baseline(array, ras_to_ijk, seed, radius_mm=0.6):
    vals = []
    for axis in np.eye(3):
        for sgn in (1.0, -1.0):
            _, hu = _profile(array, ras_to_ijk, seed, axis * sgn, radius_mm)
            vals.extend(hu.tolist())
    return float(np.median(vals)) if vals else None


def _rim_distance(array, ras_to_ijk, origin, direction, baseline,
                  prominence=RIM_PROMINENCE_HU, search_mm=RIM_SEARCH_MM):
    t, hu = _profile(array, ras_to_ijk, origin, direction, search_mm)
    peak = None
    for i in range(1, len(hu)):
        if hu[i] - baseline >= prominence and (peak is None or hu[i] > hu[peak]):
            peak = i
        elif peak is not None and hu[i] < baseline + 0.3 * prominence:
            break
    if peak is None:
        return None
    level = baseline + 0.5 * (hu[peak] - baseline)
    for i in range(peak, 0, -1):
        if hu[i - 1] < level <= hu[i]:
            frac = (level - hu[i - 1]) / (hu[i] - hu[i - 1])
            return float(t[i - 1] + frac * (t[i] - t[i - 1]))
    return None


def min_feret_aperture(array, ras_to_ijk, seed, plane="axial",
                       n_angles=FERET_N_ANGLES):
    """Narrowest rim-to-rim aperture through `seed`, over all in-plane angles.

    Takes ONE seed and no direction, so the operator's angle cannot influence
    the result. Rims are found by LOCAL prominence above the lumen baseline,
    not an absolute threshold: the SPF's posterior rim is thin enough that
    partial volume caps it near 227 HU, and a 300 HU threshold looks straight
    through it and wrongly reports the foramen as open.

    On the reference scan this read 3.63 mm where a hand-drawn chord read
    6.35 mm, against a published SPF width of 3.79 +/- 0.35 mm.
    """
    seed = np.asarray(seed, float)
    baseline = _lumen_baseline(array, ras_to_ijk, seed)
    if baseline is None or baseline >= LUMEN_HU_MAX:
        return {"status": "seed_not_in_lumen", "baseline_hu": baseline}
    best = None
    for ang in np.linspace(0.0, np.pi, int(n_angles), endpoint=False):
        if plane == "axial":
            u = np.array([np.cos(ang), np.sin(ang), 0.0])
        elif plane == "coronal":
            u = np.array([np.cos(ang), 0.0, np.sin(ang)])
        else:
            u = np.array([0.0, np.cos(ang), np.sin(ang)])
        dp = _rim_distance(array, ras_to_ijk, seed, u, baseline)
        dn = _rim_distance(array, ras_to_ijk, seed, -u, baseline)
        if dp is None or dn is None:
            continue
        total = dp + dn
        if best is None or total < best["aperture_mm"]:
            best = {"aperture_mm": total,
                    "angle_deg": float(np.degrees(ang)),
                    "p1_ras": (seed + u * dp).tolist(),
                    "p2_ras": (seed - u * dn).tolist()}
    if best is None:
        return {"status": "no_bounded_direction", "baseline_hu": baseline}
    best.update({"status": "ok", "plane": plane, "baseline_hu": baseline})
    return best


def spf_centre(array, ras_to_ijk, seed):
    """SPF centre and both apertures from one seed."""
    ax = min_feret_aperture(array, ras_to_ijk, seed, plane="axial")
    co = min_feret_aperture(array, ras_to_ijk, seed, plane="coronal")
    if ax.get("status") != "ok" or co.get("status") != "ok":
        return {"status": "aperture_failed", "axial": ax, "coronal": co}
    centre = (np.array(ax["p1_ras"]) + np.array(ax["p2_ras"])
              + np.array(co["p1_ras"]) + np.array(co["p2_ras"])) / 4.0
    return {"status": "ok", "centre_ras": centre.tolist(),
            "aperture_axial_mm": ax["aperture_mm"],
            "aperture_coronal_mm": co["aperture_mm"],
            "shift_from_seed_mm": float(np.linalg.norm(centre - np.asarray(seed, float))),
            "axial": ax, "coronal": co}


# ------------------------------------------------------- lumen cross section

def lumen_at_coronal(array, ijk_to_ras, ras_to_ijk, seed_ras, spacing,
                     target_ras=None, midline_x=None, side_sign=None,
                     prev_mask=None, prev_axis_ki=None,
                     axis_search_mm=AXIS_STEP_SEARCH_MM):
    """Mucosal and bony cross section on the coronal plane through `seed_ras`.

    Mucosal: the enclosed airway component containing the seed. Components
    touching the slice border are dropped, which is what stops the region
    flooding out of the nostril and reporting the whole head as lumen.

    Bony: rays cast from the lumen centroid to the first bone. This is a
    star-shaped approximation about the centroid, stated plainly because it
    under-reports a re-entrant recess. It is used rather than a flood fill
    because a flood fill of everything below 300 HU leaks into the orbit and
    the soft tissues of the face.
    """
    ijk = (ras_to_ijk @ np.array([*seed_ras, 1.0]))[:3]
    j = int(round(ijk[1]))
    K, J, I = array.shape
    if not (0 <= j < J):
        return {"status": "plane_outside_volume"}
    sl = array[:, j, :]                       # [k, i] = [S, R]
    dz, dx = spacing[2], spacing[0]

    air = sl < AIR_HU_MAX

    # Restrict to the ipsilateral side of the midline. Without this the two
    # nasal cavities merge, inferiorly and through the choana, into one
    # component: on the reference scan that reported 2000 mm2 and a 25 mm
    # inscribed diameter for a nasal airway, and made the profile oscillate
    # between 310 and 838 mm2 from slice to slice. This is the reason the
    # symmetry plane is computed at all.
    if midline_x is not None and side_sign is not None:
        cols = np.arange(sl.shape[1])
        xs = (ijk_to_ras @ np.stack([cols, np.full_like(cols, j),
                                     np.zeros_like(cols),
                                     np.ones_like(cols)]).astype(float))[0]
        ipsilateral = (((xs - midline_x) * side_sign) > 0) & (
            np.abs(xs - midline_x) <= NASAL_HALF_WIDTH_MM)
        air = air & ipsilateral[None, :]

    lab, n = ndimage.label(air)
    if n == 0:
        return {"status": "no_air_on_plane"}
    border = set(np.unique(np.concatenate(
        [lab[0, :], lab[-1, :], lab[:, 0], lab[:, -1]])))
    border.discard(0)
    seed_ki = (int(round(ijk[2])), int(round(ijk[0])))
    seed_ki = (min(max(seed_ki[0], 0), K - 1), min(max(seed_ki[1], 0), I - 1))
    l = lab[seed_ki]

    # Continuity along the corridor. Choosing each station's component
    # independently lets the selection jump into the maxillary sinus, which is
    # enclosed, ipsilateral, and far larger than the airway, so it looks like a
    # valid pick and silently reports a 25 mm "stent diameter". Preferring the
    # component that overlaps the previous station keeps the profile on one
    # continuous structure.
    if prev_mask is not None:
        overlaps = {}
        for cand in range(1, n + 1):
            if cand in border:
                continue
            ov = int((prev_mask & (lab == cand)).sum())
            if ov > 0:
                overlaps[cand] = ov
        if overlaps:
            l = max(overlaps, key=overlaps.get)

    if l == 0 or l in border:
        best, l = None, 0
        for cand in range(1, n + 1):
            if cand in border:
                continue
            c = np.argwhere(lab == cand)
            if len(c) < 20:
                continue
            d = np.hypot((c[:, 0] - ijk[2]) * dz, (c[:, 1] - ijk[0]) * dx).min()
            if best is None or d < best:
                best, l = d, cand
    if l == 0:
        return {"status": "no_enclosed_airway_on_plane"}

    comp = lab == l
    dt = ndimage.distance_transform_edt(comp, sampling=(dz, dx))
    cen_ki = ndimage.center_of_mass(comp)
    edge = comp & ~ndimage.binary_erosion(comp)
    ei = np.argwhere(edge)

    # The stent dimension is the largest circle that fits at this station,
    # centred on a TRACKED centreline rather than on a straight line between
    # the two seeds. Two failures on the reference scan forced this:
    #   - a global maximum over the component reported 25 mm, because the
    #     maxillary sinus is genuinely continuous with the nasal cavity through
    #     the maxillary ostium, so the component legitimately contains it;
    #   - a straight seed-to-seed axis reported 1.0 mm almost everywhere,
    #     because the nasal corridor curves and the straight line runs along
    #     the lumen edge.
    # Tracking the local distance-transform maximum near the previous station's
    # centre follows the corridor without leaving it.
    if prev_axis_ki is not None:
        r_k = max(1, int(round(axis_search_mm / dz)))
        r_i = max(1, int(round(axis_search_mm / dx)))
        k0, i0 = prev_axis_ki
        win = np.zeros_like(dt, dtype=bool)
        win[max(0, k0 - r_k):k0 + r_k + 1, max(0, i0 - r_i):i0 + r_i + 1] = True
        search = dt * (win & comp)
    else:
        # First station: stay near the operator's seed, so the corridor starts
        # where they pointed rather than in the largest air space in frame.
        r_k = max(1, int(round(12.0 / dz)))
        r_i = max(1, int(round(12.0 / dx)))
        k0, i0 = seed_ki
        win = np.zeros_like(dt, dtype=bool)
        win[max(0, k0 - r_k):k0 + r_k + 1, max(0, i0 - r_i):i0 + r_i + 1] = True
        search = dt * (win & comp)
    if search.max() <= 0:
        return {"status": "no_lumen_near_axis"}
    ax_ki = np.unravel_index(int(np.argmax(search)), search.shape)
    out = {"status": "ok",
           "component_area_mm2": float(comp.sum() * dz * dx),
           "mucosal_inscribed_diameter_mm": float(2 * dt[ax_ki]),
           "component_max_diameter_mm": float(2 * dt.max()),
           "axis_ki": (int(ax_ki[0]), int(ax_ki[1]))}

    # bony contour by ray casting from the lumen centroid
    cen_ras = (ijk_to_ras @ np.array([ax_ki[1], j, ax_ki[0], 1.0]))[:3]
    rays = []
    for ang in np.linspace(0, 2 * np.pi, BONY_N_RAYS, endpoint=False):
        u = np.array([np.cos(ang), 0.0, np.sin(ang)])
        t, hu = _profile(array, ras_to_ijk, cen_ras, u, 40.0)
        hit = np.where(hu >= BONE_HU_MIN)[0]
        rays.append(float(t[hit[0]]) if len(hit) else np.nan)
    rays = np.array(rays)
    if np.isfinite(rays).sum() >= BONY_N_RAYS * 0.8:
        r = np.nan_to_num(rays, nan=np.nanmedian(rays))
        ang = np.linspace(0, 2 * np.pi, BONY_N_RAYS, endpoint=False)
        out["bony_area_mm2"] = float(0.5 * np.sum(r ** 2) * (2 * np.pi / BONY_N_RAYS))
        out["bony_max_inscribed_diameter_mm"] = float(2 * r.min())
        out["bony_contour_method"] = "star-shaped rays from lumen centroid"
    else:
        out["bony_status"] = "too_many_rays_unbounded"

    if target_ras is not None and len(ei):
        d = np.hypot((ei[:, 0] - (ras_to_ijk @ np.array([*target_ras, 1.0]))[2]) * dz,
                     (ei[:, 1] - (ras_to_ijk @ np.array([*target_ras, 1.0]))[0]) * dx)
        out["surface_to_target_mm"] = float(d.min())
    out["centroid_ras"] = cen_ras.tolist()
    out["_mask"] = comp
    return out


def corridor_profile(array, ijk_to_ras, ras_to_ijk, p_ant, p_post, spacing,
                     step_mm=PROFILE_STEP_MM, midline_x=None, side_sign=None):
    """Lumen cross section every step_mm from the piriform seed to the choana seed.

    Sampling the corridor continuously rather than at named anatomical stations
    is deliberate. A stent is a continuous device, so taper and binding point
    are what matter, and it sidesteps the basal lamella ambiguity entirely
    (two CT studies disagree by 4.3 mm on where that landmark is).
    """
    p_ant = np.asarray(p_ant, float); p_post = np.asarray(p_post, float)
    L = float(np.linalg.norm(p_post - p_ant))
    if L <= 0:
        return {"status": "degenerate_corridor"}
    u = (p_post - p_ant) / L
    stations = []
    prev_mask = None
    prev_axis = None
    for t in np.arange(0.0, L + 1e-9, step_mm):
        p = p_ant + u * t
        r = lumen_at_coronal(array, ijk_to_ras, ras_to_ijk, p, spacing,
                             midline_x=midline_x, side_sign=side_sign,
                             prev_mask=prev_mask, prev_axis_ki=prev_axis)
        if r.get("status") == "ok":
            prev_mask = r.pop("_mask", None)
            prev_axis = r.get("axis_ki")
            stations.append({
                "t_mm": round(float(t), 2),
                "component_area_mm2": round(r["component_area_mm2"], 1),
                "mucosal_diam_mm": round(r["mucosal_inscribed_diameter_mm"], 3),
                "bony_diam_mm": (round(r["bony_max_inscribed_diameter_mm"], 3)
                                 if "bony_max_inscribed_diameter_mm" in r else None)})
    if not stations:
        return {"status": "no_valid_stations", "corridor_length_mm": L}

    # Continuity. A real airway changes smoothly along its axis, so a station
    # whose area leaps away from its neighbours means component selection
    # jumped to a different structure. Flag rather than silently average: a
    # jumpy profile is the signature of a merged or wrong lumen.
    areas = np.array([s["component_area_mm2"] for s in stations])
    jumps = 0
    for i in range(1, len(areas)):
        lo, hi = sorted((areas[i - 1], areas[i]))
        if lo > 0 and hi / lo > 2.0:
            stations[i]["discontinuity"] = True
            jumps += 1
    clean = [s for s in stations if not s.get("discontinuity")]
    md = [s["mucosal_diam_mm"] for s in clean] or [s["mucosal_diam_mm"] for s in stations]
    bd = [s["bony_diam_mm"] for s in clean if s["bony_diam_mm"] is not None]
    binding = min(clean or stations, key=lambda s: s["mucosal_diam_mm"])
    return {"status": "ok",
            "corridor_length_mm": L,
            "n_stations": len(stations),
            "n_discontinuities": jumps,
            "profile_is_smooth": bool(jumps <= max(1, len(stations) // 20)),
            "mucosal_min_diameter_mm": float(min(md)),
            "mucosal_max_diameter_mm": float(max(md)),
            "binding_point_t_mm": binding["t_mm"],
            "binding_point_diameter_mm": binding["mucosal_diam_mm"],
            "bony_min_diameter_mm": float(min(bd)) if bd else None,
            "compressible_margin_mm": (round(float(min(bd) - min(md)), 3)
                                       if bd else None),
            "stations": stations}


def derive_corridor_endpoints(array, ijk_to_ras, ras_to_ijk, spf_centre,
                              midline_x, yaw_deg, side_sign):
    """Corridor endpoints from the two structures that sit in the midline.

    Per Dr. Iloreta: the posterior edge of the septum and the anterior edge of
    the maxilla lie in the midline almost always, so both are findable on the
    symmetry plane, which is the one frame this tool derives reliably.

    Searching all midline bone does NOT work: an unbounded argmax put the
    anterior point at the chin and the posterior point behind the sphenoid,
    because the sphenoid sinus also puts air in the midline and the clivus is
    also midline bone. Bounding the search to a height band around the SPF and
    to a limited distance either side of it selects the vomer and the maxilla.

    On the reference scan this gave a corridor of 55.5 mm against a published
    nasal cavity length of about 45 to 60 mm, verified by eye. n=1, unvalidated,
    which is why the operator's seeds win whenever they are present.
    """
    sa, ss = float(spf_centre[1]), float(spf_centre[2])
    step = 0.5
    ty = np.tan(np.radians(yaw_deg))

    def midline_bone(a_mm, s_lo, s_hi):
        """Bone anywhere in a thin midline slab at this A, over an S band.

        Sampled over a +/- 1 mm slab rather than a single line: the vomer is
        thin and a one-voxel-wide probe drops in and out of it, which truncated
        the derived corridor from 55.5 mm to 33.5 mm on the reference scan.
        """
        svals = np.arange(s_lo, s_hi, step)
        for dx in (-1.0, -0.5, 0.0, 0.5, 1.0):
            xs = midline_x + ty * (a_mm - sa) + dx
            pts = np.stack([np.full_like(svals, xs),
                            np.full_like(svals, a_mm), svals], axis=-1)
            if (_sample(array, ras_to_ijk, pts) > BONE_HU_MIN).any():
                return True
        return False

    post = None
    for a in np.arange(sa, sa - CORRIDOR_POST_SEARCH_MM, -step):
        if midline_bone(a, ss - 16.0, ss + 4.0):
            post = a
    ant = None
    for a in np.arange(sa, sa + CORRIDOR_ANT_SEARCH_MM, step):
        if midline_bone(a, ss - 20.0, ss - 4.0):
            ant = a
    if post is None or ant is None:
        return None
    # Place the endpoints inside the airway on the correct side, halfway
    # between the midline and the SPF, and let the centreline tracker do the
    # rest from there.
    x = midline_x + side_sign * abs(float(spf_centre[0]) - midline_x) * 0.5
    return {"p_ant": np.array([x, float(ant), ss]),
            "p_post": np.array([x, float(post), ss]),
            "derived_length_mm": float(ant - post),
            "source": "derived from midline septal and maxillary edges"}


# ------------------------------------------------------------------- driver

def _seeds_in_scene():
    """Resolve seed fiducials by name. Never guesses a missing seed."""
    found = {}
    for node in slicer.util.getNodesByClass("vtkMRMLMarkupsFiducialNode"):
        name = node.GetName()
        for role in SEED_ROLES:
            for side in ("R", "L"):
                if name == f"{role}_{side}":
                    if node.GetNumberOfControlPoints() < 1:
                        continue
                    p = [0.0, 0.0, 0.0]
                    node.GetNthControlPointPositionWorld(0, p)
                    found[(role, side)] = np.array(p)
    return found


def pick_volume(volume_node_name=None):
    if volume_node_name:
        return slicer.mrmlScene.GetFirstNodeByName(volume_node_name)
    vols = [v for v in slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
            if v.GetImageData() is not None]
    if not vols:
        return None
    return max(vols, key=lambda v: np.prod(v.GetImageData().GetDimensions()))


def run(volume_node_name=None, write_json_dir=None, skip_symmetry=False):
    t0 = time.time()
    params = collect_params()
    report = {"script_version": SCRIPT_VERSION,
              "run_datetime": datetime.datetime.now().isoformat(timespec="seconds"),
              "config_hash": config_hash(params),
              "params": params,
              "conduction_budget_mm": CONDUCTION_BUDGET_MM,
              "checks": [], "sides": {}}

    def check(cid, severity, detail):
        report["checks"].append({"id": cid, "severity": severity, "detail": detail})

    volume = pick_volume(volume_node_name)
    if volume is None:
        check("no_volume", "error", "no scalar volume in the scene")
        report["status"] = "failed"
        return report
    report["volume"] = {"name": volume.GetName(),
                        "spacing": list(volume.GetSpacing()),
                        "dimensions": list(volume.GetImageData().GetDimensions())}
    array, M, R2I = _volume_arrays(volume)
    spacing = volume.GetSpacing()

    if not skip_symmetry:
        sym = find_midsagittal(array, M)
        report["midsagittal"] = sym
        if sym is None:
            check("symmetry_failed", "warn", "no bone found for symmetry search")
        elif not sym["sharp_enough"]:
            check("symmetry_flat", "warn",
                  f"the symmetry optimum is flat (Dice drops only "
                  f"{sym['dice_drop_at_1mm']:.4f} at 1 mm off-plane), so the "
                  f"midline estimate is not trustworthy on this scan")

    sym = report.get("midsagittal")
    midline_x = (sym.get("x_ras_mm")
                 if sym and sym.get("sharp_enough") else None)

    seeds = _seeds_in_scene()
    report["seeds_found"] = sorted(f"{r}_{s}" for (r, s) in seeds)
    if not seeds:
        check("no_seeds", "error",
              "no seed fiducials found. Expected names like SPF_R, piriform_R, "
              "choana_R. See the module docstring.")
        report["status"] = "failed"
        return report

    for side in ("R", "L"):
        if not any(s == side for (_, s) in seeds):
            continue
        rec = {}
        spf_seed = seeds.get(("SPF", side))
        if spf_seed is None:
            check(f"no_spf_{side}", "error",
                  f"side {side} has seeds but no SPF_{side}; the headline "
                  f"conduction gap cannot be computed")
        else:
            sc = spf_centre(array, R2I, spf_seed)
            rec["spf"] = sc
            if sc["status"] != "ok":
                check(f"spf_aperture_{side}", "error",
                      f"SPF_{side}: {sc['status']}. Re-place the seed well "
                      f"inside the foramen lumen.")
            else:
                if sc["shift_from_seed_mm"] > 3.0:
                    check(f"spf_centre_moved_{side}", "warn",
                          f"SPF_{side} centre sits {sc['shift_from_seed_mm']:.2f} mm "
                          f"from the seed; check the seed is in the foramen and "
                          f"not an adjacent air space")
                lum = lumen_at_coronal(array, M, R2I, sc["centre_ras"], spacing,
                                       target_ras=sc["centre_ras"],
                                       midline_x=midline_x,
                                       side_sign=(1.0 if side == "R" else -1.0))
                lum.pop("_mask", None)
                rec["spf_station"] = lum
                if lum.get("status") == "ok" and "surface_to_target_mm" in lum:
                    gap = lum["surface_to_target_mm"]
                    rec["conduction_gap_mm"] = gap
                    rec["within_budget"] = bool(gap <= CONDUCTION_BUDGET_MM)
                    if gap > CONDUCTION_BUDGET_MM:
                        check(f"budget_exceeded_{side}", "warn",
                              f"side {side}: mucosal surface to SPF is "
                              f"{gap:.2f} mm, over the "
                              f"{CONDUCTION_BUDGET_MM:.1f} mm budget")
                else:
                    check(f"spf_station_{side}", "warn",
                          f"side {side}: {lum.get('status')}")

        pir, cho = seeds.get(("piriform", side)), seeds.get(("choana", side))
        corridor_source = "seeded"
        if (pir is None or cho is None) and midline_x is not None \
                and rec.get("spf", {}).get("status") == "ok":
            der = derive_corridor_endpoints(
                array, M, R2I, rec["spf"]["centre_ras"], midline_x,
                (sym or {}).get("yaw_deg", 0.0),
                1.0 if side == "R" else -1.0)
            if der is not None:
                pir = pir if pir is not None else der["p_ant"]
                cho = cho if cho is not None else der["p_post"]
                corridor_source = "derived"
                rec["corridor_endpoints"] = der["source"]
                check(f"corridor_derived_{side}", "warn",
                      f"side {side}: corridor endpoints DERIVED from the midline "
                      f"septal and maxillary edges, not seeded. Length "
                      f"{der['derived_length_mm']:.1f} mm. This derivation is "
                      f"validated on one scan by eye only; place piriform_{side} "
                      f"and choana_{side} to override it.")
        if pir is not None and cho is not None:
            if midline_x is None:
                check(f"corridor_no_midline_{side}", "error",
                      f"side {side}: no trustworthy midsagittal plane, so the "
                      f"two nasal cavities cannot be separated. Without the "
                      f"midline the lumen component merges across the septum "
                      f"and reports roughly double the true area. Corridor "
                      f"profile refused rather than reported wrong.")
            else:
                rec["corridor_source"] = corridor_source
                rec["corridor"] = corridor_profile(
                    array, M, R2I, pir, cho, spacing,
                    midline_x=midline_x,
                    side_sign=(1.0 if side == "R" else -1.0))
                cor = rec["corridor"]
                if cor.get("status") != "ok":
                    check(f"corridor_{side}", "warn",
                          f"side {side}: {cor.get('status')}")
                elif not cor.get("profile_is_smooth"):
                    check(f"corridor_jumpy_{side}", "warn",
                          f"side {side}: {cor['n_discontinuities']} of "
                          f"{cor['n_stations']} stations jump by more than 2x "
                          f"in area against their neighbour. The lumen "
                          f"selection is switching structures; inspect the "
                          f"profile before using these dimensions.")
        else:
            missing = [n for n, s in (("piriform", pir), ("choana", cho)) if s is None]
            check(f"corridor_seeds_{side}", "warn",
                  f"side {side}: no corridor profile, missing {missing}")

        for role in ("vidian", "FR"):
            seed = seeds.get((role, side))
            if seed is None:
                rec[f"{role}_clearance_mm"] = None
                continue
            if rec.get("spf", {}).get("status") == "ok":
                d = float(np.linalg.norm(np.array(rec["spf"]["centre_ras"]) - seed))
                rec[f"{role}_clearance_mm"] = round(d, 3)
        report["sides"][side] = rec

    report["status"] = "ok" if not any(
        c["severity"] == "error" for c in report["checks"]) else "failed"
    report["elapsed_s"] = round(time.time() - t0, 1)

    if write_json_dir:
        os.makedirs(write_json_dir, exist_ok=True)
        path = os.path.join(
            write_json_dir,
            f"stent_constraints_{report['run_datetime'].replace(':', '')}.json")
        with open(path, "w") as fh:
            json.dump(report, fh, indent=2, default=str)
        report["written_to"] = path
    return report


def print_report(report):
    B = report.get("conduction_budget_mm", CONDUCTION_BUDGET_MM)
    print("=" * 72)
    print(f" NINS stent design constraints  v{report['script_version']}"
          f"   [{report.get('status')}]")
    print("=" * 72)
    print(f"  volume:  {report.get('volume', {}).get('name')}")
    print(f"  config:  {report.get('config_hash')}   "
          f"{report.get('elapsed_s')} s")
    print(f"  seeds:   {', '.join(report.get('seeds_found', [])) or 'none'}")
    sym = report.get("midsagittal")
    if sym:
        print(f"  midline: x = {sym['x_ras_mm']:.2f} mm, yaw {sym['yaw_deg']:+.2f} deg, "
              f"Dice {sym['dice']:.3f} "
              f"({'sharp' if sym['sharp_enough'] else 'FLAT, do not trust'})")
    for side, rec in report.get("sides", {}).items():
        print()
        print(f"  --- side {side} " + "-" * 52)
        spf = rec.get("spf", {})
        if spf.get("status") == "ok":
            print(f"    SPF aperture      axial {spf['aperture_axial_mm']:.2f} mm"
                  f"   coronal {spf['aperture_coronal_mm']:.2f} mm")
        gap = rec.get("conduction_gap_mm")
        if gap is not None:
            verdict = "WITHIN" if rec.get("within_budget") else "OVER"
            print(f"    CONDUCTION GAP    {gap:.2f} mm   "
                  f"(budget {B:.1f} mm -> {verdict})")
        st = rec.get("spf_station", {})
        if st.get("status") == "ok":
            print(f"    at SPF station    stent dia {st['mucosal_inscribed_diameter_mm']:.2f} mm"
                  f"   (airway component {st['component_area_mm2']:.0f} mm2)")
        cor = rec.get("corridor", {})
        if cor.get("status") == "ok":
            flag = ("" if cor.get("profile_is_smooth")
                    else f"  [JUMPY: {cor['n_discontinuities']} discontinuities]")
            src = rec.get("corridor_source", "seeded")
            print(f"    corridor length   {cor['corridor_length_mm']:.2f} mm"
                  f"   ({cor['n_stations']} stations, {src}){flag}")
            print(f"    stent diameter    min {cor['mucosal_min_diameter_mm']:.2f}"
                  f"  max {cor['mucosal_max_diameter_mm']:.2f} mm")
            print(f"    binding point     {cor['binding_point_diameter_mm']:.2f} mm"
                  f" at {cor['binding_point_t_mm']:.1f} mm from the piriform seed")
            if cor.get("compressible_margin_mm") is not None:
                print(f"    compressible      {cor['compressible_margin_mm']:.2f} mm"
                      f" (bony minus mucosal)")
        for role in ("vidian", "FR"):
            d = rec.get(f"{role}_clearance_mm")
            if d is not None:
                print(f"    {role:<8} clearance {d:.2f} mm from the SPF centre")
    if report.get("checks"):
        print()
        print("  checks:")
        for c in report["checks"]:
            print(f"    [{c['severity']}] {c['id']}: {c['detail']}")
    print("=" * 72)


print(f"(C) NINS_Stent_Constraints_0811.py v{SCRIPT_VERSION} loaded. "
      f"Nothing has run yet.")
print(f"  Seeds expected: " + ", ".join(f"{r}_R/{r}_L" for r in SEED_ROLES))
print(f"  Run:  report = run(); print_report(report)")
