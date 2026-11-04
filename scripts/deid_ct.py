#!/usr/bin/env python3
"""De-identify a CT DICOM export and independently verify the result.

Usage:
    python deid_ct.py inventory SRC
    python deid_ct.py add SRC [--library DIR] [--key DIR] [--modality CT MR] [--series N ...]

Library layout (no original identifiers anywhere, including folder names):
    DEID_DICOM/index.csv
    DEID_DICOM/P001/manifest.json
    DEID_DICOM/P001/<shifted-date>_<study>/S0003_CT_<desc>_<id>/00001_<id>.dcm

Pipeline per kept file:
  1. dicognito: replaces names/IDs, shifts dates consistently, remaps every UID
     (including nested ones) consistently, so volumes still load.
  2. Strip private tags, overlays/curves, and free-text comment fields.
  3. Remove PatientBirthDate (age is kept).

Series are kept only if they are CT or MR image storage objects (per --modality),
ImageType without LOCALIZER/SCREEN SAVE/SECONDARY/DOSE, and BurnedInAnnotation is
not YES. Scouts, dose screens, SR/PR/secondary captures are excluded.

Verification walks every element of every output file (and file meta) and checks
that no original name part, ID, institution string, date or non-standard UID
survives. Nothing here prints header values: only counts, tag keywords, indexes.

Key: --key holds a random salt and a table {salted hash of PatientID+BirthDate -> P###}.
The anonymizer is seeded per patient from it, so re-running or adding a new study for
an existing patient gives the same P-number, fake UIDs and date shift. It is a
re-identification key: keep it off the library drive and never share it.
"""

from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import re
import secrets
import sys
from pathlib import Path

import pydicom
from pydicom.dataset import Dataset
from dicognito.anonymizer import Anonymizer

IMAGE_SOP = {
    "CT": {"1.2.840.10008.5.1.4.1.1.2", "1.2.840.10008.5.1.4.1.1.2.1"},
    "MR": {"1.2.840.10008.5.1.4.1.1.4", "1.2.840.10008.5.1.4.1.1.4.1"},
}
MODALITIES = {"CT"}  # set from --modality
BAD_IMAGE_TYPE = {"LOCALIZER", "SCREEN SAVE", "SECONDARY", "DOSE_INFO", "DOSE INFO", "DERIVED_SCOUT"}

TEXT_KEYWORDS_TO_DROP = {
    "ImageComments", "StudyComments", "PatientComments", "AdditionalPatientHistory",
    "RequestAttributesSequence", "RequestedProcedureDescription", "ReasonForStudy",
    "ReasonForTheRequestedProcedure", "RequestedProcedureComments", "InterpretationText",
    "AcquisitionComments", "ImagingServiceRequestComments", "PerformedProcedureStepDescription",
    "ScheduledProcedureStepDescription", "AdmittingDiagnosesDescription", "Allergies",
    "MedicalAlerts", "PatientBirthDate", "PatientBirthTime", "PatientBirthName",
    "PatientMotherBirthName", "OtherPatientNames", "OtherPatientIDsSequence",
    "IssuerOfPatientID", "PatientAddress", "CountryOfResidence", "RegionOfResidence",
    "DeviceSerialNumber", "OperatorsName", "PerformingPhysicianName",
}

# Values we harvest from originals as "canaries" to look for in outputs.
CANARY_KEYWORDS = {
    "PatientID", "OtherPatientIDs", "AccessionNumber", "StudyID", "InstitutionName",
    "InstitutionAddress", "InstitutionalDepartmentName", "StationName", "DeviceSerialNumber",
    "PatientAddress", "PatientTelephoneNumbers", "IssuerOfPatientID", "MedicalRecordLocator",
    "RequestedProcedureID", "PerformedProcedureStepID", "ScheduledProcedureStepID",
    "PlacerOrderNumberImagingServiceRequest", "FillerOrderNumberImagingServiceRequest",
}
STD_UID_PREFIX = "1.2.840.10008."


# ---------------------------------------------------------------- helpers
SKIP_EXT = {".exe", ".dll", ".zip", ".xml", ".txt", ".def", ".chk", ".msg", ".rtf", ".inf",
            ".jpg", ".png", ".pdf", ".mp4", ".mov", ".html", ".htm", ".js", ".css", ".ini", ".json"}


def dicom_files(root: Path):
    for p in sorted(root.rglob("*")):
        if (p.is_file() and not p.name.startswith("._") and p.name != "DICOMDIR"
                and p.suffix.lower() not in SKIP_EXT):
            yield p


def read_header(p: Path):
    try:
        return pydicom.dcmread(p, stop_before_pixels=True, force=True)
    except Exception:
        return None


def series_key(ds) -> str:
    return str(ds.get("SeriesInstanceUID", "NO_SERIES"))


def classify(ds) -> tuple[bool, str]:
    sop = str(ds.get("SOPClassUID", ""))
    mod = str(ds.get("Modality", ""))
    if mod not in MODALITIES:
        return False, f"modality={mod or '?'}"
    if sop not in IMAGE_SOP.get(mod, set()):
        return False, f"not {mod} image storage"
    it = {str(v).upper() for v in (ds.get("ImageType") or [])}
    bad = it & BAD_IMAGE_TYPE
    if bad:
        return False, f"ImageType {sorted(bad)}"
    if str(ds.get("BurnedInAnnotation", "")).upper() == "YES":
        return False, "BurnedInAnnotation=YES"
    return True, f"{mod} image"


def walk_strings(ds: Dataset):
    """Yield (keyword, VR, str_value) for every element, recursing sequences."""
    for elem in ds:
        if elem.VR == "SQ":
            for item in elem.value or []:
                yield from walk_strings(item)
            continue
        if elem.tag.group == 0x7FE0:
            continue
        val = elem.value
        if val is None or isinstance(val, (bytes, bytearray)):
            continue
        vals = val if isinstance(val, (list, tuple, pydicom.multival.MultiValue)) else [val]
        for v in vals:
            yield elem.keyword or str(elem.tag), elem.VR, str(v)


def harvest_canaries(ds: Dataset, c: dict[str, set[str]]) -> None:
    for src in (ds, getattr(ds, "file_meta", Dataset())):
        for kw, vr, v in walk_strings(src):
            v = v.strip()
            if not v:
                continue
            if vr == "PN":
                for part in re.split(r"[\^=\s,]+", v):
                    if len(part) >= 3:
                        c["name"].add(part.upper())
            elif kw in CANARY_KEYWORDS and len(v) >= 3:
                c["id"].add(v.upper())
            elif vr == "DA" and re.fullmatch(r"\d{8}", v):
                c["date"].add(v)
            elif vr == "UI" and not v.startswith(STD_UID_PREFIX):
                c["uid"].add(v)


def preserve_diffusion(ds: Dataset) -> None:
    """Copy Siemens private diffusion values (numbers, not PHI) into standard tags
    so stripping private tags doesn't lose DWI b-values/gradient directions."""
    if "DiffusionBValue" not in ds and (0x0019, 0x100C) in ds:
        try:
            ds.DiffusionBValue = float(ds[0x0019, 0x100C].value)
        except (TypeError, ValueError):
            pass
    if "DiffusionGradientOrientation" not in ds and (0x0019, 0x100E) in ds:
        v = ds[0x0019, 0x100E].value
        try:
            vals = [float(x) for x in v] if not isinstance(v, bytes) else list(
                __import__("struct").unpack("<3d", v[:24]))
            if len(vals) == 3:
                ds.DiffusionGradientOrientation = vals
        except Exception:
            pass


def strip_extra(ds: Dataset) -> None:
    ds.remove_private_tags()
    for elem in list(ds):
        g = elem.tag.group
        if 0x5000 <= g <= 0x50FF or 0x6000 <= g <= 0x60FF:  # curves, overlays
            del ds[elem.tag]
        elif elem.keyword in TEXT_KEYWORDS_TO_DROP:
            del ds[elem.tag]
        elif elem.VR == "SQ":
            for item in elem.value or []:
                strip_extra(item)


def leaks(ds: Dataset, c: dict[str, set[str]]) -> list[str]:
    """Return 'kind:keyword' for each leaked canary (never the value)."""
    found = []
    for src in (ds, ds.file_meta):
        for kw, vr, v in walk_strings(src):
            u = v.upper()
            if v in c["uid"]:
                found.append(f"uid:{kw}")
            if vr in ("DA", "DT") and v[:8] in c["date"]:
                found.append(f"date:{kw}")
            if vr in ("PN", "LO", "SH", "ST", "LT", "UT", "UC", "CS", "AE"):
                if any(n in u for n in c["name"]):
                    found.append(f"name:{kw}")
                if any(i in u for i in c["id"]):
                    found.append(f"id:{kw}")
    return found


# ---------------------------------------------------------------- commands
def inventory(src: Path) -> dict:
    series = collections.OrderedDict()
    unreadable = 0
    for p in dicom_files(src):
        ds = read_header(p)
        if ds is None or "SOPClassUID" not in ds:
            unreadable += 1
            continue
        k = series_key(ds)
        if k not in series:
            keep, why = classify(ds)
            series[k] = {"n": 0, "number": ds.get("SeriesNumber"), "keep": keep, "why": why,
                         "mod": str(ds.get("Modality", "NA")),
                         "rows": ds.get("Rows"), "thick": ds.get("SliceThickness"), "files": []}
        series[k]["n"] += 1
        series[k]["files"].append(p)
    return {"series": series, "unreadable": unreadable}


def print_inventory(inv: dict) -> None:
    print(f"{'#':>3} {'SeriesNo':>8} {'files':>6} {'rows':>5} {'thick':>6}  keep  reason")
    for i, s in enumerate(inv["series"].values()):
        print(f"{i:>3} {str(s['number']):>8} {s['n']:>6} {str(s['rows']):>5} {str(s['thick']):>6}  "
              f"{'YES ' if s['keep'] else 'no  '}  {s['why']}")
    print(f"non-DICOM/unreadable files: {inv['unreadable']}")




# ---------------------------------------------------------------- patient library
def norm(v) -> str:
    return re.sub(r"\s+", " ", str(v or "")).strip().upper()


def load_key(keydir: Path) -> tuple[str, dict]:
    keydir.mkdir(parents=True, exist_ok=True)
    keydir.chmod(0o700)
    salt_f, reg_f = keydir / "salt.txt", keydir / "registry.json"
    if not salt_f.exists():
        salt_f.write_text(secrets.token_hex(32))
        salt_f.chmod(0o600)
    reg = json.loads(reg_f.read_text()) if reg_f.exists() else {}
    return salt_f.read_text().strip(), reg


def save_registry(keydir: Path, reg: dict) -> None:
    f = keydir / "registry.json"
    f.write_text(json.dumps(reg, indent=1, sort_keys=True))
    f.chmod(0o600)


def patient_hash(ds, salt: str) -> str:
    """Salted hash of PatientID + BirthDate; the raw values are never stored."""
    raw = f"{norm(ds.get('PatientID'))}|{norm(ds.get('PatientBirthDate'))}"
    return hashlib.sha256((salt + raw).encode()).hexdigest()


def safe_label(text: str, c: dict[str, set[str]], maxlen: int = 40) -> str:
    """Folder-name fragment from a de-identified description, dropped if it hits a canary."""
    u = norm(text)
    if not u or any(n in u for n in c["name"]) or any(i in u for i in c["id"]):
        return ""
    return re.sub(r"[^A-Za-z0-9]+", "_", u).strip("_")[:maxlen]


def first_date(ds) -> str:
    for k in ("StudyDate", "SeriesDate", "AcquisitionDate", "ContentDate"):
        v = str(ds.get(k, "") or "")
        if re.fullmatch(r"\d{8}", v):
            return v
    return "nodate"


def short(uid: str, n: int = 6) -> str:
    return hashlib.sha1(uid.encode()).hexdigest()[:n]


def add(src: Path, lib: Path, keydir: Path, only: list[int] | None) -> int:
    inv = inventory(src)
    print_inventory(inv)
    chosen = [s for i, s in enumerate(inv["series"].values())
              if s["keep"] and (only is None or i in only)]
    if not chosen:
        print("nothing selected")
        return 1

    salt, reg = load_key(keydir)
    canaries: dict[str, set[str]] = {"name": set(), "id": set(), "date": set(), "uid": set()}
    anonymizers: dict[str, Anonymizer] = {}
    written, new_patients = [], []
    lib.mkdir(parents=True, exist_ok=True)

    for s in chosen:
        for p in s["files"]:
            ds = pydicom.dcmread(p, force=True)
            harvest_canaries(ds, canaries)
            h = patient_hash(ds, salt)
            if h not in reg:
                reg[h] = {"pid": f"P{len(reg) + 1:03d}", "added": dt.date.today().isoformat()}
                new_patients.append(reg[h]["pid"])
            pid = reg[h]["pid"]
            # Seeded per patient: same patient => same fake UIDs/date shift on every run.
            anon = anonymizers.setdefault(h, Anonymizer(seed=salt + h))
            preserve_diffusion(ds)
            anon.anonymize(ds)
            strip_extra(ds)
            ds.PatientID = pid
            ds.PatientName = pid
            ds.PatientIdentityRemoved = "YES"
            fm = ds.file_meta
            fm.ImplementationClassUID = pydicom.uid.PYDICOM_IMPLEMENTATION_UID
            fm.ImplementationVersionName = "DEID_CT"
            if "SourceApplicationEntityTitle" in fm:
                del fm.SourceApplicationEntityTitle

            study_dir = lib / pid / f"{first_date(ds)}_{ds.Modality}_{short(str(ds.StudyInstanceUID))}"
            desc = safe_label(ds.get("SeriesDescription", ""), canaries)
            sn = ds.get("SeriesNumber")
            series_dir = study_dir / "_".join(x for x in (
                f"S{int(sn):04d}" if sn is not None else "S----", str(ds.Modality), desc,
                short(str(ds.SeriesInstanceUID), 4)) if x)
            series_dir.mkdir(parents=True, exist_ok=True)
            inst = ds.get("InstanceNumber")
            name = f"{int(inst):05d}_" if inst is not None else ""
            out = series_dir / f"{name}{short(str(ds.SOPInstanceUID), 10)}.dcm"
            ds.save_as(out, enforce_file_format=True)
            written.append(out)

    save_registry(keydir, reg)

    # Independent verification of exactly the files written in this run.
    failures, leak_files, sops = collections.Counter(), [], set()
    series_uids = collections.defaultdict(set)
    for f in written:
        ds = pydicom.dcmread(f, stop_before_pixels=True)
        sops.add(str(ds.SOPInstanceUID))
        series_uids[f.parent].add(str(ds.SeriesInstanceUID))
        if any(e.tag.is_private for e in ds.iterall()):
            failures["private tag present"] += 1
        if str(ds.get("PatientIdentityRemoved")) != "YES":
            failures["PatientIdentityRemoved not YES"] += 1
        hit = leaks(ds, canaries)
        if hit:
            leak_files.append((str(f.relative_to(lib)), sorted(set(hit))))
    # Folder names are checked too.
    for d in {f.parent for f in written} | {f.parent.parent for f in written}:
        u = norm(d.name)
        if any(n in u for n in canaries["name"]) or any(i in u for i in canaries["id"]):
            failures["identifier in folder name"] += 1

    report = {
        "files_written": len(written),
        "unique_sop_uids": len(sops) == len(written),
        "one_series_uid_per_folder": all(len(v) == 1 for v in series_uids.values()),
        "patients_in_this_run": sorted({f.relative_to(lib).parts[0] for f in written}),
        "new_patients": new_patients,
        "canaries": {k: len(v) for k, v in canaries.items()},
        "other_failures": dict(failures),
        "files_with_leaks": len(leak_files),
        "leaks_by_tag": dict(collections.Counter(h for _, hits in leak_files for h in hits)),
        "leak_example_files": [f for f, _ in leak_files[:5]],
    }
    print(json.dumps(report, indent=2))
    build_index(lib)
    # exFAT/FAT drives get macOS "._" AppleDouble sidecars for every file we write.
    for junk in lib.rglob("._*"):
        junk.unlink(missing_ok=True)
    ok = not leak_files and not failures and report["unique_sop_uids"] and report["one_series_uid_per_folder"]
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 2


def build_index(lib: Path) -> None:
    """Write LIB/index.csv and LIB/P###/manifest.json from de-identified headers only."""
    rows = []
    for pdir in sorted(d for d in lib.iterdir() if d.is_dir() and re.fullmatch(r"P\d{3,}", d.name)):
        studies = []
        for sdir in sorted(d for d in pdir.iterdir() if d.is_dir()):
            series = []
            for se in sorted(d for d in sdir.iterdir() if d.is_dir()):
                files = sorted(f for f in se.glob("*.dcm") if not f.name.startswith("._"))
                if not files:
                    continue
                h = pydicom.dcmread(files[0], stop_before_pixels=True)
                series.append({"folder": se.name, "modality": str(h.get("Modality", "")),
                               "files": len(files), "rows": h.get("Rows"), "cols": h.get("Columns"),
                               "slice_thickness": str(h.get("SliceThickness", "")),
                               "pixel_spacing": [str(x) for x in (h.get("PixelSpacing") or [])]})
            studies.append({"folder": sdir.name, "series": series})
        (pdir / "manifest.json").write_text(json.dumps({"patient": pdir.name, "studies": studies}, indent=1))
        mods = sorted({se["modality"] for st in studies for se in st["series"]})
        rows.append([pdir.name, len(studies), sum(len(st["series"]) for st in studies),
                     "+".join(mods), sum(se["files"] for st in studies for se in st["series"])])
    with open(lib / "index.csv", "w") as f:
        f.write("patient,studies,series,modalities,files\n")
        for r in rows:
            f.write(",".join(map(str, r)) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("inventory"); a.add_argument("src", type=Path)
    b = sub.add_parser("add", help="de-identify SRC into the patient library")
    b.add_argument("src", type=Path)
    b.add_argument("--library", type=Path, default=Path("/Volumes/Extreme SSD/DEID_DICOM"))
    b.add_argument("--key", type=Path, default=Path.home() / "PHI/.key",
                   help="private key dir (salt + patient registry); keep it OFF the library drive")
    b.add_argument("--series", type=int, nargs="*", help="inventory row numbers to include")
    for p in (a, b):
        p.add_argument("--modality", nargs="+", default=["CT", "MR"], choices=sorted(IMAGE_SOP))
    args = ap.parse_args()
    MODALITIES.clear(); MODALITIES.update(args.modality)
    if args.cmd == "inventory":
        print_inventory(inventory(args.src))
        return 0
    return add(args.src, args.library, args.key, args.series)


if __name__ == "__main__":
    sys.exit(main())
