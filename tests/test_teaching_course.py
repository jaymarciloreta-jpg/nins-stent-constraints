"""The teaching course quotes the NINS Pilot panel. Fail when the two drift apart."""
import re
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
COURSE = (ROOT/"docs/teaching/reading-the-mask.html").read_text(encoding="utf-8")
PANEL = (ROOT/"NINSPilot/NINSPilot.py").read_text(encoding="utf-8")
LOGIC = (ROOT/"NINSPilot/pilot_logic.py").read_text(encoding="utf-8")

STEP_TITLES = [
    "1. Source and session", "2. Independent reference annotations",
    "3. Candidate segmentation and correction", "4. Review geometry and calculate",
    "5. Save and export"]
BUTTONS = [
    "Start session", "Load synthetic demonstration", "Start / resume timer", "Pause timer",
    "Build right contact patch", "Build left contact patch", "Accept current geometry",
    "Calculate measurements", "Save progress / refusal record",
    "Export measurements and slice evidence",
    "I reviewed the masks, rims and any contact patches in the source CT",
    "HU threshold proposals", "TotalSegmentator nasal initialization"]
PANEL_MESSAGES = [
    "Start a session or load the synthetic demonstration.",
    "Synthetic fixture ready. Review the masks and patches, then accept and calculate. This is not a pilot patient.",
    "Start a session first", "Confirm anatomical review before accepting",
    "Current geometry accepted. Later changes require new acceptance.",
    "Research geometry; not a fabrication specification."]
LOGIC_MESSAGES = [
    "Start or resume a session",
    "Generate and review a segmentation on at least one side",
    "Air mask, contact point or patch settings changed; rebuild patch ",
    "Review and accept the current inputs before measuring",
    "Measurements are missing or inputs changed; review and recalculate"]
DEFAULTS = [("air_max_hu", "Air upper threshold (HU)"), ("bone_min_hu", "Bone lower threshold (HU)"),
            ("geometric_tolerance_mm", "Exploratory gap tolerance (mm)")]


def shown(value):
    """The way the course prints a panel default: true minus sign, no trailing zeros."""
    text = f"{value:.2f}".rstrip("0").rstrip(".") if value % 1 else str(int(value))
    return text.replace("-", "−")


class TeachingCourseMatchesPanel(unittest.TestCase):
    def assertQuoted(self, text, source, where):
        self.assertTrue(text in source, f"{where} no longer contains {text!r}; update the course and this list")
        self.assertTrue(text in COURSE, f"course no longer quotes {text!r}")

    def test_step_titles_and_controls(self):
        for text in STEP_TITLES+BUTTONS:
            with self.subTest(text=text):
                self.assertQuoted(text, PANEL, "NINSPilot.py")

    def test_status_and_refusal_messages(self):
        for text in PANEL_MESSAGES:
            with self.subTest(text=text):
                self.assertQuoted(text, PANEL, "NINSPilot.py")
        for text in LOGIC_MESSAGES:
            with self.subTest(text=text):
                self.assertQuoted(text, LOGIC, "pilot_logic.py")

    def test_default_settings(self):
        for key, label in DEFAULTS:
            with self.subTest(key=key):
                found = re.search(r'\("%s", _\("%s"\), (-?[\d.]+)' % (re.escape(key), re.escape(label)), PANEL)
                self.assertIsNotNone(found, f"{key} default not found in NINSPilot.py")
                row = f'<span>{label}</span><span class="sl-field">{shown(float(found.group(1)))}</span>'
                self.assertTrue(row in COURSE, f"course does not show {shown(float(found.group(1)))} for {label}")

    def test_segment_editor_effect_order(self):
        found = re.search(r"setEffectNameOrder\(\[(.*?)\]\)", PANEL)
        self.assertIsNotNone(found)
        effects = re.findall(r'"([^"]+)"', found.group(1))
        strip = COURSE[COURSE.index('id="mk-effects"'):]
        strip = strip[:strip.index("</div>")]
        self.assertEqual(re.findall(r'<button type="button">([^<]+)</button>', strip), effects)

    def test_course_embeds_no_external_images(self):
        sources = re.findall(r'<img[^>]+src="([^"]{0,40})', COURSE)
        self.assertTrue(sources)
        for src in sources:
            self.assertTrue(src.startswith("data:image/"), src)


if __name__ == "__main__":
    unittest.main()
