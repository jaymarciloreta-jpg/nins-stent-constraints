# Reading the Mask: teaching course

`reading-the-mask.html` is a self-contained interactive course for attending physicians and surgical residents on how CT and MRI autosegmentation works, tied step by step to the [NINS Pilot panel](../PILOT_PANEL.md). Open the file in any browser. It needs a network connection only for its web fonts.

The published copy is at <https://claude.ai/artifact/7t6sVfnKiTJPZVwbKcjam4>.

- Eight lessons, eleven live demos, flashcards, knowledge checks, a seven-case "spot the problem" quiz and a final review.
- Every image is generated in the browser from a synthetic phantom. The only embedded picture is the panel's synthetic demonstration screenshot from `docs/images`. The course contains no patient data, and none may be added.
- Lesson 8 includes a simplified practice copy of the panel. It quotes the panel's step titles, button labels, default thresholds and refusal messages.

## Keeping it in step with the panel

`tests/test_teaching_course.py` compares those quoted strings and defaults with `NINSPilot/NINSPilot.py` and `NINSPilot/pilot_logic.py`. If you rename a button, change a default or reword a refusal message in the panel, that test fails and names the string. Update the course text to match, then republish the page.

```bash
python -m unittest tests.test_teaching_course -v
```

The course is educational material. It is not validated clinical guidance, and its phantoms are not calibrated to any scanner.
