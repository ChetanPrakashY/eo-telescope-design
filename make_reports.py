"""Fill the report templates with numbers from the result files, so the reports always match the analysis.
    python make_reports.py      # optics/REPORT.template.md -> optics/REPORT.md, same for mechanical/
"""
import json
import pathlib
import re

HERE = pathlib.Path(__file__).resolve().parent
ctx = {p.stem.replace("results_", ""): json.loads(p.read_text()) for p in HERE.glob("*/results_*.json")}


def fill(text):
    def sub(m):
        v = ctx
        for part in m.group(1).split("."):
            v = v[int(part)] if isinstance(v, list) else v[part]
        return format(v, m.group(2) or "")
    return re.sub(r"\{\{([\w.\-]+?)(?::([^}]+))?\}\}", sub, text)


for t in HERE.glob("*/REPORT.template.md"):
    (t.parent / "REPORT.md").write_text(fill(t.read_text()))
    print("wrote", t.parent.name + "/REPORT.md")
