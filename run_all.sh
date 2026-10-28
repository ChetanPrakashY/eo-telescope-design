#!/usr/bin/env bash
# Re-run the whole study and rebuild both reports (about a minute).
set -euo pipefail
cd "$(dirname "$0")/optics"
python3 design.py && python3 performance.py && python3 tolerance.py && python3 radiometry.py
cd ../mechanical
python3 thermal_focus.py && python3 structure.py && python3 stack.py
cd .. && python3 make_reports.py
