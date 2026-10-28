"""Tolerance stack for the lateral position of the secondary mirror, compared with the decentre the optics
can accept (optics/tolerance.py). Part tolerances are assumed machining values (TBC with the workshop)."""
import json
import pathlib
import numpy as np

CHAIN = [("bore in the primary bulkhead", 20), ("tube concentricity", 50), ("tube-to-bulkhead fit", 25),
         ("spider hub position", 25), ("mirror cell bore", 20), ("mirror in its cell", 25)]   # +/- um
budget = json.load(open(pathlib.Path(__file__).resolve().parents[1] / "optics" / "results_tolerance.json"))["tolerances"]["decentre_um"]
sag = {r["material"]: r["sag_um"] for r in json.load(open("results_structure.json"))["rows"]}
wc, rss = sum(t for _, t in CHAIN), np.sqrt(sum(t ** 2 for _, t in CHAIN))
for name, t in CHAIN:
    print(f"  {name:30} +/-{t:3d} um")
print(f"worst case +/-{wc:.0f} um, root-sum-square +/-{rss:.0f} um; optics allows {budget:.0f} um of decentre")
print(f"ground-to-orbit change (gravity sag released), CFRP tube: {sag['CFRP quasi-isotropic']:.1f} um")
json.dump(dict(chain=CHAIN, worst_case_um=wc, rss_um=rss, optics_allowance_um=budget), open("results_stack.json", "w"), indent=1)
