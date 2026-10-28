"""How far the image moves out of focus when the telescope warms or cools uniformly, for each candidate
metering-structure material, with glass-ceramic mirrors (CTE ~0, so the mirrors keep their shape). The
spacing change is traced through the optical model, and the allowed focus error comes from the optics too."""
import json
import pathlib
import sys
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "optics"))
import requirements as R                                                   # noqa: E402
from telescope import design, wavefront, psf_mtf, mtf_at, pixel_mtf, NYQ  # noqa: E402
from raytrace import trace                                                 # noqa: E402
from materials import MATERIALS                                            # noqa: E402

T = design("rc")
N = R.FOCAL_LENGTH_M / R.APERTURE_M
perf = json.load(open(pathlib.Path(__file__).resolve().parents[1] / "optics" / "results_performance.json"))
Z0 = T["b"] + perf["detector_offset_um"] * 1e-6


def focus(dz2):
    X, D, _ = trace(T, 0.0, np.array([1e-4]), np.array([0.0]), dz2=dz2)
    return X[0, 2] - X[0, 0] * D[0, 2] / D[0, 0]


def worst_mtf(defocus):
    v = []
    for fd in (0.0, 0.1, perf["max_field_deg"]):
        nu, mx, my, _ = psf_mtf(wavefront(T, fd, Z0 + defocus, n=64))
        v.append(min(mtf_at(nu, mx, NYQ), mtf_at(nu, my, NYQ)) * pixel_mtf(NYQ))
    return min(v)


nominal = worst_mtf(0.0)
allowed = min(brentq(lambda dz: worst_mtf(s * dz) - (nominal - 0.02), 0, 500e-6) for s in (1, -1))
print(f"depth of focus, classical +/-2 lambda N^2: {2 * R.WAVELENGTH_M * N**2 * 1e6:.0f} um; "
      f"defocus costing 0.02 MTF in this design: {allowed*1e6:.0f} um")
rows = []
print(f"{'metering structure':22} {'CTE ppm/K':>9} {'focus shift um/K':>17} {'allowed uniform dT':>19}")
for name, m in MATERIALS.items():
    # warming lengthens the mirror spacing (secondary moves away from the primary: dz2 < 0) and the
    # detector support behind the primary (detector moves away by cte * b)
    dz2 = -m["cte"] * T["d"]
    shift = (focus(dz2) - focus(0.0)) - m["cte"] * Z0                  # image motion relative to the detector, per K
    dT = allowed / abs(shift)
    rows.append(dict(material=name, cte_ppm=m["cte"] * 1e6, shift_um_per_K=shift * 1e6, allowed_dT_K=dT))
    print(f"{name:22} {m['cte']*1e6:9.1f} {shift*1e6:17.2f} {dT:16.1f} K")
json.dump(dict(allowed_defocus_um=allowed * 1e6, dof_classical_um=2 * R.WAVELENGTH_M * N ** 2 * 1e6, rows=rows),
          open("results_thermal.json", "w"), indent=1)
