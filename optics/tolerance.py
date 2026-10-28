"""Step 3: tolerancing the secondary mirror. For despace, decentre and tilt: how large an error costs 0.02 of
system MTF at Nyquist (worst over the usable field, after refocusing the detector), then a Monte Carlo of
builds with all three errors at once."""
import json
import numpy as np
from scipy.optimize import brentq
import requirements as R
from telescope import design, wavefront, psf_mtf, mtf_at, pixel_mtf, NYQ
from raytrace import trace

T = design("rc")
perf = json.load(open("results_performance.json"))
Z0, FIELDS = T["b"] + perf["detector_offset_um"] * 1e-6, [0.0, 0.1, perf["max_field_deg"]]
PIX = pixel_mtf(NYQ)
FOCUS_GAIN = T["m"] ** 2 + 1                  # focus moves (m^2 + 1) x the secondary despace


def worst_mtf(dz2=0.0, dx2=0.0, tilt2=0.0, z=None):
    z = Z0 + FOCUS_GAIN * dz2 if z is None else z
    vals = []
    for fd in FIELDS:
        nu, mx, my, _ = psf_mtf(wavefront(T, fd, z, n=64, dz2=dz2, dx2=dx2, tilt2=tilt2))
        vals.append(min(mtf_at(nu, mx, NYQ), mtf_at(nu, my, NYQ)) * PIX)
    return min(vals)


def paraxial_focus(**pert):
    X, D, _ = trace(T, 0.0, np.array([1e-4]), np.array([0.0]), **pert)
    return X[0, 2] - X[0, 0] * D[0, 2] / D[0, 0]


def refocused(**pert):
    """Detector refocused: start from the focus shift the despace causes (keeping the nominal offset),
    try a few positions around it, keep the best."""
    zc = Z0 + paraxial_focus(dz2=pert.get("dz2", 0.0)) - paraxial_focus()   # decentre/tilt barely move focus
    return max(worst_mtf(z=zc + dz, **pert) for dz in np.linspace(-60e-6, 60e-6, 7))


dz = 10e-6
print(f"focus moves {(paraxial_focus(dz2=dz) - paraxial_focus())/dz:.2f} x the secondary despace (m^2 + 1 = {FOCUS_GAIN:.2f})")
nominal = refocused()
print(f"nominal worst-case MTF at Nyquist over 0-{FIELDS[-1]:.2f} deg: {nominal:.3f} (requirement {R.MTF_NYQUIST_MIN})")
alloc = {}
for name, key, scale, unit in [("despace", "dz2", 1e-6, "um"), ("decentre", "dx2", 1e-6, "um"),
                               ("tilt", "tilt2", np.deg2rad(1 / 3600), "arcsec")]:
    f = lambda v: refocused(**{key: v * scale}) - (nominal - 0.02)
    hi = 1.0
    while f(hi) > 0 and hi < 1e4:
        hi *= 2
    alloc[key] = brentq(f, 0, hi, xtol=0.05) * scale
    print(f"  secondary {name:9}: {alloc[key]/scale:7.1f} {unit} costs 0.02 of MTF")

rng = np.random.default_rng(3)
builds = []
for _ in range(150):
    pert = {k: rng.uniform(-v, v) for k, v in alloc.items()}     # each error anywhere within its tolerance
    builds.append(refocused(**pert))
builds = np.array(builds)
yield_ = float(np.mean(builds >= R.MTF_NYQUIST_MIN))
print(f"Monte Carlo, 150 builds, all three errors uniform within those tolerances, detector refocused:")
print(f"  worst-case MTF median {np.median(builds):.3f}, 10th percentile {np.percentile(builds, 10):.3f}; "
      f"meet the requirement: {100*yield_:.0f} %")
json.dump(dict(nominal=nominal, tolerances=dict(despace_um=alloc["dz2"] * 1e6, decentre_um=alloc["dx2"] * 1e6,
          tilt_arcsec=float(np.degrees(alloc["tilt2"]) * 3600)), focus_gain=FOCUS_GAIN,
          mc_median=float(np.median(builds)), mc_p10=float(np.percentile(builds, 10)), mc_yield=yield_),
          open("results_tolerance.json", "w"), indent=1)
