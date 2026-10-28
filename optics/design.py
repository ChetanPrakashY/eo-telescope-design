"""Step 1: first-order design from the requirements, the Ritchey-Chretien prescription, and two checks of
the ray tracer: the paraxial focus and focal length, and the coma of a classical Cassegrain against the
third-order formula."""
import json
import numpy as np
import requirements as R
from telescope import design, spot, F1, M, NYQ
from raytrace import trace, to_plane

T = design("rc")
N = R.FOCAL_LENGTH_M / R.APERTURE_M
airy = 1.22 * R.WAVELENGTH_M * N
first_order = dict(focal_length_m=R.FOCAL_LENGTH_M, f_number=N, primary_focal_length_m=F1, magnification=M,
                   plate_scale_urad_per_px=R.PIXEL_M / R.FOCAL_LENGTH_M * 1e6, nyquist_cyc_mm=NYQ,
                   optical_cutoff_cyc_mm=1 / (R.WAVELENGTH_M * N * 1e3), airy_radius_um=airy * 1e6,
                   q_factor=R.WAVELENGTH_M * N / R.PIXEL_M)
print("First-order design")
for k, v in first_order.items():
    print(f"  {k:28} {v:10.4f}")
print("\nPrescription (Ritchey-Chretien), metres")
print(f"  {'surface':10} {'shape':8} {'radius':>8} {'conic':>9} {'to next':>8} {'clear dia.':>10}")
sec_dia = R.APERTURE_M * (T["f1"] - T["d"]) / T["f1"]        # cone from the primary, cut at the secondary
print(f"  {'primary':10} {'concave':8} {T['R1']:8.4f} {T['K1']:9.4f} {T['d']:8.4f} {R.APERTURE_M:10.4f}")
print(f"  {'secondary':10} {'convex':8} {T['R2']:8.4f} {T['K2']:9.4f} {T['d'] + T['b']:8.4f} {sec_dia:10.4f}")
print("  (secondary diameter is the on-axis beam footprint; off-axis fields need a few mm more)")

# check 1: paraxial focus and focal length from a ray 0.1 mm off axis
X, D, _ = trace(T, 0.0, np.array([1e-4]), np.array([0.0]))
z_cross = X[0, 2] - X[0, 0] * D[0, 2] / D[0, 0]
efl = 1e-4 / (-D[0, 0] / D[0, 2])
print(f"\ncheck 1: paraxial focus {z_cross*1e3:.4f} mm behind the primary (design {T['b']*1e3:.1f}), "
      f"focal length {efl:.5f} m (design {R.FOCAL_LENGTH_M:.5f})")
assert abs(z_cross - T["b"]) < 1e-6 and abs(efl - R.FOCAL_LENGTH_M) < 1e-5

# check 2: the coma of a classical Cassegrain equals that of a paraboloid of the same focal ratio: tangential
# coma 3 theta / (16 N^2) radians (Schroeder, Astronomical Optics). Third-order theory holds as the field
# goes to zero, so trace rim rays at small fields and extrapolate the ratio traced/formula to zero field.
Tc = design("classical")
phi = np.linspace(0, 2 * np.pi, 360, endpoint=False)
px, py = R.APERTURE_M / 2 * np.cos(phi), R.APERTURE_M / 2 * np.sin(phi)
fields, ratios = np.array([0.0025, 0.005, 0.01, 0.02]), []
for fd in fields:
    X, D, _ = trace(Tc, fd, px, py)
    Xc, Dc, _ = trace(Tc, fd, np.array([0.0]), np.array([0.0]))
    tc = (to_plane(X, D, Tc["b"])[:, 0] - to_plane(Xc, Dc, Tc["b"])[0, 0]).max()
    ratios.append(tc / (3 * np.deg2rad(fd) / (16 * N ** 2) * R.FOCAL_LENGTH_M))
slope, intercept = np.polyfit(fields, ratios, 1)
print("check 2: classical Cassegrain tangential coma, ray trace / third-order formula")
for fd, r in zip(fields, ratios):
    print(f"  field {fd:.4f} deg: {r:.4f}")
print(f"  extrapolated to zero field: {intercept:.4f} (higher-order terms add {slope:.2f} per degree)")
assert abs(intercept - 1) < 0.01
json.dump(dict(first_order=first_order, prescription=dict(R1=T["R1"], R2=T["R2"], K1=T["K1"], K2=T["K2"], d=T["d"],
          b=T["b"], secondary_dia=sec_dia), checks=dict(coma_ratio_at_zero_field=float(intercept))),
          open("results_design.json", "w"), indent=1)
