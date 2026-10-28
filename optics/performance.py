"""Step 2: image quality across the field. For each field angle: RMS spot, Strehl ratio and the system MTF
at the detector Nyquist frequency (optics x pixel), on one flat detector placed at the focus that maximises
the worst-case MTF. The largest field meeting the MTF requirement sets the swath."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import requirements as R
from telescope import design, spot, rms_radius, best_focus, wavefront, psf_mtf, mtf_at, pixel_mtf, NYQ

T = design("rc")
fields = np.round(np.arange(0.0, 0.61, 0.05), 3)
pix = pixel_mtf(NYQ)

def metrics(z):
    out = []
    for fd in fields:
        nu, mx, my, strehl = psf_mtf(wavefront(T, fd, z, n=96))
        out.append(dict(field=float(fd), mtf_x=mtf_at(nu, mx, NYQ) * pix, mtf_y=mtf_at(nu, my, NYQ) * pix, strehl=float(strehl)))
    return out

def reach(rows):
    """Largest field angle up to which every field point meets the MTF requirement (-1 if the axis fails)."""
    best = -1.0
    for r in rows:
        if min(r["mtf_x"], r["mtf_y"]) < R.MTF_NYQUIST_MIN:
            break
        best = r["field"]
    return best

# place the flat detector where the usable field is widest (ties: the higher worst-case MTF)
zs = T["b"] + np.linspace(-0.4e-3, 0.1e-3, 11)
cands = []
for z in zs:
    m = metrics(z); f = reach(m)
    cands.append((f, min(min(r["mtf_x"], r["mtf_y"]) for r in m if r["field"] <= max(f, 0)), z, m))
max_field, _, z0, rows = max(cands, key=lambda c: (c[0], c[1]))
z0 = float(z0)
for r, fd in zip(rows, fields):
    r["rms_spot_um"] = float(rms_radius(spot(T, fd, z0)) * 1e6)
    r["best_focus_shift_um"] = float((best_focus(T, fd) - z0) * 1e6)
swath_km = 2 * R.ALTITUDE_M * np.tan(np.deg2rad(max_field)) / 1e3
assert max_field > 0, "requirement not met even on axis"

# upper bound: if a field flattener (or a curved detector) removed the field curvature, each field would be
# imaged at its own best focus; this is the wavefront with defocus removed
flat_rows = []
for fd in fields:
    nu, mx, my, strehl = psf_mtf(wavefront(T, fd, T["b"], n=96, remove_defocus=True))
    flat_rows.append(dict(field=float(fd), mtf_x=mtf_at(nu, mx, NYQ) * pix, mtf_y=mtf_at(nu, my, NYQ) * pix))
max_field_flat = reach(flat_rows)
swath_flat_km = 2 * R.ALTITUDE_M * np.tan(np.deg2rad(max_field_flat)) / 1e3

print(f"detector placed {1e6*(z0 - T['b']):+.0f} um from the paraxial focus; pixel MTF at Nyquist {pix:.3f}")
print(f"{'field':>6} {'RMS spot':>9} {'Strehl':>7} {'MTF x':>6} {'MTF y':>6} {'best focus':>11}")
for r in rows:
    print(f"{r['field']:5.2f}d {r['rms_spot_um']:8.2f}u {r['strehl']:7.3f} {r['mtf_x']:6.3f} {r['mtf_y']:6.3f} {r['best_focus_shift_um']:+9.0f}um")
print(f"flat detector: MTF at Nyquist >= {R.MTF_NYQUIST_MIN} out to {max_field:.2f} deg -> swath {swath_km:.1f} km at {R.ALTITUDE_M/1e3:.0f} km")
print(f"field curvature removed (upper bound with a field flattener): out to {max_field_flat:.2f} deg -> swath {swath_flat_km:.1f} km")
json.dump(dict(rows=rows, detector_offset_um=1e6 * (z0 - T["b"]), pixel_mtf=pix, max_field_deg=max_field, swath_km=swath_km,
               max_field_flattened_deg=max_field_flat, swath_flattened_km=swath_flat_km, flattened=flat_rows),
          open("results_performance.json", "w"), indent=1)

fig, ax = plt.subplots(1, 3, figsize=(11, 3.4))
f = [r["field"] for r in rows]
ax[0].plot(f, [r["mtf_x"] for r in rows], "o-", ms=3, label="x (cross-track)")
ax[0].plot(f, [r["mtf_y"] for r in rows], "s-", ms=3, label="y (along-track)")
ax[0].plot(f, [min(r["mtf_x"], r["mtf_y"]) for r in flat_rows], ":", color="k", lw=1, label="worse axis, curvature removed")
ax[0].axhline(R.MTF_NYQUIST_MIN, color="r", lw=0.8, ls="--", label="requirement")
ax[0].axvline(max_field, color="gray", lw=0.6)
ax[0].set(xlabel="field angle (deg)", ylabel="system MTF at Nyquist", title="MTF at 91 cycles/mm vs field"); ax[0].legend(fontsize=7, loc="upper right")
ax[1].plot(f, [r["best_focus_shift_um"] for r in rows], "o-", ms=3)
ax[1].set(xlabel="field angle (deg)", ylabel="best focus - detector (um)", title="field curvature")
for k, fd in enumerate([0.0, max_field / 2, max_field]):
    Q = spot(T, fd, z0) * 1e6; Q -= Q.mean(0)
    xc = (k - 1) * 22
    ax[2].plot(Q[:, 0] + xc, Q[:, 1], ".", ms=1, label=f"{fd:.2f} deg")
    ax[2].add_patch(plt.Rectangle((xc - 2.75, -2.75), 5.5, 5.5, fill=False, lw=1.0, zorder=5))
ax[2].set_aspect("equal"); ax[2].set_xlim(-36, 36); ax[2].set_ylim(-14, 14)
ax[2].set(xlabel="um (square = one 5.5 um pixel)", title="spots across the usable field"); ax[2].legend(fontsize=7, loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.45), markerscale=6)
for a in ax: a.title.set_fontsize(9)
fig.tight_layout(); fig.savefig("performance.png", dpi=150)
