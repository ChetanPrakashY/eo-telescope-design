# Optical design: 30 cm telescope for 1 m ground sampling from 500 km

Every number here is produced by the scripts in this folder (`../run_all.sh` rebuilds them).

## 1. Requirements

From `requirements.py`. These are design choices for this study, not a real mission's numbers.

| Item | Value | Reason |
|---|---|---|
| Altitude | 500 km | common small-satellite EO orbit |
| Ground sample distance | 1.0 m | target product |
| Detector pixel | 5.5 um | TDI CMOS class (TBC) |
| Aperture | 0.30 m, 30% linear obscuration | small-satellite envelope |
| Band | 450-850 nm, design wavelength 650 nm | panchromatic |
| Image quality | system MTF at Nyquist >= 0.10 | optics x pixel |
| Radiometry | SNR >= 100, 30% reflectance, sun 30 deg from zenith | |

## 2. First-order design

GSD = H p / F gives the focal length; the aperture then fixes the f-number.

| Quantity | Value |
|---|---|
| Focal length | {{design.first_order.focal_length_m:.3f}} m |
| f-number | {{design.first_order.f_number:.2f}} |
| Primary focal length / secondary magnification | {{design.first_order.primary_focal_length_m:.2f}} m / {{design.first_order.magnification:.3f}} |
| Detector Nyquist frequency | {{design.first_order.nyquist_cyc_mm:.1f}} cycles/mm |
| Optical cut-off 1/(lambda N) | {{design.first_order.optical_cutoff_cyc_mm:.1f}} cycles/mm |
| Airy radius | {{design.first_order.airy_radius_um:.2f}} um |
| Q = lambda N / p | {{design.first_order.q_factor:.2f}} (sampling close to Nyquist) |

## 3. Prescription

Ritchey-Chretien: both mirrors hyperbolic, chosen to cancel spherical aberration and coma.

| Surface | Shape | Radius (m) | Conic | Distance to next (m) | Clear diameter (m) |
|---|---|---|---|---|---|
| Primary | concave | {{design.prescription.R1:.4f}} | {{design.prescription.K1:.4f}} | {{design.prescription.d:.4f}} | 0.300 |
| Secondary | convex | {{design.prescription.R2:.4f}} | {{design.prescription.K2:.4f}} | {{design.prescription.d:.4f}} + {{design.prescription.b:.3f}} | {{design.prescription.secondary_dia:.4f}} (on axis) |

**Checking the ray tracer.** `raytrace.py` intersects each ray with each conic exactly. Two independent checks
(`design.py`): a near-axis ray focuses exactly at the design back focus with the design focal length; and the coma
of a classical Cassegrain traced at small fields, extrapolated to zero field, matches the third-order formula
3 theta / (16 N^2) to a ratio of {{design.checks.coma_ratio_at_zero_field:.4f}}.

## 4. Image quality across the field

`performance.py` traces each field angle, computes the wavefront error on the reference sphere, the PSF by FFT,
and the MTF, then multiplies by the pixel MTF ({{performance.pixel_mtf:.3f}} at Nyquist). One flat detector is
placed {{performance.detector_offset_um:+.0f}} um from the paraxial focus, the position that gives the widest field
meeting the requirement.

| Field (deg) | RMS spot (um) | Strehl | MTF at Nyquist, x | MTF at Nyquist, y | Best focus vs detector (um) |
|---|---|---|---|---|---|
| {{performance.rows.0.field:.2f}} | {{performance.rows.0.rms_spot_um:.1f}} | {{performance.rows.0.strehl:.2f}} | {{performance.rows.0.mtf_x:.3f}} | {{performance.rows.0.mtf_y:.3f}} | {{performance.rows.0.best_focus_shift_um:+.0f}} |
| {{performance.rows.2.field:.2f}} | {{performance.rows.2.rms_spot_um:.1f}} | {{performance.rows.2.strehl:.2f}} | {{performance.rows.2.mtf_x:.3f}} | {{performance.rows.2.mtf_y:.3f}} | {{performance.rows.2.best_focus_shift_um:+.0f}} |
| {{performance.rows.4.field:.2f}} | {{performance.rows.4.rms_spot_um:.1f}} | {{performance.rows.4.strehl:.2f}} | {{performance.rows.4.mtf_x:.3f}} | {{performance.rows.4.mtf_y:.3f}} | {{performance.rows.4.best_focus_shift_um:+.0f}} |
| {{performance.rows.6.field:.2f}} | {{performance.rows.6.rms_spot_um:.1f}} | {{performance.rows.6.strehl:.2f}} | {{performance.rows.6.mtf_x:.3f}} | {{performance.rows.6.mtf_y:.3f}} | {{performance.rows.6.best_focus_shift_um:+.0f}} |
| {{performance.rows.8.field:.2f}} | {{performance.rows.8.rms_spot_um:.1f}} | {{performance.rows.8.strehl:.2f}} | {{performance.rows.8.mtf_x:.3f}} | {{performance.rows.8.mtf_y:.3f}} | {{performance.rows.8.best_focus_shift_um:+.0f}} |
| {{performance.rows.12.field:.2f}} | {{performance.rows.12.rms_spot_um:.1f}} | {{performance.rows.12.strehl:.2f}} | {{performance.rows.12.mtf_x:.3f}} | {{performance.rows.12.mtf_y:.3f}} | {{performance.rows.12.best_focus_shift_um:+.0f}} |

![performance](performance.png)

**Finding.** With a flat detector the requirement holds out to {{performance.max_field_deg:.2f}} deg, a swath of
{{performance.swath_km:.1f}} km. The limit is field curvature: the best-focus surface bends away from the flat
detector by more than a millimetre at the edge of a wider field. If a field-flattener lens (or a curved detector)
removed the curvature, astigmatism would become the limit at {{performance.max_field_flattened_deg:.2f}} deg,
a {{performance.swath_flattened_km:.1f}} km swath. So the design choice is: accept a {{performance.swath_km:.1f}} km
swath with two mirrors only, or add a refractive field corrector for about {{performance.swath_flattened_km:.0f}} km
(the corrector design is the next step and is not done here).

## 5. Tolerances on the secondary mirror

`tolerance.py` asks how large each error can be before the worst-case MTF over the usable field drops by 0.02,
with the detector refocused after assembly (the focus is a compensator).

| Error | Tolerance |
|---|---|
| Despace (along the axis) | +/-{{tolerance.tolerances.despace_um:.0f}} um (absorbed by refocusing; focus moves {{tolerance.focus_gain:.2f}} x the despace) |
| Decentre (sideways) | +/-{{tolerance.tolerances.decentre_um:.0f}} um |
| Tilt | +/-{{tolerance.tolerances.tilt_arcsec:.0f}} arcsec |

Monte Carlo of 150 builds with all three errors drawn uniformly within those tolerances: median worst-case MTF
{{tolerance.mc_median:.3f}}, 10th percentile {{tolerance.mc_p10:.3f}}; {{tolerance.mc_yield:.0%}} of builds meet
the requirement without any adjustment beyond focus. The rest need the secondary adjusted during alignment (below).

## 6. Alignment plan

The sensitivities say what to adjust and what each symptom means.

1. **Primary to the mechanical axis.** Measure the primary's tilt with an autocollimator against the bulkhead reference.
2. **Secondary spacing.** Set with a gauge to well within the despace tolerance; fine focus comes later, at the detector
   (a 10 um spacing change moves the focus about {{tolerance.focus_gain:.0f}} times as far).
3. **Remove coma on axis.** In double pass against a reference flat (interferometer at focus), or with a star test from a
   collimator: on-axis coma means the secondary is decentred or tilted. Adjust secondary tip/tilt until the on-axis image
   is symmetric.
4. **Detector focus.** Through-focus MTF scans with a slanted-edge target in a collimator at the centre and two field points;
   place the detector where the worst of the three is highest.
5. **Verify.** MTF at Nyquist at each field point against the requirement; record the as-built values.
6. **Thermal check.** Repeat step 4 at two temperatures; the drift must agree with `../mechanical/thermal_focus.py`.

## 7. Radiometry

`radiometry.py`: blackbody sun (which reproduces the 1361 W/m^2 solar constant), Lambertian ground, no atmosphere,
optics transmission 0.8 and quantum efficiency 0.6 (both TBC). The line time is {{radiometry.t_line_us:.0f}} us and one
detector row collects {{radiometry.signal_per_stage_e:.0f}} electrons from the 30% target, SNR
{{radiometry.rows.0.snr_cmos:.0f}}. Time-delay integration sums several rows looking at the same ground line:
**{{radiometry.stages_needed}} TDI stages** reach SNR >= 100 even when every stage adds its own read noise.
Adding the atmosphere (transmission below 1, path radiance) will change the count, which is why the TDI sensor should
offer more stages than this minimum.

## 8. Stray light (first look)

Paths to check in the next phase: direct sky light past the secondary onto the focal plane (needs a primary baffle tube
and a secondary baffle cone sized so the detector cannot see the sky directly); scatter from mirror roughness and
particles (BSDF of the coatings); and the bright Earth just outside the field. The usual figure of merit is the point
source transmittance against off-axis angle, from a non-sequential ray trace. Not done here.
