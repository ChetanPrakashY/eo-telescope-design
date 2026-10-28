"""Step 4: signal and SNR per pixel for a pushbroom detector, and the number of TDI stages needed to meet
the SNR requirement. Sun: 5772 K blackbody (IAU 2015 nominal) at 1 au; Lambertian ground; no atmosphere
(a real budget adds atmospheric transmission and path radiance, e.g. from a radiative-transfer code)."""
import json
import numpy as np
from scipy.integrate import quad
import requirements as R

h, c, kB, SB = 6.62607015e-34, 2.99792458e8, 1.380649e-23, 5.670374419e-8
T_SUN, R_SUN, AU, GM, R_E = 5772.0, 6.957e8, 1.495978707e11, 3.986004418e14, 6.371e6
BAND = (450e-9, 850e-9)
RHO, SZA = 0.30, 30.0
TAU, QE, READ_E, FULL_WELL_E = 0.80, 0.60, 6.0, 20000.0       # TBC: optics transmission, detector QE, noise, well

S0 = SB * T_SUN ** 4 * (R_SUN / AU) ** 2
assert abs(S0 - 1361) / 1361 < 0.005                          # the blackbody gives the solar constant
planck = lambda l: 2 * h * c ** 2 / l ** 5 / np.expm1(h * c / (l * kB * T_SUN))
photons = quad(lambda l: planck(l) * np.pi * (R_SUN / AU) ** 2 * l / (h * c), *BAND)[0]   # photons/s/m^2 at 1 au
L_ph = RHO * photons * np.cos(np.deg2rad(SZA)) / np.pi
A_eff = np.pi * R.APERTURE_M ** 2 / 4 * (1 - R.OBSCURATION ** 2)
omega = (R.PIXEL_M / R.FOCAL_LENGTH_M) ** 2
v_ground = np.sqrt(GM / (R_E + R.ALTITUDE_M)) * R_E / (R_E + R.ALTITUDE_M)
t_line = R.GSD_M / v_ground
S1 = L_ph * A_eff * omega * TAU * QE * t_line
print(f"line time {t_line*1e6:.0f} us; signal per stage {S1:.0f} e-; single-stage SNR {S1/np.sqrt(S1 + READ_E**2):.1f}")
rows = []
for n in (1, 4, 8, 16, 32, 64, 96, 128):
    s = n * S1
    rows.append(dict(stages=n, signal_e=s, snr_ccd=s / np.sqrt(s + READ_E ** 2), snr_cmos=s / np.sqrt(s + n * READ_E ** 2),
                     ccd_well_full_on_white=bool(s / RHO > FULL_WELL_E)))   # a 100% reflector accumulated in one CCD well
need = next(r["stages"] for r in rows if r["snr_cmos"] >= R.SNR_MIN)
print(f"{'TDI':>4} {'signal e-':>10} {'SNR (one read)':>15} {'SNR (read per stage)':>21}")
for r in rows:
    print(f"{r['stages']:4d} {r['signal_e']:10.0f} {r['snr_ccd']:15.1f} {r['snr_cmos']:21.1f}"
          + ("  (CCD TDI: a white target overfills the well)" if r["ccd_well_full_on_white"] else ""))
print(f"per stage, a white (100%) target gives {S1 / RHO:.0f} e-, below the {FULL_WELL_E:.0f} e- well, so digital TDI never saturates")
print(f"SNR >= {R.SNR_MIN:.0f} for the 30% target needs {need} TDI stages (digital TDI, read noise per stage)")
json.dump(dict(t_line_us=t_line * 1e6, signal_per_stage_e=S1, stages_needed=need, rows=rows), open("results_radiometry.json", "w"), indent=1)
