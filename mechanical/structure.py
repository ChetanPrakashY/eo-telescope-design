"""The metering tube that holds the secondary mirror, modelled as a thin-walled cantilever (clamped at the
primary bulkhead) carrying the secondary assembly at its tip. For each material: mass, first natural
frequency (finite elements, checked against the closed-form tip-mass result), gravity sag of the secondary
on the ground, and launch random vibration load and stress by Miles' equation."""
import json
import numpy as np
from scipy.linalg import eigh
from materials import MATERIALS

L = 0.60                    # tube length to the secondary: mirror spacing plus the secondary cell
D_TUBE, WALL = 0.34, 2.0e-3 # mean diameter and wall
M_TIP = 0.50                # secondary mirror, cell, spider and hub (TBC)
F_MIN = 100.0               # first-mode requirement (TBC: set by the launcher's user manual)
ASD, Q = 0.1, 10.0          # random vibration input g^2/Hz (TBC) and amplification
A = np.pi * D_TUBE * WALL
I = np.pi * (D_TUBE / 2) ** 3 * WALL          # thin-walled tube


def beam(E, rho, ne=30):
    le = L / ne; n = 2 * (ne + 1)
    K = np.zeros((n, n)); Mm = np.zeros((n, n))
    k = E * I / le ** 3 * np.array([[12, 6 * le, -12, 6 * le], [6 * le, 4 * le ** 2, -6 * le, 2 * le ** 2],
                                    [-12, -6 * le, 12, -6 * le], [6 * le, 2 * le ** 2, -6 * le, 4 * le ** 2]])
    m = rho * A * le / 420 * np.array([[156, 22 * le, 54, -13 * le], [22 * le, 4 * le ** 2, 13 * le, -3 * le ** 2],
                                       [54, 13 * le, 156, -22 * le], [-13 * le, -3 * le ** 2, -22 * le, 4 * le ** 2]])
    for e in range(ne):
        s = slice(2 * e, 2 * e + 4); K[s, s] += k; Mm[s, s] += m
    Mm[-2, -2] += M_TIP
    return K[2:, 2:], Mm[2:, 2:]


rows = []
print(f"tube {D_TUBE*1e3:.0f} mm x {WALL*1e3:.1f} mm wall, {L*1e3:.0f} mm long, {M_TIP} kg at the tip")
print(f"{'material':22} {'tube kg':>7} {'f1 Hz':>7} {'check Hz':>8} {'1 g sag um':>10} {'Grms':>6} {'root stress MPa':>15} {'FoS':>5}")
for name, mt in MATERIALS.items():
    K, Mm = beam(mt["E"], mt["rho"])
    f1 = np.sqrt(eigh(K, Mm, eigvals_only=True)[0]) / (2 * np.pi)
    m_tube = mt["rho"] * A * L
    # closed form: tip mass plus 33/140 of the beam mass on a cantilever spring 3EI/L^3 (Rayleigh)
    f_check = np.sqrt(3 * mt["E"] * I / L ** 3 / (M_TIP + 33 / 140 * m_tube)) / (2 * np.pi)
    F = np.zeros(K.shape[0]); F[-2] = M_TIP * 9.81
    sag = np.linalg.solve(K, F)[-2] + mt["rho"] * A * 9.81 * L ** 4 / (8 * mt["E"] * I)   # tip mass + own weight
    grms = np.sqrt(np.pi / 2 * f1 * Q * ASD)
    load = 3 * grms * 9.81 * (M_TIP + 33 / 140 * m_tube)          # 3-sigma equivalent static tip load
    stress = load * L * (D_TUBE / 2) / I
    fos = mt["strength"] / stress
    rows.append(dict(material=name, tube_kg=m_tube, f1_hz=f1, f_check_hz=f_check, sag_um=sag * 1e6, grms=grms,
                     stress_mpa=stress / 1e6, fos=fos, meets_f1=bool(f1 >= F_MIN)))
    print(f"{name:22} {m_tube:7.2f} {f1:7.0f} {f_check:8.0f} {sag*1e6:10.1f} {grms:6.1f} {stress/1e6:15.1f} {fos:5.1f}")
    assert abs(f1 - f_check) / f1 < 0.02
json.dump(dict(L=L, D=D_TUBE, wall=WALL, m_tip=M_TIP, f_min=F_MIN, asd=ASD, Q=Q, rows=rows), open("results_structure.json", "w"), indent=1)
