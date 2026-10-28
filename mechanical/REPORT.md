# Mechanical design: holding the telescope in focus and in one piece

The optics (`../optics/REPORT.md`) set the mechanical problem: keep the secondary mirror within tolerance of the
primary through launch, gravity release and temperature changes, at the lowest mass. Numbers come from the scripts
in this folder.

## 1. Thermal focus stability: the design driver

When the metering structure between the mirrors warms, the mirror spacing grows and the focus moves by about
(m^2 + 1) times as much: the secondary magnification turns microns into tens of microns. `thermal_focus.py` traces
the spacing change through the optical model (glass-ceramic mirrors, which hardly expand) and compares the image
motion with the defocus the optics can accept: 25 um costs 0.02 of MTF at Nyquist
in this design. (The classical depth of focus +/-2 lambda N^2 would allow 109 um; this
design's small MTF margin is the tighter limit.)

| Metering structure | CTE (ppm/K) | Image motion (um per K) | Allowed uniform temperature change |
|---|---|---|---|
| Aluminium 6061-T6 | 23.6 | -193.4 | 0.1 K |
| Titanium Ti-6Al-4V | 8.6 | -70.5 | 0.4 K |
| Invar 36 | 1.3 | -10.7 | 2.3 K |
| CFRP quasi-isotropic | 1.0 | -8.2 | 3.0 K |

An aluminium structure with glass mirrors would need the telescope held to a tenth of a kelvin, which is not practical.
Low-expansion structures (CFRP, Invar) allow a few kelvin, which a heater and insulation can hold. The other classic
answer is an all-aluminium telescope with aluminium mirrors: then everything expands together, the design scales
uniformly and stays in focus for uniform temperature changes (it still suffers from gradients). A focus mechanism is the
third option. The CFRP CTE is layup-dependent and must be measured on a sample.

## 2. Stiffness, strength and gravity

`structure.py` models the tube carrying the secondary assembly as a thin-walled cantilever from the primary bulkhead
(0.34 m diameter, 0.0020 m wall, 0.60 m long, 0.5 kg
at the tip), with finite elements (checked against the closed-form tip-mass result within 2%). Launch random vibration
uses Miles' equation at the first mode, with an assumed flat input of 0.1 g^2/Hz and Q = 10,
and a 3-sigma equivalent static load.

| Material | Tube mass (kg) | First mode (Hz) | 1 g sag of the secondary (um) | Root stress (MPa) | Factor of safety |
|---|---|---|---|---|---|
| Aluminium 6061-T6 | 3.46 | 750 | 0.6 | 4.4 | 63 |
| Titanium Ti-6Al-4V | 5.68 | 814 | 0.5 | 6.4 | 138 |
| Invar 36 | 10.32 | 715 | 0.7 | 9.6 | 29 |
| CFRP quasi-isotropic | 2.00 | 759 | 0.6 | 3.3 | 92 |

A closed tube of this size is very stiff: every material clears the assumed 100 Hz requirement by a
wide margin, the stresses are small, and the gravity sag (which disappears in orbit) is far below the optical decentre
tolerance. Stiffness and strength therefore do not choose the material; thermal stability and mass do, and both point to
CFRP: the lightest tube and one of the two low-expansion options. The weak point in a real telescope is usually the spider
that holds the secondary inside the tube, which is not modelled here and is the next analysis.

## 3. Composite stiffness

`clt.py` computes the CFRP modulus used above from fibre and matrix properties: rule of mixtures for one ply (60% carbon
fibre in epoxy) and classical lamination theory for a quasi-isotropic [0/45/-45/90]s stack, which has the same stiffness
in every in-plane direction (about 52 GPa).

## 4. Tolerance stack for the secondary's position

`stack.py` adds the part tolerances between the primary and the secondary mirror (assumed machining values, to be agreed
with the workshop): worst case +/-165 um, root-sum-square +/-72 um, against
+/-194 um of decentre the optics accept. Even the worst case fits, so the secondary needs
tip/tilt adjustment during alignment (to remove residual coma) but no lateral adjustment mechanism.

## 5. Next steps

Spider design and its modes; a thermal model with gradients (sun on one side of the tube); measured CFRP CTE; CAD of the
bulkhead, tube, spider and mirror cells with flexure mounts for the mirrors, so the structure can expand without stressing
the glass.
