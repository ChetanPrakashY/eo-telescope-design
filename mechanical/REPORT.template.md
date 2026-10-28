# Mechanical design: holding the telescope in focus and in one piece

The optics (`../optics/REPORT.md`) set the mechanical problem: keep the secondary mirror within tolerance of the
primary through launch, gravity release and temperature changes, at the lowest mass. Numbers come from the scripts
in this folder.

## 1. Thermal focus stability: the design driver

When the metering structure between the mirrors warms, the mirror spacing grows and the focus moves by about
(m^2 + 1) times as much: the secondary magnification turns microns into tens of microns. `thermal_focus.py` traces
the spacing change through the optical model (glass-ceramic mirrors, which hardly expand) and compares the image
motion with the defocus the optics can accept: {{thermal.allowed_defocus_um:.0f}} um costs 0.02 of MTF at Nyquist
in this design. (The classical depth of focus +/-2 lambda N^2 would allow {{thermal.dof_classical_um:.0f}} um; this
design's small MTF margin is the tighter limit.)

| Metering structure | CTE (ppm/K) | Image motion (um per K) | Allowed uniform temperature change |
|---|---|---|---|
| {{thermal.rows.0.material}} | {{thermal.rows.0.cte_ppm:.1f}} | {{thermal.rows.0.shift_um_per_K:.1f}} | {{thermal.rows.0.allowed_dT_K:.1f}} K |
| {{thermal.rows.1.material}} | {{thermal.rows.1.cte_ppm:.1f}} | {{thermal.rows.1.shift_um_per_K:.1f}} | {{thermal.rows.1.allowed_dT_K:.1f}} K |
| {{thermal.rows.2.material}} | {{thermal.rows.2.cte_ppm:.1f}} | {{thermal.rows.2.shift_um_per_K:.1f}} | {{thermal.rows.2.allowed_dT_K:.1f}} K |
| {{thermal.rows.3.material}} | {{thermal.rows.3.cte_ppm:.1f}} | {{thermal.rows.3.shift_um_per_K:.1f}} | {{thermal.rows.3.allowed_dT_K:.1f}} K |

An aluminium structure with glass mirrors would need the telescope held to a tenth of a kelvin, which is not practical.
Low-expansion structures (CFRP, Invar) allow a few kelvin, which a heater and insulation can hold. The other classic
answer is an all-aluminium telescope with aluminium mirrors: then everything expands together, the design scales
uniformly and stays in focus for uniform temperature changes (it still suffers from gradients). A focus mechanism is the
third option. The CFRP CTE is layup-dependent and must be measured on a sample.

## 2. Stiffness, strength and gravity

`structure.py` models the tube carrying the secondary assembly as a thin-walled cantilever from the primary bulkhead
({{structure.D:.2f}} m diameter, {{structure.wall:.4f}} m wall, {{structure.L:.2f}} m long, {{structure.m_tip:.1f}} kg
at the tip), with finite elements (checked against the closed-form tip-mass result within 2%). Launch random vibration
uses Miles' equation at the first mode, with an assumed flat input of {{structure.asd}} g^2/Hz and Q = {{structure.Q:.0f}},
and a 3-sigma equivalent static load.

| Material | Tube mass (kg) | First mode (Hz) | 1 g sag of the secondary (um) | Root stress (MPa) | Factor of safety |
|---|---|---|---|---|---|
| {{structure.rows.0.material}} | {{structure.rows.0.tube_kg:.2f}} | {{structure.rows.0.f1_hz:.0f}} | {{structure.rows.0.sag_um:.1f}} | {{structure.rows.0.stress_mpa:.1f}} | {{structure.rows.0.fos:.0f}} |
| {{structure.rows.1.material}} | {{structure.rows.1.tube_kg:.2f}} | {{structure.rows.1.f1_hz:.0f}} | {{structure.rows.1.sag_um:.1f}} | {{structure.rows.1.stress_mpa:.1f}} | {{structure.rows.1.fos:.0f}} |
| {{structure.rows.2.material}} | {{structure.rows.2.tube_kg:.2f}} | {{structure.rows.2.f1_hz:.0f}} | {{structure.rows.2.sag_um:.1f}} | {{structure.rows.2.stress_mpa:.1f}} | {{structure.rows.2.fos:.0f}} |
| {{structure.rows.3.material}} | {{structure.rows.3.tube_kg:.2f}} | {{structure.rows.3.f1_hz:.0f}} | {{structure.rows.3.sag_um:.1f}} | {{structure.rows.3.stress_mpa:.1f}} | {{structure.rows.3.fos:.0f}} |

A closed tube of this size is very stiff: every material clears the assumed {{structure.f_min:.0f}} Hz requirement by a
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
with the workshop): worst case +/-{{stack.worst_case_um:.0f}} um, root-sum-square +/-{{stack.rss_um:.0f}} um, against
+/-{{stack.optics_allowance_um:.0f}} um of decentre the optics accept. Even the worst case fits, so the secondary needs
tip/tilt adjustment during alignment (to remove residual coma) but no lateral adjustment mechanism.

## 5. Next steps

Spider design and its modes; a thermal model with gradients (sun on one side of the tube); measured CFRP CTE; CAD of the
bulkhead, tube, spider and mirror cells with flexure mounts for the mirrors, so the structure can expand without stressing
the glass.
