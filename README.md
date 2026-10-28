# eo-telescope-design

Design study of a 30 cm two-mirror telescope for 1 m ground sampling from a 500 km orbit: optical design and
analysis in `optics/`, the structure that holds it in `mechanical/`. Everything is computed by the scripts here,
and the two reports are filled in from their results.

| | Report | What is in it |
|---|---|---|
| Optics | [optics/REPORT.md](optics/REPORT.md) | first-order design from the requirements, Ritchey-Chretien prescription, exact ray tracing (checked against third-order coma theory), wavefront error, Strehl ratio and MTF across the field, the usable swath and what limits it, secondary-mirror tolerances with a Monte Carlo of builds, an alignment plan, a radiometric budget with TDI, and a first look at stray light |
| Mechanical | [mechanical/REPORT.md](mechanical/REPORT.md) | thermal focus stability traced through the optical model (the design driver), a four-material trade (aluminium, titanium, Invar, CFRP), finite-element modes, gravity sag, launch random vibration by Miles' equation, CFRP stiffness from lamination theory, and a tolerance stack for the secondary's position |

Main findings:

- With two mirrors and a flat detector, the image-quality requirement holds over a **3.5 km** swath; field curvature is
  the limit. Removing it (a field-flattener lens) would allow about **6 km**, where astigmatism takes over.
- The secondary mirror tolerates about **194 um** of decentre and **151 arcsec** of tilt (focus refocused), and the
  machining tolerance stack fits inside that.
- An aluminium metering tube with glass mirrors moves the focus **193 um per kelvin**; CFRP or Invar bring that to under
  11 um per kelvin. Stiffness and strength have large margins for every material, so thermal stability and mass decide:
  **CFRP**.
- **4 TDI stages** give SNR 100 on a 30% reflectance target (no atmosphere included).

```bash
./run_all.sh          # about a minute; numpy, scipy, matplotlib only
```

Assumptions (detector, transmission, launch environment, machining tolerances) are marked TBC where they appear. The
spider holding the secondary, the field corrector and a non-sequential stray-light analysis are the next steps and are
not done here.
