"""Material properties used in the trade (room-temperature handbook values; check against the supplier's
data sheet before design). The CFRP modulus is computed from classical lamination theory in clt.py."""
from clt import quasi_isotropic

CFRP_E, CFRP_RHO = quasi_isotropic()
MATERIALS = {   # E (Pa), density (kg/m^3), CTE (1/K), yield or design strength (Pa)
    "Aluminium 6061-T6": dict(E=68.9e9, rho=2700, cte=23.6e-6, strength=276e6),
    "Titanium Ti-6Al-4V": dict(E=113.8e9, rho=4430, cte=8.6e-6, strength=880e6),
    "Invar 36": dict(E=141e9, rho=8050, cte=1.3e-6, strength=276e6),
    "CFRP quasi-isotropic": dict(E=CFRP_E, rho=CFRP_RHO, cte=1.0e-6, strength=300e6),   # CTE and strength TBC: layup dependent
}
