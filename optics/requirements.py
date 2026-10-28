"""Requirements and assumptions for the design study. Everything downstream reads its numbers from here.
Values are choices for this study, not taken from a real mission; each has a one-line reason."""
ALTITUDE_M = 500e3          # a common small-satellite Earth-observation orbit height
GSD_M = 1.0                 # ground sample distance at nadir
PIXEL_M = 5.5e-6            # detector pixel pitch (TDI CMOS class, TBC)
APERTURE_M = 0.30           # largest primary that fits the assumed small-satellite envelope
OBSCURATION = 0.30          # secondary + baffle diameter / primary diameter
WAVELENGTH_M = 0.65e-6      # centre of a 450-850 nm panchromatic band
MTF_NYQUIST_MIN = 0.10      # system MTF at the detector Nyquist frequency, optics x pixel
SNR_MIN = 100.0             # for a 30% reflectance target at 30 deg solar zenith angle
PRIMARY_F_NUMBER = 2.5      # fast primary keeps the telescope short (trade: tighter tolerances)
BACK_FOCUS_M = 0.15         # room behind the primary for the focal plane and its electronics

FOCAL_LENGTH_M = ALTITUDE_M * PIXEL_M / GSD_M           # GSD = H p / F
