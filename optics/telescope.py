"""The telescope model: Ritchey-Chretien design from the requirements, and image-quality metrics computed
from exact ray traces (spot size, wavefront error, Strehl ratio, MTF)."""
import numpy as np
from scipy.optimize import minimize_scalar
import requirements as R
from raytrace import cassegrain, trace, to_plane, pupil_grid, opd

F1 = R.PRIMARY_F_NUMBER * R.APERTURE_M
M = R.FOCAL_LENGTH_M / F1
NYQ = 1 / (2 * R.PIXEL_M * 1e3)                          # cycles per mm


def design(kind="rc"):
    return cassegrain(F1, M, R.BACK_FOCUS_M, kind)


def spot(T, field_deg, z, n=41, **pert):
    xx, yy, keep = pupil_grid(R.APERTURE_M, n, R.OBSCURATION)
    X, D, _ = trace(T, field_deg, xx[keep], yy[keep], **pert)
    return to_plane(X, D, z)[:, :2]


def rms_radius(Q):
    return np.sqrt(((Q - Q.mean(0)) ** 2).sum(1).mean())


def best_focus(T, field_deg, **pert):
    f = lambda z: rms_radius(spot(T, field_deg, z, n=25, **pert))
    return minimize_scalar(f, bounds=(T["b"] - 4e-3, T["b"] + 4e-3), method="bounded", options=dict(xatol=1e-8)).x


def wavefront(T, field_deg, z_image, n=96, remove_defocus=False, **pert):
    """Wavefront error on an n x n pupil grid (NaN outside the annulus), piston and tilt removed."""
    xx, yy, keep = pupil_grid(R.APERTURE_M, n, R.OBSCURATION)
    X, D, opl = trace(T, field_deg, xx[keep], yy[keep], **pert)
    Xc, Dc, _ = trace(T, field_deg, np.array([0.0]), np.array([0.0]), **pert)
    P = to_plane(Xc, Dc, z_image)[0]
    w = opd(X, D, opl, P, np.linalg.norm(P - np.array([0, 0, -T["d"]])))
    u, v = xx[keep] / (R.APERTURE_M / 2), yy[keep] / (R.APERTURE_M / 2)
    cols = [np.ones_like(u), u, v] + ([u ** 2 + v ** 2] if remove_defocus else [])
    A = np.stack(cols, 1)
    w = w - A @ np.linalg.lstsq(A, w, rcond=None)[0]
    W = np.full(xx.shape, np.nan); W[keep] = w
    return W


def psf_mtf(W, pad=4):
    """PSF and MTF from a wavefront map; returns MTF along x and y, frequency axis in cycles/mm."""
    n = W.shape[0]
    pupil = np.where(np.isnan(W), 0, np.exp(2j * np.pi * np.nan_to_num(W) / R.WAVELENGTH_M))
    big = np.zeros((n * pad, n * pad), complex); big[:n, :n] = pupil
    psf = np.abs(np.fft.fft2(big)) ** 2
    otf = np.abs(np.fft.fft2(psf)); otf /= otf[0, 0]
    dx = R.WAVELENGTH_M * R.FOCAL_LENGTH_M / (pad * R.APERTURE_M)   # image-plane sample, m
    nu = np.fft.fftfreq(n * pad, d=dx * 1e3)
    k = slice(0, n * pad // 2)
    perfect = np.where(np.isnan(W), 0, 1.0); bp = np.zeros_like(big); bp[:n, :n] = perfect
    strehl = psf.max() / (np.abs(np.fft.fft2(bp)) ** 2).max()
    return nu[k], otf[0, k], otf[k, 0], strehl


def mtf_at(nu, mtf, f):
    return float(np.interp(f, nu, mtf))


def pixel_mtf(f_cyc_mm):
    return abs(np.sinc(f_cyc_mm * R.PIXEL_M * 1e3))
