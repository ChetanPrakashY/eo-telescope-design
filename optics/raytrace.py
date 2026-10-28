"""Exact (non-paraxial) ray tracing through a two-mirror telescope with conic mirrors.
Light enters travelling +z, reflects off the primary (vertex at z=0), then off the secondary
(vertex at z=-d), and comes to focus behind the primary at z=+b. Units: metres."""
import numpy as np

def cassegrain(f1, m, b, kind):
    """Design equations: primary focal length f1, secondary magnification m, back focus b."""
    p = (f1 + b) / (m + 1)                  # secondary to primary focus (virtual object)
    d = f1 - p                              # mirror separation
    R1, R2 = 2 * f1, 2 * p * m * p / (m * p - p)   # R2 = 2 p p' / (p' - p), p' = m p
    beta = b / f1
    if kind == "classical":                 # parabolic primary, hyperbolic secondary
        K1, K2 = -1.0, -((m + 1) / (m - 1)) ** 2
    else:                                   # Ritchey-Chretien: no spherical aberration, no coma
        K1 = -1 - 2 * (1 + beta) / (m**2 * (m - beta))
        K2 = -((m + 1) / (m - 1)) ** 2 - 2 * m * (m + 1) / ((m - beta) * (m - 1) ** 3)
    return dict(f1=f1, m=m, b=b, d=d, R1=R1, R2=R2, K1=K1, K2=K2, F=m * f1, kind=kind)

def _hit(P, D, c, K):
    """Intersect rays with the conic c(x^2+y^2+(1+K)z^2) - 2z = 0 (vertex at origin); reflect."""
    a = c * (D[:, 0]**2 + D[:, 1]**2 + (1 + K) * D[:, 2]**2)
    bb = 2 * c * (P[:, 0]*D[:, 0] + P[:, 1]*D[:, 1] + (1 + K)*P[:, 2]*D[:, 2]) - 2 * D[:, 2]
    cc = c * (P[:, 0]**2 + P[:, 1]**2 + (1 + K) * P[:, 2]**2) - 2 * P[:, 2]
    disc = np.sqrt(bb**2 - 4 * a * cc)
    q = -0.5 * (bb + np.sign(bb) * disc)    # stable quadratic roots: q/a and cc/q
    t1 = np.where(np.abs(a) > 1e-15, q / np.where(a == 0, 1, a), np.inf)
    t2 = cc / q
    # a hyperboloid has two sheets: take the root that lands nearest the vertex (|z| smallest)
    z1, z2 = P[:, 2] + t1 * D[:, 2], P[:, 2] + t2 * D[:, 2]
    t = np.where(np.abs(z2) < np.abs(z1), t2, t1)
    X = P + t[:, None] * D
    n = np.stack([2*c*X[:, 0], 2*c*X[:, 1], 2*c*(1 + K)*X[:, 2] - 2], axis=1)
    n /= np.linalg.norm(n, axis=1)[:, None]
    D2 = D - 2 * (D * n).sum(1)[:, None] * n
    return X, D2, t

def _rot_y(v, ang):
    c, s = np.cos(ang), np.sin(ang)
    return np.stack([c * v[:, 0] + s * v[:, 2], v[:, 1], -s * v[:, 0] + c * v[:, 2]], axis=1)

def trace(T, field_deg, px, py, dz2=0.0, dx2=0.0, tilt2=0.0):
    """Trace rays at entrance-pupil coordinates (px, py) for a field angle in the x-z plane.
    dz2, dx2, tilt2 perturb the secondary (despace towards the primary is dz2 > 0).
    Returns the secondary hit points, outgoing directions and optical path lengths there."""
    th = np.deg2rad(field_deg)
    D = np.tile([np.sin(th), 0.0, np.cos(th)], (px.size, 1))
    P = np.stack([px, py, np.full(px.size, -0.01)], axis=1) - D * 0.0
    opl = (P * D).sum(1)                          # phase of the incoming plane wave at P
    X1, D1, t1 = _hit(P, D, -1 / T["R1"], T["K1"])
    opl = opl + t1
    # secondary local frame: vertex at z = -d + dz2, decentred by dx2, tilted by tilt2 about y.
    # It faces +z; flipping z makes the light travel +z in the local frame, as the routine expects.
    v0 = np.array([dx2, 0.0, -T["d"] + dz2])
    Pl = _rot_y(X1 - v0, -tilt2); Dl = _rot_y(D1, -tilt2)
    Pl[:, 2] *= -1; Dl[:, 2] *= -1
    X2l, D2l, t2 = _hit(Pl, Dl, 1 / T["R2"], T["K2"])  # convex towards +z: positive sag after the flip
    X2l[:, 2] *= -1; D2l[:, 2] *= -1
    X2 = _rot_y(X2l, tilt2) + v0; D2 = _rot_y(D2l, tilt2)
    return X2, D2, opl + t2

def to_plane(X, D, z):
    t = (z - X[:, 2]) / D[:, 2]
    return X + t[:, None] * D

def pupil_grid(D_ap, n, obsc):
    """Square grid over the entrance pupil, keeping the annulus between the obscuration and the rim."""
    u = (np.arange(n) - (n - 1) / 2) * D_ap / n
    xx, yy = np.meshgrid(u, u)
    r = np.hypot(xx, yy) / (D_ap / 2)
    keep = (r <= 1) & (r >= obsc)
    return xx, yy, keep

def opd(X, D, opl, Pimg, Rref):
    """Optical path difference on the reference sphere of radius Rref centred on the image point."""
    S = X - Pimg
    bq = (S * D).sum(1); cq = (S * S).sum(1) - Rref**2
    t = -bq + np.sign(bq) * np.sqrt(bq**2 - cq)   # root nearest t = 0
    t = np.where(np.abs(-bq - np.sqrt(bq**2 - cq)) < np.abs(t), -bq - np.sqrt(bq**2 - cq), t)
    w = opl + t
    return w - w.mean()
