"""In-plane modulus of a quasi-isotropic carbon/epoxy laminate [0/45/-45/90]s from fibre and matrix
properties: rule of mixtures for the ply, classical lamination theory for the stack."""
import numpy as np

def quasi_isotropic(Ef=230e9, Em=3.5e9, nuf=0.2, num=0.35, Vf=0.6, rho_f=1800, rho_m=1200):
    E1 = Vf * Ef + (1 - Vf) * Em
    E2 = 1 / (Vf / Ef + (1 - Vf) / Em)
    nu12 = Vf * nuf + (1 - Vf) * num
    G12 = 1 / (Vf / (Ef / (2 * (1 + nuf))) + (1 - Vf) / (Em / (2 * (1 + num))))
    d = 1 - nu12 ** 2 * E2 / E1
    Q = np.array([[E1 / d, nu12 * E2 / d, 0], [nu12 * E2 / d, E2 / d, 0], [0, 0, G12]])
    A = np.zeros((3, 3))
    for th in (0, 45, -45, 90, 90, -45, 45, 0):
        c, s = np.cos(np.radians(th)), np.sin(np.radians(th))
        T = np.array([[c * c, s * s, 2 * c * s], [s * s, c * c, -2 * c * s], [-c * s, c * s, c * c - s * s]])
        Rm = np.diag([1, 1, 2])
        A += np.linalg.inv(T) @ Q @ Rm @ T @ np.linalg.inv(Rm)       # equal ply thickness: A is a sum
    A /= 8
    Ex = 1 / np.linalg.inv(A)[0, 0]
    return Ex, Vf * rho_f + (1 - Vf) * rho_m

if __name__ == "__main__":
    E, rho = quasi_isotropic()
    print(f"quasi-isotropic CFRP: E = {E/1e9:.1f} GPa, density {rho:.0f} kg/m^3")
