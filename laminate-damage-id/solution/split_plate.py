"""Open-branch midplane split on the verified nine-mode plate."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from core import Linear, laminate


def sublaminate_D(angles, ply_t, mat):
    E1, E2, G, nu = mat["E1_Pa"], mat["E2_Pa"], mat["G12_Pa"], mat["nu12"]
    den = 1 - nu * nu * E2 / E1
    q11, q22, q12, q66 = E1 / den, E2 / den, nu * E2 / den, G
    D = np.zeros((3, 3))
    h = ply_t * len(angles)
    for i, ang in enumerate(angles):
        c, t = np.cos(np.deg2rad(ang)), np.sin(np.deg2rad(ang))
        c2, t2 = c * c, t * t
        q = np.zeros((3, 3))
        q[0, 0] = q11 * c**4 + 2 * (q12 + 2 * q66) * c2 * t2 + q22 * t**4
        q[1, 1] = q11 * t**4 + 2 * (q12 + 2 * q66) * c2 * t2 + q22 * c**4
        q[0, 1] = q[1, 0] = (q11 + q22 - 4 * q66) * c2 * t2 + q12 * (c**4 + t**4)
        q[0, 2] = q[2, 0] = (q11 - q12 - 2 * q66) * c**3 * t + (q12 - q22 + 2 * q66) * c * t**3
        q[1, 2] = q[2, 1] = (q11 - q12 - 2 * q66) * c * t**3 + (q12 - q22 + 2 * q66) * c**3 * t
        q[2, 2] = (q11 + q22 - 2 * q12 - 2 * q66) * c2 * t2 + q66 * (c**4 + t**4)
        z0 = -h / 2 + i * ply_t
        z1 = z0 + ply_t
        D += q * (z1**3 - z0**3) / 3.0
    return D


class SplitPlate:
    """Bending operator in Linear.stiffness units, including edge springs."""

    def __init__(self, spec, sensors):
        self.s = spec
        self.lin = Linear(spec, sensors)
        lay = spec["layup_deg"]
        t = spec["ply_thickness_m"]
        mat = spec["material"]
        self.D_split = sublaminate_D(lay[:4], t, mat) + sublaminate_D(lay[4:], t, mat)
        groups = laminate(spec)
        self.D_intact = groups[0] + groups[1]
        x, y = self.lin.xy.T
        al, be = self.lin.al, self.lin.be
        sn = np.sin(x[:, None] * al) * np.sin(y[:, None] * be)
        cs = np.cos(x[:, None] * al) * np.cos(y[:, None] * be)
        self.B = np.stack([sn * al**2, sn * be**2, -2 * cs * al * be], axis=1)
        nx, ny = spec["model"]["gauss_orders"]
        _, wx = leggauss(nx)
        _, wy = leggauss(ny)
        self.wt = np.outer(wx, wy).ravel() * self.lin.L * self.lin.W / 4 / self.lin.mass
        self.delta = self.D_split - self.D_intact

    def separated(self, p):
        """Separated fraction. Each zone is a C1 ellipse that is identically zero outside its support."""
        remain = np.ones(len(self.lin.xy))
        xy = self.lin.xy
        n_zones = (len(p) - 2) // 6
        for j in range(n_zones):
            xc, yc, a, b, theta, severity = p[6 * j:6 * j + 6]
            c, t = np.cos(np.deg2rad(theta)), np.sin(np.deg2rad(theta))
            dx, dy = xy[:, 0] - xc, xy[:, 1] - yc
            u, v = c * dx + t * dy, -t * dx + c * dy
            rho = (u / a) ** 2 + (v / b) ** 2
            profile = np.where(rho < 1.0, (1.0 - rho) ** 2, 0.0)
            remain *= 1.0 - severity * profile
        return 1.0 - remain

    def bending(self, p, open_branch):
        retain = 1.0 if open_branch else float(self.s["damage_model"]["closed_bending_retain"])
        sep = self.separated(p) * retain
        T = np.einsum("qci,cd,qdj,q->ij", self.B, self.D_intact, self.B, self.wt)
        T = T + np.einsum("qci,cd,qdj,q,q->ij", self.B, self.delta, self.B, sep, self.wt)
        fx, fy = p[-2], p[-1]
        springs = self.lin.edge[0] * fx + self.lin.edge[1] + self.lin.edge[2] + self.lin.edge[3] * fy
        return T + springs
