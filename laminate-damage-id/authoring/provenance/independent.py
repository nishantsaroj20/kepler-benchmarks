"""Direct-quadrature forward model for the graded snap and forecast.

Newton continuation on strains integrated at the Gauss points. No import of the
reference polynomial condensation or its trust-region equilibrium solver.
"""
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh, solve


def reduced_q(E1, E2, G, nu, angle_deg):
    den = 1.0 - nu * nu * E2 / E1
    q11, q22, q12, q66 = E1 / den, E2 / den, nu * E2 / den, G
    c, t = np.cos(np.deg2rad(angle_deg)), np.sin(np.deg2rad(angle_deg))
    c2, t2 = c * c, t * t
    q = np.zeros((3, 3))
    q[0, 0] = q11 * c**4 + 2 * (q12 + 2 * q66) * c2 * t2 + q22 * t**4
    q[1, 1] = q11 * t**4 + 2 * (q12 + 2 * q66) * c2 * t2 + q22 * c**4
    q[0, 1] = q[1, 0] = (q11 + q22 - 4 * q66) * c2 * t2 + q12 * (c**4 + t**4)
    q[0, 2] = q[2, 0] = (q11 - q12 - 2 * q66) * c**3 * t + (q12 - q22 + 2 * q66) * c * t**3
    q[1, 2] = q[2, 1] = (q11 - q12 - 2 * q66) * c * t**3 + (q12 - q22 + 2 * q66) * c**3 * t
    q[2, 2] = (q11 + q22 - 2 * q12 - 2 * q66) * c2 * t2 + q66 * (c**4 + t**4)
    return q


def bending_about_centroid(angles, ply_t, mat):
    """Bending stiffness of a ply stack about that stack's own midplane."""
    h = ply_t * len(angles)
    d = np.zeros((3, 3))
    for i, ang in enumerate(angles):
        z0 = -0.5 * h + i * ply_t
        z1 = z0 + ply_t
        d += reduced_q(mat["E1_Pa"], mat["E2_Pa"], mat["G12_Pa"], mat["nu12"], ang) * (z1**3 - z0**3) / 3.0
    return d


class DirectPlate:
    """Nine-mode von Karman plate with a breathing midplane split."""

    def __init__(self, spec, quad=(30, 26)):
        self.spec = spec
        self.L = float(spec["Lx_m"])
        self.W = float(spec["Ly_m"])
        lay = spec["layup_deg"]
        t = float(spec["ply_thickness_m"])
        mat = spec["material"]
        self.h = t * len(lay)
        self.A = np.zeros((3, 3))
        for ang in lay:
            self.A += reduced_q(mat["E1_Pa"], mat["E2_Pa"], mat["G12_Pa"], mat["nu12"], ang) * t
        self.D_intact = bending_about_centroid(lay, t, mat)
        self.D_split = bending_about_centroid(lay[:4], t, mat) + bending_about_centroid(lay[4:], t, mat)
        self.delta = self.D_split - self.D_intact
        self.mn = np.array([(m, n) for m in range(1, 4) for n in range(1, 4)])
        self.al = self.mn[:, 0] * np.pi / self.L
        self.be = self.mn[:, 1] * np.pi / self.W
        nx, ny = quad
        xi, wx = leggauss(nx)
        eta, wy = leggauss(ny)
        x = (xi + 1.0) * self.L / 2.0
        y = (eta + 1.0) * self.W / 2.0
        X, Y = np.meshgrid(x, y, indexing="ij")
        self.wt = (np.outer(wx, wy) * self.L * self.W / 4.0).ravel()
        X, Y = X.ravel(), Y.ravel()
        sn = np.sin(X[:, None] * self.al) * np.sin(Y[:, None] * self.be)
        cs = np.cos(X[:, None] * self.al) * np.cos(Y[:, None] * self.be)
        self.curv = np.stack([sn * self.al**2, sn * self.be**2, -2.0 * cs * self.al * self.be], axis=1)
        self.wx = self.h * self.al * np.cos(X[:, None] * self.al) * np.sin(Y[:, None] * self.be)
        self.wy = self.h * self.be * np.sin(X[:, None] * self.al) * np.cos(Y[:, None] * self.be)
        phi_x = self.al * np.cos(X[:, None] * self.al) * np.sin(Y[:, None] * self.be)
        phi_y = self.be * np.sin(X[:, None] * self.al) * np.cos(Y[:, None] * self.be)
        self.Bm = np.zeros((len(X), 3, 18))
        self.Bm[:, 0, :9] = phi_x
        self.Bm[:, 1, 9:] = phi_y
        self.Bm[:, 2, :9] = phi_y
        self.Bm[:, 2, 9:] = phi_x
        self.C = np.einsum("qai,ab,qbj,q->ij", self.Bm, self.A, self.Bm, self.wt, optimize=True)
        self.xy = np.column_stack([X, Y])
        self.Kg = np.einsum("qi,qj,q->ij", self.wx, self.wx, self.wt) + spec["precompression_y_to_x"] * np.einsum(
            "qi,qj,q->ij", self.wy, self.wy, self.wt
        )
        self.Nx = float(spec["precompression_N_m"])
        eta = np.zeros(9)
        eta[: len(spec["imperfection_h_coefficients"])] = spec["imperfection_h_coefficients"]
        self.w0x = self.wx @ eta
        self.w0y = self.wy @ eta
        self.lin_geo = self.Nx * (self.wx.T @ (self.wt * self.w0x)) + spec["precompression_y_to_x"] * self.Nx * (
            self.wy.T @ (self.wt * self.w0y)
        )
        load = np.array([4.0 / (m * n * np.pi**2) if m % 2 and n % 2 else 0.0 for m, n in self.mn])
        self.load = self.h * load
        rho = mat["rho_kg_m3"]
        self.mass = rho * self.h * self.L * self.W / 4.0
        self.omega = 1000.0
        self.E0 = self.mass * self.h**2 * self.omega**2
        self.p = None
        self.open_branch = True

    def separated(self, p):
        remain = np.ones(len(self.xy))
        n_zones = (len(p) - 2) // 6
        for j in range(n_zones):
            xc, yc, a, b, theta, severity = p[6 * j : 6 * j + 6]
            c, t = np.cos(np.deg2rad(theta)), np.sin(np.deg2rad(theta))
            dx, dy = self.xy[:, 0] - xc, self.xy[:, 1] - yc
            u, v = c * dx + t * dy, -t * dx + c * dy
            rho = (u / a) ** 2 + (v / b) ** 2
            remain *= 1.0 - severity * np.where(rho < 1.0, (1.0 - rho) ** 2, 0.0)
        return 1.0 - remain

    def prepare(self, p, open_branch):
        self.p = np.asarray(p, float)
        self.open_branch = bool(open_branch)
        damage = self.spec["damage_model"]
        retain = 1.0 if open_branch else float(damage["closed_bending_retain"])
        membrane_loss = float(damage["open_membrane_loss"]) if open_branch else 0.0
        raw = self.separated(self.p)
        s = raw * retain
        d = self.D_intact + s[:, None, None] * self.delta
        self.A_q = self.A * (1.0 - membrane_loss * raw)[:, None, None]
        self.C = np.einsum("qai,qab,qbj,q->ij", self.Bm, self.A_q, self.Bm, self.wt, optimize=True)
        self.Kb = self.h**2 * np.einsum("qci,qcd,qdj,q->ij", self.curv, d, self.curv, self.wt, optimize=True)
        fx, fy = self.p[-2], self.p[-1]
        scale = {"x0": fx, "xL": 1.0, "y0": 1.0, "yL": fy}
        # Rebuild the scaled spring matrix from the unit-multiplier operator stored in Ks parts.
        self.Kspring = np.zeros((9, 9))
        springs = self.spec["rotational_springs_N"]
        for edge, key, coeff, match, span in (
            (0, "x0", self.al, self.mn[:, 1], self.W),
            (1, "xL", self.al * ((-1.0) ** self.mn[:, 0]), self.mn[:, 1], self.W),
            (2, "y0", self.be, self.mn[:, 0], self.L),
            (3, "yL", self.be * ((-1.0) ** self.mn[:, 1]), self.mn[:, 0], self.L),
        ):
            self.Kspring += scale[key] * springs[key] * self.h**2 * np.outer(coeff, coeff) * (
                match[:, None] == match[None, :]
            ) * span / 2.0

    def gradient(self, q, force):
        """Condensed residual. In-plane coefficients are eliminated at the current q."""
        wx, wy = self.wx @ q, self.wy @ q
        en = np.column_stack(
            [
                self.w0x * wx + 0.5 * wx**2,
                self.w0y * wy + 0.5 * wy**2,
                self.w0x * wy + self.w0y * wx + wx * wy,
            ]
        )
        den = np.empty((len(self.wt), 3, 9))
        den[:, 0, :] = (self.w0x + wx)[:, None] * self.wx
        den[:, 1, :] = (self.w0y + wy)[:, None] * self.wy
        den[:, 2, :] = (self.w0y + wy)[:, None] * self.wx + (self.w0x + wx)[:, None] * self.wy
        membrane = np.einsum("qai,qab,qb,q->i", self.Bm, self.A_q, en, self.wt, optimize=True)
        relaxed = -solve(self.C, membrane, assume_a="pos")
        strain = en + np.einsum("qai,i->qa", self.Bm, relaxed)
        stress = np.einsum("qab,qb->qa", self.A_q, strain)
        membrane_grad = np.einsum("qci,qc,q->i", den, stress, self.wt, optimize=True)
        return (self.Kb @ q + self.Kspring @ q - self.Nx * self.Kg @ q - self.lin_geo + membrane_grad - force * self.load) / self.E0

    def hessian(self, q, force):
        """Central differences of the condensed residual. Step 1e-7 in q."""
        step = 1e-7
        columns = []
        for i in range(9):
            dq = np.zeros(9)
            dq[i] = step
            columns.append((self.gradient(q + dq, force) - self.gradient(q - dq, force)) / (2.0 * step))
        h = np.column_stack(columns)
        return 0.5 * (h + h.T)

    def equilibrium(self, force, branch, start=None):
        q = np.zeros(9) if start is None else np.array(start, float)
        if start is None:
            q[0] = float(branch) * 1.2
        for _ in range(20):
            g = self.gradient(q, force)
            if np.linalg.norm(g) < 1e-11:
                break
            h = self.hessian(q, force)
            step = solve(h, g)
            alpha = 1.0
            base = np.linalg.norm(g)
            for _ls in range(8):
                trial = np.linalg.norm(self.gradient(q - alpha * step, force))
                if trial < base:
                    break
                alpha *= 0.5
            q = q - alpha * step
        g = self.gradient(q, force)
        h = self.hessian(q, force)
        lam = float(np.linalg.eigvalsh(h)[0])
        if np.linalg.norm(g) > 1e-8 or lam <= 0.0 or int(np.sign(q[0])) != int(branch):
            raise RuntimeError("direct equilibrium missed the requested preparation")
        return q, h

    def frequencies(self, force, p, branch):
        self.prepare(p, branch > 0)
        _q, h = self.equilibrium(force, branch)
        val = eigh(h, eigvals_only=True)
        return np.sqrt(val) * self.omega / (2.0 * np.pi), float(val[0])

    def snap(self, p, branch=-1):
        """Brent root of the continued smallest eigenvalue, warm-started in force."""
        self.prepare(p, branch > 0)
        q, _ = self.equilibrium(40.0, branch)
        last_f, last_q = 40.0, q
        hi = None
        for force in np.linspace(48.0, 120.0, 19):
            try:
                q, h = self.equilibrium(float(force), branch, start=q)
            except RuntimeError:
                hi = float(force)
                break
            if float(np.linalg.eigvalsh(h)[0]) <= 0.0:
                hi = float(force)
                break
            last_f, last_q = float(force), q
        if hi is None:
            raise RuntimeError("closed branch stayed stable through 120 N")

        def lam(force):
            nonlocal last_q
            q, h = self.equilibrium(float(force), branch, start=last_q)
            last_q = q
            return float(np.linalg.eigvalsh(h)[0])

        lo = last_f
        for _ in range(28):
            mid = 0.5 * (lo + hi)
            try:
                value = lam(mid)
            except RuntimeError:
                hi = mid
                continue
            if value > 0.0:
                lo = mid
            else:
                hi = mid
        return float(lo)

    def forecast(self, p):
        rows = []
        sign = int(self.spec["forecast"]["preparation_sign"])
        count = int(self.spec["forecast"]["ascending_modes"])
        for force in self.spec["forecast"]["forces_N"]:
            freq, _lam = self.frequencies(float(force), p, sign)
            rows.append({"force_N": float(force), "frequencies_Hz": [float(v) for v in freq[:count]]})
        return rows


def authored_parameters(root):
    spec = json.loads((root / "environment/data/plate_specs.json").read_text())
    plate = json.loads((root / "authoring/provenance/parameters.json").read_text())["plate"]
    keys = ("xc", "yc", "a", "b", "theta_deg", "d")
    p = np.array([z[k] for z in plate["zones"] for k in keys] + list(plate["fixture"].values()), float)
    return spec, plate, p


def graded_forward(spec, p, quad=(30, 26)):
    model = DirectPlate(spec, quad=quad)
    return float(model.snap(p)), model.forecast(p)


if __name__ == "__main__":
    import sys

    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(root / "solution"))
    from solve_damage import Campaign

    spec, plate, p = authored_parameters(root)
    direct = DirectPlate(spec, quad=tuple(spec["model"]["gauss_orders"]))
    refined = DirectPlate(spec, quad=(40, 36))
    camp = Campaign(spec, np.array([[0.2, 0.2]]))
    camp.apply(p, 1)
    rows = []
    for force, branch in ((80.0, 1), (-90.0, 1), (-60.0, 1), (45.0, -1)):
        fd, _ = direct.frequencies(force, p, branch)
        fr = camp.observe(force, p, branch)[1]
        ff, _ = refined.frequencies(force, p, branch)
        rows.append(
            {
                "force_N": force,
                "branch": branch,
                "direct_vs_reference_max_relative": float(np.max(np.abs(fd[:7] / fr[:7] - 1.0))),
                "refined_quadrature_max_relative": float(np.max(np.abs(ff[:7] / fd[:7] - 1.0))),
            }
        )
    snap_direct = direct.snap(p)
    snap_reference = camp.snap_load(p, branch=-1)
    forecast = direct.forecast(p)
    # A correct inverse from the stored reference artifact, if present, against this forward model.
    artifact = Path("/tmp/split-damage.json")
    inverse = None
    if artifact.is_file():
        got = json.loads(artifact.read_text())["plates"][0]
        keys = ("xc", "yc", "a", "b", "theta_deg", "d")
        ph = np.array([z[k] for z in got["zones"] for k in keys] + list(got["fixture"].values()), float)
        pred = direct.forecast(ph)
        freq_err = max(
            abs(a / b - 1.0)
            for row, src in zip(pred, forecast)
            for a, b in zip(row["frequencies_Hz"], src["frequencies_Hz"])
        )
        inverse = {
            "artifact": str(artifact),
            "forecast_max_relative_vs_direct": float(freq_err),
            "snap_relative_vs_direct": abs(got["snap_load_N"] / snap_direct - 1.0),
        }
    evidence = {
        "snap_direct_N": snap_direct,
        "snap_reference_N": snap_reference,
        "snap_relative_difference": abs(snap_direct / snap_reference - 1.0),
        "frequency_comparisons": rows,
        "reference_inverse": inverse,
    }
    out = root / "authoring/evidence/independent-forward.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps(evidence, indent=2))
