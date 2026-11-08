"""Independent total-Lagrangian solver used only as a cross-check.

State is Cp^{-1}. Stress is the first Piola-Kirchhoff stress. The return map is
bisection on the logarithmic principal strains from numpy.linalg.eigh. F-bar
uses the average of the four Gauss-point Jacobians. The mesh is uniform in z
with a quadratic radial bias. This module does not import coupled_current.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from problem_spec import (
    CHI,
    E_MPA,
    EDOT0,
    L_HALF,
    NECK_Z,
    NU,
    R0,
    RHO_C,
    SHOULDER,
    T_REF,
    TM,
    Z_GAUGE,
)

_G = 1.0 / math.sqrt(3.0)
_PTS = ((-_G, -_G), (_G, -_G), (_G, _G), (-_G, _G))


def _shape(xi, eta):
    n = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta), (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
    dn = 0.25 * np.array([
        [-(1 - eta), -(1 - xi)],
        [(1 - eta), -(1 + xi)],
        [(1 + eta), (1 + xi)],
        [-(1 + eta), (1 - xi)],
    ])
    return n, dn


class RefCard:
    def __init__(self, K, eps0, n, alpha, Q, beta, C, m):
        self.K, self.eps0, self.n = float(K), float(eps0), float(n)
        self.alpha, self.Q, self.beta = float(alpha), float(Q), float(beta)
        self.C, self.m = float(C), float(m)
        self.mu = E_MPA / (2.0 * (1.0 + NU))
        self.kappa = E_MPA / (3.0 * (1.0 - 2.0 * NU))
        self.sy0 = self.K * self.eps0 ** self.n
        self.chi = CHI
        self.rho_c = RHO_C

    def y(self, ep, rate, T):
        swift = self.K * (self.eps0 + ep) ** self.n
        voce = self.sy0 + self.Q * (1.0 - np.exp(-self.beta * ep))
        base = self.alpha * swift + (1.0 - self.alpha) * voce
        fac = 1.0 + self.C * np.log(np.maximum(np.asarray(rate, float) / EDOT0, 1.0))
        theta = np.clip((np.asarray(T, float) - T_REF) / (TM - T_REF), 0.0, 1.0)
        return base * fac * (1.0 - theta ** self.m)


def mesh_uniform(dR, nr, nz, R0=R0, L=L_HALF):
    z = np.linspace(0.0, L, nz + 1)
    xy = np.zeros(((nr + 1) * (nz + 1), 2))
    for j, zz in enumerate(z):
        rad = R0 - dR * math.cos(math.pi * float(zz) / (2.0 * L)) ** 2
        for i in range(nr + 1):
            s = i / nr
            xy[j * (nr + 1) + i] = (rad * (0.85 * s + 0.15 * s * s), zz)
    conn = []
    w = nr + 1
    for j in range(nz):
        for i in range(nr):
            a = j * w + i
            conn.append((a, a + 1, a + w + 1, a + w))
    return xy, np.asarray(conn, dtype=int)


def _yield_T(card, ep, dg, T0, dt, adiabatic):
    rate = dg / dt if dt > 0 else np.zeros_like(dg)
    if not adiabatic:
        return card.y(ep + dg, rate, T0), T0
    temp = np.array(T0, copy=True)
    yv = card.y(ep + dg, rate, temp)
    for _ in range(8):
        yv = card.y(ep + dg, rate, temp)
        temp = T0 + card.chi * yv * dg / card.rho_c
    return yv, temp


def _bisect(q, ep, T, dt, card, adiabatic):
    dg = np.zeros_like(q)
    tout = np.array(T, copy=True)
    y0, _ = _yield_T(card, ep, np.zeros_like(q), T, dt, False)
    on = q > y0 * (1.0 + 1e-12)
    if not np.any(on):
        return dg, tout
    qa, epa, Ta = q[on], ep[on], T[on]
    lo, hi = np.zeros_like(qa), qa / (3.0 * card.mu)
    for _ in range(46):
        mid = 0.5 * (lo + hi)
        yv, _ = _yield_T(card, epa, mid, Ta, dt, adiabatic)
        f = qa - 3.0 * card.mu * mid - yv
        lo = np.where(f > 0, mid, lo)
        hi = np.where(f <= 0, mid, hi)
    dg[on] = 0.5 * (lo + hi)
    _, tn = _yield_T(card, epa, dg[on], Ta, dt, adiabatic)
    tout[on] = tn
    return dg, tout


class PK1Model:
    """Isothermal (or frozen-temperature) cross-check. Conduction is not repeated here;
    the thermal comparison uses coupled_current's adiabatic and k=0 modes, which are
    different algorithms inside that file, plus this mechanical cross-check.
    """

    def __init__(self, card, xy, conn, t_init, adiabatic=False):
        self.adiabatic = adiabatic
        self.card, self.xy, self.conn = card, xy, conn
        self.ne, self.nn = len(conn), len(xy)
        self.xe = xy[conn]
        self.gp = []
        for xi, eta in _PTS:
            n, dn = _shape(xi, eta)
            jm = np.einsum("ai,eaj->eij", dn, self.xe)
            det = jm[:, 0, 0] * jm[:, 1, 1] - jm[:, 0, 1] * jm[:, 1, 0]
            dndx = np.einsum("eij,aj->eai", np.linalg.inv(jm), dn)
            self.gp.append({"N": n, "dN": dndx, "R": np.einsum("a,ea->e", n, self.xe[:, :, 0]), "det": det})
        self.cp = np.broadcast_to(np.eye(3), (self.ne, 4, 3, 3)).copy()
        self.ep = np.zeros((self.ne, 4))
        self.u = np.zeros_like(xy)
        self.t_gp = np.full((self.ne, 4), float(t_init))
        zmax = float(xy[:, 1].max())
        self.end = np.where(np.abs(xy[:, 1] - zmax) < 1e-10)[0]
        self.axis = np.where(xy[:, 0] < 1e-12)[0]
        self.mid = np.where(xy[:, 1] < 1e-12)[0]
        fixed = np.zeros(2 * self.nn, dtype=bool)
        fixed[2 * self.axis] = True
        fixed[2 * self.mid + 1] = True
        fixed[2 * self.end + 1] = True
        self.free = ~fixed
        self.dofs = (2 * conn[:, :, None] + np.arange(2)).reshape(self.ne, 8)
        nr1 = int(np.sum(xy[:, 1] < 1e-12))
        self.surf = np.array([j * nr1 + nr1 - 1 for j in range(self.nn // nr1)])

    def element_force(self, ue, dt):
        ne = self.ne
        js, kin = [], []
        for geo in self.gp:
            ur = np.einsum("a,ea->e", geo["N"], ue[:, :, 0])
            fth = (geo["R"] + ur) / geo["R"]
            f2 = np.einsum("eai,eaj->eij", ue, geo["dN"]) + np.eye(2)
            j = (f2[:, 0, 0] * f2[:, 1, 1] - f2[:, 0, 1] * f2[:, 1, 0]) * fth
            js.append(j)
            kin.append((f2, fth, j, geo))
        jbar = np.mean(np.stack(js, 0), 0)
        fe = np.zeros((ne, 8))
        ep_new = np.zeros_like(self.ep)
        cp_new = np.zeros_like(self.cp)
        t_new = np.array(self.t_gp, copy=True)
        for g, (f2, fth, j, geo) in enumerate(kin):
            scale = np.cbrt(jbar / np.maximum(j, 1e-30))
            fb = np.zeros((ne, 3, 3))
            fb[:, 0, 0], fb[:, 0, 1] = f2[:, 0, 0] * scale, f2[:, 0, 1] * scale
            fb[:, 1, 0], fb[:, 1, 1] = f2[:, 1, 0] * scale, f2[:, 1, 1] * scale
            fb[:, 2, 2] = fth * scale
            be = np.einsum("eij,ejk,elk->eil", fb, self.cp[:, g], fb)
            be = 0.5 * (be + np.swapaxes(be, -1, -2))
            lam, vec = np.linalg.eigh(be)
            lam = np.maximum(lam, 1e-18)
            strain = 0.5 * np.log(lam)
            tr = strain.sum(-1)
            dev = strain - tr[:, None] / 3.0
            q = 2.0 * self.card.mu * np.sqrt(1.5 * np.sum(dev * dev, -1))
            dg, tnew = _bisect(q, self.ep[:, g], self.t_gp[:, g], dt, self.card, self.adiabatic)
            fac = np.ones_like(q)
            hit = q > 1e-12
            fac[hit] = 1.0 - 3.0 * self.card.mu * dg[hit] / q[hit]
            dev = dev * fac[:, None]
            ec = dev + tr[:, None] / 3.0
            tau_p = self.card.kappa * tr[:, None] + 2.0 * self.card.mu * dev
            tau = np.einsum("eij,ej,ekj->eik", vec, tau_p, vec)
            ftrue = np.zeros_like(fb)
            ftrue[:, 0, 0], ftrue[:, 0, 1] = f2[:, 0, 0], f2[:, 0, 1]
            ftrue[:, 1, 0], ftrue[:, 1, 1] = f2[:, 1, 0], f2[:, 1, 1]
            ftrue[:, 2, 2] = fth
            p = ((j / np.maximum(jbar, 1e-30))[:, None, None]) * tau @ np.linalg.inv(np.swapaxes(ftrue, -1, -2))
            w0 = 2.0 * math.pi * geo["R"] * geo["det"]
            dN = geo["dN"]
            fe[:, 0::2] += (
                np.einsum("e,ea->ea", w0 * p[:, 0, 0], dN[:, :, 0])
                + np.einsum("e,ea->ea", w0 * p[:, 0, 1], dN[:, :, 1])
                + np.einsum("e,a->ea", w0 * p[:, 2, 2] / geo["R"], geo["N"])
            )
            fe[:, 1::2] += (
                np.einsum("e,ea->ea", w0 * p[:, 1, 0], dN[:, :, 0])
                + np.einsum("e,ea->ea", w0 * p[:, 1, 1], dN[:, :, 1])
            )
            be_n = np.einsum("eij,ej,ekj->eik", vec, np.exp(2.0 * ec), vec)
            fbi = np.linalg.inv(fb)
            cp_new[:, g] = np.einsum("eij,ejk,elk->eil", fbi, be_n, fbi)
            ep_new[:, g] = self.ep[:, g] + dg
            t_new[:, g] = tnew
        return fe, ep_new, cp_new, t_new

    def residual(self, u, dt, need_k):
        ue = u[self.conn]
        fe, ep_new, cp_new, t_new = self.element_force(ue, dt)
        fint = np.zeros(2 * self.nn)
        np.add.at(fint, self.dofs.ravel(), fe.ravel())
        kg = None
        if need_k:
            h = 2.0e-6
            flat = ue.reshape(self.ne, 8)
            ke = np.zeros((self.ne, 8, 8))
            for k in range(8):
                pert = flat.copy()
                pert[:, k] += h
                fe2, _, _, _ = self.element_force(pert.reshape(self.ne, 4, 2), dt)
                ke[:, :, k] = (fe2 - fe) / h
            rows = np.repeat(self.dofs, 8, 1).ravel()
            cols = np.tile(self.dofs, (1, 8)).ravel()
            kg = sparse.csr_matrix((ke.ravel(), (rows, cols)), shape=(2 * self.nn, 2 * self.nn))
        return fint, kg, ep_new, cp_new, t_new

    def step(self, u_end, dt):
        u = self.u.copy()
        u[self.end, 1] = u_end
        u[self.axis, 0] = 0.0
        u[self.mid, 1] = 0.0
        for it in range(12):
            fint, kg, ep_new, cp_new, t_new = self.residual(u, dt, True)
            if not np.all(np.isfinite(fint)):
                return False
            scale = max(1.0, float(np.sum(np.abs(fint[2 * self.end + 1]))))
            rhs = -fint[self.free]
            nr = float(np.linalg.norm(rhs))
            if nr < 1e-5 * scale:
                self.u, self.ep, self.cp, self.t_gp = u, ep_new, cp_new, t_new
                return True
            du = spsolve(kg[self.free][:, self.free].tocsc(), rhs)
            if not np.all(np.isfinite(du)):
                return False
            step = np.zeros_like(fint)
            step[self.free] = du
            ok = False
            for _ls in range(6):
                trial = u.copy()
                trial.ravel()[:] += step
                trial[self.end, 1] = u_end
                trial[self.axis, 0] = 0.0
                trial[self.mid, 1] = 0.0
                f2, _, ep2, cp2, t2 = self.residual(trial, dt, False)
                if float(np.linalg.norm(f2[self.free])) < nr:
                    u, ep_new, cp_new, t_new = trial, ep2, cp2, t2
                    ok = True
                    break
                step *= 0.5
            if not ok:
                return False
        return False

    def observe(self):
        gauge = 2.0 * float(np.interp(Z_GAUGE, self.xy[self.surf, 1], self.u[self.surf, 1]))
        rmin = float(np.min(self.xy[self.surf, 0] + self.u[self.surf, 0]))
        fint, _, _, _, _ = self.residual(self.u, 1.0, False)
        neck_n = neck_d = 0.0
        for g, geo in enumerate(self.gp):
            z = self.xy[self.conn, 1] @ geo["N"]
            w0 = 2.0 * math.pi * geo["R"] * geo["det"]
            neck = z <= NECK_Z
            neck_n += float(np.sum(w0[neck] * self.t_gp[neck, g]))
            neck_d += float(np.sum(w0[neck]))
        return {
            "force_N": float(np.sum(fint[2 * self.end + 1])),
            "gauge_mm": gauge,
            "r_min_mm": rmin,
            "ep_max": float(self.ep.max()),
            "T_neck_K": neck_n / max(neck_d, 1e-30),
        }


def run_frozen(card, schedule, dR, t_init, nr, nz, nsub, adiabatic=False):
    xy, conn = mesh_uniform(dR, nr, nz)
    model = PK1Model(card, xy, conn, t_init, adiabatic=adiabatic)
    hist = []
    cursor = 0.0
    for kind, t0, dur, _a, _b in schedule.stages():
        n = nsub if kind == "ramp" else max(2, nsub // 2)
        if dur <= 0.0:
            continue
        for i in range(1, n + 1):
            target = t0 + dur * i / n
            while cursor < target - 1e-12:
                adv = target - cursor
                ok = False
                step = adv
                while step >= adv / 32.0 - 1e-15:
                    if model.step(schedule.u_end(cursor + step), step):
                        cursor += step
                        obs = model.observe()
                        obs["t_s"] = cursor
                        hist.append(obs)
                        ok = True
                        break
                    step *= 0.5
                if not ok:
                    raise RuntimeError(f"reference newton failed at t={cursor}")
    return hist
