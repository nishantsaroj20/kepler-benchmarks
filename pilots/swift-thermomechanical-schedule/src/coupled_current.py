"""Coupled axisymmetric thermomechanics, current-configuration residual.

Independent implementation of multiplicative finite-strain J2 plasticity with
Hencky elasticity, Swift-Voce hardening and a multiplicative rate/temperature
factor. Heat conduction is a staggered backward-Euler solve on the reference
configuration. This file does not import the v1 task solvers.

State at each Gauss point: elastic Finger tensor b_e (pushed forward
incrementally), equivalent plastic strain, and the deformation gradient at the
end of the last accepted step. The return map is a vectorised safeguarded
Newton solve on the logarithmic principal strains. The element is bilinear Q4
with F-bar (centroid Jacobian). The tangent is a forward difference of the
internal force. Conduction uses the true deformation gradient; plastic heating
uses Kirchhoff plastic work per reference volume.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from problem_spec import (
    A_ACT,
    CHI,
    E_MPA,
    EDOT0,
    K_COND,
    L_HALF,
    NECK_Z,
    NU,
    R0,
    RHO_C,
    SHOULDER,
    T_GRIP,
    T_REF,
    TM,
    V_ACT,
    Z_GAUGE,
    Schedule,
)

_G = 1.0 / math.sqrt(3.0)
_PTS = np.array([[-_G, -_G], [_G, -_G], [_G, _G], [-_G, _G], [0.0, 0.0]])


def _shape(xi, eta):
    n = 0.25 * np.array([
        (1.0 - xi) * (1.0 - eta),
        (1.0 + xi) * (1.0 - eta),
        (1.0 + xi) * (1.0 + eta),
        (1.0 - xi) * (1.0 + eta),
    ])
    dn = 0.25 * np.array([
        [-(1.0 - eta), -(1.0 - xi)],
        [(1.0 - eta), -(1.0 + xi)],
        [(1.0 + eta), (1.0 + xi)],
        [-(1.0 + eta), (1.0 - xi)],
    ])
    return n, dn


_N = np.stack([_shape(*p)[0] for p in _PTS])       # (5, 4)
_DN = np.stack([_shape(*p)[1] for p in _PTS])      # (5, 4, 2) d/d(xi,eta)


class Card:
    """Eight-parameter Swift-Voce card. Known elastic and thermal constants."""

    def __init__(self, K, eps0, n, alpha, Q, beta, C, m,
                 E=E_MPA, nu=NU, edot0=EDOT0, chi=CHI, rho_c=RHO_C, T0=T_REF, Tm=TM):
        self.K = float(K)
        self.eps0 = float(eps0)
        self.n = float(n)
        self.alpha = float(alpha)
        self.Q = float(Q)
        self.beta = float(beta)
        self.C = float(C)
        self.m = float(m)
        self.E = float(E)
        self.nu = float(nu)
        self.edot0 = float(edot0)
        self.chi = float(chi)
        self.rho_c = float(rho_c)
        self.T0 = float(T0)
        self.Tm = float(Tm)
        self.mu = self.E / (2.0 * (1.0 + self.nu))
        self.lam = self.E * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        self.kappa = self.lam + 2.0 * self.mu / 3.0
        self.sy0 = self.K * self.eps0 ** self.n

    def sig_h(self, ep):
        a = self.alpha
        swift = self.K * (self.eps0 + ep) ** self.n
        voce = self.sy0 + self.Q * (1.0 - np.exp(-self.beta * ep))
        return a * swift + (1.0 - a) * voce

    def flow(self, ep, rate, T):
        rate_fac = 1.0 + self.C * np.log(np.maximum(np.asarray(rate, dtype=float) / self.edot0, 1.0))
        theta = np.clip((np.asarray(T, dtype=float) - self.T0) / (self.Tm - self.T0), 0.0, 1.0)
        return self.sig_h(ep) * rate_fac * (1.0 - theta ** self.m)

    def as_dict(self):
        return {k: getattr(self, k) for k in ("K", "eps0", "n", "alpha", "Q", "beta", "C", "m", "sy0")}


def _eig2(a, b, c):
    half = 0.5 * (a - b)
    disc = np.sqrt(half * half + c * c)
    mean = 0.5 * (a + b)
    th = 0.5 * np.arctan2(2.0 * c, a - b)
    return mean + disc, mean - disc, np.cos(th), np.sin(th)


def _return_map(q, ep_n, T, dt, card, adiabatic):
    """Solve q - 3 mu dg - Y(ep+dg, dg/dt, T) = 0. T is frozen unless adiabatic."""
    dg = np.zeros_like(q)
    T_out = np.array(T, dtype=float, copy=True)
    y_out = card.flow(ep_n, np.zeros_like(q), T)
    active = q > y_out * (1.0 + 1e-12)
    if not np.any(active):
        return dg, T_out, y_out
    qa = q[active]
    epa = ep_n[active]
    Ta = T[active]
    hi = qa / (3.0 * card.mu)
    lo = np.zeros_like(qa)
    slope = 3.0 * card.mu + 1.0e4
    d = np.clip((qa - y_out[active]) / slope, 0.0, 0.5 * hi)

    def y_of(dloc, tloc):
        rate = dloc / dt if dt > 0.0 else np.zeros_like(dloc)
        if not adiabatic:
            return card.flow(epa + dloc, rate, tloc), tloc
        temp = np.array(tloc, copy=True)
        yv = card.flow(epa + dloc, rate, temp)
        for _ in range(8):
            yv = card.flow(epa + dloc, rate, temp)
            temp = Ta + card.chi * yv * dloc / card.rho_c
        return yv, temp

    for _ in range(30):
        yv, tn = y_of(d, Ta)
        f = qa - 3.0 * card.mu * d - yv
        lo = np.where(f > 0.0, d, lo)
        hi = np.where(f <= 0.0, d, hi)
        if np.max(np.abs(f)) < 1e-8 * max(card.sy0, 1.0):
            T_out[active] = tn
            y_out[active] = yv
            dg[active] = d
            return dg, T_out, y_out
        eps = np.maximum(1e-10, 1e-6 * np.maximum(d, 1e-8))
        y2, _ = y_of(d + eps, Ta)
        df = -3.0 * card.mu - (y2 - yv) / eps
        step = np.where(np.abs(df) < 1e-14, 0.0, f / df)
        dn = d - step
        bad = (~np.isfinite(dn)) | (dn <= lo) | (dn >= hi)
        d = np.where(bad, 0.5 * (lo + hi), dn)
        d = np.clip(d, lo, hi)
    yv, tn = y_of(d, Ta)
    T_out[active] = tn
    y_out[active] = yv
    dg[active] = d
    return dg, T_out, y_out


def mesh_bar(R0, dR, L, nr, nz, z_exp=1.45):
    """Power-biased axial mesh, uniform in the normalised radius. Outer profile is the cosine imperfection."""
    s = np.linspace(0.0, 1.0, nz + 1)
    zs = L * s ** z_exp
    xy = np.zeros(((nz + 1) * (nr + 1), 2))
    for j, z in enumerate(zs):
        rad = R0 - dR * math.cos(math.pi * z / (2.0 * L)) ** 2
        for i in range(nr + 1):
            xy[j * (nr + 1) + i, 0] = rad * (i / nr)
            xy[j * (nr + 1) + i, 1] = z
    conn = []
    w = nr + 1
    for j in range(nz):
        for i in range(nr):
            a = j * w + i
            conn.append((a, a + 1, a + w + 1, a + w))
    return xy, np.asarray(conn, dtype=int), zs


def _ref_geom(xy, conn):
    xe = xy[conn]
    out = []
    for g in range(5):
        ng, dng = _N[g], _DN[g]
        jmat = np.einsum("ai,eaj->eij", dng, xe)
        det = jmat[:, 0, 0] * jmat[:, 1, 1] - jmat[:, 0, 1] * jmat[:, 1, 0]
        inv = np.linalg.inv(jmat)
        dndx = np.einsum("eij,aj->eai", inv, dng)
        rg = np.einsum("a,ea->e", ng, xe[:, :, 0])
        out.append({"N": ng, "dNdX": dndx, "R": rg, "detJ0": det})
    return xe, out


def _kinematics(xe, ue, ref):
    """True F at 4 Gauss points and at the centroid (index 4)."""
    packs = []
    for g in range(5):
        geo = ref[g]
        ur = np.einsum("a,ea->e", geo["N"], ue[:, :, 0])
        r = geo["R"] + ur
        fth = r / geo["R"]
        h = np.einsum("eai,eaj->eij", ue, geo["dNdX"])
        f2 = h + np.eye(2)
        det2 = f2[:, 0, 0] * f2[:, 1, 1] - f2[:, 0, 1] * f2[:, 1, 0]
        packs.append({"F2": f2, "Fth": fth, "J": det2 * fth, "r": r, "dNdX": geo["dNdX"], "N": geo["N"], "R": geo["R"], "detJ0": geo["detJ0"]})
    return packs


def _constitutive(packs, be_n, fn2, fnth, ep_n, t_gp, dt, card, adiabatic):
    """Return Cauchy stress at the 4 Gauss points and candidate state. Centroid sets Jbar."""
    jbar = packs[4]["J"]
    ne = be_n.shape[0]
    sig = np.zeros((ne, 4, 4))
    be_new = np.zeros_like(be_n)
    fn2_new = np.zeros_like(fn2)
    fnth_new = np.zeros_like(fnth)
    ep_new = np.zeros_like(ep_n)
    t_new = np.zeros_like(t_gp)
    q_gp = np.zeros_like(ep_n)
    for g in range(4):
        pk = packs[g]
        j = pk["J"]
        # F-bar: scale toward the centroid Jacobian. a^3 * J = Jbar.
        scale = np.cbrt(jbar / np.maximum(j, 1e-30))
        f2b = pk["F2"] * scale[:, None, None]
        fthb = pk["Fth"] * scale
        jb = jbar
        f2n_inv = np.linalg.inv(fn2[:, g])
        frel = np.einsum("eij,ejk->eik", f2b, f2n_inv)
        b = np.zeros((ne, 2, 2))
        b[:, 0, 0] = be_n[:, g, 0]
        b[:, 1, 1] = be_n[:, g, 1]
        b[:, 0, 1] = b[:, 1, 0] = be_n[:, g, 2]
        bt = np.einsum("eij,ejk,elk->eil", frel, b, frel)
        bth = (fthb / fnth[:, g]) ** 2 * be_n[:, g, 3]
        l1, l2, c, s = _eig2(bt[:, 0, 0], bt[:, 1, 1], bt[:, 0, 1])
        l1 = np.maximum(l1, 1e-18)
        l2 = np.maximum(l2, 1e-18)
        bth = np.maximum(bth, 1e-18)
        e = np.stack((0.5 * np.log(l1), 0.5 * np.log(l2), 0.5 * np.log(bth)), axis=-1)
        tr = e.sum(axis=-1)
        dev = e - tr[:, None] / 3.0
        sd = 2.0 * card.mu * dev
        q = np.sqrt(1.5 * np.sum(sd * sd, axis=-1))
        dg, tout, _ = _return_map(q, ep_n[:, g], t_gp[:, g], dt, card, adiabatic)
        fac = np.ones_like(q)
        hit = q > 0.0
        fac[hit] = 1.0 - 3.0 * card.mu * dg[hit] / q[hit]
        dev = dev * fac[:, None]
        e_c = dev + tr[:, None] / 3.0
        tau_p = card.kappa * tr[:, None] + 2.0 * card.mu * dev
        t1, t2, t3 = tau_p[:, 0], tau_p[:, 1], tau_p[:, 2]
        trr = c * c * t1 + s * s * t2
        tzz = s * s * t1 + c * c * t2
        trz = c * s * (t1 - t2)
        sig[g] if False else None
        sig[:, g, 0] = trr / jb
        sig[:, g, 1] = tzz / jb
        sig[:, g, 2] = trz / jb
        sig[:, g, 3] = t3 / jb
        x1 = np.exp(2.0 * e_c[:, 0])
        x2 = np.exp(2.0 * e_c[:, 1])
        be_new[:, g, 0] = c * c * x1 + s * s * x2
        be_new[:, g, 1] = s * s * x1 + c * c * x2
        be_new[:, g, 2] = c * s * (x1 - x2)
        be_new[:, g, 3] = np.exp(2.0 * e_c[:, 2])
        fn2_new[:, g] = f2b
        fnth_new[:, g] = fthb
        ep_new[:, g] = ep_n[:, g] + dg
        t_new[:, g] = tout
        q_gp[:, g] = np.sqrt(1.5 * np.sum((2.0 * card.mu * dev) * (2.0 * card.mu * dev), axis=-1))
        # q at the end of the step equals Y for active points; recompute from the corrected stress
        # so heating uses the Kirchhoff von Mises actually carried by the stress.
        q_end = np.sqrt(1.5 * (trr - (trr + tzz + t3) / 3.0) ** 2 * 0.0)  # placeholder replaced below
        dev_tau_norm = np.sqrt(
            (trr - (trr + tzz + t3) / 3.0) ** 2
            + (tzz - (trr + tzz + t3) / 3.0) ** 2
            + (t3 - (trr + tzz + t3) / 3.0) ** 2
            + 2.0 * trz ** 2
        )
        q_gp[:, g] = np.sqrt(1.5) * dev_tau_norm
    return sig, be_new, fn2_new, fnth_new, ep_new, t_new, q_gp, jbar


def _internal_force(xe, ue, ref, state, t_gp, dt, card, adiabatic):
    packs = _kinematics(xe, ue, ref)
    sig, be, fn2, fnth, ep, tout, qgp, jbar = _constitutive(
        packs, state["be"], state["Fn2"], state["Fnth"], state["ep"], t_gp, dt, card, adiabatic
    )
    ne = xe.shape[0]
    fe = np.zeros((ne, 4, 2))
    for g in range(4):
        pk = packs[g]
        # Current volume weight uses the true Jacobian. Stress is tau/Jbar.
        w = 2.0 * math.pi * pk["R"] * pk["detJ0"] * pk["J"]
        invf = np.linalg.inv(pk["F2"])
        dndx = np.einsum("eai,eij->eaj", pk["dNdX"], invf)
        srr, szz, srz, sth = sig[:, g, 0], sig[:, g, 1], sig[:, g, 2], sig[:, g, 3]
        ng = pk["N"]
        r = pk["r"]
        fe[:, :, 0] += (
            np.einsum("e,ea->ea", w * srr, dndx[:, :, 0])
            + np.einsum("e,ea->ea", w * srz, dndx[:, :, 1])
            + np.einsum("e,a->ea", w * sth / r, ng)
        )
        fe[:, :, 1] += (
            np.einsum("e,ea->ea", w * srz, dndx[:, :, 0])
            + np.einsum("e,ea->ea", w * szz, dndx[:, :, 1])
        )
    aux = {
        "be": be, "Fn2": fn2, "Fnth": fnth, "ep": ep, "T": tout, "q": qgp,
        "packs": packs, "sig": sig, "jbar": jbar,
    }
    return fe.reshape(ne, 8), aux


def _heat_step(xy, conn, ref, aux, t_nodes, dt, card, k_cond, t_grip, end_nodes):
    """Backward Euler reference conduction. Returns new nodal T and balance terms."""
    ne = conn.shape[0]
    nn = xy.shape[0]
    ke_all = np.zeros((ne, 4, 4))
    me_all = np.zeros((ne, 4, 4))
    se_all = np.zeros((ne, 4))
    plastic = 0.0
    rate = 1.0 / dt if dt > 0.0 else 0.0
    for g in range(4):
        geo = ref[g]
        pk = aux["packs"][g]
        f2 = pk["F2"]
        # In-plane C = F2^T F2. The hoop temperature gradient is zero.
        c2 = np.einsum("eji,ejk->eik", f2, f2)
        invc = np.linalg.inv(c2)
        kappa = (k_cond * pk["J"])[:, None, None] * invc
        dndx = geo["dNdX"]
        w0 = 2.0 * math.pi * geo["R"] * geo["detJ0"]
        dn_k = np.einsum("eai,eij->eaj", dndx, kappa)
        ke_all += np.einsum("eai,ebi,e->eab", dn_k, dndx, w0)
        ng = geo["N"]
        me_all += card.rho_c * np.einsum("a,b,e->eab", ng, ng, w0)
        dg = aux["dg"][:, g]
        q = aux["q"][:, g]
        se_all += np.einsum("a,e->ea", ng, w0 * card.chi * q * dg * rate)
        plastic += float(np.sum(w0 * card.chi * q * dg))
    nodes = conn
    rr = np.repeat(nodes, 4, axis=1).ravel()
    cc = np.tile(nodes, (1, 4)).ravel()
    kmat = sparse.coo_matrix((ke_all.ravel(), (rr, cc)), shape=(nn, nn)).tocsr()
    mcons = sparse.coo_matrix((me_all.ravel(), (rr, cc)), shape=(nn, nn)).tocsr()
    # Row-sum lumping. A consistent capacity plus a large plastic source oscillates on coarse steps.
    lump = np.asarray(mcons.sum(axis=1)).ravel()
    mmat = sparse.diags(lump, format="csr")
    source = np.zeros(nn)
    np.add.at(source, nodes.ravel(), se_all.ravel())
    a = (mmat / dt + kmat) if dt > 0.0 else kmat
    rhs = (mmat / dt) @ t_nodes + source if dt > 0.0 else source.copy()
    free = np.ones(nn, dtype=bool)
    free[end_nodes] = False
    # Dirichlet contribution
    if np.any(~free):
        rhs = rhs - a[:, ~free] @ np.full((~free).sum(), t_grip)
    t_new = np.full(nn, t_grip)
    if np.any(free):
        t_free = spsolve(a[free][:, free].tocsc(), rhs[free])
        if not np.all(np.isfinite(t_free)):
            raise RuntimeError("heat solve produced non-finite temperatures")
        t_new[free] = t_free
    # Energy: boundary power into the body over the step equals stored - plastic,
    # with stored computed from the solved nodal field.
    t_gp = np.zeros_like(aux["T"])
    stored_quad = 0.0
    for g in range(4):
        ng = ref[g]["N"]
        w0 = 2.0 * math.pi * ref[g]["R"] * ref[g]["detJ0"]
        tgp = t_new[conn] @ ng
        t_gp[:, g] = tgp
        t_old = t_nodes[conn] @ ng
        stored_quad += float(np.sum(w0 * card.rho_c * (tgp - t_old)))
    full_res = a @ t_new - ((mmat / dt) @ t_nodes + source if dt > 0.0 else source)
    # Lumped capacity change. With a symmetric conductivity matrix of zero row sums,
    # lump·(T-Tn) - plastic - boundary_in = 0 when the free residual vanishes.
    stored_lump = float(lump @ (t_new - t_nodes))
    boundary_in = float(np.sum(full_res[~free])) * (dt if dt > 0.0 else 0.0)
    balance = {
        "stored_Nmm": stored_lump,
        "stored_quad_Nmm": stored_quad,
        "plastic_Nmm": plastic,
        "boundary_out_Nmm": -boundary_in,
        "imbalance_Nmm": stored_lump - plastic - boundary_in,
    }
    return t_new, t_gp, balance


class CoupledSolver:
    def __init__(self, card, xy, conn, t_init, k_cond=K_COND):
        self.card = card
        self.xy = xy
        self.conn = conn
        self.k_cond = float(k_cond)
        self.ne = conn.shape[0]
        self.nn = xy.shape[0]
        self.xe, self.ref = _ref_geom(xy, conn)
        self.u = np.zeros_like(xy)
        self.be = np.zeros((self.ne, 4, 4))
        self.be[..., 0] = 1.0
        self.be[..., 1] = 1.0
        self.be[..., 3] = 1.0
        self.fn2 = np.broadcast_to(np.eye(2), (self.ne, 4, 2, 2)).copy()
        self.fnth = np.ones((self.ne, 4))
        self.ep = np.zeros((self.ne, 4))
        self.t_nodes = np.full(self.nn, float(t_init))
        self.tgp = np.full((self.ne, 4), float(t_init))
        zmax = xy[:, 1].max()
        self.end = np.where(np.abs(xy[:, 1] - zmax) < 1e-10)[0]
        self.axis = np.where(xy[:, 0] < 1e-12)[0]
        self.mid = np.where(xy[:, 1] < 1e-12)[0]
        fixed = np.zeros(2 * self.nn, dtype=bool)
        fixed[2 * self.axis] = True
        fixed[2 * self.mid + 1] = True
        fixed[2 * self.end + 1] = True
        self.fixed = fixed
        self.free = ~fixed
        self.dofs = (2 * conn[:, :, None] + np.arange(2)).reshape(self.ne, 8)
        # surface nodes, one per axial station, for gauge and radius
        w = int(round(math.sqrt(self.nn)))  # not used
        nr1 = int(np.sum(np.abs(xy[:, 1]) < 1e-12))
        self.nr1 = nr1
        nz1 = self.nn // nr1
        self.surf = np.array([j * nr1 + (nr1 - 1) for j in range(nz1)])

    def _state(self):
        return {"be": self.be, "Fn2": self.fn2, "Fnth": self.fnth, "ep": self.ep}

    def _tgp_from_nodes(self, t_nodes):
        out = np.zeros((self.ne, 4))
        for g in range(4):
            out[:, g] = t_nodes[self.conn] @ self.ref[g]["N"]
        return out

    def _assemble(self, u, t_gp, dt, adiabatic, need_k):
        ue = u[self.conn]
        fe, aux = _internal_force(self.xe, ue, self.ref, self._state(), t_gp, dt, self.card, adiabatic)
        fint = np.zeros(2 * self.nn)
        np.add.at(fint, self.dofs.ravel(), fe.ravel())
        if not need_k:
            return fint, None, aux
        ke = np.zeros((self.ne, 8, 8))
        h = 1e-6
        flat = ue.reshape(self.ne, 8)
        for k in range(8):
            pert = flat.copy()
            pert[:, k] += h
            f2, _ = _internal_force(self.xe, pert.reshape(self.ne, 4, 2), self.ref, self._state(), t_gp, dt, self.card, adiabatic)
            ke[:, :, k] = (f2 - fe) / h
        rows = np.repeat(self.dofs, 8, axis=1).ravel()
        cols = np.tile(self.dofs, (1, 8)).ravel()
        kglob = sparse.csr_matrix((ke.ravel(), (rows, cols)), shape=(2 * self.nn, 2 * self.nn))
        return fint, kglob, aux

    def _newton(self, u, u_target_end, t_gp, dt, adiabatic):
        u = u.copy()
        u[self.end, 1] = u_target_end
        u[self.axis, 0] = 0.0
        u[self.mid, 1] = 0.0
        last = None
        for it in range(16):
            fint, kglob, aux = self._assemble(u, t_gp, dt, adiabatic, True)
            if not np.all(np.isfinite(fint)):
                return None, aux, it, np.inf
            rhs = -fint[self.free]
            scale = max(1.0, float(np.sum(np.abs(fint[2 * self.end + 1]))))
            nr = float(np.linalg.norm(rhs))
            last = (fint, aux, nr, scale)
            if nr < 1e-6 * scale:
                return u, aux, it, nr
            kff = kglob[self.free][:, self.free].tocsc()
            try:
                du = spsolve(kff, rhs)
            except Exception:
                return None, aux, it, nr
            if not np.all(np.isfinite(du)):
                return None, aux, it, nr
            step = np.zeros_like(fint)
            step[self.free] = du
            # Reject the iteration if even a short step does not decrease the residual.
            accepted = False
            trial = u
            nr2 = nr
            for _ls in range(8):
                trial = u.copy()
                trial.ravel()[:] += step
                trial[self.end, 1] = u_target_end
                trial[self.axis, 0] = 0.0
                trial[self.mid, 1] = 0.0
                f2, _, aux2 = self._assemble(trial, t_gp, dt, adiabatic, False)
                if not np.all(np.isfinite(f2)):
                    step *= 0.5
                    continue
                nr2 = float(np.linalg.norm(f2[self.free]))
                if nr2 < nr * (1.0 - 1e-4) or nr2 < 1e-6 * scale:
                    u, aux, nr = trial, aux2, nr2
                    accepted = True
                    break
                step *= 0.5
            if not accepted:
                return None, aux, it, nr
            if nr < 1e-6 * scale:
                return u, aux, it, nr
        return None, last[1] if last else None, 16, last[2] if last else np.inf

    def snapshot(self):
        """Copy every committed field. A failed attempt must restore this exactly."""
        return {
            "u": self.u.copy(),
            "be": self.be.copy(),
            "fn2": self.fn2.copy(),
            "fnth": self.fnth.copy(),
            "ep": self.ep.copy(),
            "t_nodes": self.t_nodes.copy(),
            "tgp": self.tgp.copy(),
        }

    def restore(self, snap):
        self.u = snap["u"].copy()
        self.be = snap["be"].copy()
        self.fn2 = snap["fn2"].copy()
        self.fnth = snap["fnth"].copy()
        self.ep = snap["ep"].copy()
        self.t_nodes = snap["t_nodes"].copy()
        self.tgp = snap["tgp"].copy()

    def _observe(self, u, t_gp, ep):
        # Gauge: material point initially at z = Z_GAUGE on the outer surface.
        zs = self.xy[self.surf, 1]
        uz = u[self.surf, 1]
        gauge = 2.0 * float(np.interp(Z_GAUGE, zs, uz))
        r_now = self.xy[self.surf, 0] + u[self.surf, 0]
        r_min = float(np.min(r_now))
        force = None  # filled by caller from reaction
        neck_num = neck_den = sh_num = sh_den = 0.0
        ep_neck_num = 0.0
        ep_max = float(np.max(ep))
        for g in range(4):
            z = self.xy[self.conn, 1] @ self.ref[g]["N"]
            w0 = 2.0 * math.pi * self.ref[g]["R"] * self.ref[g]["detJ0"]
            tn = t_gp[:, g]
            neck = z <= NECK_Z
            sh = (z >= SHOULDER[0]) & (z <= SHOULDER[1])
            neck_num += float(np.sum(w0[neck] * tn[neck]))
            neck_den += float(np.sum(w0[neck]))
            sh_num += float(np.sum(w0[sh] * tn[sh]))
            sh_den += float(np.sum(w0[sh]))
            ep_neck_num += float(np.sum(w0[neck] * ep[:, g][neck]))
        if neck_den <= 0.0 or sh_den <= 0.0:
            raise RuntimeError(f"monitoring region empty neck={neck_den} shoulder={sh_den}")
        return {
            "gauge_mm": gauge,
            "r_min_mm": r_min,
            "T_neck_K": neck_num / neck_den,
            "T_shoulder_K": sh_num / sh_den,
            "ep_max": ep_max,
            "ep_neck": ep_neck_num / neck_den,
        }

    def _commit_step(self, u_new, aux, t_nodes, t_gp, adiabatic):
        self.u = np.array(u_new, copy=True)
        self.be = np.array(aux["be"], copy=True)
        self.fn2 = np.array(aux["Fn2"], copy=True)
        self.fnth = np.array(aux["Fnth"], copy=True)
        self.ep = np.array(aux["ep"], copy=True)
        if adiabatic:
            self.tgp = np.array(aux["T"], copy=True)
        else:
            self.t_nodes = np.array(t_nodes, copy=True)
            self.tgp = np.array(t_gp, copy=True)

    def _attempt_step(self, u_guess, u_end, dt, mode, stagger_tol, max_stagger):
        """Return a commit payload only when mechanics and, if coupled, heat both converge.

        This function does not write committed state. ``k`` enters only through ``self.k_cond``
        inside the same heat operator used for every positive conductivity, including the
        discrete k=0 limit. End-face temperature remains the Dirichlet value T_grip for every k.
        """
        adiabatic = mode == "adiabatic"
        if mode == "adiabatic":
            t_gp = self.tgp.copy()
        else:
            t_gp = self._tgp_from_nodes(self.t_nodes)
        if mode == "isothermal":
            t_gp[:] = float(self.t_nodes[0])
        t_gp_start = t_gp.copy()
        u_new, aux, it, nr = self._newton(u_guess, u_end, t_gp, dt, adiabatic)
        if u_new is None or not np.isfinite(nr):
            return None
        balance = {"stored_Nmm": 0.0, "plastic_Nmm": 0.0, "boundary_out_Nmm": 0.0, "imbalance_Nmm": 0.0}
        ep_n = self.ep
        stagger_passes = 0
        stagger_delta = 0.0
        t_nodes = self.t_nodes
        t_mech = t_gp
        if mode == "coupled":
            coupled_ok = False
            for _s in range(max_stagger):
                aux["dg"] = np.maximum(aux["ep"] - ep_n, 0.0)
                t_nodes_new, t_gp_new, balance = _heat_step(
                    self.xy, self.conn, self.ref, aux, self.t_nodes, dt, self.card, self.k_cond, T_GRIP, self.end
                )
                stagger_delta = float(np.max(np.abs(t_gp_new - t_gp)))
                stagger_passes = _s + 1
                if stagger_delta <= stagger_tol:
                    t_mech = np.array(t_gp, copy=True)
                    t_nodes, t_gp = t_nodes_new, t_gp_new
                    coupled_ok = True
                    break
                u2, aux2, it2, nr2 = self._newton(u_new, u_end, t_gp_new, dt, False)
                if u2 is None or not np.isfinite(nr2):
                    return None
                u_new, aux, it, nr = u2, aux2, it2, nr2
                t_gp = t_gp_new
                t_nodes = t_nodes_new
            if not coupled_ok:
                return None
        elif mode == "adiabatic":
            dg = np.maximum(aux["ep"] - ep_n, 0.0)
            plastic = 0.0
            stored = 0.0
            for g in range(4):
                w0 = 2.0 * math.pi * self.ref[g]["R"] * self.ref[g]["detJ0"]
                plastic += float(np.sum(w0 * self.card.chi * aux["q"][:, g] * dg[:, g]))
                stored += float(np.sum(w0 * self.card.rho_c * (aux["T"][:, g] - t_gp_start[:, g])))
            balance = {"stored_Nmm": stored, "plastic_Nmm": plastic, "boundary_out_Nmm": 0.0, "imbalance_Nmm": stored - plastic}
            t_gp = aux["T"]
        fint, _, _ = self._assemble(np.asarray(u_new), t_mech, dt, adiabatic, False)
        force = float(np.sum(fint[2 * self.end + 1]))
        return {
            "u": u_new, "aux": aux, "t_nodes": t_nodes, "t_gp": t_gp, "it": int(it), "nr": float(nr),
            "balance": balance, "stagger_passes": stagger_passes, "stagger_delta_K": stagger_delta,
            "adiabatic": adiabatic, "t_mech": t_mech, "force_N": force,
        }

    def run(self, schedule: Schedule, nsub=8, mode="coupled", max_stagger=4, stagger_tol=0.25, verbose=False):
        """mode: coupled | adiabatic | isothermal.

        A step is committed only after the Newton residual meets its tolerance and, in coupled
        mode, the staggered temperature change meets ``stagger_tol``. Hitting ``max_stagger``
        without that tolerance rejects the step. Rejected attempts restore the previous
        displacement, plastic strain, elastic Finger tensor, deformation gradient and nodal temperature.
        """
        hist = []
        extrema = {
            "force_N": 0.0, "force_t_s": 0.0,
            "T_neck_K": -1e9, "T_neck_t_s": 0.0,
            "T_shoulder_K": -1e9, "T_shoulder_t_s": 0.0,
            "r_min_mm": 1e9, "r_min_t_s": 0.0,
            "ep_max": 0.0, "ep_neck": 0.0,
        }
        heat_stored = heat_plastic = heat_boundary = 0.0
        u_prev = self.u.copy()
        t_cursor = 0.0
        # t = 0 is the unloaded, uniformly hot (or cold) solid, before the grip Dirichlet is applied.
        obs0 = self._observe(self.u, self.tgp, self.ep)
        obs0.update({
            "t_s": 0.0, "u_end_mm": 0.0, "force_N": 0.0, "newton_iter": 0, "residual_N": 0.0,
            "dt_s": 0.0, "stagger_passes": 0, "stagger_delta_K": 0.0, "coupling_converged": True,
            "heat_imbalance_step_Nmm": 0.0, "accepted": True,
        })
        hist.append(obs0)
        self._track_extrema(extrema, obs0)
        if mode == "coupled":
            self.t_nodes[self.end] = T_GRIP
        plan = []
        for kind, t0, dur, _u0, _delta in schedule.stages():
            n = max(int(nsub if kind == "ramp" else max(4, nsub // 2)), 1)
            if dur <= 0.0:
                continue
            for i in range(1, n + 1):
                plan.append(t0 + dur * i / n)
        failures = []
        for t_target in plan:
            remaining = t_target - t_cursor
            guard = 0
            while remaining > 1e-12:
                guard += 1
                if guard > 12:
                    failures.append({"t": t_cursor, "reason": "step halved out"})
                    raise RuntimeError(f"time step collapsed at t={t_cursor}")
                dt = remaining
                ok = False
                adv = dt
                while adv >= dt / 64.0 - 1e-15:
                    snap = self.snapshot()
                    u_prev_snap = u_prev.copy()
                    u_end = schedule.u_end(t_cursor + adv)
                    dend = u_end - float(np.mean(self.u[self.end, 1]))
                    prev = float(np.mean(self.u[self.end, 1] - u_prev[self.end, 1]))
                    if abs(prev) > 1e-8 and abs(dend) <= 1.5 * abs(prev) + 0.05:
                        u_guess = self.u + (self.u - u_prev) * (dend / prev)
                    else:
                        u_guess = self.u.copy()
                        u_guess[:, 1] += (self.xy[:, 1] / max(float(self.xy[:, 1].max()), 1e-9)) * dend
                    payload = self._attempt_step(u_guess, u_end, adv, mode, stagger_tol, max_stagger)
                    if payload is None:
                        self.restore(snap)
                        u_prev = u_prev_snap
                        adv *= 0.5
                        failures.append({"t": t_cursor, "reason": "rejected_unconverged", "dt_next": adv})
                        continue
                    self._commit_step(payload["u"], payload["aux"], payload["t_nodes"], payload["t_gp"], payload["adiabatic"])
                    # The Newton residual is the equilibrium residual at the start-of-step state.
                    # Reassembling after the Finger tensor is updated would apply the step twice.
                    force = float(payload["force_N"])
                    obs = self._observe(self.u, self.tgp, self.ep)
                    obs.update({
                        "t_s": t_cursor + adv,
                        "u_end_mm": u_end,
                        "force_N": force,
                        "newton_iter": payload["it"],
                        "residual_N": payload["nr"],
                        "post_residual_N": payload["nr"],
                        "dt_s": adv,
                        "stagger_passes": payload["stagger_passes"],
                        "stagger_delta_K": payload["stagger_delta_K"],
                        "coupling_converged": True,
                        "heat_imbalance_step_Nmm": payload["balance"]["imbalance_Nmm"],
                        "accepted": True,
                    })
                    heat_stored += payload["balance"]["stored_Nmm"]
                    heat_plastic += payload["balance"]["plastic_Nmm"]
                    heat_boundary += payload["balance"]["boundary_out_Nmm"]
                    hist.append(obs)
                    self._track_extrema(extrema, obs)
                    u_prev = snap["u"]
                    t_cursor += adv
                    remaining = t_target - t_cursor
                    ok = True
                    break
                if not ok:
                    raise RuntimeError(
                        f"newton failed at t={t_cursor:.4f} u_end={schedule.u_end(t_cursor):.4f} "
                        f"T=[{self.t_nodes.min():.1f},{self.t_nodes.max():.1f}] ep={self.ep.max():.3f}"
                    )
                if verbose:
                    print(
                        f"t={t_cursor:.4f} F={force:.1f} g={obs['gauge_mm']:.4f} Tn={obs['T_neck_K']:.2f} "
                        f"r={obs['r_min_mm']:.4f} dT={obs['stagger_delta_K']:.3f}",
                        flush=True,
                    )
        return {
            "history": hist,
            "extrema": extrema,
            "final": hist[-1],
            "heat": {
                "stored_Nmm": heat_stored,
                "plastic_Nmm": heat_plastic,
                "boundary_out_Nmm": heat_boundary,
                "imbalance_Nmm": heat_stored - heat_plastic + heat_boundary,
                "balance_is": "internal discrete residual of the same heat operator, not an external flux audit",
            },
            "failures": failures,
            "ep_field_max": float(self.ep.max()),
            "settings": {
                "mode": mode, "nsub": nsub, "max_stagger": max_stagger, "stagger_tol_K": stagger_tol,
                "k_cond": self.k_cond, "newton_residual_factor": 1e-6, "newton_iteration_cap": 16,
                "end_dirichlet": "T_grip on Z=L for every k, including k=0",
            },
        }

    @staticmethod
    def _track_extrema(extrema, obs):
        if obs["force_N"] >= extrema["force_N"]:
            extrema["force_N"] = obs["force_N"]
            extrema["force_t_s"] = obs["t_s"]
        if obs["T_neck_K"] >= extrema["T_neck_K"]:
            extrema["T_neck_K"] = obs["T_neck_K"]
            extrema["T_neck_t_s"] = obs["t_s"]
        if obs["T_shoulder_K"] >= extrema["T_shoulder_K"]:
            extrema["T_shoulder_K"] = obs["T_shoulder_K"]
            extrema["T_shoulder_t_s"] = obs["t_s"]
        if obs["r_min_mm"] <= extrema["r_min_mm"]:
            extrema["r_min_mm"] = obs["r_min_mm"]
            extrema["r_min_t_s"] = obs["t_s"]
        extrema["ep_max"] = max(extrema["ep_max"], obs["ep_max"])
        extrema["ep_neck"] = max(extrema["ep_neck"], obs["ep_neck"])


def make_solver(card, dR, t_init, nr, nz, k_cond=K_COND, L=L_HALF, R0=R0):
    xy, conn, _ = mesh_bar(R0, dR, L, nr, nz)
    return CoupledSolver(card, xy, conn, t_init, k_cond=k_cond)


def quick_card(**kwargs):
    base = dict(K=680.0, eps0=0.015, n=0.18, alpha=0.5, Q=220.0, beta=14.0, C=0.012, m=0.8)
    base.update(kwargs)
    return Card(**base)
