"""Classical p-type silicon MOS capacitor: nonlinear Poisson + capacitance.

Units inside the solver: m, m^-3, V, C/m^2, F/m^2. Input doping is cm^-3,
oxide thickness nm, oxide charge elementary positive charges/cm^2.
No transport or frequency-dependent dynamics are solved.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq

Q = 1.602176634e-19
KB = 1.380649e-23
EPS0 = 8.8541878128e-12
EPS_SI = 11.7 * EPS0
EPS_OX = 3.9 * EPS0
ROOT = Path(__file__).resolve().parents[1]


def eg_ev(t: float) -> float:
    return 1.17 - 4.73e-4 * t*t / (t + 636.0)


@dataclass(frozen=True)
class Device:
    na_cm3: float = 1e16
    tox_nm: float = 10.0
    temperature: float = 300.0
    phi_ms_v: float = 0.0
    qox_cm2: float = 0.0

    def __post_init__(self) -> None:
        values = (self.na_cm3, self.tox_nm, self.temperature, self.phi_ms_v, self.qox_cm2)
        if not all(np.isfinite(values)):
            raise ValueError("Device parameters must be finite")
        if not 1e15 <= self.na_cm3 <= 1e17:
            raise ValueError("Teaching model: NA must be 1e15..1e17 cm^-3")
        if not 5 <= self.tox_nm <= 40 or not 250 <= self.temperature <= 400:
            raise ValueError("Teaching model: tox=5..40 nm, T=250..400 K")
        if abs(self.phi_ms_v) > 1 or abs(self.qox_cm2) > 5e11:
            raise ValueError("Use |phi_ms|<=1 V and |Qox/q|<=5e11 cm^-2")

    @property
    def na(self) -> float:
        return self.na_cm3 * 1e6

    @property
    def vt(self) -> float:
        return KB * self.temperature / Q

    @property
    def ni(self) -> float:
        t = self.temperature
        return float(1e16 * (t / 300)**1.5 * np.exp(
            eg_ev(300) / (2*KB/Q*300) - eg_ev(t) / (2*self.vt)))

    @property
    def p0(self) -> float:
        return (self.na + np.sqrt(self.na**2 + 4*self.ni**2)) / 2

    @property
    def n0(self) -> float:
        return self.ni**2 / self.p0

    @property
    def eta(self) -> float:
        return self.n0 / self.p0

    @property
    def ld(self) -> float:
        return float(np.sqrt(EPS_SI * self.vt / (Q*self.p0)))

    @property
    def cox(self) -> float:
        return EPS_OX / (self.tox_nm * 1e-9)

    @property
    def qox(self) -> float:
        return Q * self.qox_cm2 * 1e4

    @property
    def vfb(self) -> float:
        return self.phi_ms_v - self.qox / self.cox

    @property
    def phi_f(self) -> float:
        return float(self.vt * np.log(self.p0 / self.ni))

    @property
    def wmax(self) -> float:
        return float(np.sqrt(4 * EPS_SI * self.phi_f / (Q*self.na)))

    @property
    def threshold_depletion(self) -> float:
        return self.vfb + 2*self.phi_f + np.sqrt(4*EPS_SI*Q*self.na*self.phi_f)/self.cox

    def summary(self) -> dict:
        return {"NA_cm-3": self.na_cm3, "tox_nm": self.tox_nm,
                "temperature_K": self.temperature, "phi_ms_V": self.phi_ms_v,
                "Qox_over_q_cm-2": self.qox_cm2, "ni_cm-3": self.ni/1e6,
                "Cox_F_m2": self.cox, "Vfb_V": self.vfb,
                "phi_F_V": self.phi_f, "Debye_length_nm": self.ld*1e9,
                "depletion_threshold_V": self.threshold_depletion}


def charge_functions(d: Device, u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """rho/(q*p0) and its derivative with respect to u=psi/VT.

    expm1 preserves bulk neutrality at u=0. Large trial iterates are rejected;
    density clipping would silently change both the model and its Jacobian.
    """
    if not np.all(np.isfinite(u)) or np.max(np.abs(u)) > 150:
        raise ValueError("Potential outside safe classical exponential range")
    r = np.expm1(-u) - d.eta * np.expm1(u)
    derivative = -np.exp(-u) - d.eta * np.exp(u)
    return r, derivative


def make_grid(d: Device, nodes: int = 401, depth_factor: float = 1.0) -> np.ndarray:
    if not isinstance(nodes, int) or nodes < 41:
        raise ValueError("Use at least 41 integer nodes")
    if not np.isfinite(depth_factor) or depth_factor < .5:
        raise ValueError("depth_factor must be finite and >=0.5")
    length = max(14*d.ld, 3*d.wmax) * depth_factor
    t = np.linspace(0, 1, nodes)
    return length * np.expm1(4*t) / np.expm1(4.0)


def system(d: Device, y: np.ndarray, u: np.ndarray, vg: float,
           freeze_electrons: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Node-centered finite volumes on a stretched grid, Robin gate condition.

    Residual is scaled by eps*VT/LD. Interface half-cell charge is included,
    ensuring charge conservation with the same nodal integration weights.
    The far bulk is grounded. The returned Jacobian is in scipy band format.
    """
    h = np.diff(y)
    r, dr = charge_functions(d, u)
    if freeze_electrons:
        dr = -np.exp(-u)  # Only the AC electron response is frozen, not DC charge.
    volume = (h[:-1] + h[1:]) / 2
    beta = d.cox * d.ld / EPS_SI
    f = np.empty_like(u)
    f[0] = (u[0]-u[1])/h[0] + beta*(u[0]-(vg-d.vfb)/d.vt) - r[0]*h[0]/2
    f[1:-1] = (u[1:-1]-u[:-2])/h[:-1] + (u[1:-1]-u[2:])/h[1:] - r[1:-1]*volume
    f[-1] = u[-1]
    band = np.zeros((3, len(u)))
    band[1, 0] = 1/h[0] + beta - dr[0]*h[0]/2
    band[1, 1:-1] = 1/h[:-1] + 1/h[1:] - dr[1:-1]*volume
    band[1, -1] = 1
    band[0, 1:] = -1/h
    band[2, :-2] = -1/h[:-1]
    return f, band


def solve_bias(d: Device, vg: float, nodes: int = 401,
               initial: np.ndarray | None = None, depth_factor: float = 1.0,
               tolerance: float = 2e-10, max_iterations: int = 100) -> dict:
    if not np.isfinite(vg) or abs(vg-d.vfb) > 2:
        raise ValueError("Use finite |Vg-Vfb|<=2 V; inspect classical model limits")
    x = make_grid(d, nodes, depth_factor)
    y = x / d.ld
    u = np.zeros(nodes) if initial is None else np.array(initial, dtype=float, copy=True)
    if u.shape != x.shape or not np.all(np.isfinite(u)):
        raise ValueError("Initial potential must have one finite value per grid node")
    trace = []
    for iteration in range(max_iterations):
        f, band = system(d, y, u, vg)
        norm = float(np.max(np.abs(f)))
        if norm < tolerance:
            trace.append({"iteration": iteration, "residual_inf": norm, "damping": 0.0})
            break
        delta = solve_banded((1, 1), band, -f, check_finite=False)
        alpha = min(1.0, 3.0 / max(float(np.max(np.abs(delta))), 1e-30))
        for _ in range(35):
            trial = u + alpha*delta
            if np.max(np.abs(trial)) <= 150:
                trial_norm = np.max(np.abs(system(d, y, trial, vg)[0]))
                if trial_norm <= (1-1e-4*alpha)*norm:
                    break
            alpha /= 2
        else:
            raise RuntimeError(f"Newton line search failed at Vg={vg:g}, residual={norm:.3e}")
        trace.append({"iteration": iteration, "residual_inf": norm, "damping": alpha})
        u = trial
    else:
        raise RuntimeError(f"Newton did not converge at Vg={vg:g}")

    h = np.diff(x)
    weights = np.r_[h[0]/2, (h[:-1]+h[1:])/2, h[-1]/2]
    r, _ = charge_functions(d, u)
    rho = Q*d.p0*r
    psi = d.vt*u
    p, n = d.p0*np.exp(-u), d.n0*np.exp(u)
    qg = d.cox*(vg-d.phi_ms_v-psi[0])
    qs = float(weights @ rho)
    d_bulk = EPS_SI*(psi[-2]-psi[-1])/h[-1]
    balance = qg+d.qox+qs-d_bulk
    c_values = []
    for freeze in (False, True):
        jac = system(d, y, u, vg, freeze_electrons=freeze)[1]
        rhs = np.zeros(nodes)
        rhs[0] = d.cox*d.ld/EPS_SI/d.vt
        sensitivity = solve_banded((1, 1), jac, rhs, check_finite=False)
        c_values.append(float(d.cox*(1-d.vt*sensitivity[0])))
    classical_ratio = max(float(n.max())/(2.8e25*(d.temperature/300)**1.5),
                          float(p.max())/(1.04e25*(d.temperature/300)**1.5))
    return {"x": x, "u": u, "psi": psi, "rho": rho, "n": n, "p": p,
            "field": -np.gradient(psi, x, edge_order=2), "Vg_V": float(vg),
            "psi_s_V": float(psi[0]), "Qg_C_m2": float(qg), "Qs_C_m2": qs,
            "Cqs_F_m2": c_values[0], "Chf_frozen_F_m2": c_values[1],
            "charge_balance_C_m2": float(balance), "bulk_displacement_C_m2": float(d_bulk),
            "residual_inf": trace[-1]["residual_inf"], "iterations": len(trace)-1,
            "classical_density_to_DOS_max": classical_ratio, "trace": trace}


def sweep(d: Device, voltages: np.ndarray, nodes: int = 401,
          depth_factor: float = 1.0) -> list[dict]:
    voltages = np.asarray(voltages, dtype=float)
    if voltages.ndim != 1 or len(voltages) == 0:
        raise ValueError("Voltages must be a nonempty 1D array")
    results, guess = [], None
    for vg in voltages:
        result = solve_bias(d, float(vg), nodes, guess, depth_factor)
        results.append(result)
        guess = result["u"]  # Bias continuation, never an analytical potential.
    return results


def reference_charge(d: Device, psi_s: float) -> tuple[float, float]:
    """Independent semi-infinite first-integral reference: Qs and d(-Qs)/dpsi.

    This scalar analytical relation is never used to initialize Newton.
    """
    u = psi_s / d.vt
    if abs(u) < 1e-5:
        # Stable series of expm1(+/-u) -/+ u through sixth order.
        a = (1+d.eta)/2 + (-1+d.eta)*u/6 + (1+d.eta)*u*u/24
        a += (-1+d.eta)*u**3/120 + (1+d.eta)*u**4/720
        da = (-1+d.eta)/6 + (1+d.eta)*u/12 + (-1+d.eta)*u*u/40 + (1+d.eta)*u**3/180
        scale = np.sqrt(2*EPS_SI*Q*d.p0*d.vt)
        return float(-scale*u*np.sqrt(a)), float(scale/d.vt*(np.sqrt(a)+u*da/(2*np.sqrt(a))))
    h = np.expm1(-u)+u + d.eta*(np.expm1(u)-u)
    hp = -np.expm1(-u)+d.eta*np.expm1(u)
    scale = np.sqrt(2*EPS_SI*Q*d.p0*d.vt)
    qs = -np.sign(u)*scale*np.sqrt(h)
    cs = np.sign(u)*scale*hp/(2*d.vt*np.sqrt(h))
    return float(qs), float(cs)


def reference_bias(d: Device, vg: float) -> dict:
    if not np.isfinite(vg) or abs(vg-d.vfb) > 2:
        raise ValueError("Use finite |Vg-Vfb|<=2 V")
    def mismatch(psi):
        return d.vfb + psi - reference_charge(d, psi)[0]/d.cox - vg
    psi = brentq(mismatch, -2, 2, xtol=5e-14)
    qs, cs = reference_charge(d, psi)
    return {"Vg_V": float(vg), "psi_s_V": float(psi), "Qs_C_m2": qs,
            "Qg_C_m2": -qs-d.qox, "Cqs_F_m2": float(d.cox*cs/(d.cox+cs))}
