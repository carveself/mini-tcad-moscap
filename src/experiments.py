"""Reproducible plots, physical checks, mesh studies and parameter sweeps."""
from __future__ import annotations

import csv
from dataclasses import replace
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from moscap import Device, EPS_SI, eg_ev, reference_bias, solve_bias, sweep

COLORS = ("#146c94", "#e07a24", "#6d4c91")
SCALARS = ("Vg_V", "psi_s_V", "Qg_C_m2", "Qs_C_m2", "Cqs_F_m2",
           "Chf_frozen_F_m2", "iterations", "residual_inf", "charge_balance_C_m2",
           "bulk_displacement_C_m2", "classical_density_to_DOS_max")


def savefig(fig, out: Path, name: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out/name, dpi=160, bbox_inches="tight")
    plt.close(fig)


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def generate(d: Device, out: Path, nodes: int = 401) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    threshold_relative = d.threshold_depletion-d.vfb
    if threshold_relative+.25 > 2:
        raise ValueError("This three-regime demo needs estimated Vth-Vfb<=1.75 V; reduce tox or NA")
    volts = np.linspace(-.3, max(1.1, threshold_relative+.25), 85) + d.vfb
    states = sweep(d, volts, nodes)
    refs = [reference_bias(d, float(v)) for v in volts]
    cnum = np.array([s["Cqs_F_m2"] for s in states])
    cref = np.array([s["Cqs_F_m2"] for s in refs])
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    axes[0].plot(volts, cnum/d.cox, lw=2, label="Quasistatic: nonlinear Poisson")
    axes[0].plot(volts[::4], cref[::4]/d.cox, "o", ms=4, fillstyle="none", label="Independent first-integral reference")
    axes[0].plot(volts, [s["Chf_frozen_F_m2"]/d.cox for s in states], "--", label="HF approximation: frozen electron AC response")
    for ax in axes:
        ax.axvline(d.vfb, color="gray", ls=":", label="Flat band")
        ax.axvline(d.threshold_depletion, color="orange", ls=":", label="Depletion threshold estimate")
        ax.set_xlabel("Gate voltage (V)")
        ax.grid(alpha=.25)
    axes[0].set(ylabel="Capacitance per area / Cox", title="P-type silicon MOS capacitor")
    axes[0].legend(fontsize=8)
    axes[1].plot(volts, [s["psi_s_V"] for s in states], lw=2, label="Surface potential")
    axes[1].axhline(2*d.phi_f, color="purple", ls="--", label="2 phi_F: strong-inversion criterion")
    axes[1].set(ylabel="Surface potential (V)", title="Gate control of the surface")
    axes[1].legend(fontsize=8)
    savefig(fig, out, "cv_surface_potential.png")
    rows = []
    for s, ref in zip(states, refs):
        rows.append({**{key: s[key] for key in SCALARS},
                     "reference_Cqs_F_m2": ref["Cqs_F_m2"], "reference_psi_s_V": ref["psi_s_V"]})
    write_csv(out/"cv.csv", rows)

    biases = (-.25+d.vfb, .45*threshold_relative+d.vfb, d.threshold_depletion+.2)
    labels = ("Accumulation", "Depletion", "Inversion")
    profiles = [solve_bias(d, v, nodes) for v in biases]
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), constrained_layout=True)
    for col, (s, name, color) in enumerate(zip(profiles, labels, COLORS)):
        x = s["x"]*1e6
        ei = d.phi_f-s["psi"]
        axes[0, col].plot(x, ei+eg_ev(d.temperature)/2, color=color, label="Ec")
        axes[0, col].plot(x, ei, "--", color=color, label="Ei (midgap approximation)")
        axes[0, col].plot(x, ei-eg_ev(d.temperature)/2, color=color, alpha=.65, label="Ev")
        axes[0, col].axhline(0, color="black", ls=":", label="Si EF=0")
        axes[0, col].set(title=f"{name}, Vg={s['Vg_V']:.2f} V", ylabel="Si energy relative to EF (eV)")
        axes[1, col].semilogy(x, s["p"]/1e6, label="Holes p")
        axes[1, col].semilogy(x, s["n"]/1e6, label="Electrons n")
        axes[1, col].axhline(d.na_cm3, color="gray", ls=":", label="NA")
        axes[1, col].set(ylabel="Carrier density (cm^-3)")
        for ax in axes[:, col]:
            ax.set_xlim(0, min(s["x"][-1]*1e6, .6))
            ax.set_xlabel("Depth into silicon (um); interface at x=0")
            ax.grid(alpha=.25)
            ax.legend(fontsize=7)
        np.savetxt(out/f"profile_{name.lower()}.csv", np.column_stack(
            (s["x"], s["psi"], s["rho"], s["field"], s["n"], s["p"])), delimiter=",",
            header="x_m,psi_V,rho_C_m3,E_gradient_V_m,n_m-3,p_m-3", comments="")
    savefig(fig, out, "bands_carriers.png")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
    for s, name in zip(profiles, labels):
        axes[0].plot(s["x"]*1e6, s["rho"], label=name)
        axes[1].plot(s["x"]*1e6, s["field"], label=name)
        trace = s["trace"]
        axes[2].semilogy([r["iteration"] for r in trace], [max(r["residual_inf"], 1e-16) for r in trace], "o-", label=name)
    axes[0].set(xlabel="Depth into Si (um)", ylabel="Charge density (C/m^3)", title="Self-consistent charge")
    axes[1].set(xlabel="Depth into Si (um)", ylabel="Electric field (V/m)", title="Field from potential gradient")
    axes[2].set(xlabel="Newton iteration", ylabel="Scaled residual infinity norm", title="Damped Newton from zero initial guess")
    for ax in axes[:2]:
        ax.set_xlim(0, .6)
    for ax in axes:
        ax.legend(fontsize=8)
        ax.grid(alpha=.25)
    savefig(fig, out, "charge_field_newton.png")
    (out/"newton_traces.json").write_text(json.dumps(
        {name: s["trace"] for s, name in zip(profiles, labels)}, indent=2)+"\n")

    groups = (("tox_nm", (8., 10., 20.), "Oxide thickness (nm)"),
              ("na_cm3", (3e15, 1e16, 3e16), "Acceptor density (cm^-3)"),
              ("qox_cm2", (-1e11, 0., 1e11), "Fixed oxide charge / q (cm^-2)"))
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
    sweep_rows = []
    scan_v = np.linspace(-.3, 1.15, 65)+d.phi_ms_v
    for ax, (param, values, title) in zip(axes, groups):
        for value in values:
            dev = replace(d, **{param: value})
            scan = sweep(dev, scan_v, nodes)
            ax.plot(scan_v, np.array([r["Cqs_F_m2"] for r in scan])*100, label=f"{value:g}")
            for s in scan:
                sweep_rows.append({"parameter": param, "value": value, "Vg_V": s["Vg_V"],
                                   "Cqs_uF_cm2": s["Cqs_F_m2"]*100, "psi_s_V": s["psi_s_V"],
                                   "Vfb_V": dev.vfb, "residual_inf": s["residual_inf"],
                                   "classical_density_to_DOS_max": s["classical_density_to_DOS_max"]})
        ax.set(title=title, xlabel="Gate voltage (V)", ylabel="QS capacitance per area (uF/cm^2)")
        ax.legend(fontsize=8)
        ax.grid(alpha=.25)
    savefig(fig, out, "parameter_sweeps.png")
    write_csv(out/"parameter_sweeps.csv", sweep_rows)

    refinement = []
    check_v = (-.25+d.vfb, .4+d.vfb, 1.05+d.vfb)
    for count in (101, 201, 401, 801):
        solutions = sweep(d, np.array(check_v), count)
        reference = [reference_bias(d, v) for v in check_v]
        refinement.append({"nodes": count,
            "max_surface_error_V": max(abs(s["psi_s_V"]-r["psi_s_V"]) for s, r in zip(solutions, reference)),
            "max_capacitance_error_over_Cox": max(abs(s["Cqs_F_m2"]-r["Cqs_F_m2"])/d.cox for s, r in zip(solutions, reference))})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    counts = np.array([r["nodes"]-1 for r in refinement])
    for ax, key, ylabel in zip(axes, ("max_surface_error_V", "max_capacitance_error_over_Cox"),
                               ("Max surface-potential error (V)", "Max capacitance error / Cox")):
        errors = np.array([r[key] for r in refinement])
        ax.loglog(counts, errors, "o-", label="Independent reference comparison")
        ax.loglog(counts, errors[0]*(counts[0]/counts)**2, "--", label="Second-order guide")
        ax.set(xlabel="Number of mesh intervals", ylabel=ylabel)
        ax.legend(fontsize=8)
        ax.grid(True, which="both", alpha=.25)
    savefig(fig, out, "mesh_convergence.png")
    write_csv(out/"mesh_convergence.csv", refinement)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)
    axes[0].semilogy(volts, np.maximum([s["residual_inf"] for s in states], 1e-16))
    axes[0].axhline(2e-10, ls="--", color="gray", label="Newton tolerance")
    axes[0].set(ylabel="Scaled residual infinity norm", title="Nonlinear equation convergence")
    axes[0].legend(fontsize=8)
    axes[1].plot(volts, [s["charge_balance_C_m2"] for s in states])
    axes[1].set(ylabel="Qg + Qox + Qs - D_bulk (C/m^2)", title="Discrete Gauss-law balance")
    axes[2].plot(volts, (cnum-cref)/d.cox*100)
    axes[2].set(ylabel="(C_numeric - C_reference) / Cox (%)", title="Independent capacitance check")
    for ax in axes:
        ax.set_xlabel("Gate voltage (V)")
        ax.grid(alpha=.25)
    savefig(fig, out, "validation_dashboard.png")
    depths = [solve_bias(d, 1.05+d.vfb, 801, depth_factor=f) for f in (1., 1.5)]
    numerical_threshold = brentq(lambda v: solve_bias(d, v, nodes)["psi_s_V"]-2*d.phi_f,
                                d.vfb, d.threshold_depletion+.2, xtol=1e-9)
    metrics = {"device": d.summary(), "nodes": nodes, "bias_points": len(volts),
        "parameter_sweep_points": len(sweep_rows),
        "max_cv_error_over_Cox": float(np.max(np.abs(cnum-cref))/d.cox),
        "max_surface_error_V": max(abs(s["psi_s_V"]-r["psi_s_V"]) for s, r in zip(states, refs)),
        "max_residual_inf": max(s["residual_inf"] for s in states),
        "max_charge_balance_C_m2": max(abs(s["charge_balance_C_m2"]) for s in states),
        "max_bulk_displacement_C_m2": max(abs(s["bulk_displacement_C_m2"]) for s in states),
        "max_classical_density_to_DOS": max(s["classical_density_to_DOS_max"] for s in states),
        "sweep_max_classical_density_to_DOS": max(s["classical_density_to_DOS_max"] for s in sweep_rows),
        "classical_caution_threshold": .05,
        "profile_biases_V": dict(zip(labels, biases)),
        "threshold_numeric_V": numerical_threshold,
        "threshold_first_integral_V": d.vfb+2*d.phi_f-reference_bias_at_threshold(d)/d.cox,
        "depth_and_mesh_sensitivity_surface_V": abs(depths[0]["psi_s_V"]-depths[1]["psi_s_V"]),
        "depth_study_note": "Fixed 801 nodes, depth factors 1 and 1.5; spacing also changes, so this is a combined sensitivity check",
        "refinement": refinement,
        "HF_model": "DC inversion charge retained; AC electron density perturbation frozen; no frequency axis"}
    (out/"validation.json").write_text(json.dumps(metrics, indent=2)+"\n")
    return metrics


def reference_bias_at_threshold(d: Device) -> float:
    from moscap import reference_charge
    return reference_charge(d, 2*d.phi_f)[0]
