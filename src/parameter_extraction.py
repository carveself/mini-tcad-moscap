"""Fit tox and NA from explicitly synthetic noisy C-V data.

Data use the independent semi-infinite first-integral reference. The inverse
model uses finite-volume Poisson, so forward/inverse discretizations differ.
Work-function difference, fixed oxide charge, temperature, and capacitance
normalization per area are treated as known, not additional fitted variables.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from moscap import Device, reference_bias, sweep
from experiments import plt, savefig, write_csv


def fit_cv(voltages: np.ndarray, capacitance: np.ndarray, sigma: float,
           known: Device, nodes: int = 301, initial: tuple = (14., 2e16)) -> dict:
    voltages, capacitance = np.asarray(voltages, float), np.asarray(capacitance, float)
    if voltages.shape != capacitance.shape or voltages.ndim != 1 or len(voltages) < 8:
        raise ValueError("Need at least 8 paired 1D voltage and capacitance points")
    if not np.isfinite(sigma) or sigma <= 0 or not np.all(np.isfinite(capacitance)):
        raise ValueError("Capacitances and positive noise sigma must be finite")
    evaluations = []
    def residual(log_params):
        tox, na = np.exp(log_params)
        dev = Device(na_cm3=na, tox_nm=tox, temperature=known.temperature,
                     phi_ms_v=known.phi_ms_v, qox_cm2=known.qox_cm2)
        prediction = np.array([s["Cqs_F_m2"] for s in sweep(dev, voltages, nodes)])
        residuals = (prediction-capacitance)/sigma
        evaluations.append({"evaluation": len(evaluations), "tox_nm": float(tox),
                            "NA_cm-3": float(na), "weighted_sum_squared": float(residuals@residuals)})
        return residuals
    result = least_squares(residual, np.log(initial), bounds=(np.log([5., 1.1e15]), np.log([30., 8e16])),
                           diff_step=1e-4, ftol=1e-8, xtol=1e-8, gtol=1e-7, max_nfev=45)
    if not result.success:
        raise RuntimeError(f"Parameter fit failed: {result.message}")
    tox, na = np.exp(result.x)
    fitted = Device(na_cm3=na, tox_nm=tox, temperature=known.temperature,
                    phi_ms_v=known.phi_ms_v, qox_cm2=known.qox_cm2)
    prediction = np.array([s["Cqs_F_m2"] for s in sweep(fitted, voltages, nodes)])
    reduced = float(result.fun@result.fun/(len(voltages)-2))
    covariance = np.linalg.pinv(result.jac.T@result.jac)*reduced
    sd = np.sqrt(np.maximum(np.diag(covariance), 0))
    correlation = float(covariance[0, 1]/max(sd[0]*sd[1], 1e-30))
    return {"device": fitted, "prediction": prediction, "evaluations": evaluations,
            "metrics": {"success": bool(result.success), "message": result.message,
                        "tox_nm": float(tox), "NA_cm-3": float(na), "nodes": nodes,
                        "solver_evaluations": len(evaluations), "reduced_chi_squared": reduced,
                        "RMSE_F_m2": float(np.sqrt(np.mean((prediction-capacitance)**2))),
                        "approx_relative_1sigma_tox": float(sd[0]),
                        "approx_relative_1sigma_NA": float(sd[1]),
                        "log_parameter_correlation": correlation,
                        "jacobian_condition_number": float(np.linalg.cond(result.jac))}}


def run_extraction(out: Path, data_out: Path, nodes: int = 301) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    data_out.mkdir(parents=True, exist_ok=True)
    truth = Device(na_cm3=1.2e16, tox_nm=12.)
    voltages = np.linspace(-.3, 1.1, 61)
    clean = np.array([reference_bias(truth, float(v))["Cqs_F_m2"] for v in voltages])
    sigma = .004*truth.cox
    seed = 20261019
    noisy = clean+np.random.default_rng(seed).normal(0, sigma, size=clean.shape)
    result = fit_cv(voltages, noisy, sigma, truth, nodes)
    metrics = {"data_origin": "SYNTHETIC, not measured; independent first-integral model",
               "seed": seed, "noise_sigma_F_m2": sigma,
               "known_parameters": {"T_K": truth.temperature, "phi_ms_V": truth.phi_ms_v,
                                    "Qox_over_q_cm-2": truth.qox_cm2, "area_normalization": "per unit area"},
               "truth": truth.summary(), "fit": result["metrics"],
               "relative_error_tox": (result["device"].tox_nm/truth.tox_nm-1),
               "relative_error_NA": (result["device"].na_cm3/truth.na_cm3-1),
               "uncertainty_note": "Local linearization, conditional on assumed model and known parameters; not a global confidence guarantee"}
    rows = [{"Vg_V": float(v), "clean_reference_C_F_m2": float(c), "synthetic_noisy_C_F_m2": float(n),
             "sigma_F_m2": sigma, "fitted_Poisson_C_F_m2": float(p)}
            for v, c, n, p in zip(voltages, clean, noisy, result["prediction"])]
    write_csv(data_out/"synthetic_cv.csv", rows)
    write_csv(out/"fit_evaluations.csv", result["evaluations"])
    (out/"extraction.json").write_text(json.dumps(metrics, indent=2)+"\n")
    (data_out/"synthetic_metadata.json").write_text(json.dumps(metrics, indent=2)+"\n")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    axes[0].plot(voltages, clean*100, label="Synthetic truth: first-integral reference")
    axes[0].plot(voltages, noisy*100, ".", label="Synthetic observations + Gaussian noise")
    axes[0].plot(voltages, result["prediction"]*100, "--", label="Fitted nonlinear-Poisson model")
    axes[0].set(xlabel="Gate voltage (V)", ylabel="QS capacitance per area (uF/cm^2)",
                title="Synthetic parameter extraction (not experimental data)")
    axes[0].legend(fontsize=8)
    axes[1].plot(voltages, (result["prediction"]-noisy)/sigma, "o-", ms=3)
    axes[1].axhline(0, color="gray", ls=":")
    axes[1].set(xlabel="Gate voltage (V)", ylabel="Fit residual / noise sigma",
                title=f"tox error={metrics['relative_error_tox']:.2%}; NA error={metrics['relative_error_NA']:.2%}")
    for ax in axes:
        ax.grid(alpha=.25)
    savefig(fig, out, "parameter_extraction.png")
    return metrics
