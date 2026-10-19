# Validation record

Baseline: NA=1e16 cm^-3, tox=10 nm, T=300 K, phi_ms=0 V, Qox=0.
Main mesh: 401 stretched nodes. Inverse mesh: 301 nodes. Random seed: 20261019.

Tested environment: Python 3.12, NumPy 2.3.5, SciPy 1.18.1, Matplotlib 3.11.2.

Executed:

```bash
python run_all.py
python -m unittest discover -s tests -v
```

13 tests passed: flatband and neutrality; independent reference; charge and
far-boundary checks; surface regimes; mesh refinement; capacitance derivative;
reference derivative near flatband; frozen-electron HF behavior; oxide-charge
shift; Jacobian directional derivative; continuation versus direct solve;
noise-free inverse recovery; input guards.

The complete run produces 7 PNG figures, 7 CSV files at the figures root,
3 JSON diagnostics, and 2 synthetic-data files in figures/data. The main sweep
uses 85 biases. Parameter sweeps contain 585 solved biases in total.

Representative results:

| Quantity | Value |
| --- | ---: |
| Max surface potential error against reference | 9.73e-6 V |
| Max C error / Cox | 5.98e-5 |
| Max scaled nonlinear residual | 1.76e-10 |
| Max discrete Gauss-law balance error | 3.53e-13 C/m^2 |
| Depletion threshold | 0.855333575 V |
| Numerical threshold | 0.855334700 V |
| Synthetic fitted tox | 11.9942498 nm |
| Synthetic fitted NA | 1.2036894e16 cm^-3 |

Authoritative machine-readable values are figures/validation.json and
figures/extraction.json. Tiny differences between library versions are normal.
This record describes the default case; custom outputs require their own record.

The finite-depth check changes both depth and spacing at fixed node count;
it is explicitly a combined sensitivity check. Nondegenerate density screening
is recorded separately from numerical error. The HF curve is an approximation,
not a solution at a specified physical frequency. Synthetic inverse results are
not measured data and do not establish general extraction accuracy.
