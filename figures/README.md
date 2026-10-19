# Generated results

Regenerate with `python run_all.py`. Outputs use the baseline documented in
the root README; custom runs should use a separate output directory.

Seven plots cover C-V/surface potential, silicon bands and carriers, charge
and Newton histories, parameter sweeps, mesh convergence, validation diagnostics,
and synthetic parameter extraction. CSV headers include units. Capacitance and
charge are per unit area, not total device capacitance or total device charge.

`data/` contains explicitly synthetic noisy C-V data and its generation metadata.
The truth model uses a first integral; the inverse uses finite-volume Poisson.
The high-frequency plot is a frozen-electron AC-response approximation with
no physical frequency axis. See `../docs/MODEL.md` for the assumptions.
