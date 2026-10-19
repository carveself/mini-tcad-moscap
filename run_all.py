"""Generate the complete MOS capacitor learning project results."""
from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/"src"))
from moscap import Device
from experiments import generate
from parameter_extraction import run_extraction


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--na", type=float, default=1e16, help="Acceptor density [cm^-3]")
    parser.add_argument("--tox", type=float, default=10., help="Oxide thickness [nm]")
    parser.add_argument("--temperature", type=float, default=300., help="Temperature [K]")
    parser.add_argument("--nodes", type=int, default=401, help="Stretched grid nodes, >=41")
    parser.add_argument("--fit-nodes", type=int, default=301, help="Independent inverse-model grid")
    parser.add_argument("--output", type=Path, default=ROOT/"figures")
    parser.add_argument("--data-output", type=Path, default=None)
    args = parser.parse_args()
    data_out = args.data_output if args.data_output is not None else args.output/"data"
    try:
        d = Device(args.na, args.tox, args.temperature)
        validation = generate(d, args.output, args.nodes)
        fit = run_extraction(args.output, data_out, args.fit_nodes)
    except ValueError as exc:
        parser.error(str(exc))
    print(json.dumps({"device": d.summary(),
          "max_cv_error_over_Cox": validation["max_cv_error_over_Cox"],
          "max_residual_inf": validation["max_residual_inf"],
          "max_charge_balance_C_m2": validation["max_charge_balance_C_m2"],
          "fit_relative_error_tox": fit["relative_error_tox"],
          "fit_relative_error_NA": fit["relative_error_NA"]}, indent=2))
    print(f"Saved figures to {args.output.resolve()}; synthetic data to {data_out.resolve()}")
    print("Extraction uses a fixed, documented synthetic case, independent of baseline CLI parameters.")


if __name__ == "__main__":
    main()
