"""Small checkpoints for following docs/REPRODUCE.md without running everything."""
import argparse
import json
from pathlib import Path

import numpy as np

from moscap import Device, ROOT, reference_bias, solve_bias
from experiments import generate, write_csv
from parameter_extraction import run_extraction


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("reference", "profiles", "cv", "extract"))
    parser.add_argument("--output", type=Path, default=ROOT/"reproduction")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    d = Device()
    if args.stage == "reference":
        rows = [reference_bias(d, float(v)) for v in np.linspace(-.3, 1.1, 85)]
        write_csv(args.output/"reference.csv", rows)
        print(json.dumps(d.summary(), indent=2))
    elif args.stage == "profiles":
        for name, v in (("accumulation", -.25), ("depletion", .4), ("inversion", 1.05)):
            s = solve_bias(d, v)
            np.savetxt(args.output/f"{name}.csv", np.column_stack((s["x"], s["psi"], s["p"], s["n"])),
                       delimiter=",", header="x_m,psi_V,p_m-3,n_m-3", comments="")
            print(name, json.dumps({key: s[key] for key in ("Vg_V", "psi_s_V", "iterations", "residual_inf", "charge_balance_C_m2")}))
    elif args.stage == "cv":
        print(json.dumps(generate(d, args.output), indent=2))
    else:
        print(json.dumps(run_extraction(args.output, args.output/"data"), indent=2))
    print(f"Checkpoint outputs: {args.output.resolve()}")


if __name__ == "__main__":
    main()
