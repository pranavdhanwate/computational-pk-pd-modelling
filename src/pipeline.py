import argparse
import json
from pathlib import Path
from typing import Dict, Any, Tuple
import numpy as np
import pandas as pd

from .one_compartment import OneCompartmentModel
from .two_compartment import TwoCompartmentModel
from .pharmacodynamics import PharmacodynamicModel
from .pk_pd_link import PKPDLink
from .multiple_dosing import MultipleDosingSimulator
from .sensitivity import ParameterSensitivityAnalyzer
from .visualize import plot_single_dose, plot_hill_curve, plot_multidose_steady_state, plot_sensitivity_spider


def run_pkpd_simulation_pipeline(
    params_json: str | Path = "data/clinical_drug_parameters.json",
    selected_drug: str = "Theophylline",
    output_dir: str | Path = "results",
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Executes single-dose, pharmacodynamic, multi-dose steady-state, and sensitivity simulations."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    with open(params_json, "r", encoding="utf-8") as f:
        drugs_data = json.load(f)

    if selected_drug not in drugs_data:
        raise ValueError(f"Drug '{selected_drug}' not found in {params_json}")

    drug_info = drugs_data[selected_drug]

    # 1. One-Compartment PK Simulation
    pk_model = OneCompartmentModel(
        dose=drug_info["typical_dose_mg"],
        vd=drug_info["volume_of_distribution_l"],
        cl=drug_info["clearance_l_h"],
        ka=drug_info["absorption_rate_constant_h"],
        f=drug_info["bioavailability_f"],
    )
    pk_summary = pk_model.get_summary_parameters()
    pd.DataFrame([pk_summary]).to_csv(out_dir / "single_dose_summary.csv", index=False)

    t_single = np.linspace(0.0, 48.0, 300)
    c_iv = pk_model.simulate_iv_analytical(t_single)
    c_oral = pk_model.simulate_oral_analytical(t_single)
    tmax, cmax = pk_model.calculate_tmax_cmax()

    plot_single_dose(t_single, c_iv, c_oral, tmax, cmax, fig_dir / "single_dose_profile.png", drug_name=selected_drug)

    # 2. Pharmacodynamics (Hill Curve)
    pd_model = PharmacodynamicModel(
        e0=0.0,
        emax=drug_info["pd_emax"],
        ec50=drug_info["pd_ec50_mg_l"],
        gamma=drug_info["pd_hill_gamma"],
    )
    c_range = np.logspace(-1, 2, 200)
    effects = pd_model.calculate_effect(c_range)
    plot_hill_curve(c_range, effects, drug_info["pd_ec50_mg_l"], drug_info["pd_emax"], fig_dir / "sigmoidal_hill_curve.png")

    # 3. Multiple-Dose Steady-State Accumulation
    md_simulator = MultipleDosingSimulator(pk_model, tau=drug_info["dosing_interval_h"])
    md_summary = md_simulator.calculate_steady_state_metrics()
    pd.DataFrame([md_summary]).to_csv(out_dir / "multidose_steady_state.csv", index=False)

    t_multi, c_multi = md_simulator.simulate_multidose_profile(n_doses=7, route="oral")
    window = tuple(drug_info["therapeutic_window_mg_l"])
    plot_multidose_steady_state(t_multi, c_multi, window, fig_dir / "multiple_dosing_steady_state.png", drug_name=selected_drug)

    # 4. Parameter Sensitivity Analysis
    sens_df = ParameterSensitivityAnalyzer.run_sensitivity_sweep(pk_model)
    sens_df.to_csv(out_dir / "sensitivity_analysis.csv", index=False)
    plot_sensitivity_spider(sens_df, fig_dir / "sensitivity_spider_plot.png")

    return pk_summary, md_summary


def main():
    parser = argparse.ArgumentParser(description="Computational PK-PD Modelling Engine")
    parser.add_argument("--params", "-p", type=str, default="data/clinical_drug_parameters.json", help="Drug parameters JSON")
    parser.add_argument("--drug", "-d", type=str, default="Theophylline", help="Drug name")
    parser.add_argument("--output", "-o", type=str, default="results", help="Output directory")
    args = parser.parse_args()

    pk_res, md_res = run_pkpd_simulation_pipeline(args.params, args.drug, output_dir=args.output)
    print(f"\nPK-PD Simulation Complete for {args.drug}!")
    print("Single-Dose Pharmacokinetics:")
    for k, v in pk_res.items():
        print(f"  {k}: {v}")
    print("\nSteady-State Kinetics:")
    for k, v in md_res.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
