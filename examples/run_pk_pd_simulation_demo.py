import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from src.pipeline import run_pkpd_simulation_pipeline

if __name__ == "__main__":
    params_path = project_root / "data" / "clinical_drug_parameters.json"
    out_dir = project_root / "results"

    print("Running Computational PK-PD Simulation Demo (Theophylline)...")
    pk_sum, md_sum = run_pkpd_simulation_pipeline(params_path, selected_drug="Theophylline", output_dir=out_dir)

    print("\nPharmacokinetic Parameters:")
    for k, v in pk_sum.items():
        print(f"  {k}: {v}")

    print("\nMultiple-Dose Steady-State Metrics:")
    for k, v in md_sum.items():
        print(f"  {k}: {v}")
