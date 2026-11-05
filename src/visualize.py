from pathlib import Path
from typing import Tuple, List, Optional
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_single_dose(
    t: np.ndarray,
    c_iv: np.ndarray,
    c_oral: np.ndarray,
    tmax: float,
    cmax: float,
    output_path: str | Path,
    drug_name: str = "Theophylline",
) -> Path:
    """Renders single-dose comparative IV bolus vs Oral concentration-time curves."""
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=300)

    ax.plot(t, c_iv, color="#C62828", linestyle="--", linewidth=2.0, label="IV Bolus Administration")
    ax.plot(t, c_oral, color="#1565C0", linewidth=2.2, label="Oral Administration (First-Order Absorption)")

    # Highlight Cmax / Tmax
    ax.scatter([tmax], [cmax], color="#E65100", s=80, zorder=5, label=f"Cmax = {cmax:.1f} mg/L (Tmax = {tmax:.1f} h)")
    ax.axvline(tmax, color="#E65100", linestyle=":", linewidth=1.0, alpha=0.7)

    ax.set_xlabel("Time Post-Dose (hours)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Plasma Drug Concentration (mg/L)", fontsize=11, fontweight="bold")
    ax.set_title(f"Pharmacokinetic Concentration-Time Profile: {drug_name}", fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="upper right", fontsize=9.5)

    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
    return out_path


def plot_hill_curve(
    concentrations: np.ndarray,
    effects: np.ndarray,
    ec50: float,
    emax: float,
    output_path: str | Path,
) -> Path:
    """Renders Sigmoidal Emax Hill pharmacodynamic dose-response relationship."""
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=300)

    ax.plot(concentrations, effects, color="#2E7D32", linewidth=2.5, label="Hill Sigmoidal Response Curve")
    ax.axvline(ec50, color="#FF8F00", linestyle="--", linewidth=1.2, label=f"EC50 = {ec50:.1f} mg/L")
    ax.axhline(emax * 0.5, color="#FF8F00", linestyle=":", linewidth=1.0, label=f"50% Emax ({emax * 0.5:.1f}%)")

    ax.set_xscale("log")
    ax.set_xlabel("Drug Concentration (mg/L, log scale)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Pharmacodynamic Effect (% Emax)", fontsize=11, fontweight="bold")
    ax.set_title("Pharmacodynamic Concentration-Response Relationship", fontsize=12, fontweight="bold")
    ax.grid(True, which="both", linestyle=":", alpha=0.6)
    ax.legend(loc="lower right", fontsize=9.5)

    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
    return out_path


def plot_multidose_steady_state(
    t: np.ndarray,
    concentrations: np.ndarray,
    window: Tuple[float, float],
    output_path: str | Path,
    drug_name: str = "Theophylline",
) -> Path:
    """Renders multi-dose accumulation profile against therapeutic index window."""
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    mec, mtc = window
    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=300)

    # Therapeutic window shading
    ax.axhspan(mec, mtc, color="#E8F5E9", alpha=0.7, label=f"Therapeutic Window ({mec} - {mtc} mg/L)")
    ax.axhline(mtc, color="#C62828", linestyle="--", linewidth=1.2, label=f"Minimum Toxic Conc (MTC = {mtc} mg/L)")
    ax.axhline(mec, color="#F57C00", linestyle="--", linewidth=1.2, label=f"Minimum Effective Conc (MEC = {mec} mg/L)")

    # Concentration profile
    ax.plot(t, concentrations, color="#0D47A1", linewidth=2.0, label="Repeated Oral Dosing (q12h)")

    ax.set_xlabel("Time (hours)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Plasma Drug Concentration (mg/L)", fontsize=11, fontweight="bold")
    ax.set_title(f"Multiple-Dose Pharmacokinetics & Steady-State Accumulation: {drug_name}", fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="upper right", fontsize=9)

    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
    return out_path


def plot_sensitivity_spider(
    sens_df: pd.DataFrame,
    output_path: str | Path,
) -> Path:
    """Renders spider sensitivity plot showing impact of physiological parameters on AUC."""
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=300)

    colors = {"CL": "#C62828", "VD": "#1565C0", "KA": "#2E7D32", "F": "#6A1B9A"}

    for param, group in sens_df.groupby("parameter"):
        sorted_g = group.sort_values(by="perturbation_percent")
        ax.plot(
            sorted_g["perturbation_percent"],
            sorted_g["auc_change_percent"],
            marker="o",
            linewidth=2.0,
            color=colors.get(param, "black"),
            label=f"Parameter {param}",
        )

    ax.axhline(0.0, color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(0.0, color="gray", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Parameter Variation (%)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Change in Systemic Exposure AUC (%)", fontsize=11, fontweight="bold")
    ax.set_title("Pharmacokinetic Parameter Sensitivity Analysis (AUC Impact)", fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="upper left", fontsize=9.5)

    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
    return out_path
