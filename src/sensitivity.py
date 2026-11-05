from typing import Dict, Any, List
import pandas as pd
from .one_compartment import OneCompartmentModel


class ParameterSensitivityAnalyzer:
    """Evaluates local sensitivity of PK outcomes (Cmax, Tmax, AUC, t1/2) to parameter perturbations."""

    @staticmethod
    def run_sensitivity_sweep(
        base_model: OneCompartmentModel,
        perturbations: List[float] = [-0.30, -0.15, 0.0, 0.15, 0.30],
    ) -> pd.DataFrame:
        params_to_vary = ["cl", "vd", "ka", "f"]
        records = []

        base_summary = base_model.get_summary_parameters()
        base_cmax = base_summary["cmax_mg_l"]
        base_tmax = base_summary["tmax_h"]
        base_auc = base_summary["auc_inf_mg_h_l"]
        base_thalf = base_summary["elimination_half_life_h"]

        for param in params_to_vary:
            base_val = getattr(base_model, param)
            for delta in perturbations:
                new_val = base_val * (1.0 + delta)

                # Instantiate perturbed model kwargs
                kwargs = {
                    "dose": base_model.dose,
                    "vd": base_model.vd,
                    "cl": base_model.cl,
                    "ka": base_model.ka,
                    "f": base_model.f,
                }
                kwargs[param] = new_val

                try:
                    perturbed = OneCompartmentModel(**kwargs)
                    p_summary = perturbed.get_summary_parameters()

                    records.append({
                        "parameter": param.upper(),
                        "perturbation_percent": round(delta * 100.0, 1),
                        "parameter_value": round(new_val, 4),
                        "cmax_mg_l": p_summary["cmax_mg_l"],
                        "tmax_h": p_summary["tmax_h"],
                        "auc_inf": p_summary["auc_inf_mg_h_l"],
                        "half_life_h": p_summary["elimination_half_life_h"],
                        "cmax_change_percent": round(((p_summary["cmax_mg_l"] - base_cmax) / base_cmax) * 100.0, 2),
                        "auc_change_percent": round(((p_summary["auc_inf_mg_h_l"] - base_auc) / base_auc) * 100.0, 2),
                    })
                except Exception:
                    pass

        return pd.DataFrame(records)
