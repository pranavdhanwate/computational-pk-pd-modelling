import math
from typing import Dict, Any, Tuple
import numpy as np
from .one_compartment import OneCompartmentModel


class MultipleDosingSimulator:
    """Simulates multi-dose accumulation and steady-state concentration kinetics."""

    def __init__(self, pk_model: OneCompartmentModel, tau: float = 12.0):
        self.pk = pk_model
        self.tau = float(tau)

        # Accumulation index: R = 1 / (1 - exp(-ke * tau))
        self.accumulation_factor = 1.0 / (1.0 - math.exp(-self.pk.ke * self.tau))
        # Average steady-state concentration: Cav,ss = (F * Dose) / (CL * tau)
        self.cav_ss = (self.pk.f * self.pk.dose) / (self.pk.cl * self.tau)

    def simulate_multidose_profile(
        self,
        n_doses: int = 7,
        points_per_interval: int = 100,
        route: str = "oral",
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Simulates full concentration-time curve across consecutive doses."""
        total_time = n_doses * self.tau
        t_eval = np.linspace(0.0, total_time, n_doses * points_per_interval)
        concentrations = np.zeros_like(t_eval)

        for dose_idx in range(n_doses):
            dose_time = dose_idx * self.tau
            # Mask for times after this dose
            mask = t_eval >= dose_time
            dt = t_eval[mask] - dose_time

            if route.lower() == "oral":
                c_dose = self.pk.simulate_oral_analytical(dt)
            else:
                c_dose = self.pk.simulate_iv_analytical(dt)

            concentrations[mask] += c_dose

        return t_eval, concentrations

    def calculate_steady_state_metrics(self) -> Dict[str, float]:
        tmax, cmax_single = self.pk.calculate_tmax_cmax()
        cmax_ss = cmax_single * self.accumulation_factor
        cmin_ss = cmax_ss * math.exp(-self.pk.ke * (self.tau - tmax))
        fluctuation = ((cmax_ss - cmin_ss) / self.cav_ss) * 100.0 if self.cav_ss > 0 else 0.0

        return {
            "dosing_interval_tau_h": self.tau,
            "accumulation_factor_R": round(self.accumulation_factor, 3),
            "cav_ss_mg_l": round(self.cav_ss, 2),
            "cmax_ss_mg_l": round(cmax_ss, 2),
            "cmin_ss_mg_l": round(cmin_ss, 2),
            "percent_fluctuation": round(fluctuation, 2),
            "time_to_steady_state_h": round(4.5 * self.pk.half_life, 2),
        }
