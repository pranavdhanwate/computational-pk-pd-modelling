import math
from typing import Dict, Any, Tuple
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp


class OneCompartmentModel:
    """One-compartment pharmacokinetic model for IV bolus and oral first-order absorption."""

    def __init__(
        self,
        dose: float,
        vd: float,
        cl: float,
        ka: float = 1.0,
        f: float = 1.0,
    ):
        self.dose = float(dose)
        self.vd = float(vd)
        self.cl = float(cl)
        self.ka = float(ka)
        self.f = float(f)

        # Derived parameters
        self.ke = self.cl / self.vd
        self.half_life = math.log(2.0) / self.ke
        self.auc_inf = (self.f * self.dose) / self.cl

    def calculate_tmax_cmax(self) -> Tuple[float, float]:
        """Calculates exact analytical Tmax (h) and Cmax (mg/L) for oral administration."""
        if self.ka <= self.ke:
            raise ValueError("Absorption rate constant ka must be greater than elimination ke.")
        tmax = math.log(self.ka / self.ke) / (self.ka - self.ke)
        cmax = self.simulate_oral_analytical(np.array([tmax]))[0]
        return round(float(tmax), 3), round(float(cmax), 3)

    def simulate_iv_analytical(self, t: np.ndarray) -> np.ndarray:
        """C(t) = (Dose / Vd) * exp(-ke * t)"""
        c0 = self.dose / self.vd
        return c0 * np.exp(-self.ke * t)

    def simulate_oral_analytical(self, t: np.ndarray) -> np.ndarray:
        """C(t) = [F * Dose * ka / (Vd * (ka - ke))] * [exp(-ke * t) - exp(-ka * t)]"""
        prefactor = (self.f * self.dose * self.ka) / (self.vd * (self.ka - self.ke))
        return prefactor * (np.exp(-self.ke * t) - np.exp(-self.ka * t))

    def simulate_oral_ode(self, t_eval: np.ndarray) -> np.ndarray:
        """Numerically integrates ODEs: dAd/dt = -ka*Ad, dAc/dt = ka*Ad - ke*Ac"""
        def odes(t, y):
            ad, ac = y
            dad_dt = -self.ka * ad
            dac_dt = (self.ka * ad) - (self.ke * ac)
            return [dad_dt, dac_dt]

        y0 = [self.f * self.dose, 0.0]
        sol = solve_ivp(odes, [t_eval[0], t_eval[-1]], y0, t_eval=t_eval, rtol=1e-6, atol=1e-9)
        ac = sol.y[1]
        return ac / self.vd

    def get_summary_parameters(self) -> Dict[str, Any]:
        tmax, cmax = self.calculate_tmax_cmax() if self.ka > self.ke else (0.0, self.dose / self.vd)
        return {
            "dose_mg": self.dose,
            "volume_distribution_l": self.vd,
            "clearance_l_h": self.cl,
            "elimination_rate_ke_h": round(self.ke, 4),
            "elimination_half_life_h": round(self.half_life, 2),
            "absorption_rate_ka_h": self.ka,
            "bioavailability_f": self.f,
            "tmax_h": tmax,
            "cmax_mg_l": cmax,
            "auc_inf_mg_h_l": round(self.auc_inf, 2),
        }
