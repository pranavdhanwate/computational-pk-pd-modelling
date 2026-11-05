import math
from typing import Dict, Any, Tuple
import numpy as np
from scipy.integrate import solve_ivp


class TwoCompartmentModel:
    """Two-compartment pharmacokinetic model with central and peripheral distribution."""

    def __init__(
        self,
        dose: float,
        vc: float,
        vp: float,
        cl: float,
        q: float = 2.0,  # Inter-compartmental clearance (L/h)
    ):
        self.dose = float(dose)
        self.vc = float(vc)
        self.vp = float(vp)
        self.cl = float(cl)
        self.q = float(q)

        # Micro-constants (h^-1)
        self.k10 = self.cl / self.vc
        self.k12 = self.q / self.vc
        self.k21 = self.q / self.vp

        # Hybrid macro-constants (alpha and beta)
        sum_k = self.k10 + self.k12 + self.k21
        prod_k = self.k10 * self.k21
        disc = math.sqrt(max(0.0, sum_k**2 - (4.0 * prod_k)))

        self.alpha = (sum_k + disc) / 2.0
        self.beta = (sum_k - disc) / 2.0

        # Coefficients for IV bolus
        denom = self.vc * (self.alpha - self.beta)
        self.A = (self.dose * (self.alpha - self.k21)) / denom if denom > 0 else 0.0
        self.B = (self.dose * (self.k21 - self.beta)) / denom if denom > 0 else 0.0

        self.t_half_alpha = math.log(2.0) / self.alpha if self.alpha > 0 else 0.0
        self.t_half_beta = math.log(2.0) / self.beta if self.beta > 0 else 0.0

    def simulate_iv_analytical(self, t: np.ndarray) -> np.ndarray:
        """C(t) = A * exp(-alpha * t) + B * exp(-beta * t)"""
        return (self.A * np.exp(-self.alpha * t)) + (self.B * np.exp(-self.beta * t))

    def simulate_iv_ode(self, t_eval: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Numerically integrates dAc/dt and dAp/dt using solve_ivp."""
        def odes(t, y):
            ac, ap = y
            dac_dt = -(self.k10 + self.k12) * ac + (self.k21 * ap)
            dap_dt = (self.k12 * ac) - (self.k21 * ap)
            return [dac_dt, dap_dt]

        y0 = [self.dose, 0.0]
        sol = solve_ivp(odes, [t_eval[0], t_eval[-1]], y0, t_eval=t_eval, rtol=1e-6, atol=1e-9)
        c_central = sol.y[0] / self.vc
        c_peripheral = sol.y[1] / self.vp
        return c_central, c_peripheral
