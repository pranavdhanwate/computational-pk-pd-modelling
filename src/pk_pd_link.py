from typing import Tuple
import numpy as np
from scipy.integrate import solve_ivp
from .pharmacodynamics import PharmacodynamicModel


class PKPDLink:
    """Links plasma pharmacokinetic concentrations to pharmacodynamic response via direct or effect-compartment."""

    def __init__(self, pd_model: PharmacodynamicModel, ke0: float = 0.5):
        self.pd = pd_model
        self.ke0 = float(ke0)  # Effect-site equilibration rate constant (h^-1)

    def simulate_direct_link(self, cp: np.ndarray) -> np.ndarray:
        """Simulates instantaneous pharmacological response directly driven by plasma concentration."""
        return self.pd.calculate_effect(cp)

    def simulate_effect_compartment(self, t: np.ndarray, cp: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Simulates hypothetical effect compartment concentration Ce(t) and resulting delayed PD response:
        dCe/dt = ke0 * (Cp(t) - Ce(t))
        """
        from scipy.interpolate import interp1d
        cp_func = interp1d(t, cp, kind="linear", fill_value="extrapolate")

        def ode(time, ce):
            cp_val = float(cp_func(time))
            return [self.ke0 * (cp_val - ce[0])]

        sol = solve_ivp(ode, [t[0], t[-1]], [0.0], t_eval=t, rtol=1e-5, atol=1e-8)
        ce = sol.y[0]
        effect = self.pd.calculate_effect(ce)
        return ce, effect
