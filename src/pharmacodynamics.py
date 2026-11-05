from typing import Union
import numpy as np


class PharmacodynamicModel:
    """Implements Sigmoidal Emax (Hill) and Inhibitory (Imax) pharmacodynamic models."""

    def __init__(
        self,
        e0: float = 0.0,
        emax: float = 100.0,
        ec50: float = 10.0,
        gamma: float = 1.0,  # Hill coefficient (sigmoidicity factor)
    ):
        self.e0 = float(e0)
        self.emax = float(emax)
        self.ec50 = float(ec50)
        self.gamma = float(gamma)

    def calculate_effect(self, concentration: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """E(C) = E0 + [Emax * C^gamma] / [EC50^gamma + C^gamma]"""
        c = np.maximum(0.0, np.array(concentration))
        c_pow = np.power(c, self.gamma)
        ec50_pow = np.power(self.ec50, self.gamma)

        effect = self.e0 + ((self.emax * c_pow) / (ec50_pow + c_pow))
        return float(effect) if np.isscalar(concentration) else effect

    def calculate_inhibitory_effect(self, concentration: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """E(C) = E0 - [Imax * C^gamma] / [IC50^gamma + C^gamma]"""
        c = np.maximum(0.0, np.array(concentration))
        c_pow = np.power(c, self.gamma)
        ic50_pow = np.power(self.ec50, self.gamma)

        effect = self.e0 - ((self.emax * c_pow) / (ic50_pow + c_pow))
        return float(effect) if np.isscalar(concentration) else effect
