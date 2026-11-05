import pytest
import numpy as np
from src.two_compartment import TwoCompartmentModel


def test_two_compartment_simulation():
    # Gentamicin IV: Dose = 120 mg, Vc = 12 L, Vp = 15 L, CL = 4.5 L/h, Q = 2.0 L/h
    model = TwoCompartmentModel(dose=120.0, vc=12.0, vp=15.0, cl=4.5, q=2.0)

    t_eval = np.linspace(0.0, 12.0, 50)
    c_central, c_periph = model.simulate_iv_ode(t_eval)

    # Initial central concentration at t=0 must equal Dose / Vc = 120 / 12 = 10 mg/L
    assert pytest.approx(c_central[0], abs=0.1) == 10.0
    # Peripheral concentration at t=0 must equal 0
    assert pytest.approx(c_periph[0], abs=0.1) == 0.0
    # At t=12, concentrations should decline
    assert c_central[-1] < c_central[0]
