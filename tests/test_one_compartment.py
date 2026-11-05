import pytest
import numpy as np
from src.one_compartment import OneCompartmentModel


def test_analytical_vs_ode_simulation():
    # Dose = 300 mg, Vd = 35 L, CL = 2.8 L/h, ka = 1.2 h^-1, F = 0.95
    model = OneCompartmentModel(dose=300.0, vd=35.0, cl=2.8, ka=1.2, f=0.95)

    t_eval = np.linspace(0.0, 24.0, 100)
    c_analytical = model.simulate_oral_analytical(t_eval)
    c_ode = model.simulate_oral_ode(t_eval)

    # Analytical and numerical ODE solutions must be nearly identical (within 1%)
    np.testing.assert_allclose(c_analytical, c_ode, rtol=1e-2, atol=1e-3)


def test_tmax_cmax_calculation():
    model = OneCompartmentModel(dose=300.0, vd=35.0, cl=2.8, ka=1.2, f=0.95)
    tmax, cmax = model.calculate_tmax_cmax()

    assert 2.0 < tmax < 3.0
    assert 5.0 < cmax < 8.0
    assert model.half_life > 8.0
    assert model.auc_inf > 90.0
