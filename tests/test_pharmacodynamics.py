import pytest
from src.pharmacodynamics import PharmacodynamicModel


def test_hill_emax_model():
    pd_model = PharmacodynamicModel(e0=0.0, emax=100.0, ec50=10.0, gamma=1.0)

    # At C = 0, effect must be E0 = 0
    assert pd_model.calculate_effect(0.0) == 0.0

    # At C = EC50 = 10, effect must be exactly 50% of Emax = 50.0
    assert pytest.approx(pd_model.calculate_effect(10.0), abs=1e-5) == 50.0

    # At very high concentration, effect approaches Emax = 100.0
    assert pytest.approx(pd_model.calculate_effect(10000.0), abs=0.5) == 100.0
