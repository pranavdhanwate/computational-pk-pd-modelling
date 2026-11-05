import pytest
from src.one_compartment import OneCompartmentModel
from src.multiple_dosing import MultipleDosingSimulator


def test_accumulation_factor_and_steady_state():
    model = OneCompartmentModel(dose=300.0, vd=35.0, cl=2.8, ka=1.2, f=0.95)
    simulator = MultipleDosingSimulator(model, tau=12.0)

    metrics = simulator.calculate_steady_state_metrics()

    # R must be > 1.0 (accumulation occurs)
    assert metrics["accumulation_factor_R"] > 1.0
    assert metrics["cmax_ss_mg_l"] > metrics["cmin_ss_mg_l"]
    assert metrics["cav_ss_mg_l"] > 0.0
