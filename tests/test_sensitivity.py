import pytest
from src.one_compartment import OneCompartmentModel
from src.sensitivity import ParameterSensitivityAnalyzer


def test_sensitivity_sweep():
    model = OneCompartmentModel(dose=300.0, vd=35.0, cl=2.8, ka=1.2, f=0.95)
    df = ParameterSensitivityAnalyzer.run_sensitivity_sweep(model)

    assert len(df) > 0
    assert "cmax_change_percent" in df.columns
    assert "auc_change_percent" in df.columns

    # Increasing clearance must decrease AUC
    cl_pos = df[(df["parameter"] == "CL") & (df["perturbation_percent"] > 0)]
    assert cl_pos.iloc[0]["auc_change_percent"] < 0
