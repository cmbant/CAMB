from __future__ import annotations

import math
import numbers
import os

from . import model
from .baseconfig import CAMB_Structure, CAMBValueError
from .inifile import IniFile

_initial_condition_names = [
    "initial_vector",
    "initial_adiabatic",
    "initial_iso_CDM",
    "initial_iso_baryon",
    "initial_iso_neutrino",
    "initial_iso_neutrino_vel",
]
_massive_nu_method_names = ["Nu_int", "Nu_trunc", "Nu_approx", "Nu_best"]
_roundtrip_float_tolerance_paths = {"params.Transfer.kmax"}


class CambIniFile(IniFile):
    def set(self, key: str, value) -> None:
        if key not in self.params:
            self.readOrder.append(key)
        self.params[key] = format_value(value)

    def set_sequence(self, key: str, values) -> None:
        self.set(key, " ".join(format_value(value) for value in values))

    def write_fields(
        self,
        obj,
        *,
        names: tuple[str, ...] | list[str] | None = None,
        rename: dict[str, str] | None = None,
    ) -> None:
        rename = rename or {}
        if names is None:
            names = [name for name, _ in obj.get_all_fields() if not name.startswith("_")]
        for name in names:
            self.set(rename.get(name, name), getattr(obj, name))


def format_value(value) -> str:
    if isinstance(value, bool):
        return "T" if value else "F"
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return repr(value)
    return str(value)


def _nonlinear_mode_value(mode_name: str) -> int:
    return model.NonLinear_names.index(mode_name)


def _massive_nu_method_value(mode_name: str) -> int:
    return _massive_nu_method_names.index(mode_name)


def _scalar_initial_condition_value(condition_name: str) -> int:
    return _initial_condition_names.index(condition_name)


def _update_source_terms_ini(params: model.CAMBparams, state: CambIniFile) -> None:
    state.set("limber_windows", params.SourceTerms.limber_windows)
    state.set("limber_phiphi", params.SourceTerms.limber_phi_lmin)
    state.set("Do21cm", params.Do21cm)
    state.set("counts_density", params.SourceTerms.counts_density)
    state.set("counts_redshift", params.SourceTerms.counts_redshift)
    state.set("DoRedshiftLensing", params.SourceTerms.counts_lensing)
    state.set("counts_velocity", params.SourceTerms.counts_velocity)
    state.set("counts_radial", params.SourceTerms.counts_radial)
    state.set("counts_timedelay", params.SourceTerms.counts_timedelay)
    state.set("counts_ISW", params.SourceTerms.counts_ISW)
    state.set("counts_potential", params.SourceTerms.counts_potential)
    state.set("counts_evolve", params.SourceTerms.counts_evolve)
    state.set("line_basic", params.SourceTerms.line_basic)
    state.set("line_distortions", params.SourceTerms.line_distortions)
    state.set("line_extra", params.SourceTerms.line_extra)
    state.set("line_phot_dipole", params.SourceTerms.line_phot_dipole)
    state.set("line_phot_quadrupole", params.SourceTerms.line_phot_quadrupole)
    state.set("line_reionization", params.SourceTerms.line_reionization)
    state.set("use_mK", params.SourceTerms.use_21cm_mK)
    state.set("Kmax_Boost", params.Accuracy.KmaxBoost)

    if not params.SourceWindows:
        state.set("num_redshiftwindows", 0)
        return

    state.set("num_redshiftwindows", len(params.SourceWindows))
    for index, window in enumerate(params.SourceWindows, start=1):
        window.write_ini(state, index)


def _update_ini_state_from_params(params: model.CAMBparams, state: CambIniFile) -> None:
    state.set("get_scalar_cls", params.WantCls and params.WantScalars)
    state.set("get_vector_cls", params.WantCls and params.WantVectors)
    state.set("get_tensor_cls", params.WantCls and params.WantTensors)
    state.set("want_CMB", params.Want_CMB)
    state.set("want_CMB_lensing", params.Want_CMB_lensing)
    state.set("get_transfer", params.WantTransfer)
    state.set("do_nonlinear", _nonlinear_mode_value(params.NonLinear))
    state.set("evolve_baryon_cs", params.Evolve_baryon_cs)
    state.set("evolve_delta_xe", params.Evolve_delta_xe)
    state.set("evolve_delta_ts", params.Evolve_delta_Ts)
    state.set("l_min", params.min_l)
    state.set("use_physical", True)
    state.set("hubble", params.H0)
    state.set("ombh2", params.ombh2)
    state.set("omch2", params.omch2)
    state.set("omnuh2", params.omnuh2)
    state.set("omk", params.omk)
    state.set("temp_cmb", params.TCMB)
    state.set("helium_fraction", params.YHe)
    state.set("massless_neutrinos", params.num_nu_massless)
    state.set("nu_mass_eigenstates", params.nu_mass_eigenstates)
    if params.nu_mass_eigenstates:
        state.set_sequence("massive_neutrinos", params.nu_mass_numbers[: params.nu_mass_eigenstates])
    else:
        state.set("massive_neutrinos", 0)
    state.set("share_delta_neff", params.share_delta_neff)
    if params.num_nu_massive > 0:
        if not params.share_delta_neff:
            state.set_sequence("nu_mass_degeneracies", params.nu_mass_degeneracies[: params.nu_mass_eigenstates])
        state.set_sequence("nu_mass_fractions", params.nu_mass_fractions[: params.nu_mass_eigenstates])
    state.set("Alens", params.Alens)
    state.set("derived_parameters", params.WantDerivedParameters)
    state.set("accurate_polarization", params.Accuracy.AccuratePolarization)
    state.set("accurate_reionization", params.Accuracy.AccurateReionization)
    state.set("accurate_BB", params.Accuracy.AccurateBB)
    state.set("accuracy_boost", params.Accuracy.AccuracyBoost)
    state.set("l_accuracy_boost", params.Accuracy.lAccuracyBoost)
    state.set("l_sample_boost", params.Accuracy.lSampleBoost)
    state.set("min_l_logl_sampling", params.min_l_logl_sampling)
    state.set("do_late_rad_truncation", params.DoLateRadTruncation)
    state.set("massive_nu_approx", _massive_nu_method_value(params.MassiveNuMethod))

    if params.WantCls and (params.WantScalars or params.WantVectors):
        state.set("l_max_scalar", params.max_l)
        state.set("k_eta_max_scalar", params.max_eta_k)
        state.set("lens_output_margin", params.lens_output_margin)
        if params.WantScalars:
            state.set("do_lensing", params.DoLensing)
    if params.WantCls and params.WantTensors:
        state.set("l_max_tensor", params.max_l_tensor)
        state.set("k_eta_max_tensor", params.max_eta_k_tensor)

    params.DarkEnergy.write_ini(state)
    params.Reion.write_ini(state)
    params.InitPower.write_ini(state)
    params.Recomb.write_ini(state)
    _update_source_terms_ini(params, state)

    if params.WantTransfer:
        state.set("transfer_high_precision", params.Transfer.high_precision)
        state.set("accurate_massive_neutrino_transfers", params.Transfer.accurate_massive_neutrinos)
        state.set("transfer_kmax", params.Transfer.kmax / (params.H0 / 100.0))
        state.set("transfer_k_per_logint", params.Transfer.k_per_logint)
        state.set("transfer_num_redshifts", params.Transfer.PK_num_redshifts)
        for index, redshift in enumerate(params.Transfer.PK_redshifts[: params.Transfer.PK_num_redshifts], start=1):
            state.set(f"transfer_redshift({index})", redshift)
        state.set("transfer_21cm_cl", params.transfer_21cm_cl)

    try:
        params.NonLinearModel.write_ini(state)
    except CAMBValueError:
        if params.NonLinear != model.NonLinear_none:
            raise

    state.set("initial_condition", _scalar_initial_condition_value(params.scalar_initial_condition))
    if params.scalar_initial_condition == "initial_vector":
        state.set_sequence("initial_vector", params.InitialConditionVector)
    state.set("use_cl_spline_template", params.use_cl_spline_template)


def _roundtrip_expected_params(params: model.CAMBparams) -> model.CAMBparams:
    expected = params.copy()
    if not expected.WantCls:
        defaults = model.CAMBparams()
        expected.WantScalars = False
        expected.WantVectors = False
        expected.WantTensors = False
        expected.DoLensing = False
        expected.max_l = defaults.max_l
        expected.max_eta_k = defaults.max_eta_k
        expected.lens_output_margin = defaults.lens_output_margin
        expected.InitPower.r = 0
        expected.InitPower.nt = 0
        expected.InitPower.ntrun = 0
    return expected


def _roundtrip_mismatch(actual, expected, path: str = "params") -> str | None:
    if isinstance(actual, dict) and isinstance(expected, dict):
        if actual.keys() != expected.keys():
            return f"{path} keys differ: {actual.keys()} != {expected.keys()}"
        for key in actual:
            if mismatch := _roundtrip_mismatch(actual[key], expected[key], f"{path}.{key}"):
                return mismatch
        return None
    if isinstance(actual, (list, tuple)) and isinstance(expected, (list, tuple)):
        if len(actual) != len(expected):
            return f"{path} lengths differ: {len(actual)} != {len(expected)}"
        for index, (actual_value, expected_value) in enumerate(zip(actual, expected, strict=True)):
            if mismatch := _roundtrip_mismatch(actual_value, expected_value, f"{path}[{index}]"):
                return mismatch
        return None
    if isinstance(actual, CAMB_Structure) and isinstance(expected, CAMB_Structure):
        return _roundtrip_mismatch(actual.__getstate__(), expected.__getstate__(), path)
    if (
        isinstance(actual, numbers.Real)
        and isinstance(expected, numbers.Real)
        and not isinstance(actual, (bool, numbers.Integral))
        and not isinstance(expected, (bool, numbers.Integral))
    ):
        if path in _roundtrip_float_tolerance_paths:
            if not math.isclose(float(actual), float(expected), rel_tol=1e-12, abs_tol=1e-14):
                return f"{path} differs: {actual!r} != {expected!r}"
            return None
        if actual != expected:
            return f"{path} differs: {actual!r} != {expected!r}"
        return None
    if actual != expected:
        return f"{path} differs: {actual!r} != {expected!r}"
    return None


def write_ini(
    params: model.CAMBparams,
    ini_filename,
    *,
    validate: bool = True,
) -> None:
    if not isinstance(params, model.CAMBparams):
        raise TypeError("params must be an instance of CAMBparams")

    ini_path = os.fspath(ini_filename)
    state = CambIniFile()
    _update_ini_state_from_params(params, state)
    state.saveFile(ini_path)
    if validate:
        from .camb import read_ini

        reparsed = read_ini(ini_path)
        expected = _roundtrip_expected_params(params)
        if repr(reparsed) == repr(expected):
            return
        if mismatch := _roundtrip_mismatch(reparsed, expected):
            raise CAMBValueError(f"Saved ini did not round-trip via read_ini ({ini_path}): {mismatch}")
