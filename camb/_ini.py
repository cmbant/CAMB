from __future__ import annotations

import os

from . import constants, dark_energy, initialpower, model, nonlinear, recombination, reionization, sources
from .baseconfig import CAMBValueError
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
_recfast_hswitch_offset = 1.14 - (1.105 + 0.02)


class CambIniFile(IniFile):
    pass


def format_value(value) -> str:
    if isinstance(value, bool):
        return "T" if value else "F"
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.17g}"
    return str(value)


def set_param(state: CambIniFile, key: str, value) -> None:
    if key not in state.params:
        state.readOrder.append(key)
    state.params[key] = format_value(value)


def _set_ini_sequence(state: CambIniFile, key: str, values) -> None:
    set_param(state, key, " ".join(format_value(value) for value in values))


def _nonlinear_mode_value(mode_name: str) -> int:
    return model.NonLinear_names.index(mode_name)


def _massive_nu_method_value(mode_name: str) -> int:
    return _massive_nu_method_names.index(mode_name)


def _scalar_initial_condition_value(condition_name: str) -> int:
    return _initial_condition_names.index(condition_name)


def _update_dark_energy_ini(params: model.CAMBparams, state: CambIniFile) -> None:
    dark_energy_model = params.DarkEnergy
    if isinstance(dark_energy_model, dark_energy.DarkEnergyFluid):
        set_param(state, "dark_energy_model", "fluid")
        if dark_energy_model.use_tabulated_w:
            raise CAMBValueError("write_ini cannot serialize tabulated dark energy")
        else:
            set_param(state, "w", dark_energy_model.w)
            set_param(state, "wa", dark_energy_model.wa)
        set_param(state, "cs2_lam", dark_energy_model.cs2)
        return

    if isinstance(dark_energy_model, dark_energy.DarkEnergyPPF):
        set_param(state, "dark_energy_model", "ppf")
        if dark_energy_model.use_tabulated_w:
            raise CAMBValueError("write_ini cannot serialize tabulated dark energy")
        else:
            set_param(state, "w", dark_energy_model.w)
            set_param(state, "wa", dark_energy_model.wa)
        set_param(state, "cs2_lam", dark_energy_model.cs2)
        return

    if isinstance(dark_energy_model, dark_energy.AxionEffectiveFluid):
        set_param(state, "dark_energy_model", "AxionEffectiveFluid")
        set_param(state, "AxionEffectiveFluid_w_n", dark_energy_model.w_n)
        set_param(state, "AxionEffectiveFluid_fde_zc", dark_energy_model.fde_zc)
        set_param(state, "AxionEffectiveFluid_zc", dark_energy_model.zc)
        set_param(state, "AxionEffectiveFluid_theta_i", dark_energy_model.theta_i)
        return

    if isinstance(dark_energy_model, dark_energy.EarlyQuintessence):
        raise CAMBValueError("write_ini cannot serialize EarlyQuintessence")

    raise CAMBValueError(f"write_ini does not support dark energy class {dark_energy_model.__class__.__name__}")


def _update_reionization_ini(params: model.CAMBparams, state: CambIniFile) -> None:
    reion_model = params.Reion
    set_param(state, "reionization", reion_model.Reionization)
    if isinstance(reion_model, reionization.BaseTauWithHeReionization):
        set_param(state, "re_use_optical_depth", reion_model.use_optical_depth)
        if reion_model.use_optical_depth:
            set_param(state, "re_optical_depth", reion_model.optical_depth)
        else:
            set_param(state, "re_redshift", reion_model.redshift)
        set_param(state, "re_ionization_frac", reion_model.fraction)
        set_param(state, "include_helium_fullreion", reion_model.include_helium_fullreion)
        set_param(state, "re_helium_redshift", reion_model.helium_redshift)
        set_param(state, "re_helium_delta_redshift", reion_model.helium_delta_redshift)
        set_param(state, "re_helium_redshiftstart", reion_model.helium_redshiftstart)
        set_param(state, "max_zrei", reion_model.max_redshift)
    if isinstance(reion_model, reionization.TanhReionization):
        set_param(state, "re_delta_redshift", reion_model.delta_redshift)
    elif isinstance(reion_model, reionization.ExpReionization):
        set_param(state, "reion_redshift_complete", reion_model.reion_redshift_complete)
        set_param(state, "reion_exp_power", reion_model.reion_exp_power)
        set_param(state, "reion_exp_smooth_width", reion_model.reion_exp_smooth_width)


def _update_initial_power_ini(params: model.CAMBparams, state: CambIniFile) -> None:
    initial_power_model = params.InitPower
    if isinstance(initial_power_model, initialpower.InitialPowerLaw):
        set_param(state, "pivot_scalar", initial_power_model.pivot_scalar)
        set_param(state, "pivot_tensor", initial_power_model.pivot_tensor)
        set_param(state, "initial_power_num", 1)
        set_param(state, "scalar_amp(1)", initial_power_model.As)
        set_param(state, "scalar_spectral_index(1)", initial_power_model.ns)
        set_param(state, "scalar_nrun(1)", initial_power_model.nrun)
        set_param(state, "scalar_nrunrun(1)", initial_power_model.nrunrun)
        set_param(
            state,
            "tensor_parameterization",
            initialpower.tensor_parameterization_names.index(initial_power_model.tensor_parameterization) + 1,
        )
        set_param(state, "tensor_spectral_index(1)", initial_power_model.nt)
        set_param(state, "tensor_nrun(1)", initial_power_model.ntrun)
        if initial_power_model.tensor_parameterization == "tensor_param_AT":
            set_param(state, "tensor_amp(1)", initial_power_model.At)
        else:
            set_param(state, "initial_ratio(1)", initial_power_model.r)
        return

    raise CAMBValueError(f"write_ini does not support initial power class {initial_power_model.__class__.__name__}")


def _update_recombination_ini(params: model.CAMBparams, state: CambIniFile) -> None:
    recomb_model = params.Recomb
    if isinstance(recomb_model, recombination.Recfast):
        set_param(state, "recombination_model", "Recfast")
        recfast_fudge = recomb_model.RECFAST_fudge
        if recomb_model.RECFAST_Hswitch:
            recfast_fudge += _recfast_hswitch_offset
        set_param(state, "RECFAST_fudge", recfast_fudge)
        set_param(state, "RECFAST_fudge_He", recomb_model.RECFAST_fudge_He)
        set_param(state, "RECFAST_Heswitch", recomb_model.RECFAST_Heswitch)
        set_param(state, "RECFAST_Hswitch", recomb_model.RECFAST_Hswitch)
        set_param(state, "AGauss1", recomb_model.AGauss1)
        set_param(state, "AGauss2", recomb_model.AGauss2)
        set_param(state, "zGauss1", recomb_model.zGauss1)
        set_param(state, "zGauss2", recomb_model.zGauss2)
        set_param(state, "wGauss1", recomb_model.wGauss1)
        set_param(state, "wGauss2", recomb_model.wGauss2)
        set_param(state, "RECFAST_nz", recomb_model.Nz)
    elif isinstance(recomb_model, recombination.CosmoRec):
        set_param(state, "recombination_model", "CosmoRec")
        set_param(state, "cosmorec_runmode", recomb_model.runmode)
        set_param(state, "cosmorec_accuracy", recomb_model.accuracy)
        set_param(state, "cosmorec_fdm", recomb_model.fdm)
    elif isinstance(recomb_model, recombination.HyRec):
        set_param(state, "recombination_model", "HyRec")


def _update_source_terms_ini(params: model.CAMBparams, state: CambIniFile) -> None:
    set_param(state, "limber_windows", params.SourceTerms.limber_windows)
    set_param(state, "limber_phiphi", params.SourceTerms.limber_phi_lmin)
    set_param(state, "Do21cm", params.Do21cm)
    set_param(state, "counts_density", params.SourceTerms.counts_density)
    set_param(state, "counts_redshift", params.SourceTerms.counts_redshift)
    set_param(state, "DoRedshiftLensing", params.SourceTerms.counts_lensing)
    set_param(state, "counts_velocity", params.SourceTerms.counts_velocity)
    set_param(state, "counts_radial", params.SourceTerms.counts_radial)
    set_param(state, "counts_timedelay", params.SourceTerms.counts_timedelay)
    set_param(state, "counts_ISW", params.SourceTerms.counts_ISW)
    set_param(state, "counts_potential", params.SourceTerms.counts_potential)
    set_param(state, "counts_evolve", params.SourceTerms.counts_evolve)
    set_param(state, "line_basic", params.SourceTerms.line_basic)
    set_param(state, "line_distortions", params.SourceTerms.line_distortions)
    set_param(state, "line_extra", params.SourceTerms.line_extra)
    set_param(state, "line_phot_dipole", params.SourceTerms.line_phot_dipole)
    set_param(state, "line_phot_quadrupole", params.SourceTerms.line_phot_quadrupole)
    set_param(state, "line_reionization", params.SourceTerms.line_reionization)
    set_param(state, "use_mK", params.SourceTerms.use_21cm_mK)
    set_param(state, "Kmax_Boost", params.Accuracy.KmaxBoost)

    if not params.SourceWindows:
        set_param(state, "num_redshiftwindows", 0)
        return

    if any(not isinstance(window, sources.GaussianSourceWindow) for window in params.SourceWindows):
        raise CAMBValueError("write_ini only supports GaussianSourceWindow")

    set_param(state, "num_redshiftwindows", len(params.SourceWindows))
    for index, window in enumerate(params.SourceWindows, start=1):
        set_param(state, f"redshift({index})", window.redshift)
        set_param(state, f"redshift_kind({index})", window.source_type)
        if window.source_type == "21cm":
            set_param(state, f"redshift_sigma_Mhz({index})", window.sigma * constants.f_21cm / 1e6)
        else:
            set_param(state, f"redshift_sigma({index})", window.sigma)
        if window.source_type == "counts":
            set_param(state, f"redshift_bias({index})", window.bias)
            set_param(state, f"redshift_dlog10Ndm({index})", window.dlog10Ndm)


def _update_ini_state_from_params(params: model.CAMBparams, state: CambIniFile) -> None:
    set_param(state, "get_scalar_cls", params.WantScalars)
    set_param(state, "get_vector_cls", params.WantVectors)
    set_param(state, "get_tensor_cls", params.WantTensors)
    set_param(state, "want_CMB", params.Want_CMB)
    set_param(state, "want_CMB_lensing", params.Want_CMB_lensing)
    set_param(state, "get_transfer", params.WantTransfer)
    set_param(state, "do_nonlinear", _nonlinear_mode_value(params.NonLinear))
    set_param(state, "evolve_baryon_cs", params.Evolve_baryon_cs)
    set_param(state, "evolve_delta_xe", params.Evolve_delta_xe)
    set_param(state, "evolve_delta_ts", params.Evolve_delta_Ts)
    set_param(state, "l_min", params.min_l)
    set_param(state, "use_physical", True)
    set_param(state, "hubble", params.H0)
    set_param(state, "ombh2", params.ombh2)
    set_param(state, "omch2", params.omch2)
    set_param(state, "omnuh2", params.omnuh2)
    set_param(state, "omk", params.omk)
    set_param(state, "temp_cmb", params.TCMB)
    set_param(state, "helium_fraction", params.YHe)
    set_param(state, "massless_neutrinos", params.num_nu_massless)
    set_param(state, "nu_mass_eigenstates", params.nu_mass_eigenstates)
    if params.nu_mass_eigenstates:
        _set_ini_sequence(state, "massive_neutrinos", params.nu_mass_numbers[: params.nu_mass_eigenstates])
    else:
        set_param(state, "massive_neutrinos", 0)
    set_param(state, "share_delta_neff", params.share_delta_neff)
    if params.num_nu_massive > 0:
        if not params.share_delta_neff:
            _set_ini_sequence(state, "nu_mass_degeneracies", params.nu_mass_degeneracies[: params.nu_mass_eigenstates])
        _set_ini_sequence(state, "nu_mass_fractions", params.nu_mass_fractions[: params.nu_mass_eigenstates])
    set_param(state, "Alens", params.Alens)
    set_param(state, "derived_parameters", params.WantDerivedParameters)
    set_param(state, "accurate_polarization", params.Accuracy.AccuratePolarization)
    set_param(state, "accurate_reionization", params.Accuracy.AccurateReionization)
    set_param(state, "accurate_BB", params.Accuracy.AccurateBB)
    set_param(state, "accuracy_boost", params.Accuracy.AccuracyBoost)
    set_param(state, "l_accuracy_boost", params.Accuracy.lAccuracyBoost)
    set_param(state, "l_sample_boost", params.Accuracy.lSampleBoost)
    set_param(state, "do_late_rad_truncation", params.DoLateRadTruncation)
    set_param(state, "massive_nu_approx", _massive_nu_method_value(params.MassiveNuMethod))

    if params.WantCls and (params.WantScalars or params.WantVectors):
        set_param(state, "l_max_scalar", params.max_l)
        set_param(state, "k_eta_max_scalar", params.max_eta_k)
        if params.WantScalars:
            set_param(state, "do_lensing", params.DoLensing)
    if params.WantCls and params.WantTensors:
        set_param(state, "l_max_tensor", params.max_l_tensor)
        set_param(state, "k_eta_max_tensor", params.max_eta_k_tensor)

    _update_dark_energy_ini(params, state)
    _update_reionization_ini(params, state)
    _update_initial_power_ini(params, state)
    _update_recombination_ini(params, state)
    _update_source_terms_ini(params, state)

    if params.WantTransfer:
        set_param(state, "transfer_high_precision", params.Transfer.high_precision)
        set_param(state, "accurate_massive_neutrino_transfers", params.Transfer.accurate_massive_neutrinos)
        set_param(state, "transfer_kmax", params.Transfer.kmax / (params.H0 / 100.0))
        set_param(state, "transfer_k_per_logint", params.Transfer.k_per_logint)
        set_param(state, "transfer_num_redshifts", params.Transfer.PK_num_redshifts)
        for index, redshift in enumerate(params.Transfer.PK_redshifts[: params.Transfer.PK_num_redshifts], start=1):
            set_param(state, f"transfer_redshift({index})", redshift)
        set_param(state, "transfer_21cm_cl", params.transfer_21cm_cl)

    if isinstance(params.NonLinearModel, nonlinear.Halofit):
        set_param(state, "nonlinear_model", "Halofit")
        set_param(state, "halofit_version", nonlinear.halofit_version_names[params.NonLinearModel.halofit_version])
        set_param(state, "HMCode_A_baryon", params.NonLinearModel.HMCode_A_baryon)
        set_param(state, "HMCode_eta_baryon", params.NonLinearModel.HMCode_eta_baryon)
        set_param(state, "HMcode_logT_AGN", params.NonLinearModel.HMCode_logT_AGN)
    elif isinstance(params.NonLinearModel, nonlinear.SPkNonLinear):
        base_model = params.NonLinearModel.BaseModel
        if not isinstance(base_model, nonlinear.Halofit):
            raise CAMBValueError(
                f"write_ini does not support SPkNonLinear base model class {base_model.__class__.__name__}"
            )

        set_param(state, "nonlinear_model", "SPkNonLinear")
        set_param(state, "halofit_version", nonlinear.halofit_version_names[base_model.halofit_version])
        set_param(state, "HMCode_A_baryon", base_model.HMCode_A_baryon)
        set_param(state, "HMCode_eta_baryon", base_model.HMCode_eta_baryon)
        set_param(state, "HMcode_logT_AGN", base_model.HMCode_logT_AGN)
        set_param(state, "SPk_feedback", params.NonLinearModel.SPk_feedback)
        set_param(state, "SPk_SO", params.NonLinearModel.SPk_SO)
        set_param(state, "SPk_relation_kind", params.NonLinearModel.SPk_relation_kind)
        set_param(state, "SPk_fb_a", params.NonLinearModel.SPk_fb_a)
        set_param(state, "SPk_fb_pow", params.NonLinearModel.SPk_fb_pow)
        set_param(state, "SPk_fb_pivot", params.NonLinearModel.SPk_fb_pivot)
        set_param(state, "SPk_alpha", params.NonLinearModel.SPk_alpha)
        set_param(state, "SPk_beta", params.NonLinearModel.SPk_beta)
        set_param(state, "SPk_gamma", params.NonLinearModel.SPk_gamma)
        set_param(state, "SPk_epsilon", params.NonLinearModel.SPk_epsilon)
        set_param(state, "SPk_m_pivot", params.NonLinearModel.SPk_m_pivot)
    elif params.NonLinear != model.NonLinear_none:
        raise CAMBValueError(
            f"write_ini does not support non-linear model class {params.NonLinearModel.__class__.__name__}"
        )

    set_param(state, "initial_condition", _scalar_initial_condition_value(params.scalar_initial_condition))
    if params.scalar_initial_condition == "initial_vector":
        _set_ini_sequence(state, "initial_vector", params.InitialConditionVector)
    set_param(state, "use_cl_spline_template", params.use_cl_spline_template)


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
        if repr(reparsed) != repr(params):
            raise CAMBValueError(f"Saved ini did not round-trip via read_ini ({ini_path})")
