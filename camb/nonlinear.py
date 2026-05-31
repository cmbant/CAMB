import math
from ctypes import POINTER, byref, c_bool, c_double, c_int

import numpy as np
from numpy.ctypeslib import ndpointer

from .baseconfig import AllocatableObject, CAMBValueError, F2003Class, fortran_class, numpy_1d


class NonLinearModel(F2003Class):
    """
    Abstract base class for non-linear correction models.
    """

    _fields_ = (("Min_kh_nonlinear", c_double, "minimum k/h at which to apply non-linear corrections"),)


halofit_original = "original"
halofit_bird = "bird"
halofit_peacock = "peacock"
halofit_takahashi = "takahashi"
halofit_mead = "mead"
halofit_halomodel = "halomodel"
halofit_casarini = "casarini"
halofit_mead2015 = "mead2015"
halofit_mead2016 = "mead2016"
halofit_mead2020 = "mead2020"
halofit_mead2020_feedback = "mead2020_feedback"

halofit_default = halofit_mead2020

halofit_version_names = {
    halofit_original: 1,
    halofit_bird: 2,
    halofit_peacock: 3,
    halofit_takahashi: 4,
    halofit_mead: 5,
    halofit_halomodel: 6,
    halofit_casarini: 7,
    halofit_mead2015: 8,
    halofit_mead2016: 5,
    halofit_mead2020: 9,
    halofit_mead2020_feedback: 10,
}


@fortran_class
class Halofit(NonLinearModel):
    """
    Various specific approximate non-linear correction models based on HaloFit.
    """

    _fields_ = (
        ("halofit_version", c_int, {"names": halofit_version_names}),
        ("HMCode_A_baryon", c_double, "HMcode parameter A_baryon"),
        ("HMCode_eta_baryon", c_double, "HMcode parameter eta_baryon"),
        ("HMCode_logT_AGN", c_double, "HMcode parameter log10(T_AGN/K)"),
    )

    _fortran_class_module_ = "NonLinear"
    _fortran_class_name_ = "THalofit"

    def get_halofit_version(self):
        return self.halofit_version

    def set_params(
        self, halofit_version=halofit_default, HMCode_A_baryon=3.13, HMCode_eta_baryon=0.603, HMCode_logT_AGN=7.8
    ):
        """
        Set the halofit model for non-linear corrections.

        :param halofit_version: One of

            - original: `astro-ph/0207664 <https://arxiv.org/abs/astro-ph/0207664>`_
            - bird: `arXiv:1109.4416 <https://arxiv.org/abs/1109.4416>`_
            - peacock: `Peacock fit <http://www.roe.ac.uk/~jap/haloes/>`_
            - takahashi: `arXiv:1208.2701 <https://arxiv.org/abs/1208.2701>`_
            - mead: HMCode `arXiv:1602.02154 <https://arxiv.org/abs/1602.02154>`_
            - halomodel: basic halomodel
            - casarini: PKequal `arXiv:0810.0190 <https://arxiv.org/abs/0810.0190>`_, `arXiv:1601.07230 <https://arxiv.org/abs/1601.07230>`_
            - mead2015: original 2015 version of HMCode `arXiv:1505.07833 <https://arxiv.org/abs/1505.07833>`_
            - mead2016: Alias for 'mead'.
            - mead2020: 2020 version of HMcode `arXiv:2009.01858 <https://arxiv.org/abs/2009.01858>`_
            - mead2020_feedback: 2020 version of HMcode with baryonic feedback `arXiv:2009.01858 <https://arxiv.org/abs/2009.01858>`_
        :param HMCode_A_baryon: HMcode parameter A_baryon. Default 3.13.
                                Used only in models mead2015 and mead2016 (and its alias mead).
        :param HMCode_eta_baryon: HMcode parameter eta_baryon. Default 0.603.
                                Used only in mead2015 and mead2016 (and its alias mead).
        :param HMCode_logT_AGN: HMcode parameter logT_AGN. Default 7.8. Used only in model mead2020_feedback.
        """
        self.halofit_version = halofit_version
        self.HMCode_A_baryon = HMCode_A_baryon
        self.HMCode_eta_baryon = HMCode_eta_baryon
        self.HMCode_logT_AGN = HMCode_logT_AGN


@fortran_class
class SPkNonLinear(NonLinearModel):
    """SP(k) baryon suppression model applied on top of a base non-linear model."""

    _fields_ = (
        ("BaseModel", AllocatableObject(NonLinearModel)),
        ("SPk_feedback", c_bool, "Enable SP(k) suppression"),
        ("SPk_SO", c_int, "SP(k) spherical overdensity (200 or 500)"),
        (
            "SPk_relation_kind",
            c_int,
            "SP(k) relation kind: 1=power_law, 2=cosmo_power_law, 3=double_power_law",
        ),
        ("SPk_fb_a", c_double, "Power-law relation normalization"),
        ("SPk_fb_pow", c_double, "Power-law relation exponent"),
        ("SPk_fb_pivot", c_double, "Power-law relation pivot mass [M_sun]"),
        ("SPk_alpha", c_double, "Relation alpha parameter (kinds 2/3)"),
        ("SPk_beta", c_double, "Relation beta parameter (kinds 2/3)"),
        ("SPk_gamma", c_double, "Relation gamma parameter (kinds 2/3)"),
        ("SPk_epsilon", c_double, "Relation epsilon parameter (kind 3)"),
        ("SPk_m_pivot", c_double, "Relation pivot mass [M_sun] (kind 3)"),
    )

    _fortran_class_module_ = "SPkNonLinear"
    _fortran_class_name_ = "TSPkNonLinear"

    def __init__(self, **kwargs):
        super().__init__()
        self.BaseModel = Halofit()
        self.set_params(**kwargs)

    def _validate(self):
        if self.SPk_SO not in (200, 500):
            raise CAMBValueError("SPk_SO must be 200 or 500")
        if self.SPk_relation_kind not in (1, 2, 3):
            raise CAMBValueError(
                "SPk_relation_kind must be 1 (power_law), 2 (cosmo_power_law), or 3 (double_power_law)"
            )
        if self.SPk_relation_kind == 1 and self.SPk_fb_pivot <= 0:
            raise CAMBValueError("SPk_fb_pivot must be > 0 for power_law relation")
        if self.SPk_relation_kind == 3 and self.SPk_m_pivot <= 0:
            raise CAMBValueError("SPk_m_pivot must be > 0 for double_power_law relation")

        if self.SPk_feedback and isinstance(self.BaseModel, Halofit):
            if isinstance(self.BaseModel.halofit_version, str):
                halofit_version_int = halofit_version_names[self.BaseModel.halofit_version]
            else:
                halofit_version_int = int(self.BaseModel.halofit_version)
            if halofit_version_int == halofit_version_names[halofit_mead2020_feedback]:
                raise CAMBValueError(
                    "SP(k) is not compatible with halofit_version='mead2020_feedback'. "
                    "Use halofit_version='mead2020' (or another non-feedback option) when enabling SPk_feedback."
                )

            hmcode_2015_2016_versions = {
                halofit_version_names[halofit_mead],
                halofit_version_names[halofit_mead2015],
                halofit_version_names[halofit_mead2016],
            }
            if halofit_version_int in hmcode_2015_2016_versions and (
                (not math.isclose(self.BaseModel.HMCode_A_baryon, 3.13, rel_tol=0.0, abs_tol=1e-12))
                or (not math.isclose(self.BaseModel.HMCode_eta_baryon, 0.603, rel_tol=0.0, abs_tol=1e-12))
            ):
                raise CAMBValueError(
                    "SP(k) cannot be combined with HMCode_A_baryon/HMCode_eta_baryon baryonic corrections in HMCode 2015/2016"
                )

    def set_params(
        self,
        SPk_feedback=False,
        SPk_SO=200,
        SPk_relation_kind=1,
        SPk_fb_a=1.0,
        SPk_fb_pow=0.0,
        SPk_fb_pivot=1.0,
        SPk_alpha=0.0,
        SPk_beta=0.0,
        SPk_gamma=0.0,
        SPk_epsilon=0.0,
        SPk_m_pivot=1.0,
        halofit_version=halofit_default,
        HMCode_A_baryon=3.13,
        HMCode_eta_baryon=0.603,
        HMCode_logT_AGN=7.8,
    ):
        """
        Configure the SP(k) baryon suppression model.

        References:
          - SP(k) model: `MNRAS 523, 2247 (2023) <https://doi.org/10.1093/mnras/stad1474>`_
          - pyspk: https://github.com/jemme07/pyspk

        The base model is evaluated first (Halofit by default), then SP(k)
        suppression is applied to CAMB's non-linear ratio as:

        ``sqrt(P_NL/P_L) -> sqrt(P_NL/P_L) * sqrt(SPk_suppression)``

        **SP(k) relation kinds:**

        - **kind=1** (power_law):
          ``f_b / (Omega_b/Omega_m) = SPk_fb_a * (M / SPk_fb_pivot)^SPk_fb_pow``

        - **kind=2** (cosmo_power_law):
          ``f_b / (Omega_b/Omega_m) = (exp(SPk_alpha)/100) * (M_500c/1e14)^(SPk_beta - 1) * (E(z)/E(0.3))^SPk_gamma``

        - **kind=3** (double_power_law):
          ``f_b / (Omega_b/Omega_m) = 0.5 * SPk_epsilon * ((M/SPk_m_pivot)^SPk_alpha + (M/SPk_m_pivot)^SPk_beta) * (E(z)/E(0.3))^SPk_gamma``

        :param SPk_feedback: If True, apply SP(k) suppression on top of the base model.
        :param SPk_SO: Spherical overdensity calibration (200 or 500).
        :param SPk_relation_kind: Relation type: 1 (power_law), 2 (cosmo_power_law), 3 (double_power_law).
        :param SPk_fb_a: Power-law normalization (kind=1).
        :param SPk_fb_pow: Power-law exponent (kind=1).
        :param SPk_fb_pivot: Power-law pivot mass in M_sun (kind=1).
        :param SPk_alpha: Alpha parameter (kinds 2, 3).
        :param SPk_beta: Beta parameter (kinds 2, 3).
        :param SPk_gamma: Gamma parameter (kinds 2, 3).
        :param SPk_epsilon: Epsilon parameter (kind=3).
        :param SPk_m_pivot: Pivot mass in M_sun (kind=3).
        :param halofit_version: Base Halofit version for the wrapped non-linear model.
        :param HMCode_A_baryon: HMcode A_baryon parameter (mead2015/2016 only). Default 3.13.
        :param HMCode_eta_baryon: HMcode eta_baryon parameter (mead2015/2016 only). Default 0.603.
        :param HMCode_logT_AGN: HMcode log10(T_AGN/K) parameter (mead2020_feedback only). Default 7.8.
        :return: Self, for fluent configuration.
        :raises CAMBValueError: If parameters are invalid or incompatible with the base model.

        **Cobaya usage:**

        Cobaya passes keys from ``extra_args`` directly to ``set_params()``.
        Parameters under the theory ``params:`` block are also forwarded and can be sampled.

        - **extra_args** (fixed): ``non_linear_model``, ``halofit_version``, ``SPk_feedback``,
          ``SPk_SO``, ``SPk_relation_kind``, and pivot masses.
        - **params** (sampled): continuous relation parameters (e.g. ``SPk_fb_a``, ``SPk_fb_pow``).

        Example YAML (kind=3, double_power_law)::

            params:
              SPk_epsilon:
                prior: {min: 0.24, max: 0.35}
                ref: {dist: norm, loc: 0.30, scale: 0.02}
              SPk_alpha:
                prior: {min: -0.12, max: 0.34}
              SPk_beta:
                prior: {min: -0.74, max: 0.77}
              SPk_gamma:
                prior: {min: -0.5, max: 1.20}
              log10_SPk_m_pivot:
                prior: {min: 13, max: 14}
                drop: true
              SPk_m_pivot:
                value: "lambda log10_SPk_m_pivot: 10**log10_SPk_m_pivot"

            theory:
              camb:
                extra_args:
                  non_linear_model: SPkNonLinear
                  halofit_version: mead2020
                  SPk_feedback: true
                  SPk_SO: 200
                  SPk_relation_kind: 3

        **Boundary behavior:**

        - z outside [0, 3]: no suppression applied (identity).
        - k > 12 h/Mpc: suppression clamped to value at k = 12 (saturates).
        - f_b outside calibrated limits: ``nonlin_ratio`` set to NaN.

        NaN propagates through
        :meth:`~camb.results.CAMBdata.get_matter_power_interpolator`
        (Cobaya rejects the sample). Set priors to keep f_b in range.

        .. warning::
           :meth:`~camb.results.CAMBdata.get_matter_power_spectrum` does not
           preserve NaN — use
           :meth:`~camb.results.CAMBdata.get_linear_matter_power_spectrum`
           with ``nonlinear=True`` to inspect which (z, k) cells are invalid.

        Cannot be combined with ``halofit_version='mead2020_feedback'``.
        """
        if self.BaseModel is None:
            self.BaseModel = Halofit()
        if isinstance(self.BaseModel, Halofit):
            self.BaseModel.set_params(
                halofit_version=halofit_version,
                HMCode_A_baryon=HMCode_A_baryon,
                HMCode_eta_baryon=HMCode_eta_baryon,
                HMCode_logT_AGN=HMCode_logT_AGN,
            )

        self.SPk_feedback = SPk_feedback
        self.SPk_SO = SPk_SO
        self.SPk_relation_kind = SPk_relation_kind
        self.SPk_fb_a = SPk_fb_a
        self.SPk_fb_pow = SPk_fb_pow
        self.SPk_fb_pivot = SPk_fb_pivot
        self.SPk_alpha = SPk_alpha
        self.SPk_beta = SPk_beta
        self.SPk_gamma = SPk_gamma
        self.SPk_epsilon = SPk_epsilon
        self.SPk_m_pivot = SPk_m_pivot
        self._validate()
        return self


@fortran_class
class SecondOrderPK(NonLinearModel):
    """
    Third-order Newtonian perturbation theory results for the non-linear correction.
    Only intended for use at very high redshift (z>10) where corrections are perturbative, it will not give
    sensible results at low redshift.

    See Appendix F of `astro-ph/0702600 <https://arxiv.org/abs/astro-ph/0702600>`_ for equations and references.

    Not intended for production use, it's mainly to serve as an example alternative non-linear model implementation.
    """

    _fortran_class_module_ = "SecondOrderPK"
    _fortran_class_name_ = "TSecondOrderPK"

    def set_params(self):
        pass


@fortran_class
class ExternalNonLinearRatio(NonLinearModel):
    """
    Non-linear model that applies a user-supplied ratio ``sqrt(P_NL/P_L)``
    from an external source, for example CCL or axionHMcode.

    Use :meth:`set_ratio` to provide the ratio grid, then assign the instance
    to ``params.NonLinearModel`` before calling :func:`camb.get_results`.
    This can also be used after computing time transfers, so lensed ``C_l``
    values are generated with a consistent external non-linear prescription.
    Requested ``k_h`` or ``z`` values outside the supplied grid are clamped to
    the nearest grid boundary.
    """

    _fortran_class_module_ = "ExternalNonLinearRatio"
    _fortran_class_name_ = "TExternalNonLinearRatio"

    _methods_ = (
        (
            "SetRatio",
            [
                POINTER(c_int),
                POINTER(c_int),
                numpy_1d,
                numpy_1d,
                ndpointer(c_double, flags="F_CONTIGUOUS", ndim=2),
            ],
        ),
        ("ClearRatio", []),
    )

    def set_params(self):
        pass

    def set_ratio(self, k_h, z, ratio):
        """
        Set the non-linear ratio grid sqrt(P_NL/P_L).

        :param k_h: 1D array of k values in h/Mpc units (ascending)
        :param z: 1D array of redshift values (ascending)
        :param ratio: 2D array of sqrt(P_NL/P_L), shape (len(z), len(k_h)),
                      matching the convention of CAMB's get_matter_power_spectrum.
                      Values requested outside the supplied grid are clamped to
                      the nearest tabulated boundary.
        """
        k_h = np.ascontiguousarray(k_h, dtype=np.float64)
        z = np.ascontiguousarray(z, dtype=np.float64)
        if ratio.shape != (len(z), len(k_h)):
            raise ValueError(f"ratio shape {ratio.shape} must be (len(z), len(k_h)) = ({len(z)}, {len(k_h)})")
        # Fortran expects (nk, nz) column-major; C-order (nz, nk) has the same memory layout
        ratio_f = np.asfortranarray(ratio.T, dtype=np.float64)
        self.f_SetRatio(byref(c_int(len(k_h))), byref(c_int(len(z))), k_h, z, ratio_f)

    def clear_ratio(self):
        """
        Clear the stored ratio grid and release the interpolation data.
        """
        self.f_ClearRatio()
