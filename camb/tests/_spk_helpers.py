"""Shared utilities for SP(k) tests and validation scripts.

Not collected by pytest (underscore prefix).
"""

import numpy as np

import camb
from camb.nonlinear import Halofit, SPkNonLinear  # type: ignore[attr-defined]


class CambEFunc:
    """Wrapper providing E(z) = H(z)/H0 from CAMB results, for pySPK."""

    def __init__(self, results):
        self._results = results
        self._h0 = results.Params.H0

    def efunc(self, z):
        z_arr = np.asarray(z, dtype=float)
        flat = z_arr.reshape(-1)
        vals = np.array([self._results.hubble_parameter(float(zz)) / self._h0 for zz in flat], dtype=float)
        vals = vals.reshape(z_arr.shape)
        if np.isscalar(z):
            return float(vals)
        return vals


def get_pk(model_obj, z=0.5, kmax=12.0, k_per_logint=100):
    """Run CAMB with the given NonLinearModel and return (k, pk, results)."""
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.5, ombh2=0.02237, omch2=0.12, mnu=0.06)
    pars.InitPower.set_params(As=2.1e-9, ns=0.965)
    pars.set_matter_power(redshifts=[z], kmax=kmax, k_per_logint=k_per_logint)
    pars.NonLinear = camb.model.NonLinear_both
    pars.NonLinearModel = model_obj
    results = camb.get_results(pars)
    kh, _z, pk = results.get_nonlinear_matter_power_spectrum()
    return kh, pk[0], results


def pyspk_kwargs(relation_kind, so, z, k, camb_results, params):
    """Build keyword arguments for pyspk.sup_model()."""
    base = {"SO": so, "z": z, "k_array": k, "errors": False}
    if relation_kind == 1:
        base.update({"fb_a": params["SPk_fb_a"], "fb_pow": params["SPk_fb_pow"], "fb_pivot": params["SPk_fb_pivot"]})
    elif relation_kind == 2:
        base.update(
            {
                "alpha": params["SPk_alpha"],
                "beta": params["SPk_beta"],
                "gamma": params["SPk_gamma"],
                "cosmo": CambEFunc(camb_results),
            }
        )
    elif relation_kind == 3:
        base.update(
            {
                "epsilon": params["SPk_epsilon"],
                "alpha": params["SPk_alpha"],
                "beta": params["SPk_beta"],
                "gamma": params["SPk_gamma"],
                "m_pivot": params["SPk_m_pivot"],
                "cosmo": CambEFunc(camb_results),
            }
        )
    else:
        raise ValueError(f"Unknown relation kind {relation_kind}")
    return base


def make_spk_model(relation_kind, so, params, halofit_version="mead2020"):
    """Create a configured SPkNonLinear instance."""
    spk = SPkNonLinear()
    spk.set_params(
        halofit_version=halofit_version,
        SPk_feedback=True,
        SPk_SO=so,
        SPk_relation_kind=relation_kind,
        SPk_fb_a=params.get("SPk_fb_a", 1.0),
        SPk_fb_pow=params.get("SPk_fb_pow", 0.0),
        SPk_fb_pivot=params.get("SPk_fb_pivot", 1.0),
        SPk_alpha=params.get("SPk_alpha", 0.0),
        SPk_beta=params.get("SPk_beta", 0.0),
        SPk_gamma=params.get("SPk_gamma", 0.0),
        SPk_epsilon=params.get("SPk_epsilon", 0.0),
        SPk_m_pivot=params.get("SPk_m_pivot", 1.0),
    )
    return spk


def get_param_sets(relation_kind):
    """Return a list of parameter dicts for the given relation kind."""
    if relation_kind == 1:
        return [
            {"SPk_fb_a": 0.25, "SPk_fb_pow": 0.10, "SPk_fb_pivot": 1e14},
            {"SPk_fb_a": 0.40, "SPk_fb_pow": 0.30, "SPk_fb_pivot": 1e14},
            {"SPk_fb_a": 0.70, "SPk_fb_pow": 0.05, "SPk_fb_pivot": 5e13},
        ]
    if relation_kind == 2:
        return [
            {"SPk_alpha": 3.2, "SPk_beta": 1.0, "SPk_gamma": 0.10},
            {"SPk_alpha": 3.6, "SPk_beta": 1.2, "SPk_gamma": 0.20},
            {"SPk_alpha": 4.0, "SPk_beta": 0.9, "SPk_gamma": 0.05},
        ]
    if relation_kind == 3:
        return [
            {"SPk_epsilon": 0.45, "SPk_alpha": 0.2, "SPk_beta": 0.9, "SPk_gamma": 0.1, "SPk_m_pivot": 1e14},
            {"SPk_epsilon": 0.60, "SPk_alpha": 0.1, "SPk_beta": 1.1, "SPk_gamma": 0.2, "SPk_m_pivot": 8e13},
            {"SPk_epsilon": 0.80, "SPk_alpha": 0.3, "SPk_beta": 0.8, "SPk_gamma": 0.0, "SPk_m_pivot": 2e14},
        ]
    raise ValueError(f"Unknown relation kind {relation_kind}")


def compute_suppression_data(relation_kinds, sos, redshifts, kmax=12.0):
    """Compute SPk suppression curves for all (relation_kind, so, param_index, z) combos.

    Caches the base Halofit P(k) per redshift (only 1 CAMB run per z).

    Returns
    -------
    dict keyed by (relation_kind, so, param_index, z) with values
        (k, sup_camb, sup_ref, rel)
    where rel = sup_camb / sup_ref - 1.
    """
    import pyspk as _pyspk

    # One base CAMB run per redshift
    base_cache = {}
    for z in redshifts:
        base = Halofit()
        base.set_params(halofit_version="mead2020")
        k_full, pk_full, results = get_pk(base, z=z, kmax=kmax)
        mask = k_full < 12.0
        base_cache[z] = (k_full[mask], pk_full[mask], results)

    data = {}
    for relation_kind in relation_kinds:
        for ip, params in enumerate(get_param_sets(relation_kind)):
            for so in sos:
                for z in redshifts:
                    k, pk_base, results = base_cache[z]

                    spk_model = make_spk_model(relation_kind, so, params)
                    k_spk, pk_spk_full, _ = get_pk(spk_model, z=z, kmax=kmax)
                    pk_spk = pk_spk_full[k_spk < 12.0]

                    py_kw = pyspk_kwargs(relation_kind, so, z, k, results, params)
                    _, sup_ref = _pyspk.sup_model(**py_kw)
                    sup_camb = pk_spk / pk_base

                    rel = np.full_like(sup_camb, np.nan)
                    valid = np.isfinite(sup_camb) & np.isfinite(sup_ref) & np.not_equal(sup_ref, 0.0)
                    rel[valid] = sup_camb[valid] / sup_ref[valid] - 1.0

                    data[(relation_kind, so, ip, z)] = (k, sup_camb, sup_ref, rel)
    return data
