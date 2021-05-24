from .baseconfig import F2003Class, fortran_class
from ctypes import c_int, c_double


class NonLinearModel(F2003Class):
    """
    Abstract base class for non-linear correction models
    """
    _fields_ = [("Min_kh_nonlinear", c_double, "minimum k/h at which to apply non-linear corrections")]


halofit_original = 'original'
halofit_bird = 'bird'
halofit_peacock = 'peacock'
halofit_takahashi = 'takahashi'
halofit_mead = 'mead'
halofit_halomodel = 'halomodel'
halofit_casarini = 'casarini'
halofit_mead2015 = 'mead2015'
halofit_mead2016 = 'mead2016'
halofit_mead2020 = 'mead2020'
halofit_mead2020_feedback = 'mead2020_feedback'

halofit_default = halofit_mead2020

halofit_version_names = {halofit_original: 1,
                         halofit_bird: 2,
                         halofit_peacock: 3,
                         halofit_takahashi: 4,
                         halofit_mead: 5,
                         halofit_halomodel: 6,
                         halofit_casarini: 7,
                         halofit_mead2015: 8,
                         halofit_mead2016: 5,
                         halofit_mead2020: 9,
                         halofit_mead2020_feedback: 10}


@fortran_class
class Halofit(NonLinearModel):
    """
    Various specific approximate non-linear correction models based on HaloFit.
    """
    _fields_ = [
        ("halofit_version", c_int, {"names": halofit_version_names}),
        ("HMCode_A_baryon", c_double, "HMcode parameter A_baryon"),
        ("HMCode_eta_baryon", c_double, "HMcode parameter eta_baryon"),
        ("HMCode_logT_AGN", c_double, "HMcode parameter log10(T_AGN/K)")
    ]

    _fortran_class_module_ = 'NonLinear'
    _fortran_class_name_ = 'THalofit'

    def get_halofit_version(self):
        return self.halofit_version

    def set_params(self, halofit_version=halofit_default, HMCode_A_baryon=3.13, HMCode_eta_baryon=0.603,
                   HMCode_logT_AGN=7.8):
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
            - mead2020: 2020 version of HMcode `arXiv:2009.01858 <https://arxiv.org/abs/2009.01858>`_
            - mead2020_feedback: 2020 version of HMcode with baryonic feedback `arXiv:2009.01858 <https://arxiv.org/abs/2009.01858>`_
        :param HMCode_A_baryon: HMcode parameter A_baryon. Default 3.13.
        :param HMCode_eta_baryon: HMcode parameter eta_baryon. Default 0.603.
        :param HMCode_logT_AGN: HMcode parameter logT_AGN. Default 7.8.
        """
        self.halofit_version = halofit_version
        self.HMCode_A_baryon = HMCode_A_baryon
        self.HMCode_eta_baryon = HMCode_eta_baryon
        self.HMCode_logT_AGN = HMCode_logT_AGN


@fortran_class
class SecondOrderPK(NonLinearModel):
    """
    Third-order Newtonian perturbation theory results for the non-linear correction.
    Only intended for use at very high redshift (z>10) where corrections are perturbative, it will not give
    sensible results at low redshift.

    See Appendix F of `astro-ph/0702600 <https://arxiv.org/abs/astro-ph/0702600>`_ for equations and references.

    Not intended for production use, it's mainly to serve as an example alternative non-linear model implementation.
    """

    _fortran_class_module_ = 'SecondOrderPK'
    _fortran_class_name_ = 'TSecondOrderPK'

    def set_params(self):
        pass
