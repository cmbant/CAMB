from .baseconfig import F2003Class, fortran_class
from ctypes import c_int, c_double

# ---Parameters in halofit_ppf.f90

halofit_original = 'original'
halofit_bird = 'bird'
halofit_peacock = 'peacock'
halofit_takahashi = 'takahashi'
halofit_mead = 'mead'
halofit_halomodel = 'halomodel'
halofit_casarini = 'casarini'
halofit_mead2015 = 'mead2015'

halofit_default = halofit_takahashi

halofit_version_names = [halofit_original, halofit_bird, halofit_peacock, halofit_takahashi, halofit_mead,
                         halofit_halomodel, halofit_casarini, halofit_mead2015]


class NonLinearModel(F2003Class):
    """
    Abstract base class for non-linear correction models
    """
    pass


@fortran_class
class Halofit(NonLinearModel):
    """
    Various specific approximate non-linear correction models based on HaloFit.
    """
    _fields_ = [
        ("Min_kh_nonlinear", c_double),
        ("halofit_version", c_int, {"names": halofit_version_names, "start": 1})
    ]

    _fortran_class_module_ = 'NonLinear'
    _fortran_class_name_ = 'THalofit'

    def get_halofit_version(self):
        return self.halofit_version

    def set_params(self, halofit_version=halofit_default):
        """
        Set the halofit model for non-linear corrections.

        :param version: One of

            - original: `astro-ph/0207664 <http://arxiv.org/abs/astro-ph/0207664>`_
            - bird: `arXiv:1109.4416 <http://arxiv.org/abs/1109.4416>`_
            - peacock: `Peacock fit <http://www.roe.ac.uk/~jap/haloes/>`_
            - takahashi: `arXiv:1208.2701 <http://arxiv.org/abs/1208.2701>`_
            - mead: HMCode `arXiv:1602.02154 <http://arxiv.org/abs/1602.02154>`_
            - halomodel: basic halomodel
            - casarini: PKequal `arXiv:0810.0190 <http://arxiv.org/abs/0810.0190>`_, `arXiv:1601.07230 <http://arxiv.org/abs/1601.07230>`_
            - mead2015: original 2015 version of HMCode `arXiv:1505.07833 <http://arxiv.org/abs/1505.07833>`_

        """
        self.halofit_version = halofit_version
