from .baseconfig import CAMB_Structure
from ctypes import c_int

# ---Parameters in halofit_ppf.f90

halofit_original = 1
halofit_bird = 2
halofit_peacock = 3
halofit_takahashi = 4
halofit_mead = 5
halofit_halomodel = 6
halofit_casarini = 7
halofit_mead2015 = 8

halofit_default = halofit_takahashi

halofit_version_names = ['original', 'bird', 'peacock', 'takahashi', 'mead', 'halomodel', 'casarini', 'mead2015']


class NonLinearModel(CAMB_Structure):
    _fields_ = [
        ("Min_kh_nonlinear", c_int),
        ("halofit_version", c_int)
    ]

    def get_halofit_version(self):
        return halofit_version_names[self.halofit_version - 1]

    def set_halofit_version(self, version='takahashi'):
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
        self.halofit_version = halofit_version_names.index(version) + 1
