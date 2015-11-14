# Initial power spectrum parameters

from ctypes import *
from baseconfig import CAMB_Structure, CAMBError

nnmax = 5

tensor_param_indeptilt = 1
tensor_param_rpivot = 2
tensor_param_AT = 3


class InitialPowerParams(CAMB_Structure):
    _fields_ = [
        ("tensor_parameterization", c_int),
        ("nn", c_int),
        ("an", c_double * nnmax),
        ("n_run", c_double * nnmax),
        ("n_runrun", c_double * nnmax),
        ("ant", c_double * nnmax),
        ("nt_run", c_double * nnmax),
        ("rat", c_double * nnmax),
        ("k_0_scalar", c_double),
        ("k_0_tensor", c_double),
        ("ScalarPowerAmp", c_double * nnmax),
        ("TensorPowerAmp", c_double * nnmax)
    ]

    def set_params(self, As=2e-9, ns=0.96, nrun=0, nrunrun=0, r=0, nt=None, ntrun=0,
                   pivot_scalar=0.05, pivot_tensor=0.05, parameterization=tensor_param_rpivot):

        if not parameterization in [tensor_param_rpivot, tensor_param_indeptilt]:
            raise CAMBError('Initial power paramterization not supported here')
        self.tensor_parameterization = parameterization
        self.ScalarPowerAmp[0] = As
        self.nn = 1  # number of different power spectra
        self.an[0] = ns
        self.n_run[0] = nrun
        self.n_runrun[0] = nrunrun
        if nt is None:
            # set from inflationary consistency
            if ntrun: raise CAMBError('ntrun set but using inflation consistency (nt=None)')
            if tensor_param_rpivot != tensor_param_rpivot:
                raise CAMBError('tensor parameterization not tensor_param_rpivot with inflation consistency')
            self.ant[0] = - r / 8.0 * (2.0 - ns - r / 8.0)
            self.nt_run[0] = r / 8.0 * (r / 8.0 + ns - 1)
        else:
            self.ant[0] = nt
            self.nt_run[0] = ntrun
        self.rat[0] = r
        self.k_0_scalar = pivot_scalar
        self.k_0_tensor = pivot_tensor

    def has_tensors(self):
        return self.rat[0]
