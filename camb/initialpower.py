# Initial power spectrum parameters

from .baseconfig import F2003Class, CAMBError, fortran_class, \
    c_int, c_double, POINTER, byref, numpy_1d, np

tensor_parameterization_names = ["tensor_param_indeptilt", "tensor_param_rpivot", "tensor_param_AT"]
tensor_param_indeptilt = 1
tensor_param_rpivot = 2
tensor_param_AT = 3


class InitialPower(F2003Class):
    """
    Abstract base class for initial power spectrum classes
    """
    _fortran_class_module_ = 'InitialPower'

    def set_params(self):
        pass


@fortran_class
class SplinedInitialPower(InitialPower):
    """
    Object to store a generic primordial spectrum set from a set of sampled k_i, P(k_i) values
    """
    _fortran_class_name_ = 'TSplinedInitialPower'

    _fields_ = [
        ('effective_ns_for_nonlinear', c_double, "Effective n_s to use for approximate non-linear correction models")]

    _methods_ = [('HasTensors', [], c_int),
                 ('SetScalarTable', [POINTER(c_int), numpy_1d, numpy_1d]),
                 ('SetTensorTable', [POINTER(c_int), numpy_1d, numpy_1d]),
                 ('SetScalarLogRegular', [POINTER(c_double), POINTER(c_double), POINTER(c_int), numpy_1d]),
                 ('SetTensorLogRegular', [POINTER(c_double), POINTER(c_double), POINTER(c_int), numpy_1d])]

    def __init__(self, **kwargs):
        if kwargs.get('PK', None) is not None:
            self.set_scalar_table(kwargs['ks'], kwargs['PK'])

    def has_tensors(self):
        """
        Is the tensor spectrum set?

        :return: True if tensors
        """
        return self.f_HasTensors() != 0

    def set_scalar_table(self, k, PK):
        """
        Set arrays of k and P(k) values for cublic spline interpolation.
        Note that using :meth:`set_scalar_log_regular` may be better
        (faster, and easier to get fine enough spacing a low k)

        :param k: array of k values (Mpc^{-1})
        :param PK: array of scalar power spectrum values
        """
        self.f_SetScalarTable(byref(c_int(len(k))), np.asarray(k), np.asarray(PK))

    def set_tensor_table(self, k, PK):
        """
        Set arrays of k and P_t(k) values for cublic spline interpolation

        :param k: array of k values (Mpc^{-1})
        :param PK: array of tensor power spectrum values
        """
        self.f_SetTensorTable(byref(c_int(len(k))), np.asarray(k), np.asarray(PK))

    def set_scalar_log_regular(self, kmin, kmax, PK):
        """
        Set log-regular cublic spline interpolation for P(k)

        :param kmin: minimum k value (not minimum log(k))
        :param kmax: maximum k value (inclusive)
        :param PK: array of scalar power spectrum values, with PK[0]=P(kmin) and PK[-1]=P(kmax)
        """
        self.f_SetScalarLogRegular(byref(c_double(kmin)), byref(c_double(kmax)), byref(c_int(len(PK))),
                                   np.asarray(PK))

    def set_tensor_log_regular(self, kmin, kmax, PK):
        """
        Set log-regular cublic spline interpolation for tensor spectrum P_t(k)

        :param kmin: minimum k value (not minimum log(k))
        :param kmax: maximum k value (inclusive)
        :param PK: array of scalar power spectrum values, with PK[0]=P_t(kmin) and PK[-1]=P_t(kmax)
        """

        self.f_SetTensorLogRegular(byref(c_double(kmin)), byref(c_double(kmax)), byref(c_int(len(PK))),
                                   np.asarray(PK))


@fortran_class
class InitialPowerLaw(InitialPower):
    """
    Object to store parameters for the primordial power spectrum in the standard power law expansion.

    """
    _fields_ = [
        ("tensor_parameterization", c_int, {"names": tensor_parameterization_names, "start": 1}),
        ("ns", c_double),
        ("nrun", c_double),
        ("nrunrun", c_double),
        ("nt", c_double),
        ("ntrun", c_double),
        ("r", c_double),
        ("pivot_scalar", c_double),
        ("pivot_tensor", c_double),
        ("As", c_double),
        ("At", c_double)
    ]

    _fortran_class_name_ = 'TInitialPowerLaw'

    def __init__(self, **kwargs):
        self.set_params(**kwargs)

    def set_params(self, As=2e-9, ns=0.96, nrun=0, nrunrun=0.0, r=0.0, nt=None, ntrun=0.0,
                   pivot_scalar=0.05, pivot_tensor=0.05, parameterization="tensor_param_rpivot"):
        r"""
        Set parameters using standard power law parameterization. If nt=None, uses inflation consistency relation.

        :param As: comoving curvature power at k=pivot_scalar (:math:`A_s`)
        :param ns: scalar spectral index :math:`n_s`
        :param nrun: running of scalar spectral index :math:`d n_s/d \log k`
        :param nrunrun: running of running of spectral index, :math:`d^2 n_s/d (\log k)^2`
        :param r: tensor to scalar ratio at pivot
        :param nt: tensor spectral index :math:`n_t`. If None, set using inflation consistency
        :param ntrun: running of tensor spectral index
        :param pivot_scalar: pivot scale for scalar spectrum
        :param pivot_tensor:  pivot scale for tensor spectrum
        :param parameterization: See CAMB notes. One of
            - tensor_param_indeptilt = 1
            - tensor_param_rpivot = 2
            - tensor_param_AT = 3
        :return: self
        """

        if parameterization not in [tensor_param_rpivot, tensor_param_indeptilt, "tensor_param_rpivot",
                                    "tensor_param_indeptilt"]:
            raise CAMBError('Initial power parameterization not supported here')
        self.tensor_parameterization = parameterization
        self.As = As
        self.ns = ns
        self.nrun = nrun
        self.nrunrun = nrunrun
        if nt is None:
            # set from inflationary consistency
            if ntrun:
                raise CAMBError('ntrun set but using inflation consistency (nt=None)')
            if tensor_param_rpivot != tensor_param_rpivot:
                raise CAMBError('tensor parameterization not tensor_param_rpivot with inflation consistency')
            self.nt = - r / 8.0 * (2.0 - ns - r / 8.0)
            self.ntrun = r / 8.0 * (r / 8.0 + ns - 1)
        else:
            self.nt = nt
            self.ntrun = ntrun
        self.r = r
        self.pivot_scalar = pivot_scalar
        self.pivot_tensor = pivot_tensor
        return self

    def has_tensors(self):
        """
        Do these settings have non-zero tensors?

        :return: True if non-zero tensor amplitude
        """
        return self.r > 0
