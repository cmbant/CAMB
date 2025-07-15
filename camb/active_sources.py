# Defines a custom python class and uses the pre-existing CAMB architecture to relate it to a fortran class.
# This class is used to make the active source correlator eigenvector table calculated in python available for use in fortran.

from .baseconfig import F2003Class, fortran_class, numpy_1d, np
from ctypes import c_int, c_double, byref, POINTER

class BaseClass(F2003Class):
    pass

@fortran_class
class ActiveSources(BaseClass):
    _fortran_class_module_ = 'ActiveSources'
    _fortran_class_name_ = 'TActiveSources'
    _fields_ = []

    _methods_ = [

        ('SetCorrelatorTable', [
            numpy_1d, POINTER(c_int),          # k_grid_in, nk_in
            numpy_1d, POINTER(c_int),          # tau_grid_in, ntau_in
            numpy_1d, POINTER(c_int),          # eigenfuncs_flat_in, num_eigen_types_in
            POINTER(c_int),                    # nmodes_in
            numpy_1d,                          # evals_S_flat_in
            numpy_1d,                          # evals_00_flat_in

            numpy_1d,                          # evals_V_flat_in
            numpy_1d,                          # evals_T_flat_in
            numpy_1d,                          # eigenfunc_derivs_logkt_flat_in
            POINTER(c_double),                 # mu_in
            POINTER(c_double)                  # weighting_param_in
        ]),
        ('SetActiveEigenmode', [POINTER(c_int)]),
    ]

    def set_correlator_table(self, k_grid, tau_grid, eigenfunctions,
                              eigenfunctions_d_dlogkt,
                              eigenvalues_S,eigenvalues_00, eigenvalues_V, eigenvalues_T,
                              string_params_mu, nmodes_param, weighting_param):
        print("Python: Preparing eigenvector table for Fortran...")
        # NOTE! For some reason only passing flat contiguous arrays seems to be allowed (at least only these are allowed easily)
        # We pass arrays as flat and restructure them in fortran
        k_grid_c = np.ascontiguousarray(k_grid, dtype=np.float64)
        tau_grid_c = np.ascontiguousarray(tau_grid, dtype=np.float64)
        eigenfunctions_flat_c = np.ascontiguousarray(eigenfunctions, dtype=np.float64).ravel(order='C')

        eigenfunctions_derivs_flat_c = np.ascontiguousarray(eigenfunctions_d_dlogkt, dtype=np.float64).ravel(order='C')

        eigenvalues_S_flat_c = np.ascontiguousarray(eigenvalues_S, dtype=np.float64).ravel(order='C')
        eigenvalues_00_flat_c = np.ascontiguousarray(eigenvalues_00, dtype=np.float64).ravel(order='C')

        eigenvalues_V_flat_c = np.ascontiguousarray(eigenvalues_V, dtype=np.float64).ravel(order='C')
        eigenvalues_T_flat_c = np.ascontiguousarray(eigenvalues_T, dtype=np.float64).ravel(order='C')

        nk = k_grid_c.shape[0]
        ntau = tau_grid_c.shape[0]
        num_eigen_types = eigenfunctions.shape[1]

        print(f"  Python Sending: nk={nk}, ntau={ntau}, Ntypes={num_eigen_types}, Nmodes={nmodes_param}")

        self.f_SetCorrelatorTable(
            k_grid_c, byref(c_int(nk)),
            tau_grid_c, byref(c_int(ntau)),
            eigenfunctions_flat_c, byref(c_int(num_eigen_types)),
            byref(c_int(nmodes_param)),
            eigenvalues_S_flat_c,
            eigenvalues_00_flat_c,

            eigenvalues_V_flat_c,
            eigenvalues_T_flat_c,
            eigenfunctions_derivs_flat_c,
            byref(c_double(string_params_mu)),
            byref(c_double(weighting_param))
        )
        print("Python: Fortran call to SetcorrelatorTable completed.")
        return self

    def set_active_eigenmode(self, mode_idx):
            if not isinstance(mode_idx, int):
                raise TypeError("mode_idx must be an integer.")

            print(f"Python: Setting active correlator eigenmode to: {mode_idx}")
            self.f_SetActiveEigenmode(byref(c_int(mode_idx)))
            return self
