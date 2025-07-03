# Simple test to check if eigenvector data can be read by fortran.

import numpy as np
import os
import sys
# Get the absolute path to the CAMB root directory
CAMB_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if CAMB_dir not in sys.path:
    sys.path.insert(0, CAMB_dir)

from camb.active_eigenvectors import activesources

# --- Eigenvector data ---
NPZ_FILENAME = "correlator_table.npz"

# ------------------------------------------------------------------------------
# Load Correlator Data and Initialize Custom Class
# ------------------------------------------------------------------------------
print(f"Loading eigenvector data from: {NPZ_FILENAME}...")
if not os.path.exists(NPZ_FILENAME):
    print(f"ERROR: Eigenvector data file '{NPZ_FILENAME}' not found. Please generate it first.")
    exit()

correlator_data = np.load(NPZ_FILENAME)
k_grid = correlator_data['k_grid']
tau_grid = correlator_data['ktau_grid']
all_eigenfunctions = correlator_data['eigenfunctions']
all_eigenfunctions_d_dlogkt = correlator_data['eigenfunctions_d_dlogkt']
all_eigenvalues_S = correlator_data['eigenvalues_S']
all_eigenvalues_00 = correlator_data['eigenvalues_00']
all_eigenvalues_V = correlator_data['eigenvalues_V']
all_eigenvalues_T = correlator_data['eigenvalues_T']
string_p_mu = correlator_data['string_params_mu'].item()
nmodes_from_file = correlator_data['nmodes'].item()
weighting_from_file = correlator_data['weighting_gamma'].item()
print("Correlator data loaded.")

print("Initializing custom Fortran object...")
my_custom_obj = activesources()
my_custom_obj.set_correlator_table(
    k_grid=k_grid,
    tau_grid=tau_grid,
    eigenfunctions=all_eigenfunctions,
    eigenfunctions_d_dlogkt=all_eigenfunctions_d_dlogkt,
    eigenvalues_S=all_eigenvalues_S,
    eigenvalues_00=all_eigenvalues_00,
    eigenvalues_V=all_eigenvalues_V,
    eigenvalues_T=all_eigenvalues_T,
    string_params_mu=string_p_mu,
    nmodes_param=nmodes_from_file,
    weighting_param=weighting_from_file
)
print("Fortran data transfer successful")