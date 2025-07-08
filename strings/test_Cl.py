#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import sys
import os
# Get the absolute path to the CAMB root directory
CAMB_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if CAMB_dir not in sys.path:
    sys.path.insert(0, CAMB_dir)

from camb.active_sources import activesources
import camb
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
NPZ_FILENAME = os.path.join(script_dir, "correlator_table.npz")

N_MODES_TO_SUM = 32  # Number of UETC eigenmodes to sum for the string signal
LMAX_PLOT = 4000    # Max multipole for plotting
CMB_UNIT_OUTPUT = 'muK' # 'muK' for muK^2 units, 'K' for K^2 units
pol_mode_idx=0;# !indices: TT, EE, BB, TE

# --- Plotting Style (Optional) ---
plt.style.use('seaborn-v0_8-colorblind') # Or any other style you prefer

# ------------------------------------------------------------------------------
# 1. CAMB Parameter Setup
# ------------------------------------------------------------------------------
print("Setting up CAMB parameters...")
pars = camb.CAMBparams()
# Have turned off massive neutrinos to avoid issues with the vector modes
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.0, omk=0, tau=0.06)
pars.max_l_tensor = 1500

pars.WantScalars = False
pars.WantVectors = False
pars.WantTensors = True
pars.DoLensing = False

# if pars.WantTensors:
#     pars.InitPower.set_params(As=2e-9, ns=0.965, r=0.1) 
# else:
#     pars.InitPower.set_params(As=2e-9, ns=0.965, r=0.0) 

# ------------------------------------------------------------------------------
# 2. Load UETC Data and Initialize Custom Object
# ------------------------------------------------------------------------------
print(f"Loading UETC data from: {NPZ_FILENAME}...")
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
pars.ActiveSources = my_custom_obj 

# ------------------------------------------------------------------------------
# 3. Calculate Baseline C_l^{EE} (No UETC Sources)
# ------------------------------------------------------------------------------
print("Calculating baseline C_l^{EE} (UETC sources OFF)...")
my_custom_obj.set_active_eigenmode(0) 
results_baseline = camb.get_results(pars)

if pars.WantVectors:
    power_spectra_baseline = results_baseline.get_vector_cls(CMB_unit=CMB_UNIT_OUTPUT, raw_cl=False)
    cl_baseline_dl = power_spectra_baseline[:,pol_mode_idx]
else:
    power_spectra_baseline = results_baseline.get_cmb_power_spectra(pars, CMB_unit=CMB_UNIT_OUTPUT, raw_cl=False)
    cl_baseline_dl = power_spectra_baseline['total'][:,pol_mode_idx]

lmax_calc = cl_baseline_dl.shape[0] - 1
ls_calc = np.arange(lmax_calc + 1)

print(f"Baseline C_l^{{EE}} calculated up to LMAX={lmax_calc}.")

# ------------------------------------------------------------------------------
# 4. Calculate UETC C_l^{EE} by Summing Modes
# ------------------------------------------------------------------------------

if pars.WantTensors and not pars.WantVectors and not pars.WantScalars:
    # pars.InitPower.set_params(As=1, ns=4, r=1, nt=3, pivot_scalar=1.0, pivot_tensor=1.0) 
    # scale_factor = 16/(2*np.pi**2)
    scale_factor=1
elif pars.WantVectors and not pars.WantTensors and not pars.WantScalars:
    pars.InitPower.set_params(As=1, ns=4, r=0, nt=3, pivot_scalar=1.0, pivot_tensor=1.0) 
    # Fudge factor of 2
    # scale_factor = 2 * 8/(2*np.pi**2)
    scale_factor=1
elif pars.WantScalars and not pars.WantTensors and not pars.WantVectors:
    # pars.InitPower.set_params(As=1, ns=4, r=0, nt=3, pivot_scalar=1.0, pivot_tensor=1.0) 
    pars.scalar_initial_condition = 0
    scale_factor=1
    # scale_factor = 1/(2*np.pi**2)
else:
    print("Error: Invalid combination of modes.")
    exit()

actual_n_modes_to_sum = min(N_MODES_TO_SUM, nmodes_from_file)
print(f"Calculating UETC C_l^{{EE}} by summing {actual_n_modes_to_sum} eigenmodes...")

cl_strings_sum_dl = np.zeros_like(cl_baseline_dl)

# pars.InitPower.set_params(As=1, ns=1, r=1) 
for i_mode in range(1, actual_n_modes_to_sum + 1):
    print(f"  Processing eigenmode {i_mode}/{actual_n_modes_to_sum}...")
    my_custom_obj.set_active_eigenmode(i_mode)
    
    pars.ActiveSources = my_custom_obj # Re-assign to ensure CAMB sees the updated active_mode

    results_mode_i = camb.get_results(pars)

    if pars.WantVectors:
        power_spectra_mode_i = results_mode_i.get_vector_cls(CMB_unit=CMB_UNIT_OUTPUT, raw_cl=False)
        cl_mode_i_dl = power_spectra_mode_i[:,pol_mode_idx]
    else:
        power_spectra_mode_i = results_mode_i.get_cmb_power_spectra(pars, CMB_unit=CMB_UNIT_OUTPUT, raw_cl=False)
        cl_mode_i_dl = power_spectra_mode_i['total'][:,pol_mode_idx] # !indices: TT, EE, BB, TE
        
    min_len = min(len(cl_strings_sum_dl), len(cl_mode_i_dl))
    
    cl_strings_sum_dl[:min_len] += cl_mode_i_dl[:min_len]

my_custom_obj.set_active_eigenmode(0) 
print("UETC C_l^{EE} calculation finished.")

# Assumes linear addition of power spectra.
Cl_strings = scale_factor * cl_strings_sum_dl

# ------------------------------------------------------------------------------
# 5. Plotting
# ------------------------------------------------------------------------------
print("Plotting results...")
fig, axs = plt.subplots(1, 2, figsize=(16, 6))  # 1 row, 2 columns

plot_mask = (ls_calc >= 1) & (ls_calc <= LMAX_PLOT)
ls_plot = ls_calc[plot_mask]

# First plot: Strings EE
axs[0].plot(ls_plot, Cl_strings[plot_mask]*(2e-7)**2, label='Strings', color='C1', linestyle='-')
axs[0].set_xlabel(r'$\ell$')
axs[0].set_ylabel(r'$\ell(\ell+1)C_\ell/2\pi \, [\mu K^2]$')
axs[0].set_title('Strings Power Spectrum')
axs[0].set_xlim([2, LMAX_PLOT])
axs[0].set_xscale('log')
axs[0].legend()
axs[0].grid(True, which="both", ls="-", alpha=0.5)


# Second plot: Baseline EE
axs[1].plot(ls_plot, cl_baseline_dl[plot_mask], label='Baseline', color='C0', linestyle='-')
axs[1].set_xlabel(r'$\ell$')
axs[1].set_ylabel(r'$\ell(\ell+1)C_\ell/2\pi \, [\mu K^2]$')
axs[1].set_title('Baseline Power Spectrum')
axs[1].set_xlim([2, LMAX_PLOT])
axs[1].set_ylim(ymin=0)
axs[1].set_xscale('log')
axs[1].legend()
axs[1].grid(True, which="both", ls="-", alpha=0.5)

plt.tight_layout()
# plt.savefig("cl_comparison_side_by_side.png")
# print("Plot saved to cl_comparison_side_by_side.png")

plt.show()

print("--- Script Finished ---")