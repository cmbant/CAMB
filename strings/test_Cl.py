#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import sys
import os
import argparse
# Get the absolute path to the CAMB root directory
CAMB_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if CAMB_dir not in sys.path:
    sys.path.insert(0, CAMB_dir)

from camb.active_sources import ActiveSources
import camb
import matplotlib.pyplot as plt

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Calculate CMB power spectra with UETC string sources',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Combined mode selection for both UETC source type and CMB calculation type
    parser.add_argument(
        '--mode', '-m',
        type=str,
        choices=['scalar', 'vector', 'tensor'],
        default='scalar',
        help='Type of perturbation mode (scalar, vector, or tensor) - determines both UETC source type and CMB calculation type'
    )

    # Number of eigenmodes
    parser.add_argument(
        '--nmodes', '-n',
        type=int,
        default=32,
        help='Number of UETC eigenmodes to sum'
    )

    # Polarization component
    parser.add_argument(
        '--pol', '-p',
        type=str,
        choices=['TT', 'EE', 'BB', 'TE'],
        default='TT',
        help='Polarization component to plot'
    )

    # Maximum multipole for plotting
    parser.add_argument(
        '--lmax', '-l',
        type=int,
        default=4000,
        help='Maximum multipole for plotting'
    )

    # CMB units
    parser.add_argument(
        '--units', '-u',
        type=str,
        choices=['muK', 'K'],
        default='muK',
        help='CMB units for output'
    )

    # Data file path
    parser.add_argument(
        '--datafile', '-d',
        type=str,
        default=None,
        help='Path to correlator data file (default: correlator_table.npz in script directory)'
    )

    # Output file
    parser.add_argument(
        '--output', '-o',
        type=str,
        default=None,
        help='Output filename for plot (if not specified, plot is shown)'
    )

    # Disable plotting
    parser.add_argument(
        '--no-plot',
        action='store_true',
        help='Disable plotting'
    )

    return parser.parse_args()

def setup_camb_params(args):
    """Setup CAMB parameters based on arguments"""
    print("Setting up CAMB parameters...")
    pars = camb.CAMBparams()

    # Have turned off massive neutrinos to avoid issues with the vector modes
    pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.0, omk=0, tau=0.06)
    pars.max_l_tensor = 1500

    # Set mode types based on argument
    if args.mode == 'scalar':
        pars.WantScalars = True
        pars.WantVectors = False
        pars.WantTensors = False
    elif args.mode == 'vector':
        pars.WantScalars = False
        pars.WantVectors = True
        pars.WantTensors = False
    elif args.mode == 'tensor':
        pars.WantScalars = False
        pars.WantVectors = False
        pars.WantTensors = True

    pars.DoLensing = False

    return pars

def load_correlator_data(args):
    """Load UETC correlator data"""
    if args.datafile is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        npz_filename = os.path.join(script_dir, "correlator_table.npz")
    else:
        npz_filename = args.datafile

    print(f"Loading UETC data from: {npz_filename}...")
    correlator_data = np.load(npz_filename)

    return {
        'k_grid': correlator_data['k_grid'],
        'tau_grid': correlator_data['ktau_grid'],
        'eigenfunctions': correlator_data['eigenfunctions'],
        'eigenfunctions_d_dlogkt': correlator_data['eigenfunctions_d_dlogkt'],
        'eigenvalues_S': correlator_data['eigenvalues_S'],
        'eigenvalues_00': correlator_data['eigenvalues_00'],
        'eigenvalues_V': correlator_data['eigenvalues_V'],
        'eigenvalues_T': correlator_data['eigenvalues_T'],
        'string_p_mu': correlator_data['string_params_mu'].item(),
        'nmodes_from_file': correlator_data['nmodes'].item(),
        'weighting_from_file': correlator_data['weighting_gamma'].item()
    }

def setup_active_sources(correlator_data, args):
    """Setup ActiveSources object with correlator data"""
    print("Initializing custom Fortran object...")
    active_sources = ActiveSources()
    active_sources.set_correlator_table(
        k_grid=correlator_data['k_grid'],
        tau_grid=correlator_data['tau_grid'],
        eigenfunctions=correlator_data['eigenfunctions'],
        eigenfunctions_d_dlogkt=correlator_data['eigenfunctions_d_dlogkt'],
        eigenvalues_S=correlator_data['eigenvalues_S'],
        eigenvalues_00=correlator_data['eigenvalues_00'],
        eigenvalues_V=correlator_data['eigenvalues_V'],
        eigenvalues_T=correlator_data['eigenvalues_T'],
        string_params_mu=correlator_data['string_p_mu'],
        nmodes_param=correlator_data['nmodes_from_file'],
        weighting_param=correlator_data['weighting_from_file']
    )
    print("Fortran data transfer successful")

    return active_sources

def get_polarization_index(pol_component):
    """Get index for polarization component"""
    pol_map = {'TT': 0, 'EE': 1, 'BB': 2, 'TE': 3}
    return pol_map[pol_component]

def get_eigenvalues_for_mode(correlator_data, mode):
    """Get the appropriate eigenvalues for the specified mode"""
    if mode == 'scalar':
        return correlator_data['eigenvalues_S']
    elif mode == 'vector':
        return correlator_data['eigenvalues_V']
    elif mode == 'tensor':
        return correlator_data['eigenvalues_T']
    else:
        raise ValueError(f"Unknown mode: {mode}")

def calculate_power_spectra(pars, args):
    """Calculate power spectra for baseline and UETC sources"""
    pol_mode_idx = get_polarization_index(args.pol)

    # Calculate baseline (no UETC sources)
    print("Calculating baseline C_l (UETC sources OFF)...")
    pars.ActiveSources.set_active_eigenmode(0)
    results_baseline = camb.get_results(pars)

    if pars.WantVectors:
        power_spectra_baseline = results_baseline.get_cmb_power_spectra(pars, CMB_unit=args.units, raw_cl=False, spectra=("vector",))
        cl_baseline_dl = power_spectra_baseline['vector'][:,pol_mode_idx]
    else:
        power_spectra_baseline = results_baseline.get_cmb_power_spectra(pars, CMB_unit=args.units, raw_cl=False)
        cl_baseline_dl = power_spectra_baseline['total'][:,pol_mode_idx]

    lmax_calc = cl_baseline_dl.shape[0] - 1
    ls_calc = np.arange(lmax_calc + 1)

    print(f"Baseline C_l^{{{args.pol}}} calculated up to LMAX={lmax_calc}.")

    # Calculate UETC sources by summing modes
    correlator_data = load_correlator_data(args)
    actual_n_modes_to_sum = min(args.nmodes, correlator_data['nmodes_from_file'])
    print(f"Calculating UETC C_l^{{{args.pol}}} by summing {actual_n_modes_to_sum} eigenmodes...")
    print(f"Using {args.mode} mode")

    cl_strings_sum_dl = np.zeros_like(cl_baseline_dl)

    for i_mode in range(1, actual_n_modes_to_sum + 1):
        print(f"  Processing eigenmode {i_mode}/{actual_n_modes_to_sum}...")
        pars.ActiveSources.set_active_eigenmode(i_mode)

        results_mode_i = camb.get_results(pars)

        if pars.WantVectors:
            power_spectra_mode_i = results_mode_i.get_cmb_power_spectra(pars, CMB_unit=args.units, raw_cl=False, spectra=("vector",))
            cl_mode_i_dl = power_spectra_mode_i['vector'][:,pol_mode_idx]
        else:
            power_spectra_mode_i = results_mode_i.get_cmb_power_spectra(pars, CMB_unit=args.units, raw_cl=False)
            cl_mode_i_dl = power_spectra_mode_i['total'][:,pol_mode_idx]

        min_len = min(len(cl_strings_sum_dl), len(cl_mode_i_dl))
        cl_strings_sum_dl[:min_len] += cl_mode_i_dl[:min_len]

    pars.ActiveSources.set_active_eigenmode(0)
    print("UETC C_l calculation finished.")

    # Assumes linear addition of power spectra.
    Cl_strings = cl_strings_sum_dl

    return ls_calc, cl_baseline_dl, Cl_strings

def plot_results(ls_calc, cl_baseline_dl, Cl_strings, args):
    """Plot the results"""
    if args.no_plot:
        return

    print("Plotting results...")

    # Use a more compatible matplotlib style
    try:
        plt.style.use('seaborn-v0_8-colorblind')
    except:
        try:
            plt.style.use('seaborn-colorblind')
        except:
            pass  # Use default style if seaborn not available

    fig, axs = plt.subplots(1, 2, figsize=(16, 6))

    plot_mask = (ls_calc >= 1) & (ls_calc <= args.lmax)
    ls_plot = ls_calc[plot_mask]

    # First plot: Strings
    axs[0].plot(ls_plot, Cl_strings[plot_mask]*(2e-7)**2, label=f'Strings ({args.mode})', color='C1', linestyle='-')
    axs[0].set_xlabel(r'$\ell$')
    units_label = r'\mu K^2' if args.units == 'muK' else r'K^2'
    axs[0].set_ylabel(rf'$\ell(\ell+1)C_\ell/2\pi \, [{units_label}]$')
    axs[0].set_title(f'Strings Power Spectrum ({args.pol})')
    axs[0].set_xlim([2, args.lmax])
    axs[0].set_xscale('log')
    axs[0].legend()
    axs[0].grid(True, which="both", ls="-", alpha=0.5)

    # Second plot: Baseline
    axs[1].plot(ls_plot, cl_baseline_dl[plot_mask], label='Baseline', color='C0', linestyle='-')
    axs[1].set_xlabel(r'$\ell$')
    axs[1].set_ylabel(rf'$\ell(\ell+1)C_\ell/2\pi \, [{units_label}]$')
    axs[1].set_title(f'Baseline Power Spectrum ({args.pol})')
    axs[1].set_xlim([2, args.lmax])
    axs[1].set_ylim(ymin=0)
    axs[1].set_xscale('log')
    axs[1].legend()
    axs[1].grid(True, which="both", ls="-", alpha=0.5)

    plt.tight_layout()

    if args.output:
        plt.savefig(args.output)
        print(f"Plot saved to {args.output}")
    else:
        plt.show()

def main():
    """Main function"""
    args = parse_arguments()

    print(f"Configuration:")
    print(f"  Mode: {args.mode}")
    print(f"  Number of modes: {args.nmodes}")
    print(f"  Polarization: {args.pol}")
    print(f"  Max multipole: {args.lmax}")
    print(f"  Units: {args.units}")
    print()

    # Setup CAMB parameters
    pars = setup_camb_params(args)

    # Load correlator data
    correlator_data = load_correlator_data(args)

    # Setup ActiveSources object
    active_sources = setup_active_sources(correlator_data, args)
    pars.ActiveSources = active_sources

    # Calculate power spectra
    ls_calc, cl_baseline_dl, Cl_strings = calculate_power_spectra(pars, args)

    # Plot results
    plot_results(ls_calc, cl_baseline_dl, Cl_strings, args)

    print("--- Script Finished ---")

if __name__ == "__main__":
    main()