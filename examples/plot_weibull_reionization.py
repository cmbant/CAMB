#!/usr/bin/env python
"""
Example script to plot ionization history using the Weibull reionization model.

This script demonstrates how to use the WeibullReionization model in CAMB
and compares it with the standard TanhReionization model.

The Weibull reionization model is based on arXiv:2505.15899v1 and provides
a flexible parameterization of the reionization history with three parameters:
- reion_redshift_complete: redshift where reionization is complete (5% neutral)
- reion_duration: duration parameter (Delta z_90)
- reion_asymmetry: asymmetry parameter (A_z)

Author: CAMB team
"""

import numpy as np
import matplotlib.pyplot as plt
import camb
from camb.reionization import WeibullReionization, TanhReionization


def get_ionization_history(reion_model, z_array):
    """
    Get ionization history for a given reionization model.
    
    Parameters:
    -----------
    reion_model : camb.reionization object
        The reionization model to use
    z_array : array_like
        Array of redshifts to evaluate
        
    Returns:
    --------
    xe_values : array
        Ionization fraction at each redshift
    results : camb.CAMBdata
        CAMB results object
    """
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
    pars.Reion = reion_model
    
    # Get background evolution
    results = camb.get_background(pars)
    
    # Get ionization fraction for each redshift
    xe_values = []
    for z in z_array:
        try:
            # Get background variables at this redshift
            bg_vars = results.get_background_time_evolution([results.conformal_time(z)], ['x_e'])
            xe_values.append(bg_vars[0, 0])
        except:
            # Fallback for edge cases
            xe_values.append(1.0 if z < 6 else 0.0)
    
    return np.array(xe_values), results


def plot_reionization_comparison():
    """
    Create a comparison plot of Weibull vs Tanh reionization models.
    """
    # Create redshift array
    z_range = np.linspace(0, 15, 200)
    
    # Configure Weibull reionization with typical parameters
    weibull_reion = WeibullReionization()
    weibull_reion.set_extra_params(
        reion_redshift_complete=5.8,  # z where reionization completes
        reion_duration=1.0,           # Duration parameter
        reion_asymmetry=1.5           # Asymmetry parameter
    )
    weibull_reion.set_tau(0.065)
    
    # Configure Tanh reionization for comparison
    tanh_reion = TanhReionization()
    tanh_reion.set_tau(0.065)
    
    # Get ionization histories
    xe_weibull, results_weibull = get_ionization_history(weibull_reion, z_range)
    xe_tanh, results_tanh = get_ionization_history(tanh_reion, z_range)
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    
    # Plot Weibull model
    plt.plot(z_range, xe_weibull, 'b-', linewidth=2, 
             label=f'Weibull (τ={results_weibull.Params.Reion.optical_depth:.3f}, z_re={results_weibull.Params.Reion.redshift:.1f})')
    
    # Plot Tanh model for comparison
    plt.plot(z_range, xe_tanh, 'r--', linewidth=2,
             label=f'Tanh (τ={results_tanh.Params.Reion.optical_depth:.3f}, z_re={results_tanh.Params.Reion.redshift:.1f})')
    
    # Add parameter annotations
    plt.text(0.05, 0.95, 'Weibull Parameters:\n' + 
             f'z_complete = {weibull_reion.reion_redshift_complete}\n' +
             f'duration = {weibull_reion.reion_duration}\n' +
             f'asymmetry = {weibull_reion.reion_asymmetry}',
             transform=plt.gca().transAxes, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.xlabel('Redshift z', fontsize=12)
    plt.ylabel('Ionization fraction x_e', fontsize=12)
    plt.title('Reionization History: Weibull vs Tanh Models', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    plt.xlim(0, 15)
    plt.ylim(0, 1.1)
    
    # Invert x-axis to show time evolution (higher z = earlier time)
    plt.gca().invert_xaxis()
    
    plt.tight_layout()
    return plt.gcf()


if __name__ == "__main__":
    # Create comparison plot
    print("Creating Weibull vs Tanh comparison plot...")
    fig1 = plot_reionization_comparison()
    fig1.savefig('weibull_reionization_history.png', dpi=150, bbox_inches='tight')
    print("Saved: weibull_reionization_history.png")
    
    print("Plot generation complete!")
