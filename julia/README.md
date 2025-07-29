# CAMB.jl

This directory contains a minimal Julia implementation of basic CAMB functionality.
The goal is to provide a pure Julia port that demonstrates high-performance
computations without relying on Fortran code. Currently implemented features
include:

- Definition of cosmological parameters
  (including dark energy parameters ``w0`` and ``wa``,
  neutrino effective number ``Neff``, and total neutrino mass ``mnu``)
 - Eisenstein & Hu no-wiggle transfer function
 - Full Eisenstein & Hu transfer function with baryon oscillations
 - Linear growth factor via numerical integration
 - Growth rate ``f(a)`` computed with automatic differentiation
 - Linear matter power spectrum using the Eisenstein & Hu transfer function
- Non-linear matter power spectrum via a simplified Halofit fit
- Basic distance measures (comoving radial, angular diameter, luminosity)
- Distance modulus for convenient magnitude calculations
- Hubble parameter as a function of redshift
- Lookback and age of the universe calculations
 - CPL dark energy evolution via parameter ``wa``
 - Radiation density from photons, massless neutrinos, and a massive
   neutrino component via ``mnu``
- Crude recombination history calculation
- Approximate CMB temperature power spectrum
- Comoving sound horizon at decoupling
- Drag epoch redshift and sound horizon
- RMS density fluctuations `sigma_R` and `sigma8`
- Critical density calculation (optionally at redshift)
- Hubble distance helper
- Basic math utilities: Romberg integration, root finding via `Roots.jl`,
  and Gauss-Legendre quadrature
 - Convenience helpers for density fractions (`Omega_de`, `Omega_m_z`,
   `Omega_de_z`, and `Omega_nu`)
- Physical constants translated from the CAMB Fortran code
- Basic configuration module with global settings and error handling
- Simple initial power spectrum module providing scalar and tensor
  power spectra
- Parameter structs (`TransferParams`, `AccuracyParams`, `CAMBparams`)
  with a helper to construct defaults
- Results container type for storing background outputs using in-place
   vectorized loops and convenience wrappers for distance modulus and
   fluctuation calculations
- Simple reionization module with a hyperbolic tangent model and
  optical depth calculation
- Spherical Bessel functions and derivatives via `SpecialFunctions`
- Background caching with cubic-spline interpolation for fast distance
  and time queries

This is a proof-of-concept and does not yet match the full capabilities of the
original CAMB code.

## Running the tests

A simple test suite is included under `test/`. From this directory run

```bash
julia --project=. test/runtests.jl
```

to exercise the implemented functions.

## Generating example plots

An example script is provided under `examples/` that computes several
cosmological quantities and saves plots to the current directory.
Run it with:

```bash
# from the repository root
julia --project=julia julia/examples/plots.jl

# or from within the `julia` directory
julia --project=. examples/plots.jl
```

This will produce PNG files of background distances, the distance modulus,
lookback time, the growth factor, linear and nonlinear matter power spectra,
a simple recombination history, and an approximate CMB temperature power
spectrum.
