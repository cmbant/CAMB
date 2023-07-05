.. _transfer-variables:

Matter power spectrum and matter transfer function variables
=============================================================

The various matter power spectrum functions, e.g. :func:`~camb.get_matter_power_interpolator`, can calculate power
spectra for various quantities. Each variable used to form the power spectrum has a name as follows:

=====================  ======  =====================================================================
name                   number  description
=====================  ======  =====================================================================
k/h                      1     :math:`k/h`
delta_cdm                2     :math:`\Delta_c`, CDM density
delta_baryon             3     :math:`\Delta_b`, baryon density
delta_photon             4     :math:`\Delta_\gamma`, photon density
delta_neutrino           5     :math:`\Delta_r`, for massless neutrinos
delta_nu                 6     :math:`\Delta_\nu` for massive neutrinos
delta_tot                7     :math:`\frac{\rho_c\Delta_c+\rho_b\Delta_b+\rho_\nu\Delta_\nu}{\rho_c+\rho_b+\rho_\nu}`,
                               CDM+baryons+massive neutrino density
delta_nonu               8     :math:`\frac{\rho_c\Delta_c+\rho_b\Delta_b}{\rho_b+\rho_c}`, CDM+baryon  density
delta_tot_de             9     :math:`\frac{\rho_c\Delta_c+\rho_b\Delta_b+\rho_\nu\Delta_\nu +\rho_{ de}\Delta_{de}}{\rho_c+\rho_b+\rho_\nu}`,
                               CDM+baryons+massive neutrinos+ dark energy (numerator only)  density
Weyl                    10     :math:`k^2\Psi\equiv k^2(\phi+\psi)/2`,
                               the Weyl potential scaled by :math:`k^2` to scale in :math:`k` like a density.
v_newtonian_cdm         11     :math:`-v_{N,c}\, k/{\cal H}` (where :math:`v_{N,c}` is the
                               Newtonian-gauge CDM velocity)
v_newtonian_baryon      12     :math:`-v_{N,b}\,k/{\cal H}` (Newtonian-gauge baryon velocity :math:`v_{N,b}`)
v_baryon_cdm            13     :math:`v_b-v_c`, relative baryon-CDM velocity
=====================  ======  =====================================================================

The number here corresponds to a corresponding numerical index, in Fortran these are the same as *model.name*,
where *name* are the Transfer_xxx variable names:
Transfer_kh=1,Transfer_cdm=2, Transfer_b=3, Transfer_g=4, Transfer_r=5, Transfer_nu=6, Transfer_tot=7,
Transfer_nonu=8, Transfer_tot_de=9, Transfer_Weyl=10, Transfer_Newt_vel_cdm=11, Transfer_Newt_vel_baryon=12,
Transfer_vel_baryon_cdm = 13.

So for example, requesting var1='delta_b', var2='Weyl' or alternatively var1=model.Transfer_b, var2=model.Transfer_Weyl
would get the power spectrum for the cross-correlation of the baryon density with the Weyl potential.
All density variables :math:`\Delta_i` here are synchronous gauge.

For transfer function variables (rather than matter power spectra), the variables are normalized corresponding to
unit primordial curvature perturbation on super-horizon scales. The
:meth:`~camb.results.CAMBdata.get_matter_transfer_data` function returns the above quantities
divided by :math:`k^2` (so they are roughly constant at low :math:`k` on super-horizon scales).

The  `example notebook <https://camb.readthedocs.io/en/latest/CAMBdemo.html>`_  has various examples of getting the
matter power spectrum, relating the Weyl-potential spectrum to lensing, and calculating the
baryon-dark matter relative velocity spectra. There is also an explicit example of how to calculate the matter
power spectrum manually from the matter transfer functions.

When generating dark-age 21cm power spectra (do21cm is set) the transfer functions are instead the *model.name*
variables (see equations 20 and 25 of `astro-ph/0702600 <https://arxiv.org/abs/astro-ph/0702600>`_)

========================  ======  =====================================================================
name                      number  description
========================  ======  =====================================================================
Transfer_kh                 1     :math:`k/h`
Transfer_cdm                2     :math:`\Delta_c`, CDM density
Transfer_b                  3     :math:`\Delta_b`, baryon density
Transfer_monopole           4     :math:`\Delta_s+(r_\tau-1)(\Delta_{b}-\Delta_{T_s})`, 21cm monopole source
Transfer_vnewt              5     :math:`r_\tau kv_{N,b}/\mathcal{H}`, 21cm Newtonian-gauge velocity source
Transfer_Tmat               6     :math:`\Delta_{T_m}`, matter temperature perturbation
Transfer_tot                7     :math:`\frac{\rho_c\Delta_c+\rho_b\Delta_b+\rho_\nu\Delta_\nu}{\rho_c+\rho_b+\rho_\nu}`,
                                  CDM+baryons+massive neutrino density
Transfer_nonu               8     :math:`\frac{\rho_c\Delta_c+\rho_b\Delta_b}{\rho_b+\rho_c}`, CDM+baryon  density
Transfer_tot_de             9     :math:`\frac{\rho_c\Delta_c+\rho_b\Delta_b+\rho_\nu\Delta_\nu +\rho_{ de}\Delta_{de}}{\rho_c+\rho_b+\rho_\nu}`,
                                  CDM+baryons+massive neutrinos+ dark energy (numerator only)  density
Transfer_Weyl              10     :math:`k^2\Psi\equiv k^2(\phi+\psi)/2`,
                                  the Weyl potential scaled by :math:`k^2` to scale in :math:`k` like a density.
Transfer_Newt_vel_cdm      11     :math:`-v_{N,c}\, k/{\cal H}` (where :math:`v_{N,c}` is the
                                  Newtonian-gauge CDM velocity)
Transfer_Newt_vel_baryon   12     :math:`-v_{N,b}\,k/{\cal H}` (Newtonian-gauge baryon velocity :math:`v_{N,b}`)
Transfer_vel_baryon_cdm    13     :math:`v_b-v_c`, relative baryon-CDM velocity
========================  ======  =====================================================================

If use_21cm_mK is set the 21cm results are multiplied by :math:`T_b` to give results in mK units.
