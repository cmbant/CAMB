from ctypes import *
from baseconfig import camblib, CAMB_Structure, CAMBError
import reionization as ion
import recombination as recomb
import initialpower as ipow
import constants
import numpy as np
import bbn
import math

# ---Parameters

max_nu = 5
max_transfer_redshifts = 150
nthermo_derived = 13
Transfer_kh = 1
Transfer_cdm = 2
Transfer_b = 3
Transfer_g = 4
Transfer_r = 5
Transfer_nu = 6
Transfer_tot = 7
Transfer_nonu = 8
Transfer_tot_de = 9
Transfer_Weyl = 10
Transfer_Newt_vel_cdm = 11
Transfer_Newt_vel_baryon = 12
Transfer_vel_baryon_cdm = 13
Transfer_max = Transfer_vel_baryon_cdm

NonLinear_none = 0
NonLinear_pk = 1
NonLinear_lens = 2
NonLinear_both = 3

derived_names = ['age', 'zstar', 'rstar', 'thetastar', 'DAstar', 'zdrag',
                 'rdrag', 'kd', 'thetad', 'zeq', 'keq', 'thetaeq', 'thetarseq']

# ---Variables in modules.f90
# To set the value please just put 
# variable_name.value = new_value

DebugParam = c_double.in_dll(camblib, "__modelparams_MOD_debugparam")
# DebugParam.value = 1000000*2

# logical
do_bispectrum = c_bool.in_dll(camblib, "__modelparams_MOD_do_bispectrum")
# do_bispectrum.value = False

max_bessels_l_index = c_int.in_dll(camblib, "__modelparams_MOD_max_bessels_l_index")
# max_bessels_l_index.value = 1000000

max_bessels_etak = c_double.in_dll(camblib, "__modelparams_MOD_max_bessels_etak")
# max_bessels_etak.value = 1000000*2

# logical
call_again = c_bool.in_dll(camblib, "__modelparams_MOD_call_again")
# call_again.value = False

grhom = c_double.in_dll(camblib, "__modelparams_MOD_grhom")
grhog = c_double.in_dll(camblib, "__modelparams_MOD_grhog")
grhor = c_double.in_dll(camblib, "__modelparams_MOD_grhor")
grhob = c_double.in_dll(camblib, "__modelparams_MOD_grhob")
grhoc = c_double.in_dll(camblib, "__modelparams_MOD_grhoc")
grhov = c_double.in_dll(camblib, "__modelparams_MOD_grhov")
grhornomass = c_double.in_dll(camblib, "__modelparams_MOD_grhornomass")
grhok = c_double.in_dll(camblib, "__modelparams_MOD_grhok")

taurst = c_double.in_dll(camblib, "__modelparams_MOD_taurst")
dtaurec = c_double.in_dll(camblib, "__modelparams_MOD_dtaurec")
taurend = c_double.in_dll(camblib, "__modelparams_MOD_taurend")
tau_maxvis = c_double.in_dll(camblib, "__modelparams_MOD_tau_maxvis")
adotrad = c_double.in_dll(camblib, "__modelparams_MOD_adotrad")

c_array = c_double * max_nu
grhormass = c_array.in_dll(camblib, "__modelparams_MOD_grhormass")
nu_masses = c_array.in_dll(camblib, "__modelparams_MOD_nu_masses")

akthom = c_double.in_dll(camblib, "__modelparams_MOD_akthom")
fHe = c_double.in_dll(camblib, "__modelparams_MOD_fhe")
Nnow = c_double.in_dll(camblib, "__modelparams_MOD_nnow")

limber_phiphi = c_int.in_dll(camblib, "__modelparams_MOD_limber_phiphi")
# limber_phiphi.value = 0

num_extra_redshiftwindows = c_int.in_dll(camblib, "__modelparams_MOD_num_extra_redshiftwindows")
# num_extra_redshiftwindows.value = 0

num_redshiftwindows = c_int.in_dll(camblib, "__modelparams_MOD_num_redshiftwindows")

# logical
use_spline_template = c_bool.in_dll(camblib, "__modelparams_MOD_use_spline_template")
# use_spline_template.value = True

c_array2 = c_double * nthermo_derived
ThermoDerivedParams = c_array2.in_dll(camblib, "__modelparams_MOD_thermoderivedparams")
# ThermoDerivedParams.value = 1.

# logical
Log_lvalues = c_bool.in_dll(camblib, "__lvalues_MOD_log_lvalues")
# Log_lvalues.value = False

# Variables from module ModelData

# logical
has_cl_2D_array = c_bool.in_dll(camblib, "__modeldata_MOD_has_cl_2d_array")
# has_cl_2D_array.value = False


lmax_lensed = c_int.in_dll(camblib, "__modeldata_MOD_lmax_lensed")

# Variable from module Transfer
# logical
transfer_interp_matterpower = c_bool.in_dll(camblib, "__transfer_MOD_transfer_interp_matterpower")
# transfer_interp_matterpower.value = False

transfer_power_var = c_int.in_dll(camblib, "__transfer_MOD_transfer_power_var")
# transfer_power_var.value = Transfer_tot

# logical
get_growth_sigma8 = c_bool.in_dll(camblib, "__transfer_MOD_get_growth_sigma8")
# get_growth_sigma8.value = True

CAMB_validateparams = camblib.__camb_MOD_camb_validateparams
CAMB_validateparams.restype = c_bool

CAMB_setinitialpower = camblib.__handles_MOD_camb_setinitialpower


# ---Derived Types in modules.f90

class TransferParams(CAMB_Structure):
    _fields_ = [
        ("high_precision", c_int),  # logical
        ("num_redshifts", c_int),
        ("kmax", c_double),
        ("k_per_logint", c_int),
        ("redshifts", c_double * max_transfer_redshifts),
        ("PK_redshifts", c_double * max_transfer_redshifts),
        ("NLL_redshifts", c_double * max_transfer_redshifts),
        ("PK_redshifts_index", c_int * max_transfer_redshifts),
        ("NLL_redshifts_index", c_int * max_transfer_redshifts),
        ("PK_num_redshifts", c_int),
        ("NLL_num_redshifts", c_int)
    ]


class CAMBparams(CAMB_Structure):
    def __init__(self):
        getattr(camblib, '__camb_MOD_camb_setdefparams')(byref(self))

    _fields_ = [
        ("WantCls", c_int),  # logical
        ("WantTransfer", c_int),  # logical
        ("WantScalars", c_int),  # logical
        ("WantTensors", c_int),  # logical
        ("WantVectors", c_int),  # logical
        ("DoLensing", c_int),  # logical
        ("want_zstar", c_int),  # logical
        ("want_zdrag", c_int),  # logical
        ("PK_WantTransfer", c_int),  # logical
        ("NonLinear", c_int),
        ("Want_CMB", c_int),  # logical
        ("max_l", c_int),
        ("max_l_tensor", c_int),
        ("max_eta_k", c_double),
        ("max_eta_k_tensor", c_double),
        ("omegab", c_double),
        ("omegac", c_double),
        ("omegav", c_double),
        ("omegan", c_double),
        ("H0", c_double),
        ("TCMB", c_double),
        ("YHe", c_double),
        ("num_nu_massless", c_double),
        ("num_nu_massive", c_int),
        ("nu_mass_eigenstates", c_int),
        ("share_delta_neff", c_int),  # logical
        ("nu_mass_degeneracies", c_double * max_nu),
        ("nu_mass_fractions", c_double * max_nu),
        ("nu_mass_numbers", c_int * max_nu),
        ("scalar_initial_condition", c_int),
        ("OutputNormalization", c_int),
        ("AccuratePolarization", c_int),  # logical
        ("AccurateBB", c_int),  # logical
        ("AccurateReionization", c_int),  # logical
        ("MassiveNuMethod", c_int),
        ("InitPower", ipow.InitialPowerParams),
        ("Reion", ion.ReionizationParams),
        ("Recomb", recomb.RecombinationParams),
        ("Transfer", TransferParams),
        ("InitialConditionVector", c_double * 10),
        ("OnlyTransfers", c_int),  # logical
        ("DerivedParameters", c_int),  # logical
        ("ReionHist", ion.ReionizationHistory),
        ("flat", c_int),  # logical
        ("closed", c_int),  # logical
        ("open", c_int),  # logical
        ("omegak", c_double),
        ("curv", c_double),
        ("r", c_double),
        ("Ksign", c_double),
        ("tau0", c_double),
        ("chi0", c_double)
    ]

    def validate(self):
        return CAMB_validateparams(byref(self))

    def set_initial_power(self, initial_power_params):
        assert (isinstance(initial_power_params, ipow.InitialPowerParams))
        CAMB_setinitialpower(byref(self), byref(initial_power_params))

    def set_bbn_helium(self, ombh2, delta_nnu, tau_neutron):
        Yp = bbn.yhe_fit(ombh2, delta_nnu, tau_neutron)
        self.YHe = bbn.ypBBN_to_yhe(Yp)

    def set_cosmology(self, H0=67, ombh2=0.022, omch2=0.12, omk=0.0, num_massive_neutrinos=1, mnu=0.06, nnu=3.046,
                      YHe=None, meffsterile=0, standard_neutrino_neff=3.046, TCMB=constants.COBE_CMBTemp,tau = None,
                      tau_neutron=bbn.tau_n):

        if YHe is None:
            # use BBN prediction
            self.set_bbn_helium(ombh2, nnu - standard_neutrino_neff, tau_neutron)
        else:
            self.YHe = YHe
        self.TCMB = TCMB
        self.H0 = H0
        fac = (self.H0 / 100.0) ** 2
        self.omegab = ombh2 / fac
        self.omegac = omch2 / fac

        neutrino_mass_fac = 94.07
        # conversion factor for thermal with Neff=3 TCMB=2.7255

        omnuh2 = mnu / neutrino_mass_fac * (standard_neutrino_neff / 3.0) ** 0.75
        omnuh2_sterile = meffsterile / neutrino_mass_fac
        if omnuh2_sterile > 0 and nnu < standard_neutrino_neff:
            raise CAMBError('sterile neutrino mass required Neff>3.046')

        omnuh2 = omnuh2 + omnuh2_sterile
        self.omegan = omnuh2 / fac
        self.omegav = 1 - omk - self.omegab - self.omegac - self.omegan
        self.share_delta_neff = False
        self.nu_mass_eigenstates = 0
        self.num_nu_massless = nnu
        self.nu_mass_numbers[0] = 0
        self.num_nu_massive = 0
        if omnuh2 > omnuh2_sterile:
            neff_massive_standard = num_massive_neutrinos * standard_neutrino_neff / 3.0
            self.num_nu_massive = num_massive_neutrinos
            self.nu_mass_eigenstates = self.nu_mass_eigenstates + 1
            if nnu > neff_massive_standard:
                self.num_nu_massless = nnu - neff_massive_standard
            else:
                self.num_nu_massless = 0
                neff_massive_standard = nnu

            self.nu_mass_numbers[self.nu_mass_eigenstates - 1] = num_massive_neutrinos
            self.nu_mass_degeneracies[self.nu_mass_eigenstates - 1] = neff_massive_standard
            self.nu_mass_fractions[self.nu_mass_eigenstates - 1] = (omnuh2 - omnuh2_sterile) / omnuh2
        else:
            neff_massive_standard = 0
        if omnuh2_sterile > 0:
            if nnu < standard_neutrino_neff:
                raise CAMBError('nnu < 3.046 with massive sterile')
            self.num_nu_massless = standard_neutrino_neff - neff_massive_standard
            self.num_nu_massive = self.num_nu_massive + 1
            self.nu_mass_eigenstates = self.nu_mass_eigenstates + 1
            self.nu_mass_numbers[self.nu_mass_eigenstates - 1] = 1
            self.nu_mass_degeneracies[self.nu_mass_eigenstates - 1] = max(1e-6, nnu - standard_neutrino_neff)
            self.nu_mass_fractions[self.nu_mass_eigenstates - 1] = omnuh2_sterile / omnuh2

        if tau is not None:
            self.Reion.set_tau(tau)

    def set_dark_energy(self, w=-1.0, sound_speed=1.0, dark_energy_model='fluid'):
        # Variables from module LambdaGeneral
        if dark_energy_model != 'fluid':
            raise CAMBError('This version only supports the fluid energy model')
        if w != -1 or sound_speed != 1:
            print('Warning: currently dark energy parameters are changed globally, not per parameter set')
        w_lam = c_double.in_dll(camblib, "__lambdageneral_MOD_w_lam")
        w_lam.value = w
        cs2_lam = c_double.in_dll(camblib, "__lambdageneral_MOD_cs2_lam")
        cs2_lam.value = sound_speed

    def get_omega_k(self):
        return 1 - self.omegab - self.omegac - self.omegan - self.omegav

    def set_matter_power(self, redshifts=[0.], kmax=1.2, k_per_logint=None):

        self.WantTransfer = True
        self.Transfer.high_precision = True
        self.Transfer.kmax = kmax
        if not k_per_logint:
            self.Transfer.k_per_logint = 0
        else:
            self.Transfer.k_per_logint = k_per_logint
        zs = sorted(redshifts, reverse=True)
        if np.any(np.array(zs) - np.array(redshifts) != 0):
            print "Note: redshifts have been re-sorted (earliest first)"
        self.Transfer.PK_num_redshifts = len(redshifts)
        for i, z in enumerate(zs):
            self.Transfer.PK_redshifts[i] = z

    def set_for_lmax(self, lmax, max_eta_k=None, lens_potential_accuracy=0,
                     lens_margin=150, k_eta_fac=2.5, lens_k_eta_reference=18000.0):
        if self.DoLensing:
            self.max_l = lmax + lens_margin
        else:
            self.max_l = lmax
        self.max_eta_k = max_eta_k or min(self.max_l, 3000) * k_eta_fac
        if lens_potential_accuracy:
            if self.NonLinear == NonLinear_none:
                self.NonLinear = NonLinear_lens
            else:
                self.NonLinear = NonLinear_both
            self.max_eta_k = max(self.max_eta_k, lens_k_eta_reference * lens_potential_accuracy)


def Transfer_SetForNonlinearLensing(P):
    camblib.__transfer_MOD_transfer_setfornonlinearlensing(byref(P))


def Transfer_SortAndIndexRedshifts(P):
    camblib.__transfer_MOD_transfer_sortandindexredshifts(byref(P))
