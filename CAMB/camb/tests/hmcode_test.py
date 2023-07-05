import os
import sys
import unittest
import numpy as np

try:
    import camb
except ImportError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
    import camb


class HMcodeTest(unittest.TestCase):

    def testHMcode(self):

        # Parameters
        # TODO: Not clear that the neutrino setup here agrees exactly with my HMx setup
        # TODO: Larger errors than I would like with high-massive-nu and high-baryon cosmologies
        # TODO: Normalisation via As not agreeing exactly with sigma_8 values
        kmax_calc = 2e2  # Maximum wavenumber for the CAMB calculation [h/Mpc]
        pi = np.pi  # Lovely pi
        twopi = 2. * pi  # Twice as lovely pi
        neutrino_number = 94.14  # Neutrino closure mass [eV] TODO: VERY LAZY (should it be 93.14, 94.07eV?)
        eps_k = 1e-6  # Fractional error tolerated in k
        eps_a = 1e-6  # Fractional error tolerated in a
        eps_Pk = 5e-3  # Fractional error tolerated in Pk # TODO: This should agree to 1e-3 as in HMx
        kmin_test = 1e-2  # Minimum wavenumber for test [h/Mpc]
        kmax_test = 1e1  # Maximum wavenumber for test [h/Mpc]
        amin_test = 0.333  # Minimum scale factor for test
        amax_test = 1.000  # Maximum scale factor for test
        Y_He = 0.24  # Helium fraction
        T_CMB = 2.725  # CMB temperature
        nnu = 3  # Number of massive neutrinos
        neff = 3.046  # Effective number of radiative neutrinos
        verbose = False  # Verbose tests
        norm_sig8 = False  # Normalisation by sigma_8 or As?
        a_in = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])  # Scale factors for tests

        # Dictionary for HMcode versions
        HMcode_version_ihm = {
            'mead': 51,
            'mead2016': 51,
            'mead2015': 66,
            'mead2020': 123,
            'mead2020_feedback': 124
        }

        # Read in and sort Mead benchmark datawi
        def read_Mead_benchmark(infile):

            # Read data and split into k, a and Pk
            data = np.loadtxt(infile)
            k = data[:, 0]
            a = a_in
            D2 = data[:, 1:]
            Pk = (D2.T / (4. * pi * (k / twopi) ** 3)).T  # Convert Delta^2(k) -> P(k)

            # Return results
            return k, a, Pk

        def HMcode_test_cosmologies(icos):

            # Common test parameters
            Om_m = 0.30
            Om_b = 0.05
            h = 0.7
            ns = 0.96
            sig8 = 0.8
            As = 1.97269e-9
            w0 = -1.
            wa = 0.
            mnu = 0.
            logT = 7.8

            # Uncommon test parameters
            if icos == 56:
                Om_m = 0.3158
                Om_b = 0.04939
                h = 0.6732
                ns = 0.96605
                sig8 = 0.8120
                As = 2.08924e-9
                w0 = -1.
                wa = 0.
                mnu = 0.06
            elif icos == 241:
                w0 = -0.7
                As = 2.46981e-9
            elif icos == 242:
                w0 = -1.3
                As = 1.75070e-9
            elif icos == 243:
                mnu = 0.3
                As = 2.35868e-9
            elif icos == 244:
                mnu = 0.9
                As = 3.43507e-9
            elif icos == 245:
                ns = 0.7
                As = 2.38515e-9
            elif icos == 246:
                ns = 1.3
                As = 1.47601e-9
            elif icos == 247:
                Om_b = 0.01
                As = 1.20028e-9
            elif icos == 248:
                Om_b = 0.1
                As = 3.88822e-9
            elif icos == 249:
                wa = 0.9
                As = 3.35538e-9
            elif icos == 250:
                logT = 7.6
            elif icos == 251:
                logT = 8.0

            # Return cosmological parameters
            return Om_m, Om_b, h, ns, sig8, As, w0, wa, mnu, logT

        # Set up a CAMB parameter set for a Mead cosmology
        def setup_HMcode_test(a, icos):

            # Get the cosmological parameters for the cosmology
            Om_m, Om_b, h, ns, sig8, As, w0, wa, mnu, logT = HMcode_test_cosmologies(icos)

            # Redshifts
            z = -1. + 1. / a

            # Derive and set CAMB parameters
            H0 = 100. * h
            wnu = mnu / neutrino_number
            Om_nu = wnu / h ** 2
            Om_c = Om_m - Om_b - Om_nu
            Om_k = 0.
            wb = Om_b * h ** 2
            wc = Om_c * h ** 2

            # Set parameters using the traditional Fortran language
            # TODO: Check carefully against my ini files
            pars = camb.CAMBparams(WantCls=False, H0=H0, ombh2=wb, omch2=wc, omnuh2=wnu, omk=Om_k, YHe=Y_He, TCMB=T_CMB,
                                   num_nu_massive=nnu,
                                   num_nu_massless=neff - nnu,
                                   nu_mass_numbers=[nnu],
                                   nu_mass_degeneracies=[neff / nnu],
                                   nu_mass_fractions=[1.],
                                   share_delta_neff=True,
                                   )
            pars.set_dark_energy(w=w0, wa=wa, dark_energy_model='ppf')
            pars.InitPower.set_params(As=As, ns=ns)

            # Set sigma_8 normalisation (wasteful)
            if norm_sig8:
                pars.set_matter_power(redshifts=[0.], kmax=kmax_calc)
                pars.NonLinear = camb.model.NonLinear_none
                results = camb.get_results(pars)
                sig8_init = results.get_sigma8()
                As = As * (sig8 / sig8_init) ** 2
                pars.InitPower.set_params(As=As, ns=ns)  # Reset As

            # Main calculation step
            pars.set_matter_power(redshifts=z, kmax=kmax_calc)
            pars.NonLinear = camb.model.NonLinear_both
            results = camb.get_results(pars)
            if verbose:
                print(pars)

            # Return results and AGN temperature (ugly, but it needs to come out)
            return results, logT

        # Get the HMcode power from CAMB
        def get_HMcode_power_from_CAMB(results, k_in, a_in, logT, HMcode_version):

            # k and z ranges for results
            kmin = k_in[0]
            kmax = k_in[-1]
            nk = len(k_in)

            # Get non-linear spectra
            results.Params.NonLinearModel.set_params(halofit_version=HMcode_version, HMCode_logT_AGN=logT)
            k, z, Pk = results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints=nk)
            Pk = Pk.T[:, ::-1]
            z = np.array(z)[::-1]
            a = 1. / (1. + z)
            if verbose:
                print('sigma_8:', results.get_sigma8()[-1])

            # Return non-linear power
            return k, a, Pk

        # Input file name
        def HMcode_benchmark_file(icos, ihm):
            return 'HMcode_test_outputs/HMcode_cos%d_hm%d.txt' % (icos, ihm)

        # Whitespace
        if verbose:
            print('')

        # Loop over cosmologies
        for icos in [26, 56, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251]:

            # Init CAMB
            results, logT = setup_HMcode_test(a_in, icos)

            # Loop over HMcode versions
            for HMcode_version in ['mead2015', 'mead2016', 'mead2020', 'mead2020_feedback']:

                # Read benchmark data
                ihm = HMcode_version_ihm[HMcode_version]
                infile = HMcode_benchmark_file(icos, ihm)
                if verbose:
                    print('Infile:', infile)
                k_in, a_in, Pk_in = read_Mead_benchmark(infile)

                # Get power from CAMB
                k_nl, a_nl, Pk_nl = get_HMcode_power_from_CAMB(results, k_in, a_in, logT, HMcode_version)

                # Compare benchmark to calculation
                for ik in range(len(k_in)):
                    self.assertAlmostEqual(k_nl[ik] / k_in[ik], 1., delta=eps_k)
                for ia in range(len(a_in)):
                    self.assertAlmostEqual(a_nl[ia] / a_in[ia], 1., delta=eps_a)
                for ia in range(len(a_in)):
                    for ik in range(len(k_in)):
                        if kmin_test <= k_in[ik] <= kmax_test and amin_test <= a_in[ia] <= amax_test:
                            self.assertAlmostEqual(Pk_nl[ik, ia] / Pk_in[ik, ia], 1., delta=eps_Pk)
