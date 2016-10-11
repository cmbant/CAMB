# Make BBN interpolation table from AlterBBN or Parthenelope
# Used AlterBBN 1.4, compiled after changing exit code to 0, and increasing sf of the text output in alter_etannutau
# Fitting formula for Parthenelope via Julien Lesourges Dec 2014

import numpy as np
import subprocess
import sys

hbar = 1.05457e-34
c = 299792458.
kB = 1.380650e-23
MeV = 1.6021764e-13
eV = MeV / 1e6
G = 6.6738e-11
TCMB = 2.7255
mP = np.sqrt(hbar * c / G)
Mpc = 3.0856776e22

m_proton = 1.672621637e-27
m_H = 1.673575e-27
not4 = 3.9715
m_He = m_H * not4

zeta3 = 1.202056903
tau_n = 880.3

n_photon = (kB * TCMB / hbar / c) ** 3 * zeta3 * 2 / np.pi ** 2

omegafac = (1e5 / Mpc) ** 2 / (8 * np.pi * G) * 3

alterBBN = False


def yhe_to_ypBBN(Yp):
    return -4 * m_H * Yp / (Yp * m_He - 4 * Yp * m_H - m_He)


def ypBBN_to_yhe(YBBN):
    return -YBBN * m_He / (-YBBN * m_He + 4 * YBBN * m_H - 4 * m_H)


# Parthenelope fits, from Julien 3 Dec 2014
def yhe_fit(omegab, dneff, taun):
    return (
               0.2311 + 0.9502 * omegab - 11.27 * omegab * omegab
               + dneff * (0.01356 + 0.008581 * omegab - 0.1810 * omegab * omegab)
               + dneff * dneff * (-0.0009795 - 0.001370 * omegab + 0.01746 * omegab * omegab)
           ) * pow(taun / 880.3, 0.728)


def dh_fit(omegab, dneff, taun):
    return (
               18.754 - 1534.4 * omegab + 48656. * omegab * omegab - 552670. * omegab * omegab * omegab
               + dneff * (2.4914 - 208.11 * omegab + 6760.9 * omegab * omegab - 78007. * omegab * omegab * omegab)
               + dneff * dneff * (
                   0.012907 - 1.3653 * omegab + 37.388 * omegab * omegab - 267.78 * omegab * omegab * omegab)
           ) * pow(taun / 880.3, 0.418)


def runBBN(eta, DeltaN=0, tau=tau_n, prog='./alter_etannutau.x'):
    return subprocess.check_output([prog, str(eta), str(3 + DeltaN), str(tau)], shell=False, stderr=subprocess.STDOUT)


# BBN predictions from Parthenope
def ypBBN_Parthenope(omegab, neff):
    # neff here is including 0.046 factor
    return 0.1817 + 0.9020 * omegab - 10.33 * omegab * omegab + (
                                                                    0.01948 + 0.02440 * omegab - 0.4404 * omegab * omegab) * neff + \
           (-0.001002 - 0.002778 * omegab + 0.04719 * omegab * omegab) * neff * neff;


def yhe_to_ypBBN(Yp, ombh2):
    return -4 * m_H * Yp / (Yp * m_He - 4 * Yp * m_H - m_He)


def runBBNplot():
    # standard BBN plot
    etas = 10 ** np.linspace(-12, -8, num=50, endpoint=True)
    lines = []
    for eta in etas:
        outtxt = runBBN(eta, prog='./main.x')
        (YBBN, D, He3, Li7, Li6, Be7) = [float(s.strip()) for s in outtxt.split('\n')[0].split()]
        line = ('%12.5e' * 7) % (eta, YBBN, D, He3, Li7, Li6, Be7)
        print(line)
        lines += [line]
    f = open('bbn_abundance_eta.dat', 'w')
    f.write("\n".join(lines))
    f.close()
    sys.exit()


if __name__ == "__main__":
    # Make interpolation table

    ombh2s = np.hstack((np.linspace(0.005, 0.02, num=int((0.02 - 0.005) / 0.001), endpoint=False),
                        np.linspace(0.020, 0.024, num=16, endpoint=False),
                        np.linspace(0.024, 0.04, num=int((0.04 - 0.024) / 0.001) + 1, endpoint=True)))

    DeltaNs = [-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5, 6, 7]

    lines = []
    for DeltaN in DeltaNs:
        # this is fiducial mass fraction. For low ombh2 has no effect, for subsequent steps use previous value
        Yp_fid = 0.2;
        for ombh2 in ombh2s:
            rho_b = ombh2 * omegafac;
            if not alterBBN:
                YBBN = yhe_fit(ombh2, DeltaN, tau_n)
                Yp_fid = ypBBN_to_yhe(YBBN)
            n_baryon = (4 * Yp_fid / m_He + (1 - Yp_fid) / m_H) * rho_b
            eta = n_baryon / n_photon
            if alterBBN:

                outtxt = runBBN(eta, DeltaN)

                (YBBN, D, He3, Li7, Li6, Be7) = [float(s.strip()) for s in outtxt.split('\n')[7].split()[1:]]
                (sigYBBN, sigdD, sigHe3, sigLi7, sigLi6, sigBe7) = [float(s.strip()) for s in
                                                                    outtxt.split('\n')[8].split()[2:]]

                actual_ombh2 = eta * n_photon * (YBBN * m_He / 4 + (1 - YBBN) * m_H) / omegafac
                Yp = (eta * n_photon / actual_ombh2 / omegafac - 1 / m_H) / (4 / m_He - 1 / m_H)
                Yp_fid = Yp
                #        print ombh2, DeltaN, eta, YBBN, actual_ombh2, Yp
                line = (('%12.5f ') * 6 + ('%12.3e %12.2e') * 3) % (
                    actual_ombh2, eta * 1e10, DeltaN, Yp, YBBN, sigYBBN, D, sigdD, He3, sigHe3, Li7, sigLi7)
                lines += [line]
            else:
                D = dh_fit(ombh2, DeltaN, tau_n) * 1e-5
                Yp = Yp_fid
                # Yp error already includes tau_n error??
                sigdD = 6e-7  # np.sqrt(4.e-07**2 + ((dh_fit(ombh2,DeltaN,tau_n+1.1)-dh_fit(ombh2,DeltaN,tau_n-1.1))/2*1e-5)**2 +  (0.002*D)**2)
                sigYBBN = 0.0003  # np.sqrt(0.0003**2 + ((yhe_fit(ombh2,DeltaN,tau_n+1.1)-yhe_fit(ombh2,DeltaN,tau_n-1.1))/2)**2 + (0.002*Yp)**2)
                line = (('%12.5f ') * 6 + ('%12.3e %12.2e')) % (ombh2, eta * 1e10, DeltaN, Yp, YBBN, sigYBBN, D, sigdD)
                lines += [line]
            print(line)
        lines += ['']

    if alterBBN:
        f = open('BBN_full_alterBBN_' + str(tau_n) + '.dat', 'w')
        f.write('''#BBN prediction of the primordial Helium abundance $Y_p$ as 
#function of the baryon density $\omega_b h^2$ and number of 
#extra radiation degrees of freedom $\Delta N$.
#Calculated with AlterBBN v1.4 [http://superiso.in2p3.fr/relic/alterbbn/] for a 
#neutron lifetime of %.1f s and CMB temperature %.4f K
#Yp^BBN is the BBN-standard nucleon number fraction, Yp is the mass fraction for CMB codes
    
#      ombh2        eta10       DeltaN           Yp       Yp^BBN     sigma_Yp          D/H     err D/H       He3/H    err He3/H         Li7      sig Li7
     
    ''' % (tau_n, TCMB))
        f.write("\n".join(lines))
        f.close()

    else:
        f = open('BBN_full_Parthenelope_' + str(tau_n) + '.dat', 'w')
        f.write('''#BBN prediction of the primordial Helium abundance $Y_p$ as 
#function of the baryon density $\omega_b h^2$ and number of 
#extra radiation degrees of freedom $\Delta N$.
#Calculated from Parthenelope fitting function Dec 2014. Errors are guesstimates. 
#neutron lifetime of %.1f s and CMB temperature %.4f K
#Yp^BBN is the BBN-standard nucleon number fraction, Yp is the mass fraction for CMB codes
    
#      ombh2        eta10       DeltaN           Yp       Yp^BBN     sigma_Yp          D/H     err D/H

''' % (tau_n, TCMB))
        f.write("\n".join(lines))
        f.close()
