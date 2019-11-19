"""

This module defines the scalar linear perturbation equations for standard LCDM cosmology, using sympy.
It uses the covariant perturbation notation, but includes functions to project into the
Newtonian or synchronous gauge, as well as constructing general gauge invariant quantities.
It uses "t" as the conformal time variable (=tau in the fortran code).

For a guide to usage and content see the `ScalEqs notebook <http://camb.readthedocs.org/en/latest/ScalEqs.html>`_

As well as defining standard quantities, and how they map to CAMB variables, there are also functions for
converting a symbolic expression to CAMB source code, and compiling custom sources for use with CAMB
(as used by :meth:`.model.CAMBparams.set_custom_scalar_sources`, :meth:`.results.CAMBdata.get_time_evolution`)

A Lewis July 2017
"""

import ctypes
import os
import six
import sympy
from sympy import diff, Eq, simplify, Function, Symbol
from sympy.abc import t, kappa, K, k

# Background variables

tau0 = Symbol('tau0', description='conformal time today')
tau_maxvis = Symbol('tau_maxvis', description='conformal time of peak visibility')

f_K = Function('f_K', description='comoving angular diameter distance')

H = Function('H', description='comoving hubble parameter', camb_var='adotoa')(t)
rho = Function('rho', description='total density', camb_var='grho', camb_sub='grho/kappa/a**2')(t)
P = Function('P', description='total pressure', camb_var='gpres', camb_sub='gpres/kappa/a**2')(t)
a = Function('a', description='scale factor')(t)

rho_b = Function('rho_b', description='baryon density', camb_var='grhob_t', camb_sub='grhob_t/kappa/a**2')(t)
rho_c = Function('rho_c', description='CDM density', camb_var='grhoc_t', camb_sub='grhoc_t/kappa/a**2')(t)
rho_g = Function('rho_g', description='photon density', camb_var='grhog_t', camb_sub='grhog_t/kappa/a**2')(t)
rho_r = Function('rho_r', description='massless neutrino density', camb_var='grhor_t', camb_sub='grhor_t/kappa/a**2')(t)
rho_nu = Function('rho_nu', description='massive neutrino density', camb_var='grhonu_t',
                  camb_sub='grhonu_t/kappa/a**2')(t)
rho_de = Function('rho_de', description='dark energy density', camb_var='grhov_t', camb_sub='grhov_t/kappa/a**2')(t)

p_b = Function('p_b', description='baryon pressure', camb_sub='0')(t)
p_nu = Function('p_nu', description='massive neutrino pressure')(t)
w_de = Function('w_de', camb_var='w_lam', description='fluid dark energy equation of state')(t)
p_g = rho_g / 3
p_r = rho_r / 3
p_c = 0
p_de = w_de * rho_de

opacity = Function('opacity', description='opacity, a n_e sigma_t')(t)
visibility = Function('visibility', description='ionization visibility')(t)
exptau = Function('exptau', description='exp(-tau)')(t)


def subs(eqs, expr):
    # generalization to act on lists of equations, and support using Eq equations.
    # lists are substituted irrespective of order, so no RHS variables are substituted by other elements
    if isinstance(expr, (list, tuple)):
        res = [subs(eqs, ex) for ex in expr]
        return [x for x in res if x is not True]
    if not isinstance(expr, sympy.Expr):
        return expr
    if isinstance(eqs, dict):
        return expr.subs(eqs)
    else:
        if not isinstance(eqs, (list, tuple)):
            return expr.subs(eqs.lhs, eqs.rhs)
        eqs = [(eq.lhs, eq.rhs) for eq in eqs]
        # must use dict in order for rhs not to be substituted by subsequent lhs
        return expr.subs(dict(eqs))


def solve(eq, x):
    res = sympy.solve(eq, x)
    if len(res) == 1:
        return res[0]
    else:
        return res


half = sympy.Rational(1, 2)
third = sympy.Rational(1, 3)
Kf = sympy.IndexedBase('Kf', shape=(sympy.oo,))
K_fac = Symbol('Kf_1')  # currently sympy bugs just using Kf[1], so treat as separate symbol
K_fac_sub = 1 - 3 * K / k ** 2
K_sub = Eq(K, solve(K_fac_sub - K_fac, K))

H_t = diff(a, t) / a
dH = -a ** 2 / 6 * kappa * (rho + 3 * P)
Friedmann = Eq(H ** 2, a ** 2 * kappa * rho / 3 - K)
Friedmann_subs = [Friedmann, Eq(diff(H, t), dH), Eq(diff(a, t), a * H)]
Friedmann_Kfac_subs = subs(K_sub, Friedmann_subs)

delta_frame = Function('uprime')(t)


def LinearPerturbation(name, species=None, camb_var=None, camb_sub=None, frame_dependence=None, description=None):
    """
    Returns as linear perturbation variable, a function of conformal time t.
    Use help(x) to quickly view all quantities defined for the result.

    :param name: sympy name for the Function
    :param species: tag for the species if relevant (not used)
    :param camb_var: relevant CAMB fortran variable
    :param camb_sub:  if not equal to camb_var, and string giving the expression in CAMB variables
    :param frame_dependence: the change in the perturbation when the frame 4-velocity u change
             from u to u + delta_frame. Should be a sumpy expression involving delta_frame.
    :param description: string describing variable
    :return: sympy Function instance (function of t), with attributes set to the arguments above.
    """
    if isinstance(camb_var, list):
        camb_var = tuple(camb_var)
    f = Function(name, species=species, camb_var=camb_var, camb_sub=camb_sub, perturbation_order=1,
                 frame_dependence=frame_dependence, description=description)
    return f(t)


def list_perturbations(expr, lst=None):
    if lst is None:
        lst = []
    if getattr(expr, 'perturbation_order', None) and expr not in lst:
        lst.append(expr)
    for arg in expr.args:
        list_perturbations(arg, lst)
    return lst


def list_frame_dependent_vars(expr, lst=None):
    if lst is None:
        lst = []
    if getattr(expr, 'frame_dependence', None) and expr not in lst:
        lst.append(expr)
    for arg in expr.args:
        list_frame_dependent_vars(arg, lst)
    return lst


def simplify_sum(expr):
    if isinstance(expr, sympy.Add):
        return sympy.Add(*[simplify(term) for term in expr.args])
    else:
        return expr


def frame_change(expr, delta_u=None, total=False):
    if getattr(expr, 'frame_dependence', None):
        res = expr.frame_dependence
        if delta_u is not None:
            res = res.subs(delta_frame, delta_u)
        if total:
            res += expr
        return res
    else:
        if isinstance(expr, (list, tuple)):
            return [frame_change(x) for x in expr]
        perts = list_frame_dependent_vars(expr)
        res = subs([Eq(pert, pert + pert.frame_dependence) for pert in perts], expr)
        if delta_u is not None and isinstance(res, sympy.Expr):
            res = res.subs(delta_frame, delta_u)
        if not total:
            res -= expr
        res = simplify(res).expand()
        return simplify_sum(res.collect(list_frame_dependent_vars(res)))


# Perturbation variables

# gauge-invariant potential
phi = LinearPerturbation('phi', description='Weyl potential')

eta = LinearPerturbation('eta', camb_var='etak', camb_sub='-2*etak/k', description='three-curvature',
                         frame_dependence=-2 * K_fac * H * delta_frame / k)

sigma = LinearPerturbation('sigma', description='shear',
                           frame_dependence=delta_frame)

A = LinearPerturbation('A', camb_sub=0, description='acceleration',
                       frame_dependence=(diff(delta_frame, t) + H * delta_frame) / k)

z = LinearPerturbation('z', description='expansion rate perturbation',
                       #                        frame_dependence = delta_frame - 3*(dH- H**2)*delta_frame/k**2)
                       frame_dependence=K_fac * delta_frame + 3 * kappa * a ** 2 * (rho + P) * delta_frame / (
                               2 * k ** 2))

hdot = LinearPerturbation('hdot', camb_sub=k / 3 * z, description='time derivative of scale factor perturbation',
                          frame_dependence=k * delta_frame / 3 - diff(delta_frame * H, t) / k)

delta = LinearPerturbation('delta', camb_var='dgrho', camb_sub='dgrho/kappa/a**2',
                           description='total density perturbation',
                           frame_dependence=3 * H * (rho + P) * delta_frame / k)

delta_P = LinearPerturbation('delta_P', camb_sub='error', description='total pressure perturbation',
                             frame_dependence=-diff(P, t) * delta_frame / k)

q = LinearPerturbation('q', camb_var='dgq', camb_sub='dgq/kappa/a**2', description='total heat flux',
                       frame_dependence=-(rho + P) * delta_frame)

Pi = LinearPerturbation('Pi', camb_var='dgpi', camb_sub='dgpi/a**2/kappa', description='total anisotropic stress')

# quadrupole source
polter = Function('polter')(t)

# Newtonian gauge variables (in general equal to gauge invariant potentials)
Phi_N = LinearPerturbation('Phi_N', description='Newtonian gauge curvature potential')
Psi_N = LinearPerturbation('Psi_N', description='Newtonian gauge acceleration potential')

# Synchronous gauge variables
eta_s = LinearPerturbation('eta_s', description='Synchronous gauge curvature variable')
hdot_s = LinearPerturbation('hdot_s', description='Synchronous gauge diff(h,t)')

# velocities, and dimensionless heat flux [(rho_i+P_i)v_i = rho_i q_i]

q_r = LinearPerturbation('q_r', species='r', camb_var='qr', description='massless neutrino heat flux',
                         frame_dependence=-4 * delta_frame / 3)
q_g = LinearPerturbation('q_g', species='g', camb_var='qg', description='photon heat flux',
                         frame_dependence=-4 * delta_frame / 3)
q_nu = LinearPerturbation('q_nu', species='nu', camb_var='qnu', description='massive neutrino heat flux',
                          frame_dependence=-(1 + p_nu / rho_nu) * delta_frame)

v_c = LinearPerturbation('v_c', species='c', description='CDM velocity',
                         frame_dependence=-delta_frame)
v_b = LinearPerturbation('v_b', species='b', camb_var='vb', description='baryon velocity',
                         frame_dependence=-delta_frame)
v_de = LinearPerturbation('v_de', species='de', camb_var=['qde', 'w_lam'], camb_sub='qde/(1+w_lam)',
                          description='dark energy velocity', frame_dependence=-delta_frame)
Delta_b = LinearPerturbation('Delta_b', species='b', camb_var='clxb',
                             description='fractional baryon density perturbation',
                             frame_dependence=3 * H * (1 + p_b / rho_b) * delta_frame / k)
Delta_P_b = LinearPerturbation('Delta_P_b', camb_sub='delta_p_b',
                               description='fractional baryon pressure perturbation',
                               frame_dependence=-diff(p_b, t) / rho_b * delta_frame / k)

Delta_c = LinearPerturbation('Delta_c', species='c', camb_var='clxc', description='fractional CDM density perturbation',
                             frame_dependence=3 * H * delta_frame / k)
Delta_r = LinearPerturbation('Delta_r', species='r', camb_var='clxr',
                             description='fractional massless neutrino density perturbation',
                             frame_dependence=4 * H * delta_frame / k)

Delta_nu = LinearPerturbation('Delta_nu', species='nu', camb_var='clxnu',
                              description='fractional massive neutrino density perturbation',
                              frame_dependence=3 * H * (1 + p_nu / rho_nu) * delta_frame / k)
Delta_P_nu = LinearPerturbation('Delta_P_nu', camb_sub='error',
                                description='fractional massive neutrino pressure perturbation',
                                frame_dependence=-diff(p_nu, t) / rho_nu * delta_frame / k)

Delta_g = LinearPerturbation('Delta_g', species='g', camb_var='clxg', description='fractional CDM density perturbation',
                             frame_dependence=4 * H * delta_frame / k)
Delta_de = LinearPerturbation('Delta_de', species='de', camb_var='clxde',
                              description='fractional dark energy density perturbation',
                              frame_dependence=3 * H * (1 + w_de) * delta_frame / k)
# Sound speeds
csq_b = LinearPerturbation('c_sb^2', species='b', camb_var=['delta_p_b', 'clxb'], camb_sub='delta_p_b/clxb',
                           description='baryon sound speed')
csq_b.frame_dependence = (csq_b * Delta_b - diff(p_b, t) / rho_b * delta_frame / k) / (
        Delta_b + 3 * H * (1 + p_b / rho_b) * delta_frame / k) - csq_b

# dark energy sound speed defined in rest frame so gauge invariant
csqhat_de = LinearPerturbation('chat_sde^2', species='de', camb_var='cs2_lam',
                               description='rest frame dark energy sound speed')

# Anisotropic stress
pi_g = LinearPerturbation('pi_g', species='g', camb_var='pig', description='photon anisotropic stress')
pi_r = LinearPerturbation('pi_r', species='r', camb_var='pir', description='massless neutrino anisotropic stress')
pi_nu = LinearPerturbation('pi_nu', species='r', camb_var='pinu', description='massive neutrino anisotropic stress')

E_2 = LinearPerturbation('E_2', species='E', description='E polarization quadrupole')
E_3 = LinearPerturbation('E_3', species='E', description='E polarization octopole')
J_3 = LinearPerturbation('J_3', species='g', description='photon octopole')

polter_t = sympy.Rational(2, 15) * (3 * pi_g / 4 + 9 * E_2 / 2)
polter_sub = Eq(polter, polter_t)

background_eqs = [
    Eq(diff(a, t), H * a),
    Eq(diff(H, t), dH),
    Eq(diff(exptau, t), visibility)
]

# constraint equations relating variables
cons1 = k ** 3 * K_fac * phi + kappa / 2 * a ** 2 * k * (delta + Pi * K_fac) + 3 * kappa / 2 * a ** 2 * H * q
cons2 = k ** 2 * eta - kappa * a ** 2 * delta + 2 * k * H * z
cons3 = 2 * third * (k / a) ** 2 * (z - K_fac * sigma) + kappa * q
cons4 = hdot - k / 3 * z + H * A
constraints = [cons1, cons2, cons3, cons4]


def constraint_subs_for_variable_set(variables=(z, sigma, phi, hdot)):
    return solve(constraints, variables)


# Various substitutions for variables in terms of other variables (using constraints)

var_subs = constraint_subs_for_variable_set()

hdot_sub, phi_sub, sigma_sub, z_sub = [Eq(variable, subs(var_subs, variable)) for variable in [hdot, phi, sigma, z]]
q_sub = Eq(q, subs(Eq(z, subs(var_subs, z)), solve(cons3, q)))

# Evolution equations

dz = -H * z - kappa * a ** 2 / k / 2 * (delta + 3 * delta_P) + (3 * kappa * a ** 2 * (
        rho + P) / 2) * A / k + k * K_fac * A
dsigma = -H * sigma + k * (phi + A) - half * kappa * a ** 2 / k * Pi
deta = -1 / k * (2 * K * z + kappa * a ** 2 * q + 2 * K_fac * k * H * A)
dphi = -H * phi + half / k ** 2 * kappa * a ** 2 * (k * (rho + P) * sigma + k * q - diff(Pi, t) - H * Pi)

# Alternative equation
deta_2 = Eq(diff(eta, t), K_fac * (2 * hdot - 2 * k / 3 * sigma))

pert_eqs = [
    Eq(diff(z, t), dz),
    Eq(diff(sigma, t), dsigma),
    Eq(diff(eta, t), deta),
    Eq(diff(phi, t), dphi)]

drag_t = opacity * (4 * v_b / 3 - q_g)

total_eqs = [
    Eq(diff(q, t), k * delta_P - 4 * H * q - 2 * third * k * K_fac * Pi - (rho + P) * k * A),
    Eq(diff(delta, t), -k * q - 3 * H * (delta + delta_P) - 3 * hdot * (rho + P))
]

tot_eqs = total_eqs + pert_eqs + background_eqs

# frame names specify a frame where a certain variable is zero
frame_names = {'CDM': v_c, 'Newtonian': sigma, 'comoving': q,
               'flat': eta, 'constant density': delta}


def make_frame_invariant(expr, frame='CDM'):
    """
    Makes the quantity gauge invariant, assuming currently evaluated in frame 'frame'.
    frame can either be a string frame name, or a variable that is zero in the current frame,

    e.g. frame = Delta_g gives the constant photon density frame.
    So make_frame_invariant(sigma, frame=Delta_g) will return the combination of sigma and Delta_g
    that is frame invariant (and equal to just sigma when Delta_g=0).
    """

    if isinstance(frame, six.string_types):
        if frame not in frame_names:
            raise ValueError('Unknown frame names: %s' % frame)
        frame = frame_names[frame]
    if isinstance(expr, Eq):
        return simplify(
            Eq(make_frame_invariant(expr.lhs, frame), make_frame_invariant(expr.rhs, frame)))
    if isinstance(expr, (list, tuple)):
        return [make_frame_invariant(x, frame) for x in expr]
    # do frame change to make frame variable zero
    # special case of frame=A, equivalent to v_c frame by evolution equation
    if frame == A:
        frame = v_c
    delta_u = solve(frame_change(frame, total=True), delta_frame)
    if delta_frame in delta_u.atoms(Function):
        raise ValueError(
            'Cannot solve for change of frame. Currently only supports algebraic frame changes: %s' % delta_u)
    return frame_change(expr, delta_u=delta_u, total=True)


# Things for relating to Newtonian (zero shear) gauge

# Newtonian gauge variables in general. Note sigma here is the shear, Pi the anisotropic stress
Newt_vars = [Eq(Psi_N, -A + (diff(sigma, t) + H * sigma) / k),
             Eq(Phi_N, -eta / 2 / K_fac - H * sigma / k)]

Newtonian_var_subs = [Eq(Phi_N, phi + half * a ** 2 * kappa * Pi / k ** 2),
                      Eq(Psi_N, phi - half * a ** 2 * kappa * Pi / k ** 2)]

Newtonian_subs = [Eq(A, -Psi_N), Eq(diff(sigma, t, t), 0), Eq(diff(sigma, t), 0),
                  Eq(sigma, 0), Eq(phi, (Phi_N + Psi_N) / 2), Eq(eta, -2 * Phi_N * K_fac), Eq(hdot, -diff(Phi_N, t))]
Newtonian_subs += [Eq(z, subs(Newtonian_subs, solve(cons4, z)))]


def newtonian_gauge(x):
    r"""
    Evaluates an expression in the Newtonian gauge (zero shear, sigma=0).
    Converts to using conventional metric perturbation variables for metric

    .. math:: ds^2 = a^2\left( (1+2\Psi_N)d t^2 - (1-2\Phi_N)\delta_{ij}dx^idx^j\right)

    :param x: expression
    :return: expression evaluated in the Newtonian gauge
    """

    if isinstance(x, (list, tuple)):
        res = [newtonian_gauge(y) for y in x]
        return [x for x in res if x is not True]
    res = subs(Newtonian_subs, x)
    if isinstance(res, sympy.Expr):
        res = simplify(res.doit())
        res2 = simplify(res.subs(Psi_N, Phi_N - a ** 2 * kappa * Pi / k ** 2))
        if len(str(res2)) < len(str(res)):
            res = res2
    return res


# Things for relating to cdm frame or traditional synchronous gauge variables

# In the "cdm frame" we just fix to zero acceleration, comoving with CDM
# the covariant hdot and eta variables differ from Ma & Bertshinger by factors (the latter have _s after the name)

cdm_subs = [Eq(diff(A, t), 0), Eq(A, 0), Eq(diff(v_c, t), 0), Eq(v_c, 0)]


def cdm_gauge(x):
    r"""
    Evaluates an expression in the CDM frame :math:`(v_c=0, A=0)`. Equivalent to the synchronous gauge
    but using the covariant variable names.

    :param x: expression
    :return: expression evaluated in CDM frame.
    """

    if isinstance(x, (list, tuple)):
        return [cdm_gauge(y) for y in x]
    return simplify(subs(cdm_subs, x))


synchronous_vars = [Eq(eta_s, -make_frame_invariant(eta, v_c) / 2 / K_fac),
                    Eq(hdot_s, 6 * (hdot + H * A + v_c / k * (k ** 2 / 3 * K_fac + kappa * a ** 2 / 2 * (rho + P))))]

synchronous_subs = [Eq(eta, -2 * K_fac * eta_s), Eq(hdot, hdot_s / 6), Eq(phi, subs(var_subs, phi))]
synchronous_subs += [subs(synchronous_subs, Eq(z, solve(cons4, z))),
                     Eq(sigma, (hdot_s + diff(6 * eta_s, t)) / 2 / k)]

synchronous_subs = cdm_gauge(synchronous_subs) + cdm_subs


def synchronous_gauge(x):
    """
    evaluates an expression in the synchronous gauge, using conventional synchronous-gauge variables.

    :param x: expression
    :return: synchronous gauge variable expression
    """
    if isinstance(x, (list, tuple)):
        res = [synchronous_gauge(y) for y in x]
        return [x for x in res if x is not True]
    res = subs(synchronous_subs, x)
    if isinstance(res, sympy.Expr):
        return simplify(res.doit())
    return simplify(res)


# Fluid components
density_eqs = [
    Eq(diff(rho_b, t), -3 * H * (rho_b + p_b)),
    Eq(diff(rho_c, t), -3 * H * rho_c),
    Eq(diff(rho_g, t), -4 * H * rho_g),
    Eq(diff(rho_r, t), -4 * H * rho_r),
    Eq(diff(rho_nu, t), -3 * H * (rho_nu + p_nu)),
    Eq(diff(rho_de, t), -3 * H * (rho_de * (1 + w_de)))
]

delta_eqs = [
    Eq(diff(Delta_r, t), -4 * hdot - k * q_r),
    Eq(diff(Delta_g, t), -4 * hdot - k * q_g),
    Eq(diff(Delta_b, t), -(1 + p_b / rho_b) * (3 * hdot + k * v_b)
       + (p_b / rho_b - csq_b) * 3 * H * Delta_b),
    Eq(diff(Delta_c, t), -3 * hdot - k * v_c),
    Eq(diff(Delta_nu, t), -3 * (1 + p_nu / rho_nu) * hdot
       - k * q_nu + 3 * H * (-Delta_P_nu + Delta_nu * p_nu / rho_nu)),
    Eq(diff(Delta_de, t),
       -3 * (1 + w_de) * hdot - (1 + w_de) * k * v_de - 3 * H * (csqhat_de - w_de) * (
               Delta_de + 3 * H * (1 + w_de) * v_de / k)
       - 3 * H * diff(w_de, t) * v_de / k)
]

vel_eqs = [
    Eq(diff(q_r, t), -2 * third * k * pi_r * K_fac + third * k * Delta_r - 4 * k / 3 * A),
    Eq(diff(q_g, t), -2 * third * k * pi_g * K_fac + third * k * Delta_g + drag_t - 4 * k / 3 * A),
    Eq(diff(v_b, t), -k * A - H * v_b - diff(p_b, t) * v_b / (rho_b + p_b)
       + 1 / (rho_b + p_b) * (rho_b * csq_b * k * Delta_b - rho_g * drag_t)),
    Eq(diff(v_c, t), -H * v_c - k * A),
    Eq(diff(q_nu, t), -H * (1 - 3 * p_nu / rho_nu) * q_nu
       - k / 3 * (2 * K_fac * pi_nu - 3 * Delta_P_nu) - (1 + p_nu / rho_nu) * k * A),
    Eq(diff(v_de, t), k * csqhat_de * Delta_de / (1 + w_de) - H * (1 - 3 * csqhat_de) * v_de - k * A)
]

component_eqs = density_eqs + delta_eqs + vel_eqs

rho_t = rho_b + rho_c + rho_r + rho_g + rho_nu + rho_de
P_t = third * (rho_r + rho_g) + p_b + p_nu + w_de * rho_de
tot_subs = [
    Eq(rho, rho_t),
    Eq(P, P_t)
]

# Note that csqhat_de is defined in the dark energy rest-frame,
# so this is general-gauge result for pressure perturbation:
Delta_P_de = (
        csqhat_de * Delta_de + 3 * H * v_de / k * (1 + w_de) * (csqhat_de - w_de + diff(w_de, t) / 3 / H / (1 + w_de)))
tot_pert_subs = [
    Eq(Pi, rho_g * pi_g + rho_r * pi_r + rho_nu * pi_nu),
    Eq(delta, rho_g * Delta_g + rho_r * Delta_r + rho_b * Delta_b + rho_c * Delta_c
       + rho_nu * Delta_nu + rho_de * Delta_de),
    Eq(q, rho_g * q_g + rho_r * q_r + rho_c * v_c + (rho_b + p_b) * v_b + rho_nu * q_nu + rho_de * (1 + w_de) * v_de),
    Eq(delta_P, rho_nu * Delta_P_nu + third * (rho_g * Delta_g + rho_r * Delta_r) + csq_b * Delta_b * rho_b
       + Delta_P_de * rho_de)

]


def define_variable(name, namespace=None, order=1):
    namespace = namespace or globals()
    if name not in namespace:
        namespace[name] = sympy.Function(name, perturbation_order=order)(t)
    return namespace[name]


def define_variables(names, namespace=None, order=1):
    return [define_variable(name, namespace, order) for name in names.split()]


def _make_index_func(name, l, namespace=None):
    name += '_' + str(l)
    return define_variable(name, namespace)


# Boltzmann hierarchies. Haven't included massive neutrinos here as yet.
def J_eq(L):
    # photons
    assert (L > 1)
    Gl = _make_index_func('J', L)
    Glp = _make_index_func('J', L + 1)
    Glm = _make_index_func('J', L - 1)
    eq = -k / (2 * L + 1) * ((L + 1) * Kf[L] * Glp - L * Glm) - opacity * Gl
    if L == 2:
        eq = eq + 8 * k / 15 * sigma + opacity * polter
    return Eq(diff(Gl, t), eq).subs({sympy.sympify('J_2(t)'): pi_g, sympy.sympify('J_1(t)'): q_g})


def G_eq(L):
    # massless neutrinos
    assert (L > 1)
    Gl = _make_index_func('G', L)
    Glp = _make_index_func('G', L + 1)
    Glm = _make_index_func('G', L - 1)
    eq = -k / (2 * L + 1) * ((L + 1) * Kf[L] * Glp - L * Glm)
    if L == 2:
        eq = eq + 8 * k / 15 * sigma
    return Eq(diff(Gl, t), eq).subs({sympy.sympify('G_2(t)'): pi_r, sympy.sympify('G_1(t)'): q_r})


def E_eq(L):
    # E polarization
    assert (L > 1)
    El = _make_index_func('E', L)
    Elp = _make_index_func('E', L + 1)
    Elm = _make_index_func('E', L - 1)
    eq = -k / (2 * L + 1) * ((L + 3) * (L - 1) * Kf[L] * Elp / (L + 1) - L * Elm) - opacity * El
    if L == 2:
        eq = eq + polter * opacity
    return Eq(diff(El, t), eq).subs(sympy.sympify('E_1(t)'), 0)


def get_hierarchies(lmax=5):
    """
    Get Boltzmann hierarchies up to lmax for photons (J), E polarization and massless neutrinos (G).

    :param lmax: maxmimum multipole
    :return: list of equations
    """

    eqs = []
    for l in range(2, lmax):
        eqs += [J_eq(l), G_eq(l), E_eq(l)]
    return eqs


hierarchies = get_hierarchies()

scalar_E_source = visibility * 15 * polter / 8 / (f_K(tau0 - t) * k) ** 2


def get_scalar_temperature_sources(checks=False):
    """
    Derives terms in line of sight source, after integration by parts so that only integrated against
    a Bessel function (no derivatives).

    :param checks:  True to do consistency checks on result
    :return: monopole_source, ISW, doppler, quadrupole_source
    """

    # Line of sight temperature anisotropy source, by integrating by parts so only integrated against j_l
    source1 = k * sigma * exptau + visibility * 15 * polter / 8
    source = source1 / 3 + diff(source1, t, t) / k ** 2 + visibility * Delta_g / 4 \
             - hdot * exptau + diff(visibility * v_b - k * exptau * A, t) / k
    src = subs(var_subs, subs(background_eqs + pert_eqs,
                              subs(background_eqs + pert_eqs, source).expand().doit()).simplify()).simplify().expand()
    ISW = src.coeff(exptau) * exptau
    if checks:
        # Check ISW is just 2*d phi/d eta
        assert (subs(Friedmann_Kfac_subs,
                     subs(Friedmann_Kfac_subs,
                          subs(var_subs, 2 * dphi * exptau - ISW).simplify()).simplify()).simplify() == 0)
    doppler1 = subs(var_subs, subs(pert_eqs, diff(visibility * (v_b + sigma), t) / k)).simplify().collect(visibility)
    remainder = src - ISW - doppler1
    remainder = remainder.simplify().collect(visibility).subs(delta, solve(phi_sub, delta)).simplify()
    ISW = 2 * diff(phi, t) * exptau
    monopole_source = (Delta_g / 4 + (2 * phi + eta / (2 * K_fac))) * visibility
    quadrupole_source = (remainder - monopole_source).simplify()
    doppler = diff((v_b + sigma) * visibility, t) / k
    if checks:
        assert (subs(var_subs, subs(pert_eqs, (doppler - doppler1))).simplify() == 0)

    return monopole_source, ISW, doppler, quadrupole_source


# for translation in source output routine (does not include all equation variables)
_camb_cache = {}


def camb_fortran(expr, name='camb_function', frame='CDM', expand=False):
    """
    Convert symbolic expression to CAMB fortran code, using CAMB variable notation.
    This is not completely general, but it will handle conversion of Newtonian gauge
    variables like Psi_N, and most derivatives up to second order.

    :param expr: symbolic sympy expression using camb.symbolic variables and functions (plus any
                 standard general functions that CAMB can convert to fortran).
    :param name: lhs variable string to assign result to
    :param frame: frame in which to interret non gauge-invariant expressions.
                  By default uses CDM frame (synchronous gauge), as used natively by CAMB.
    :param expand: do a sympy expand before generating code
    :return: fortran code snippet
    """

    camb_diff_vars = 'etakdot qgdot qrdot vbdot pigdot pirdot pinudot ' + \
                     'octg octgdot polterdot polterddot diff_rhopi sigmadot phidot  ' + \
                     'ddvisibility dvisibility dopacity ddopacity'
    camb_arr_vars = 'Edot E'

    tau = _camb_cache.setdefault('tau', Symbol('tau'))

    etakdot, qgdot, qrdot, vbdot, pigdot, pirdot, pinudot, octg, octgdot, \
    polterdot, polterddot, diff_rhopi, sigmadot, phidot, \
    ddvisibility, dvisibility, dopacity, ddopacity = \
        define_variables(camb_diff_vars, _camb_cache)
    Edot, E = sympy.symbols(camb_arr_vars, cls=sympy.IndexedBase, shape=(sympy.oo,))

    # Keep everything except baryons pressure which is very small
    camb_diff_subs = [(diff(q_g, t), qgdot), (diff(pi_g, t), pigdot), (diff(pi_r, t), pirdot),
                      (diff(q_r, t), qrdot), (diff(v_b, t), vbdot),
                      (diff(visibility, t, t), ddvisibility), (diff(visibility, t), dvisibility),
                      (diff(opacity, t, t), ddopacity), (diff(opacity, t), dopacity),
                      (diff(J_3, t), octgdot), (diff(E_2, t), Edot[2]),
                      (diff(E_3, t), Edot[3]), (diff(pi_nu, t), pinudot),
                      (diff(eta, t), -2 * etakdot / k), (diff(Pi, t), diff_rhopi / kappa / a ** 2),
                      (diff(polter, t, t), polterddot), (diff(polter, t), polterdot),
                      (diff(sigma, t), sigmadot), (diff(phi, t), phidot), (diff(p_b, t), 0)]

    if frame != 'CDM':
        expr = make_frame_invariant(expr, frame)

    # substitute for variables not available in CAMB function
    expr = cdm_gauge(subs(Newt_vars + synchronous_vars + [hdot_sub, z_sub], expr)).doit().simplify()
    res = cdm_gauge(subs(background_eqs + total_eqs, expr).subs(camb_diff_subs))

    if 'Derivative' in str(res):
        raise Exception(
            'Unknown derivatives, generally can only handle up to second.\nRemaining derivatives: ' + str(res))

    res = cdm_gauge(subs([K_sub, hdot_sub, z_sub], res))
    camb_var_subs = []
    for var in res.atoms(Function):
        camb_var = getattr(var, 'camb_var', None)
        if camb_var:
            if isinstance(camb_var, (list, tuple)):
                for x in camb_var:
                    define_variable(x)
                camb_sub = getattr(var, 'camb_sub', None)
                if not camb_sub:
                    raise Exception('must have camb_sub if camb_var has more than one variable')
            else:
                camb_var = define_variable(camb_var)
                camb_sub = getattr(var, 'camb_sub', None) or camb_var
            if camb_sub:
                if isinstance(camb_sub, six.string_types):
                    camb_sub = eval(camb_sub)
                camb_var_subs.append((var, camb_sub))

    camb_subs = camb_var_subs + [(p_b, 0), (E_2, E[2]), (E_3, E[3]), (J_3, octg), (K_fac, Kf[1])]
    res = res.subs(camb_subs).simplify()
    no_arg_funcs = [f for f in res.atoms(Function) if f.args[0] == t and f is not f_K]
    res = res.subs(zip(no_arg_funcs, [Symbol(str(x.func)) for x in no_arg_funcs]))
    res = res.subs(t, tau)
    if expand:
        res = res.expand()
    res = res.collect([Symbol(str(x.func)) for x in
                       [k, sigma, opacity, visibility, dopacity, dvisibility, ddvisibility]])
    res = sympy.fcode(res, source_format='free', standard=95, assign_to=name, contract=False)
    import textwrap

    if 'if ' not in res:
        lines = res.split('\n')
        for i, line in enumerate(lines):
            if '=' in line:
                res = '\n'.join(lines[i:])
                break
        res = ''.join([x.strip() for x in res.split('&')])
        res = ' &\n    '.join(textwrap.wrap(res))
    return res


_func_cache = {}
_source_file_count = 0

_default_compiler = None
_default_flags = None


def get_default_compiler():
    global _default_compiler, _default_flags
    if _default_compiler:
        return _default_compiler
    from .baseconfig import gfortran
    if gfortran:
        _default_compiler = 'gfortran'
        _default_flags = "-shared -fPIC -O1 -fmax-errors=4"
    else:
        import platform
        _default_compiler = 'ifort'
        if platform.system() == 'Darwin':
            _default_flags = "-dynamiclib -O1 -W0 -WB"
        else:
            _default_flags = "-shared -fpic -O1 -W0 -WB"
    # _default_flags="-shared -fPIC -g -fbounds-check -fbacktrace -ffpe-trap=invalid,overflow,zero",
    return _default_compiler


_first_compile = True


def compile_source_function_code(code_body, file_path='', compiler=None, fflags=None, cache=True):
    """
    Compile fortran code into function pointer in compiled shared library.
    The function is not intended to be called from python, but for passing back to compiled CAMB.

    :param code_body: fortran code to do calculation and assign sources(i) output array.
     Can start with declarations of temporary variables if needed.
    :param file_path: optional output path for generated f90 code
    :param compiler: compiler, usually on path
    :param fflags: options for compiler
    :param cache: whether to cache the result
    :return: function pointer for compiled code
    """

    if cache and code_body in _func_cache:
        return _func_cache[code_body].source_func_

    global _source_file_count

    template = """
    REAL*8 function source_func(sources, tau, a, adotoa, grho, gpres,w_lam, cs2_lam,  &
        grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grhonu_t, &
        k,etak, etakdot, phi, phidot, sigma, sigmadot, &
        dgrho, clxg,clxb,clxc,clxr,clxnu, clxde, delta_p_b, &
        dgq, qg, qr, qde, vb, qgdot, qrdot, vbdot, &
        dgpi, pig, pir, pigdot, pirdot, diff_rhopi, &
        polter, polterdot, polterddot, octg, octgdot, E, Edot, &
        opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau, &
        tau0, tau_maxvis, Kf, f_K)
    implicit none
    real*8, intent(out) :: sources(:)
    REAL*8, intent(in) ::  tau, a, adotoa, grho, gpres,w_lam, cs2_lam,  &
            grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grhonu_t, &
            k,  etak, etakdot, phi, phidot, sigma, sigmadot, &
            dgrho, clxg,clxb,clxc,clxr,clxnu, clxde, delta_p_b, &
            dgq, qg, qr, qde, vb, qgdot, qrdot, vbdot, &
            dgpi, pig, pir,  pigdot, pirdot, diff_rhopi, &
            polter, polterdot, polterddot, octg, octgdot, E(2:3), Edot(2:3), &
            opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau, &
            tau0, tau_maxvis, Kf(*)
    real*8, external :: f_K

    %s
    end function
    """

    import subprocess
    import tempfile
    from ._compilers import is_32_bit, is_windows, compiler_environ, check_gfortran

    compiler = compiler or get_default_compiler()
    fflags = fflags or _default_flags

    if is_32_bit:
        fflags = "-m32 " + fflags
    if is_windows:
        global _first_compile
        fflags += ' -static'
        if _first_compile:
            check_gfortran(msg=True)
    workdir = file_path or tempfile.gettempdir()
    if not os.access(workdir, os.F_OK):
        os.mkdir(workdir)

    oldwork = os.getcwd()
    source_file = None
    try:
        os.chdir(workdir)
        _source_file_count += 1
        while True:
            name_tag = 'camb_source%s' % _source_file_count
            dll_name = name_tag + '.dll'
            if not os.path.exists(dll_name):
                break
            try:
                os.remove(dll_name)
            except:
                _source_file_count += 1

        source_file = name_tag + '.f90'
        with open(source_file, 'w') as f:
            f.write(template % code_body)

        command = " ".join([compiler, fflags, source_file, "-o", dll_name])
        try:
            subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True, cwd=workdir, env=compiler_environ)
        except subprocess.CalledProcessError as E:
            print(command)
            print('Error compiling generated code:')
            print(E.output)
            print('Source is:\n %s' % code_body)
            raise
    finally:
        if not file_path and source_file:
            os.remove(source_file)
        os.chdir(oldwork)

    # Had weird crashes when LoadLibrary path was relative to current dir
    dll_name = os.path.join(workdir, dll_name)
    func_lib = ctypes.LibraryLoader(ctypes.CDLL).LoadLibrary(dll_name)

    if cache:
        _func_cache[code_body] = func_lib

    if not file_path:
        # won't work on Windows while DLL in use
        try:
            os.remove(dll_name)
        except:
            pass

    return func_lib.source_func_


def compile_sympy_to_camb_source_func(sources, code_path=None, frame='CDM'):
    code = ''
    if not isinstance(sources, (list, tuple)):
        sources = [sources]
    for i, source in enumerate(sources):
        code += camb_fortran(source, 'sources(%s)' % (i + 1), frame=frame) + '\n'
    return compile_source_function_code(code, file_path=code_path)


def internal_consistency_checks():
    print('Sympy: ', sympy.__version__)

    # All equations should be gauge invariant
    for cons in constraints:
        assert (simplify(subs(Friedmann_Kfac_subs, make_frame_invariant(cons)) - cons) == 0)

    # Check deta equations consistent
    assert (subs(K_sub, subs(var_subs, subs(q_sub, deta).simplify() -
                             subs(deta_2, diff(eta, t)))).simplify() == 0)

    # check consistency of constraint equations with evolution equations
    assert (subs(var_subs, subs(tot_eqs, diff(cons3, t))).simplify() == 0)
    assert (subs(Friedmann, subs(var_subs, subs(tot_eqs, diff(cons2, t))).simplify()).simplify() == 0)
    assert (subs(K_sub, subs(Friedmann, subs(var_subs, subs(tot_eqs, diff(cons1, t))).simplify())).simplify() == 0)

    # Check Weyl potential from Newtonian gauge vars
    assert (subs(var_subs, subs(pert_eqs, subs(Newt_vars, Phi_N + Psi_N - 2 * phi))).simplify() == 0)

    # Gauge functions checks
    assert (newtonian_gauge(dsigma) == 0)
    assert (subs(K_sub, subs(var_subs, subs(pert_eqs, subs(cdm_gauge(synchronous_vars),
                                                           subs(synchronous_subs, sigma) - sigma).doit()))
                 .simplify().collect(eta)).simplify() == 0)

    # check all equations are gauge invariant
    for eq in delta_eqs + vel_eqs:
        assert (simplify(subs(component_eqs + background_eqs, make_frame_invariant(eq.lhs).doit()
                              - make_frame_invariant(eq.rhs))) == 0)

    try:
        # Check consistency of fluid equations with equations from total stress-energy conservation
        for eq in total_eqs:
            assert (subs(component_eqs + tot_subs, subs(tot_pert_subs, eq).doit()).doit().simplify())
    except TypeError:
        print('If this test fails, you probably have an old sympy version. You need version 1 or higher.')
        raise

    get_scalar_temperature_sources(True)
    print("All symbolic relation tests OK")


if __name__ == "__main__":
    internal_consistency_checks()
