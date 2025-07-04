CAMB Variables and Gauge Conventions
====================================

This page provides a comprehensive guide to the variable naming conventions used in CAMB and their relationship to different gauge choices, particularly the synchronous gauge. For detailed mathematical derivations and symbolic manipulation, see the `ScalEqs notebook <https://camb.readthedocs.io/en/latest/ScalEqs.html>`_, the :doc:`symbolic` module documentation, and the `technical notes <https://cosmologist.info/notes/CAMB.pdf>`_.

**Units Convention:** CAMB uses natural units where :math:`c = 1` and distances are measured in Mpc. The gravitational coupling is :math:`\kappa = 8\pi G` with units of :math:`\text{Mpc}^{-2}`.

Overview
--------

CAMB uses a covariant perturbation formalism that can be expressed in different gauge choices. The code internally works in the CDM frame (equivalent to synchronous gauge) but provides functions to convert to other gauges like the Newtonian gauge. The variable naming follows specific conventions that encode both the physical quantity and the species.



Key Physical Quantities
-----------------------

The fundamental perturbation variables in CAMB represent different aspects of the perturbed spacetime and matter fields:

**Metric Perturbations:**

* :math:`\phi` **(phi)**: Weyl potential - gauge-invariant gravitational potential
* :math:`\eta` **(eta)**: Three-curvature perturbation (CAMB variable: ``etak``)
* :math:`\sigma` **(sigma)**: Shear perturbation
* :math:`z`: Expansion rate perturbation
* :math:`\dot{h}` **(hdot)**: Time derivative of scale factor perturbation
* :math:`A`: Acceleration perturbation

**Matter Perturbations:**

* :math:`\delta` **(delta)**: Density perturbation
* :math:`q`: Heat flux (momentum density)
* :math:`\Pi` **(Pi)**: Anisotropic stress

Variable Naming Conventions
---------------------------

CAMB uses systematic naming conventions for variables:

**Species Indices:**
  * ``g`` - photons
  * ``r`` - massless neutrinos
  * ``c`` - CDM (Cold Dark Matter)
  * ``b`` - baryons
  * ``nu`` - massive neutrinos
  * ``de`` - dark energy

**Variable Prefixes:**
  * ``clx`` - fractional density perturbations (:math:`\Delta_i = \delta\rho_i/\rho_i`)
  * ``q`` - heat flux variables
  * ``v`` - velocity perturbations
  * ``pi`` - anisotropic stress components
  * ``grho`` - background densities (with ``g`` prefix for ":math:`\kappa` times :math:`\rho`")
  * ``dg`` - total (summed over species) perturbation quantities

CAMB Fortran Variables
----------------------

This shows the correspondence between symbolic variables (in CDM frame) and CAMB Fortran variables:

.. list-table:: **Density Perturbations**
   :widths: 20 20 60
   :header-rows: 1

   * - Symbolic
     - CAMB Variable
     - Description
   * - :math:`\Delta_c`
     - ``clxc``
     - CDM fractional density perturbation
   * - :math:`\Delta_b`
     - ``clxb``
     - Baryon fractional density perturbation
   * - :math:`\Delta_g`
     - ``clxg``
     - Photon fractional density perturbation
   * - :math:`\Delta_r`
     - ``clxr``
     - Massless neutrino fractional density perturbation
   * - :math:`\Delta_{\nu}`
     - ``clxnu``
     - Massive neutrino fractional density perturbation
   * - :math:`\Delta_{de}`
     - ``clxde``
     - Dark energy fractional density perturbation

.. list-table:: **Velocity and Heat Flux**
   :widths: 20 20 60
   :header-rows: 1

   * - Symbolic
     - CAMB Variable
     - Description
   * - :math:`v_b`
     - ``vb``
     - Baryon velocity
   * - :math:`q_g`
     - ``qg``
     - Photon heat flux
   * - :math:`q_r`
     - ``qr``
     - Massless neutrino heat flux
   * - :math:`q_{\nu}`
     - ``qnu``
     - Massive neutrino heat flux

.. list-table:: **Anisotropic Stress**
   :widths: 20 20 60
   :header-rows: 1

   * - Symbolic
     - CAMB Variable
     - Description
   * - :math:`\pi_g`
     - ``pig``
     - Photon anisotropic stress
   * - :math:`\pi_r`
     - ``pir``
     - Massless neutrino anisotropic stress
   * - :math:`\pi_{\nu}`
     - ``pinu``
     - Massive neutrino anisotropic stress

.. list-table:: **Metric and Total Quantities**
   :widths: 20 20 60
   :header-rows: 1

   * - Symbolic
     - CAMB Variable
     - Description
   * - :math:`\eta`
     - ``etak``
     - Three-curvature (:math:`\mathrm{etak} = k\eta_s = -k\eta/2` in CDM frame)
   * - :math:`\delta` (total)
     - ``dgrho``
     - Total density perturbation (:math:`\kappa a^2 \sum_i \rho_i\Delta_i`)
   * - :math:`q` (total)
     - ``dgq``
     - Total heat flux (:math:`\kappa a^2 \sum_i \rho_i q_i`)
   * - :math:`\Pi` (total)
     - ``dgpi``
     - Total anisotropic stress



Physical Interpretation
-----------------------

**Heat Flux Relations:**
For each species i, the heat flux :math:`q_i` is related to the velocity :math:`v_i` by:

.. math::
   \rho_i q_i = (\rho_i + p_i)v_i

This gives:

* For relativistic species (photons, massless neutrinos): :math:`q_i = \frac{4}{3}v_i`
* For non-relativistic matter: :math:`q_i \approx v_i`

**Density Perturbations:**
The fractional density perturbations :math:`\Delta_i = \delta\rho_i/\rho_i` are gauge-dependent. In synchronous gauge, they represent the fractional density contrast in the frame comoving with the CDM. Note the distinction: :math:`\delta\rho_i` is the absolute density perturbation, while :math:`\Delta_i` is the fractional (relative) density perturbation.

**Anisotropic Stress:**
The anisotropic stress :math:`\pi_i` represents the traceless part of the stress tensor and is gauge-invariant. It is important for:

* Photon polarization (:math:`\pi_g`)
* Free-streaming neutrinos (:math:`\pi_r`, :math:`\pi_{\nu}`)
* Gravitational wave generation

Background Variables
--------------------

CAMB also defines background (unperturbed) quantities with specific naming:

.. list-table:: **Background Densities and Pressures**
   :widths: 20 20 60
   :header-rows: 1

   * - Symbolic
     - CAMB Variable
     - Description
   * - :math:`\rho_b`
     - ``grhob_t``
     - Baryon background density (:math:`\kappa\rho_b a^2`)
   * - :math:`\rho_c`
     - ``grhoc_t``
     - CDM background density (:math:`\kappa\rho_c a^2`)
   * - :math:`\rho_g`
     - ``grhog_t``
     - Photon background density (:math:`\kappa\rho_g a^2`)
   * - :math:`\rho_r`
     - ``grhor_t``
     - Massless neutrino background density (:math:`\kappa\rho_r a^2`)
   * - :math:`\rho_{\nu}`
     - ``grhonu_t``
     - Massive neutrino background density (:math:`\kappa\rho_{\nu} a^2`)
   * - :math:`\rho_{de}`
     - ``grhov_t``
     - Dark energy background density (:math:`\kappa\rho_{de} a^2`)
   * - :math:`H`
     - ``adotoa``
     - Hubble parameter (conformal time)

**Note:** The ``g`` prefix in CAMB variables stands for ":math:`\kappa` times" where :math:`\kappa = 8\pi G`, and densities are stored as :math:`\kappa\rho a^2` with units of :math:`\text{Mpc}^{-2}` for numerical convenience.

Synchronous Gauge Details
-------------------------

**CDM Frame**

CAMB natively works in the CDM frame where:

* CDM velocity: :math:`v_c = 0`
* Acceleration: :math:`A = 0`

This is equivalent to the synchronous gauge with the gauge choice that the CDM is at rest.

**Synchronous Gauge Metric**

In synchronous gauge, the metric takes the form:

.. math::
   ds^2 = a^2(\tau)[d\tau^2 - (\delta_{ij} + h_{ij})dx^i dx^j]

where the metric perturbation :math:`h_{ij}` can be decomposed into scalar, vector, and tensor parts.

**CAMB's Synchronous Gauge Variables:**

* :math:`\eta_s`: Related to the trace of :math:`h_{ij}` (synchronous gauge curvature)
* :math:`\dot{h}_s`: Time derivative of the trace
* **etak**: :math:`k\eta_s` (the variable actually stored in CAMB)

**Conversion Relations:**

From CAMB's covariant variables to synchronous gauge:

.. math::
   \eta_s = -\frac{\eta}{2K_{\mathrm{fac}}} \quad \text{where} \quad K_{\mathrm{fac}} = 1 - \frac{3K}{k^2}

.. math::
   \mathrm{etak} = k\eta_s = -\frac{k\eta}{2K_{\mathrm{fac}}}

In the flat case (:math:`K = 0`), this simplifies to:

.. math::
   \mathrm{etak} = -\frac{k\eta}{2}

In the CDM frame (:math:`v_c = 0`, :math:`A = 0`), the relationship of CAMB's `hdot` to the synchronous
gauge variable is given by:

.. math::
   \dot{h}_s = 6\dot{h} = 2kz

where :math:`z` is the expansion rate perturbation.

Gauge Transformation Examples
-----------------------------

**Example 1: CDM Frame to Newtonian Gauge**

To transform from CDM frame to Newtonian gauge, apply:

* Set :math:`\sigma = 0` (zero shear condition)
* :math:`\Phi_N = \phi + \frac{1}{2}\frac{a^2\kappa\Pi}{k^2}`
* :math:`\Psi_N = \phi - \frac{1}{2}\frac{a^2\kappa\Pi}{k^2}`

**Example 2: Frame-Dependent Variables**

Some variables change under gauge transformations:

* Density perturbations: :math:`\Delta_i \to \Delta_i + \frac{3H(1+w_i)\delta u}{k}`
* Velocities: :math:`v_i \to v_i - \delta u`
* Heat flux: :math:`q_i \to q_i - \frac{(\rho_i+p_i)\delta u}{\rho_i}`

where :math:`\delta u` is the frame transformation parameter.

Practical Usage Notes
---------------------

**For Transfer Functions:**

* All density variables :math:`\Delta_i` are in synchronous gauge
* Velocities depend on the specific context and gauge choice
* The Weyl potential is gauge-invariant

**For Custom Sources:**

* Use :func:`camb.symbolic.make_frame_invariant` to create gauge-invariant combinations
* The symbolic module provides automatic conversion between gauges
* See the ScalEqs notebook for practical examples

**Common Pitfalls:**

* Don't mix variables from different gauges without proper transformation
* Remember that CAMB's "synchronous gauge" is specifically the CDM frame
* Anisotropic stress components (:math:`\pi_g`, :math:`\pi_r`, :math:`\pi_{\nu}`) and total anisotropic stress :math:`\Pi` are gauge-invariant

Cross-References
----------------

* :doc:`symbolic` - Complete symbolic equation system and gauge transformations
* `ScalEqs notebook <https://camb.readthedocs.io/en/latest/ScalEqs.html>`_ - Interactive examples with variable definitions
* :doc:`transfer_variables` - Transfer function variables and their meanings
* :doc:`model` - Parameter and variable definitions for the Python interface

For the complete mathematical framework and equation derivations, see the technical notes referenced in the main documentation and the symbolic module documentation.
