.. _code-modifications:

Modifying the code
==================

Although CAMB supports some non-standard models by default (e.g. some early dark energy models), when you have a new
model you'll generally need to modify the code. Simple cases that do not need code modification are:

- Dark energy fluid models with a given equation of state but constant sound speed (see :doc:`dark_energy`)
- Different primordial power spectra (see :doc:`initialpower`)
- Different BBN mappings for the Helium abundance (which is pure Python, see :doc:`bbn`)

In these cases, you can just pass in an interpolation table from Python to encapsulate the modified physics.

Defining new classes
--------------------


For other changes to the dark energy, initial power, reionization, recombination, or non-linear correction, you can usually
define new classes that inherit from the standard base classes. The classes are defined in both Python and Fortran, so you
will need to modify both. Ensure variables are defined in the same order so that the Python interface works consistently.
The easiest way to do this is probably to look at the source code for, e.g., `AxionEffectiveFluid` in both Fortran and Python and follow the same pattern.

In Python, CAMB uses the `@fortran_class` decorator to implement the Fortran wrapping. This implements a general but custom mapping between
F2008 Fortran classes (types) and Python classes. For example, the `AxionEffectiveFluid` Python class looks like this::

    @fortran_class
    class AxionEffectiveFluid(DarkEnergyModel):

        _fields_ = [
            ("w_n", c_double, "effective equation of state parameter"),
            ("fde_zc", c_double, "energy density fraction at z=zc"),
            ("zc", c_double, "decay transition redshift (not the same as peak of energy density fraction)"),
            ("theta_i", c_double, "initial condition field value")
        ]

        _fortran_class_name_ = 'TAxionEffectiveFluid'
        _fortran_class_module_ = 'DarkEnergyFluid'

The `_fields_` and `_fortran_class_name_` class variables are special (metaclass) metadata for the Fortran mapping, so that it knows which class to map to and
what the corresponding variable names are. It is not necessary to include all Fortran variables in the `_fields_` array, but they must be in
the right order, and only missing items at the end.

In Fortran, the corresponding class is defined as::

    type, extends(TDarkEnergyModel) :: TAxionEffectiveFluid
        real(dl) :: w_n = 1._dl ! Effective equation of state when oscillating
        real(dl) :: fde_zc = 0._dl ! Energy density fraction at a_c (not the same as peak dark energy fraction)
        real(dl) :: zc  ! Transition redshift (scale factor a_c)
        real(dl) :: theta_i = const_pi/2 ! Initial value
        ! om is Omega of the early DE component today (assumed to be negligible compared to omega_lambda)
        ! omL is the lambda component of the total dark energy omega
        real(dl), private :: a_c, pow, om, omL, acpow, freq, n ! Cached internally
    contains
    procedure :: ReadParams =>  TAxionEffectiveFluid_ReadParams
    procedure, nopass :: PythonClass => TAxionEffectiveFluid_PythonClass
    procedure, nopass :: SelfPointer => TAxionEffectiveFluid_SelfPointer
    procedure :: Init => TAxionEffectiveFluid_Init
    procedure :: w_de => TAxionEffectiveFluid_w_de
    procedure :: grho_de => TAxionEffectiveFluid_grho_de
    procedure :: PerturbedStressEnergy => TAxionEffectiveFluid_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TAxionEffectiveFluid_PerturbationEvolve
    end type TAxionEffectiveFluid

Here, the `a_c`, `pow`, etc., variables are not mapped to Python because they are only used in the Fortran code.
The rest of the Fortran code defines the relevant methods (`TAxionEffectiveFluid_grho_de`, etc.) to implement the modified physics.

All Fortran classes that map to Python must have a `SelfPointer` function as defined above. This always takes exactly the same form::

    subroutine TMyClass_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TMyClass), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TMyClass_SelfPointer

and hence can be trivially modified for any new classes. The function is essential so that the Python wrapper can
map Python and strongly typed Fortran classes consistently (not something that is supported by standard `iso_c_binding`).
The `ReadParams` function is only needed if you want to also be able to load parameters from `.ini` files, rather than just
using them via Python. The `PythonClass` method is not strictly needed.

Other code changes
------------------

For quintessence models with a different potential, you will need to modify the Fortran to define the new potential function.
See :class:`~camb.dark_energy.Quintessence` and the `DarkEnergyQuintessence.f90` Fortran source. This may be a simple change,
though you may also need more complicated changes to consistently map input parameters into initial conditions for the evolution.

More generally, you will need to modify the equations at both the background and the perturbation level, usually in `equations.f90`.
The `CAMB notes <https://cosmologist.info/notes/CAMB.pdf>`_ provide some guidance on conventions and variable definitions.

Code updates, testing, and gotchas
----------------------------------


Make sure you recompile the Fortran after making any changes (see :doc:`fortran_compilers`).
Changing the version number in both Python and Fortran will give you an automatic run-time check that the Python being run matches the
intended Fortran source.

The default accuracy parameters are designed for Simons Observatory-like precision for standard models. Check your results are stable to
increasing accuracy parameters `AccuracyBoost` and `lAccuracyBoost` (in :class:`~camb.model.AccuracyParams`). If not, changing specific accuracy parameters as needed may be much more efficient that using the high-level parameter
`AccuracyBoost` (which increases the accuracy of many things at once).

There are a number of possible gotchas when using Python-wrapped Fortran types. Firstly, types derived directly from `CAMB_Structure` are intended to map directly
to Fortran types (via the standard `ctypes` interface), for example, `AccuracyParams` is inherited directly from `CAMB_Structure`. These should generally not be
instantiated directly in Python as they are only intended to be used as sub-components of larger types. For example, a new Python instance of :class:`~camb.model.AccuracyParams` will
give a zero Fortran array, which does not correspond to the default values for the accuracy parameters.

Fortran-mapped classes in Python inherit from `F2003Class`. These also map data in a Fortran class type (the `_fields_` defined above).
If they are an allocatable subcomponent of another `F2003Class`, they may be created dynamically to match the underlying structure.
This can give unexpected results if you try to add variables to only the Python class. For example, if `pars` is a :class:`~camb.model.CAMBparams` instance and `test` is not defined
then doing this::

    pars.DarkEnergy.test = 'x'
    print(pars.DarkEnergy.test)

will not give you 'x'; it will give you an undefined variable error. This is because the Python code doesn't 'know' that the Fortran code is not modifying the
DarkEnergy structure, so `pars.DarkEnergy` is generating a new instance mapped to the underlying Fortran data whenever you access it.
You can avoid this by always defining fields in both Fortran and Python, or only using Python variables in container-level classes like :class:`~camb.model.CAMBparams`.

When using dark energy models, make sure you are not setting `thetastar` in Python before setting the dark energy parameters: it needs to know the dark
energy model to map `thetastar` into `H0` consistently.

When accessing array-like members of a structure, e.g., `CAMBparams.z_outputs`, you may need to explicitly cast to a list to see the elements.

Interfacing with Cobaya
-----------------------

The `Cobaya sampler <https://cobaya.readthedocs.org>`_ can do parameter inference for your custom models. It uses introspection to determine which
variables the linked CAMB version supports, so if you add new variables e.g., to :class:`~camb.model.CAMBparams` or as arguments to :meth:`~camb.model.CAMBparams.set_cosmology` or the `set_params`
method of the dark energy, reionization, etc. classes, you should automatically be able to use them in Cobaya. For other new variables, you may need to modify :func:`~camb.get_valid_numerical_params`.

For supporting new primordial power spectra or multiple bins there are `test examples <https://github.com/CobayaSampler/cobaya/blob/master/tests/test_cosmo_multi_theory.py>`_.
This also shows how to use `get_class_options` to dynamically define multiple parameters based on an input parameter.

You can only directly sample scalar parameters, but it is also easy to `map vector parameters <https://cobaya.readthedocs.io/en/latest/params_prior.html#vector-parameters>`_.
Cobaya will automatically identify numerical arguments to the `set_params`
function of custom classes (e.g. dark energy), but for vector parameters to be picked up for sampling you need define
them with a default value of `None`.

The `CosmoCoffee <https://cosmocoffee.info/viewforum.php?f=11>`_ discussion forum can be used to ask questions and to see previous answers.