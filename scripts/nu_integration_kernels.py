"""Reproduce the massive-neutrino sampling rules used by CAMB.

Reproduces the massive-neutrino momentum sampling rules hard-coded in
fortran/massive_neutrinos.f90, following the Appendix A moment-matching
construction in arXiv:1201.3654. This script is a maintained, self-contained
replacement for the old broken Mathematica notebook reference
``NeutrinoIntegrationKernels.nb``; run with --help for the full case list,
current findings, and usage notes.
"""

import argparse
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.optimize import differential_evolution, least_squares, minimize
from scipy.special import gamma, zeta

FERMI_DIRAC_CONST = 7.0 * np.pi**4 / 120.0
REPRESENTATIVE_AMS = (0.3, 1.0, 3.0, 10.0)
BENCHMARK_AM = np.geomspace(0.03, 12.0, 15)


@dataclass(frozen=True)
class QuadratureRule:
    name: str
    q: np.ndarray
    weights: np.ndarray
    condition_labels: tuple[str, ...]
    residuals: np.ndarray


@dataclass(frozen=True)
class MatchCondition:
    label: str
    func: Callable[[np.ndarray], np.ndarray]
    target: float
    scale: float


@dataclass(frozen=True)
class BenchmarkSummary:
    max_abs_relative_error: float
    rms_relative_error: float
    worst_case_label: str


@dataclass(frozen=True)
class CompiledCambCase:
    name: str
    neutrino_q_boost: float
    expected_nqmax: int
    rule: QuadratureRule | None = None


@dataclass(frozen=True)
class RelativeErrorSummary:
    max_abs_relative_error: float
    rms_relative_error: float


@dataclass(frozen=True)
class KernelFunction:
    label: str
    formula: str
    func: Callable[[np.ndarray], np.ndarray]
    target: float


@dataclass(frozen=True)
class CompiledCambModel:
    mnu: float
    num_massive_neutrinos: int


def sigmoid(value: float) -> float:
    return 1.0 / (1.0 + np.exp(-value))


def fermi_dirac_integral(power: int) -> float:
    if power == 0:
        return float(np.log(2.0))
    return float((1.0 - 2.0 ** (-power)) * gamma(power + 1.0) * zeta(power + 1.0))


def normalized_moment(exponent: int) -> float:
    if exponent == -4:
        return 1.0 / (8.0 * FERMI_DIRAC_CONST)
    if exponent < -4:
        raise ValueError(f"Unsupported exponent {exponent}")
    return ((exponent + 4.0) / 4.0) * fermi_dirac_integral(exponent + 3) / FERMI_DIRAC_CONST


def normalized_kernel_density(q: float) -> float:
    exp_minus_q = np.exp(-q)
    return 0.25 * q**4 * exp_minus_q / (1.0 + exp_minus_q) ** 2 / FERMI_DIRAC_CONST


def exact_integral(func: Callable[[np.ndarray], np.ndarray]) -> float:
    value, _ = quad(lambda q: normalized_kernel_density(q) * float(func(np.array(q))), 0.0, np.inf, limit=200)
    return float(value)


def velocity(q: np.ndarray, am: float) -> np.ndarray:
    return 1.0 / np.sqrt(1.0 + (am / q) ** 2)


def make_moment_condition(exponent: int) -> MatchCondition:
    target = normalized_moment(exponent)
    return MatchCondition(
        label=f"n={exponent}",
        func=lambda q, exponent=exponent: q**exponent,
        target=target,
        scale=max(1.0, abs(target)),
    )


def make_velocity_function(am: float, power: int) -> Callable[[np.ndarray], np.ndarray]:
    if power == 1:

        def func(q: np.ndarray) -> np.ndarray:
            return velocity(q, am)
    elif power == -1:

        def func(q: np.ndarray) -> np.ndarray:
            return 1.0 / velocity(q, am)
    else:
        raise ValueError(f"Unsupported velocity power {power}")
    return func


def make_velocity_condition(am: float, power: int) -> MatchCondition:
    if power == 1:
        func = make_velocity_function(am, power)
        label = f"v(am={am:g})"
    elif power == -1:
        func = make_velocity_function(am, power)
        label = f"1/v(am={am:g})"
    else:
        raise ValueError(f"Unsupported velocity power {power}")
    target = exact_integral(func)
    return MatchCondition(label=label, func=func, target=target, scale=max(1.0, abs(target)))


def match_residuals(q: np.ndarray, weights: np.ndarray, conditions: tuple[MatchCondition, ...]) -> np.ndarray:
    return np.array(
        [(np.dot(weights, condition.func(q)) - condition.target) / condition.scale for condition in conditions]
    )


def make_benchmark_function(exponent: int, am: float, velocity_power: int) -> Callable[[np.ndarray], np.ndarray]:
    if velocity_power == 0:

        def benchmark_func(q: np.ndarray) -> np.ndarray:
            return q**exponent
    elif velocity_power == 1:

        def benchmark_func(q: np.ndarray) -> np.ndarray:
            return q**exponent * velocity(q, am)
    elif velocity_power == -1:

        def benchmark_func(q: np.ndarray) -> np.ndarray:
            return q**exponent / velocity(q, am)
    else:
        raise ValueError(f"Unsupported benchmark velocity power {velocity_power}")
    return benchmark_func


def make_velocity_power_function(exponent: int, am: float, velocity_power: int) -> Callable[[np.ndarray], np.ndarray]:
    def func(q: np.ndarray) -> np.ndarray:
        return q**exponent * velocity(q, am) ** velocity_power

    return func


def make_actual_kernel_functions(am_values: np.ndarray | tuple[float, ...]) -> list[KernelFunction]:
    """Return q-kernel functions that appear in the massive-neutrino perturbation integrals.

    The scalar/tensor source sums in fortran/equations.f90 use velocity factors 1/v, v and v**2.
    The perturbatively relativistic tail and pinudot terms additionally bring q^-2 corrections,
    including q^-2 v^3 from vdot. These kernels are only a proxy for the full CAMB error because
    the perturbation multipoles are not q-independent, but they are a closer proxy than the original
    sparse {1, q^-2} x {1, v, 1/v} benchmark.
    """

    kernels: list[KernelFunction] = []
    for exponent in (-4, -2, -1, 0, 1, 2):
        func = make_benchmark_function(exponent, am=1.0, velocity_power=0)
        kernels.append(KernelFunction(f"q^{exponent}", f"q^{exponent}", func, normalized_moment(exponent)))
    for am in am_values:
        for exponent in (0, -2):
            for velocity_power in (-1, 1, 2):
                formula = f"q^{exponent} v^{velocity_power}"
                if velocity_power == -1:
                    formula = f"q^{exponent} / v"
                elif velocity_power == 1:
                    formula = f"q^{exponent} v"
                func = make_velocity_power_function(exponent, float(am), velocity_power)
                kernels.append(
                    KernelFunction(
                        label=f"{formula}(am={am:.6g})",
                        formula=formula,
                        func=func,
                        target=exact_integral(func),
                    )
                )
        func = make_velocity_power_function(-2, float(am), 3)
        kernels.append(
            KernelFunction(
                label=f"q^-2 v^3(am={am:.6g}) [vdot]",
                formula="q^-2 v^3",
                func=func,
                target=exact_integral(func),
            )
        )
    return kernels


def benchmark_functions() -> list[tuple[str, Callable[[np.ndarray], np.ndarray], float]]:
    benchmarks: list[tuple[str, Callable[[np.ndarray], np.ndarray], float]] = []
    for am in BENCHMARK_AM:
        for exponent in (0, -2):
            for velocity_power in (-1, 0, 1):
                if velocity_power == 0:
                    label = f"q^{exponent}(am={am:.3g})"
                elif velocity_power == 1:
                    label = f"q^{exponent} v(am={am:.3g})"
                else:
                    label = f"q^{exponent} / v(am={am:.3g})"
                func = make_benchmark_function(exponent, am, velocity_power)
                benchmarks.append((label, func, exact_integral(func)))
    return benchmarks


BENCHMARKS = benchmark_functions()


def solve_system(
    name: str,
    conditions: tuple[MatchCondition, ...],
    initial_guess: np.ndarray,
    unpack: Callable[[np.ndarray], tuple[np.ndarray, np.ndarray]],
) -> QuadratureRule:
    solution = least_squares(lambda variables: match_residuals(*unpack(variables), conditions), initial_guess)
    if not solution.success:
        raise RuntimeError(f"{name} solve failed: {solution.message}")
    q, weights = unpack(solution.x)
    residuals = match_residuals(q, weights, conditions)
    return QuadratureRule(
        name=name,
        q=q,
        weights=weights,
        condition_labels=tuple(condition.label for condition in conditions),
        residuals=residuals,
    )


def solve_weights_for_exact_moments(q: np.ndarray, exponents: tuple[int, ...]) -> np.ndarray:
    matrix = np.array([[q_i**exponent for q_i in q] for exponent in exponents])
    rhs = np.array([normalized_moment(exponent) for exponent in exponents])
    return np.linalg.solve(matrix, rhs)


def pack_free_four_point_nodes(q: np.ndarray) -> np.ndarray:
    return np.array(
        [
            np.log(q[0]),
            np.log(q[1] - q[0]),
            np.log(q[2] - q[1]),
            np.log(q[3] - q[2]),
        ]
    )


def unpack_free_four_point_nodes(variables: np.ndarray) -> np.ndarray:
    q1 = np.exp(variables[0])
    q2 = q1 + np.exp(variables[1])
    q3 = q2 + np.exp(variables[2])
    q4 = q3 + np.exp(variables[3])
    return np.array([q1, q2, q3, q4])


def pack_free_five_point_nodes(q: np.ndarray) -> np.ndarray:
    return np.array(
        [
            np.log(q[0]),
            np.log(q[1] - q[0]),
            np.log(q[2] - q[1]),
            np.log(q[3] - q[2]),
            np.log(q[4] - q[3]),
        ]
    )


def unpack_free_five_point_nodes(variables: np.ndarray) -> np.ndarray:
    q1 = np.exp(variables[0])
    q2 = q1 + np.exp(variables[1])
    q3 = q2 + np.exp(variables[2])
    q4 = q3 + np.exp(variables[3])
    q5 = q4 + np.exp(variables[4])
    return np.array([q1, q2, q3, q4, q5])


def solve_exact_moment_rule(
    name: str,
    weight_exponents: tuple[int, ...],
    constrained_exponents: tuple[int, ...],
    fit_conditions: tuple[MatchCondition, ...],
    initial_guess: np.ndarray,
    unpack_nodes: Callable[[np.ndarray], np.ndarray],
) -> QuadratureRule:
    all_exact_exponents = weight_exponents + constrained_exponents
    exact_conditions = tuple(make_moment_condition(exponent) for exponent in all_exact_exponents)
    rng = np.random.default_rng(0)

    def rule_from_variables(variables: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        q = unpack_nodes(variables)
        weights = solve_weights_for_exact_moments(q, weight_exponents)
        return q, weights

    def objective(variables: np.ndarray) -> float:
        q, weights = rule_from_variables(variables)
        if np.any(weights <= 0.0):
            return 1e6 + float(np.sum(np.square(np.minimum(weights, 0.0))))
        residuals = match_residuals(q, weights, fit_conditions)
        return 0.5 * float(np.dot(residuals, residuals))

    constraints = [
        {
            "type": "eq",
            "fun": lambda variables, exponent=exponent: (
                np.dot(rule_from_variables(variables)[1], rule_from_variables(variables)[0] ** exponent)
                - normalized_moment(exponent)
            ),
        }
        for exponent in constrained_exponents
    ]
    best_solution = None
    start_points = [initial_guess]
    start_points.extend(initial_guess + rng.normal(scale=0.25, size=initial_guess.shape) for _ in range(23))
    for start_point in start_points:
        solution = minimize(
            objective,
            start_point,
            method="SLSQP",
            constraints=constraints,
            options={"maxiter": 2000},
        )
        if not solution.success:
            continue
        if best_solution is None or solution.fun < best_solution.fun:
            best_solution = solution
    if best_solution is None:
        raise RuntimeError(f"{name} solve failed: {solution.message}")
    q, weights = rule_from_variables(best_solution.x)
    residuals = match_residuals(q, weights, exact_conditions + fit_conditions)
    return QuadratureRule(
        name=name,
        q=q,
        weights=weights,
        condition_labels=tuple(condition.label for condition in exact_conditions + fit_conditions),
        residuals=residuals,
    )


def unpack_three_point(variables: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    q1 = np.exp(variables[0])
    q2 = q1 + np.exp(variables[1])
    q3 = q2 + np.exp(variables[2])
    weights = np.exp(variables[3:6])
    return np.array([q1, q2, q3]), weights


def solve_three_point() -> QuadratureRule:
    conditions = tuple(make_moment_condition(exponent) for exponent in (-4, -2, -1, 0, 1, 2))
    initial_guess = np.log(np.array([0.8, 2.2, 4.7, 0.08, 3.0, 2.2]))
    return solve_system("3-point", conditions, initial_guess, unpack_three_point)


def solve_code_four_point() -> QuadratureRule:
    rule = solve_four_point_physics_targeted()
    return QuadratureRule("4-point", rule.q, rule.weights, rule.condition_labels, rule.residuals)


def solve_code_five_point() -> QuadratureRule:
    rule = solve_five_point_exact_ls()
    return QuadratureRule("5-point", rule.q, rule.weights, rule.condition_labels, rule.residuals)


def solve_three_point_physics_targeted() -> QuadratureRule:
    conditions = (
        *(make_moment_condition(exponent) for exponent in (-4, -2, -1, 0)),
        *(make_velocity_condition(am, 1) for am in REPRESENTATIVE_AMS),
        *(make_velocity_condition(am, -1) for am in REPRESENTATIVE_AMS),
    )
    base_rule = solve_three_point()
    initial_guess = np.array(
        [
            np.log(base_rule.q[0]),
            np.log(base_rule.q[1] - base_rule.q[0]),
            np.log(base_rule.q[2] - base_rule.q[1]),
            *np.log(base_rule.weights),
        ]
    )
    return solve_system("3-point physics-targeted", conditions, initial_guess, unpack_three_point)


def unpack_four_point(variables: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    q1 = 0.7
    q4 = 12.0
    q2 = q1 + sigmoid(variables[0]) * (q4 - q1)
    q3 = q2 + sigmoid(variables[1]) * (q4 - q2)
    weights = np.exp(variables[2:6])
    return np.array([q1, q2, q3, q4]), weights


def solve_four_point() -> QuadratureRule:
    conditions = tuple(make_moment_condition(exponent) for exponent in (-4, -2, -1, 0, 1, 2))
    initial_guess = np.array([-1.6, -0.6, np.log(0.02), np.log(1.8), np.log(3.5), np.log(0.29)])
    return solve_system("4-point", conditions, initial_guess, unpack_four_point)


def pack_free_four_point_guess(q: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return np.array(
        [
            np.log(q[0]),
            np.log(q[1] - q[0]),
            np.log(q[2] - q[1]),
            np.log(q[3] - q[2]),
            *np.log(weights),
        ]
    )


def unpack_free_four_point(variables: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    q1 = np.exp(variables[0])
    q2 = q1 + np.exp(variables[1])
    q3 = q2 + np.exp(variables[2])
    q4 = q3 + np.exp(variables[3])
    weights = np.exp(variables[4:8])
    return np.array([q1, q2, q3, q4]), weights


def solve_four_point_power_only() -> QuadratureRule:
    conditions = tuple(make_moment_condition(exponent) for exponent in (-4, -2, -1, 0, 1, 2, 3, 4))
    base_rule = solve_four_point()
    initial_guess = pack_free_four_point_guess(base_rule.q, base_rule.weights)
    return solve_system("4-point power-only", conditions, initial_guess, unpack_free_four_point)


def solve_four_point_physics_targeted() -> QuadratureRule:
    conditions = (
        *(make_moment_condition(exponent) for exponent in (-4, -2, -1, 0, 1, 2)),
        *(make_velocity_condition(am, 1) for am in REPRESENTATIVE_AMS),
        *(make_velocity_condition(am, -1) for am in REPRESENTATIVE_AMS),
    )
    base_rule = solve_four_point()
    initial_guess = pack_free_four_point_guess(base_rule.q, base_rule.weights)
    return solve_system("4-point physics-targeted", conditions, initial_guess, unpack_free_four_point)


def solve_four_point_exact_ls() -> QuadratureRule:
    fit_conditions = (
        *(make_velocity_condition(am, 1) for am in REPRESENTATIVE_AMS),
        *(make_velocity_condition(am, -1) for am in REPRESENTATIVE_AMS),
    )
    base_rule = solve_four_point()
    initial_guess = pack_free_four_point_nodes(base_rule.q)
    return solve_exact_moment_rule(
        "4-point exact+LS",
        (-4, -2, -1, 0),
        (1, 2),
        fit_conditions,
        initial_guess,
        unpack_free_four_point_nodes,
    )


def unpack_five_point(variables: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    q2 = 2.0
    q3 = 4.0
    q5 = 13.0
    q1 = sigmoid(variables[0]) * q2
    q4 = q3 + sigmoid(variables[1]) * (q5 - q3)
    weights = np.exp(variables[2:7])
    return np.array([q1, q2, q3, q4, q5]), weights


def pack_free_five_point_guess(q: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return np.array(
        [
            np.log(q[0]),
            np.log(q[1] - q[0]),
            np.log(q[2] - q[1]),
            np.log(q[3] - q[2]),
            np.log(q[4] - q[3]),
            *np.log(weights),
        ]
    )


def unpack_free_five_point(variables: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    q1 = np.exp(variables[0])
    q2 = q1 + np.exp(variables[1])
    q3 = q2 + np.exp(variables[2])
    q4 = q3 + np.exp(variables[3])
    q5 = q4 + np.exp(variables[4])
    weights = np.exp(variables[5:10])
    return np.array([q1, q2, q3, q4, q5]), weights


def solve_five_point() -> QuadratureRule:
    conditions = tuple(make_moment_condition(exponent) for exponent in (-4, -2, -1, 0, 1, 2, 3))
    initial_guess = np.array([-0.9, -0.6, np.log(0.008), np.log(0.69), np.log(2.8), np.log(2.05), np.log(0.13)])
    return solve_system("5-point", conditions, initial_guess, unpack_five_point)


def solve_five_point_physics_targeted() -> QuadratureRule:
    conditions = (
        *(make_moment_condition(exponent) for exponent in (-4, -2, -1, 0, 1, 2, 3)),
        *(make_velocity_condition(am, 1) for am in REPRESENTATIVE_AMS),
        *(make_velocity_condition(am, -1) for am in REPRESENTATIVE_AMS),
    )
    base_rule = solve_five_point()
    initial_guess = pack_free_five_point_guess(base_rule.q, base_rule.weights)
    return solve_system("5-point physics-targeted", conditions, initial_guess, unpack_free_five_point)


def solve_five_point_exact_ls() -> QuadratureRule:
    fit_conditions = (
        *(make_velocity_condition(am, 1) for am in REPRESENTATIVE_AMS),
        *(make_velocity_condition(am, -1) for am in REPRESENTATIVE_AMS),
    )
    base_rule = solve_five_point()
    initial_guess = pack_free_five_point_nodes(base_rule.q)
    return solve_exact_moment_rule(
        "5-point exact+LS",
        (-4, -2, -1, 0, 1),
        (2,),
        fit_conditions,
        initial_guess,
        unpack_free_five_point_nodes,
    )


def make_uniform_rule(nqmax: int) -> QuadratureRule:
    # fortran/massive_neutrinos.f90 computes this%nqmax/5 with integer division
    # (truncating), not a true division, so match that here for exact reproduction.
    dq = (12 + nqmax // 5) / nqmax
    q = np.array([(index + 0.5) * dq for index in range(nqmax)])
    dlfdlq = -q / (1.0 + np.exp(-q))
    raw_kernel = dq * q**3 / (np.exp(q) + 1.0) * (-0.25 * dlfdlq)
    conditions = tuple(make_moment_condition(exponent) for exponent in (-4, -2, -1, 0, 1, 2))
    weights = raw_kernel / FERMI_DIRAC_CONST
    residuals = match_residuals(q, weights, conditions)
    return QuadratureRule(
        name=f"{nqmax}-point uniform",
        q=q,
        weights=weights,
        condition_labels=tuple(condition.label for condition in conditions),
        residuals=residuals,
    )


def format_fortran_array(values: np.ndarray, precision: int = 6) -> str:
    return "(/" + ", ".join(f"{value:.{precision}g}" for value in values) + "/)"


def print_rule(rule: QuadratureRule) -> None:
    raw_kernel = rule.weights * FERMI_DIRAC_CONST
    print(rule.name)
    print(f"  q      = {format_fortran_array(rule.q)}")
    print(f"  kernel = {format_fortran_array(raw_kernel)}")
    print(f"  normalized kernel = {format_fortran_array(rule.weights)}")
    print("  matched conditions")
    for label, residual in zip(rule.condition_labels, rule.residuals, strict=True):
        condition = next(condition for condition in all_conditions() if condition.label == label)
        matched = np.dot(rule.weights, condition.func(rule.q))
        print(
            f"    {label:<14}: exact={condition.target:.16e}  matched={matched:.16e}  scaled residual={residual:+.3e}"
        )
    print(f"  max |scaled residual| = {np.max(np.abs(rule.residuals)):.3e}")
    print()


def print_rule_summary(rule: QuadratureRule) -> None:
    raw_kernel = rule.weights * FERMI_DIRAC_CONST
    print(rule.name)
    print(f"  q      = {format_fortran_array(rule.q, precision=16)}")
    print(f"  kernel = {format_fortran_array(raw_kernel, precision=16)}")
    print(f"  normalized kernel = {format_fortran_array(rule.weights, precision=16)}")
    print()


def all_conditions() -> tuple[MatchCondition, ...]:
    conditions = [*(make_moment_condition(exponent) for exponent in (-4, -2, -1, 0, 1, 2, 3, 4))]
    for am in REPRESENTATIVE_AMS:
        conditions.extend(make_velocity_condition(am, power) for power in (1, -1))
    return tuple(conditions)


def quadrature_integral(rule: QuadratureRule, func: Callable[[np.ndarray], np.ndarray]) -> float:
    return float(np.dot(rule.weights, func(rule.q)))


def benchmark_summary(rule: QuadratureRule) -> BenchmarkSummary:
    relative_errors: list[float] = []
    worst_case_label = ""
    worst_case_error = -1.0
    for label, func, exact_value in BENCHMARKS:
        matched = quadrature_integral(rule, func)
        relative_error = abs(matched - exact_value) / abs(exact_value)
        relative_errors.append(relative_error)
        if relative_error > worst_case_error:
            worst_case_error = relative_error
            worst_case_label = label
    return BenchmarkSummary(
        max_abs_relative_error=worst_case_error,
        rms_relative_error=float(np.sqrt(np.mean(np.square(relative_errors)))),
        worst_case_label=worst_case_label,
    )


def kernel_relative_errors(rule: QuadratureRule, kernels: list[KernelFunction]) -> np.ndarray:
    return np.array([(quadrature_integral(rule, kernel.func) - kernel.target) / kernel.target for kernel in kernels])


def kernel_summary(rule: QuadratureRule, kernels: list[KernelFunction]) -> BenchmarkSummary:
    relative_errors = kernel_relative_errors(rule, kernels)
    worst_index = int(np.argmax(np.abs(relative_errors)))
    return BenchmarkSummary(
        max_abs_relative_error=float(np.max(np.abs(relative_errors))),
        rms_relative_error=float(np.sqrt(np.mean(relative_errors**2))),
        worst_case_label=kernels[worst_index].label,
    )


def print_kernel_functions(kernels: list[KernelFunction]) -> None:
    print("Actual-kernel proxy family")
    print("  Integral form: integral dq normalized_kernel_density(q) f(q; am)")
    print("  normalized_kernel_density(q) = q^4 exp(-q) / [4 (1 + exp(-q))^2 FermiDiracConst]")
    print("  v(q, am) = 1 / sqrt(1 + (am / q)^2)")
    for kernel in kernels:
        print(f"  {kernel.label:<34} formula={kernel.formula:<10} target={kernel.target:.16e}")
    print()


def pack_free_rule_guess(q: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return np.array([np.log(q[0]), *np.log(np.diff(q)), *np.log(weights)])


def unpack_free_rule(variables: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    node_count = variables.size // 2
    q = np.cumsum(np.exp(variables[:node_count]))
    weights = np.exp(variables[node_count:])
    return q, weights


def soft_max_abs(values: np.ndarray, alpha: float = 80.0) -> float:
    scaled = alpha * np.abs(values)
    max_scaled = float(np.max(scaled))
    return (max_scaled + float(np.log(np.sum(np.exp(scaled - max_scaled))))) / alpha


def minimax_objective(variables: np.ndarray, kernels: list[KernelFunction], max_q: float) -> float:
    q, weights = unpack_free_rule(variables)
    if q[-1] > max_q:
        return 1e3 + q[-1]
    return soft_max_abs(kernel_relative_errors(QuadratureRule("", q, weights, (), np.array([])), kernels))


def minimax_start_rules(node_count: int) -> list[QuadratureRule]:
    if node_count == 4:
        return [solve_four_point(), solve_code_four_point(), solve_four_point_exact_ls(), make_uniform_rule(4)]
    if node_count == 5:
        return [solve_five_point(), solve_five_point_physics_targeted(), solve_code_five_point(), make_uniform_rule(5)]
    raise ValueError(f"Unsupported minimax node count {node_count}")


def polish_minimax(start: np.ndarray, kernels: list[KernelFunction], max_q: float) -> tuple[np.ndarray, float]:
    result = minimize(
        lambda variables: minimax_objective(variables, kernels, max_q),
        start,
        method="Nelder-Mead",
        options={"maxiter": 20000, "xatol": 1e-12, "fatol": 1e-12},
    )
    result2 = minimize(
        lambda variables: minimax_objective(variables, kernels, max_q),
        result.x,
        method="BFGS",
        options={"maxiter": 5000, "gtol": 1e-12},
    )
    best = result2 if result2.fun < result.fun else result
    return best.x, float(best.fun)


def solve_minimax_rule(
    node_count: int,
    kernels: list[KernelFunction],
    *,
    global_search: bool,
    seed: int,
    max_q: float,
) -> QuadratureRule:
    best_variables: np.ndarray | None = None
    best_value = np.inf
    for start_rule in minimax_start_rules(node_count):
        variables, value = polish_minimax(pack_free_rule_guess(start_rule.q, start_rule.weights), kernels, max_q)
        if value < best_value:
            best_variables = variables
            best_value = value

    if global_search:
        bounds = [(np.log(0.2), np.log(12.0))] * node_count
        bounds.extend((np.log(1e-7), np.log(1.5)) for _ in range(node_count))
        solution = differential_evolution(
            lambda variables: minimax_objective(variables, kernels, max_q),
            bounds,
            seed=seed,
            maxiter=450,
            popsize=25,
            tol=1e-7,
            polish=False,
            updating="immediate",
            workers=1,
        )
        variables, value = polish_minimax(solution.x, kernels, max_q)
        if value < best_value:
            best_variables = variables

    if best_variables is None:
        raise RuntimeError(f"{node_count}-point minimax solve failed")
    q, weights = unpack_free_rule(best_variables)
    residuals = kernel_relative_errors(QuadratureRule("", q, weights, (), np.array([])), kernels)
    return QuadratureRule(f"{node_count}-point dense minimax", q, weights, tuple(k.label for k in kernels), residuals)


def print_kernel_comparison(rules: list[QuadratureRule], kernels: list[KernelFunction], *, top_errors: int = 6) -> None:
    print("Actual-kernel proxy benchmark")
    for rule in rules:
        summary = kernel_summary(rule, kernels)
        errors = kernel_relative_errors(rule, kernels)
        print(rule.name)
        print(f"  max |relative error| = {summary.max_abs_relative_error:.3e}")
        print(f"  rms relative error   = {summary.rms_relative_error:.3e}")
        print(f"  worst case           = {summary.worst_case_label}")
        for index in np.argsort(np.abs(errors))[-top_errors:][::-1]:
            print(f"    {kernels[index].label:<34} {errors[index]:+.3e}")
    print()


def print_benchmark_comparison(rules: list[QuadratureRule]) -> None:
    print("Benchmark family: {1, q^-2} x {1, v, 1/v} over am in [0.03, 12]")
    for rule in rules:
        summary = benchmark_summary(rule)
        print(rule.name)
        print(f"  max |relative error| = {summary.max_abs_relative_error:.3e}")
        print(f"  rms relative error   = {summary.rms_relative_error:.3e}")
        print(f"  worst case           = {summary.worst_case_label}")
    print()


def compiled_camb_cases(accuracy_boost: float, reference_nqmax: int) -> list[CompiledCambCase]:
    return [
        CompiledCambCase(
            name="3-point default",
            neutrino_q_boost=0.8,
            expected_nqmax=3,
            rule=solve_three_point(),
        ),
        CompiledCambCase(
            name="4-point default",
            neutrino_q_boost=1.0,
            expected_nqmax=4,
            rule=solve_code_four_point(),
        ),
        CompiledCambCase(
            name="5-point default",
            neutrino_q_boost=2.1,
            expected_nqmax=5,
            rule=solve_code_five_point(),
        ),
        CompiledCambCase(
            name=f"{reference_nqmax}-point reference",
            neutrino_q_boost=reference_nqmax / 10.0 / accuracy_boost,
            expected_nqmax=reference_nqmax,
            rule=make_uniform_rule(reference_nqmax),
        ),
    ]


def compiled_camb_models() -> list[CompiledCambModel]:
    return [
        CompiledCambModel(mnu=0.06, num_massive_neutrinos=1),
        CompiledCambModel(mnu=0.12, num_massive_neutrinos=1),
        CompiledCambModel(mnu=0.3, num_massive_neutrinos=3),
    ]


def build_camb_params(
    case: CompiledCambCase,
    model: CompiledCambModel,
    accuracy_boost: float,
    l_sample_boost: float,
    lmax: int,
    pk_kmax: float,
):
    import camb

    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=67.5,
        ombh2=0.022,
        omch2=0.122,
        mnu=model.mnu,
        omk=0.0,
        tau=0.06,
        num_massive_neutrinos=model.num_massive_neutrinos,
    )
    pars.InitPower.set_params(As=2e-9, ns=0.965)
    pars.set_for_lmax(lmax, lens_potential_accuracy=1)
    pars.set_accuracy(AccuracyBoost=accuracy_boost, lSampleBoost=l_sample_boost, lAccuracyBoost=1.0)
    pars.Accuracy.neutrino_q_boost = case.neutrino_q_boost
    pars.set_matter_power(redshifts=[0.0], kmax=pk_kmax, accurate_massive_neutrino_transfers=True)
    return pars


def stable_relative_error(values: np.ndarray, reference: np.ndarray) -> np.ndarray:
    scale = np.maximum(np.abs(reference), 1e-12 * np.max(np.abs(reference), axis=0, keepdims=True))
    return np.abs(values - reference) / scale


def summarize_relative_error(values: np.ndarray, reference: np.ndarray) -> RelativeErrorSummary:
    relative_error = stable_relative_error(values, reference)
    return RelativeErrorSummary(
        max_abs_relative_error=float(np.max(relative_error)),
        rms_relative_error=float(np.sqrt(np.mean(np.square(relative_error)))),
    )


def print_compiled_case_samples(case: CompiledCambCase) -> None:
    print(f"{case.name} (expected nqmax={case.expected_nqmax})")
    if case.rule is None:
        print("  sample points unavailable")
    elif case.expected_nqmax <= 5:
        print(f"  q      = {format_fortran_array(case.rule.q)}")
        print(f"  kernel = {format_fortran_array(case.rule.weights * FERMI_DIRAC_CONST)}")
    else:
        print(f"  q[0:6] = {format_fortran_array(case.rule.q[:6])} ...")
        print(f"  q[-3:] = {format_fortran_array(case.rule.q[-3:])}")
    print()


def run_compiled_camb_comparison(
    models: list[CompiledCambModel],
    accuracy_boost: float,
    l_sample_boost: float,
    lmax: int,
    pk_kmax: float,
    reference_nqmax: int,
) -> None:
    import camb

    cases = compiled_camb_cases(accuracy_boost, reference_nqmax)
    for case in cases:
        print_compiled_case_samples(case)
    print()
    for model in models:
        results_by_case: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for case in cases:
            pars = build_camb_params(case, model, accuracy_boost, l_sample_boost, lmax, pk_kmax)
            results = camb.get_results(pars)
            cls = results.get_cmb_power_spectra(pars, spectra=["lensed_scalar"], CMB_unit="muK")["lensed_scalar"]
            kh, _, pk = results.get_matter_power_spectrum(
                minkh=1e-4,
                maxkh=pk_kmax,
                npoints=250,
                have_power_spectra=True,
            )
            results_by_case[case.name] = (cls, kh, pk)

        reference_name = f"{reference_nqmax}-point reference"
        reference_cls, reference_kh, reference_pk = results_by_case[reference_name]
        print(f"Compiled CAMB comparison against nqmax={reference_nqmax} reference")
        print("  Note: this compares the current compiled CAMB defaults only.")
        print(
            "  Testing alternative candidate rules in compiled CAMB requires editing fortran/massive_neutrinos.f90 and rebuilding."
        )
        print(
            f"  mnu={model.mnu:g} eV, num_massive_neutrinos={model.num_massive_neutrinos}, "
            f"AccuracyBoost={accuracy_boost:g}, lSampleBoost={l_sample_boost:g}, lmax={lmax}, pk_kmax={pk_kmax:g}"
        )
        print()
        for case in cases:
            if case.name == reference_name:
                continue
            cls, kh, pk = results_by_case[case.name]
            if not np.allclose(kh, reference_kh):
                raise RuntimeError(f"Matter power k-grid mismatch for case {case.name} at mnu={model.mnu}")
            cls_summary = summarize_relative_error(cls[2:, :3], reference_cls[2:, :3])
            pk_summary = summarize_relative_error(pk, reference_pk)
            print(case.name)
            print(f"  lensed Cl max |relative error| = {cls_summary.max_abs_relative_error:.3e}")
            print(f"  lensed Cl rms relative error   = {cls_summary.rms_relative_error:.3e}")
            print(f"  matter P(k) max |relative error| = {pk_summary.max_abs_relative_error:.3e}")
            print(f"  matter P(k) rms relative error   = {pk_summary.rms_relative_error:.3e}")
            print()
        print()


EPILOG = """\
What it does:
  - By default (--case all), reproduces the node values and kernels currently
    hard-coded in fortran/massive_neutrinos.f90: prints the current 3-point
    rule, the adopted 4-point rule, and the adopted 5-point rule.
  - Can still solve the original and alternative 4-point and 5-point candidate
    rules for comparison (--case 4-compare / 5-compare).
  - Compares candidate rules on a perturbation-motivated benchmark family
    {1, q^-2} x {1, v, 1/v} over a range of am.
  - Can dump a closer actual-kernel proxy family including v^2 and q^-2 v^3
    terms used by the perturbation equations (--case dump-kernels).
  - Can run dense minimax searches over this proxy family for diagnostic
    candidate rules (--case 4-minimax / 5-minimax / minimax-compare).
  - Can run compiled CAMB against an nqmax reference (default 80) to validate
    the currently compiled defaults (--case camb-compare).

Current compiled-CAMB results with AccuracyBoost=1.2, lSampleBoost=3, and
nqmax=80 as reference indicate:
  - the 4-point physics-targeted rule remains the best 4-point option overall.
  - the 4-point exact+LS rule improves on the fixed exact 4-point rule, but
    not enough to beat the 4-point physics-targeted rule.
  - the 5-point exact+LS rule is the strongest 5-point option for mnu=0.06
    and 0.12, and is roughly tied with the 5-point physics-targeted rule at
    mnu=0.3.
  - the 3-point physics-targeted rule is not a robust replacement for the
    default 3-point rule.
  - cross-checks against nqmax=160 show nqmax=80 is stable enough for these
    rankings, with remaining reference-level differences typically at the
    1e-5 level.
  - dense minimax rules can improve the proxy-kernel errors, but compiled-CAMB
    testing shows these proxy improvements do not translate into better
    spectra; the current hard-coded 4- and 5-point rules remain preferred.

Notes:
  - The script fits the non-relativistic matching part using representative
    am values (0.3, 1, 3, 10).
  - The dense minimax proxy is useful as a diagnostic, but the evolved
    perturbation multipoles are not q-independent, so compiled-CAMB
    validation is required before adopting any minimax candidate.
  - --case camb-compare checks the currently compiled CAMB defaults only.
  - Testing alternative candidate rules in compiled CAMB requires temporarily
    editing fortran/massive_neutrinos.f90 and rebuilding.

Example invocations:
  python scripts/nu_integration_kernels.py
  python scripts/nu_integration_kernels.py --case 3
  python scripts/nu_integration_kernels.py --case 4-compare
  python scripts/nu_integration_kernels.py --case 5-compare
  python scripts/nu_integration_kernels.py --case camb-compare
  python scripts/nu_integration_kernels.py --case camb-compare --reference-nqmax 160
  python scripts/nu_integration_kernels.py --case dump-kernels
  python scripts/nu_integration_kernels.py --case minimax-compare --minimax-global-search
  python scripts/nu_integration_kernels.py --case 4
  python scripts/nu_integration_kernels.py --case 5
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reproduce the massive-neutrino momentum sampling values used in CAMB's "
            "TNuPerturbations_init routine from the Appendix A moment conditions of "
            "arXiv:1201.3654."
        ),
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--case",
        choices=[
            "3",
            "4",
            "4-physics",
            "4-exact-ls",
            "4-compare",
            "5",
            "5-physics",
            "5-exact-ls",
            "5-compare",
            "dump-kernels",
            "4-minimax",
            "5-minimax",
            "minimax-compare",
            "all",
            "camb-compare",
        ],
        default="all",
        help=(
            "Which rule to solve. By default this prints the 3-point rule and the adopted 4-point and 5-point rules used in the code."
        ),
    )
    parser.add_argument("--accuracy-boost", type=float, default=1.2, help="AccuracyBoost for compiled CAMB runs")
    parser.add_argument("--l-sample-boost", type=float, default=3.0, help="lSampleBoost for compiled CAMB runs")
    parser.add_argument("--lmax", type=int, default=3000, help="Maximum multipole for the compiled CAMB comparison")
    parser.add_argument("--pk-kmax", type=float, default=5.0, help="Maximum k/h for matter power comparison")
    parser.add_argument("--reference-nqmax", type=int, default=80, help="nqmax to use for compiled CAMB reference runs")
    parser.add_argument("--minimax-am-min", type=float, default=0.03, help="Minimum am for dense minimax kernels")
    parser.add_argument("--minimax-am-max", type=float, default=20.0, help="Maximum am for dense minimax kernels")
    parser.add_argument("--minimax-am-count", type=int, default=41, help="Number of log-spaced am values")
    parser.add_argument("--minimax-max-q", type=float, default=25.0, help="Maximum allowed largest q node")
    parser.add_argument("--minimax-seed", type=int, default=12345, help="Random seed for differential evolution")
    parser.add_argument(
        "--minimax-global-search",
        action="store_true",
        help="Run differential evolution before local polishing. Slower, but better for checking global minima.",
    )
    parser.add_argument("--top-errors", type=int, default=6, help="Number of worst kernel errors to print")
    parser.add_argument(
        "--mnu",
        type=float,
        nargs="+",
        default=None,
        help="Override the default neutrino masses in eV for compiled CAMB comparison runs",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    minimax_ams = np.geomspace(args.minimax_am_min, args.minimax_am_max, args.minimax_am_count)
    actual_kernels = make_actual_kernel_functions(minimax_ams)
    if args.case == "camb-compare":
        models = compiled_camb_models()
        if args.mnu is not None:
            models = [CompiledCambModel(mnu=mnu, num_massive_neutrinos=1) for mnu in args.mnu]
        run_compiled_camb_comparison(
            models,
            args.accuracy_boost,
            args.l_sample_boost,
            args.lmax,
            args.pk_kmax,
            args.reference_nqmax,
        )
        return
    if args.case == "dump-kernels":
        print_kernel_functions(actual_kernels)
        return
    if args.case in ("4-minimax", "5-minimax"):
        node_count = int(args.case[0])
        rule = solve_minimax_rule(
            node_count,
            actual_kernels,
            global_search=args.minimax_global_search,
            seed=args.minimax_seed,
            max_q=args.minimax_max_q,
        )
        print_rule_summary(rule)
        print_kernel_comparison([rule], actual_kernels, top_errors=args.top_errors)
        return
    if args.case == "minimax-compare":
        rules = [
            solve_four_point(),
            solve_code_four_point(),
            solve_four_point_exact_ls(),
            solve_minimax_rule(
                4,
                actual_kernels,
                global_search=args.minimax_global_search,
                seed=args.minimax_seed,
                max_q=args.minimax_max_q,
            ),
            solve_five_point(),
            solve_five_point_physics_targeted(),
            solve_code_five_point(),
            solve_minimax_rule(
                5,
                actual_kernels,
                global_search=args.minimax_global_search,
                seed=args.minimax_seed + 1,
                max_q=args.minimax_max_q,
            ),
        ]
        print_kernel_comparison(rules, actual_kernels, top_errors=args.top_errors)
        return

    solvers = {
        "3": [solve_three_point],
        "4": [solve_code_four_point],
        "4-physics": [solve_four_point_physics_targeted],
        "4-exact-ls": [solve_four_point_exact_ls],
        "4-compare": [
            solve_four_point,
            solve_code_four_point,
            solve_four_point_exact_ls,
        ],
        "5": [solve_code_five_point],
        "5-physics": [solve_five_point_physics_targeted],
        "5-exact-ls": [solve_five_point_exact_ls],
        "5-compare": [solve_five_point, solve_five_point_physics_targeted, solve_code_five_point],
        "all": [
            solve_three_point,
            solve_code_four_point,
            solve_code_five_point,
        ],
    }
    rules = [solver() for solver in solvers[args.case]]
    for rule in rules:
        print_rule(rule)
    if args.case == "4-compare":
        print_benchmark_comparison([rule for rule in rules if rule.name.startswith("4-point")])
    if args.case == "5-compare":
        print_benchmark_comparison([rule for rule in rules if rule.name.startswith("5-point")])


if __name__ == "__main__":
    main()
