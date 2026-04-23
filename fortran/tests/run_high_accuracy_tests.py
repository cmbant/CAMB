"""
High-accuracy CAMB regression helper.

Example workflow:
1) Generate reference outputs with default high-accuracy settings
   python fortran/tests/run_high_accuracy_tests.py generate_reference /tmp/camb_hiacc_work \\
       --reference_dir /tmp/camb_hiacc_ref --max_tests 3 --no_sources --clean
2) Compare default-accuracy runs to that reference
   python fortran/tests/run_high_accuracy_tests.py compare /tmp/camb_hiacc_work \\
       --reference_dir /tmp/camb_hiacc_ref --max_tests 3 --no_sources --clean

The compare command is expected to return non-zero whenever default-accuracy outputs differ from
the high-accuracy reference, and it always prints RMS/worst-case summaries.
"""

import argparse
import math
import os
import shutil
import stat
import sys
import tempfile
import time
from pathlib import Path

PASSTHROUGH_KEYS = {"CMB_outputscale"}

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
FORTRAN_DIR = os.path.abspath(os.path.join(TESTS_DIR, ".."))
REPO_ROOT = os.path.abspath(os.path.join(FORTRAN_DIR, ".."))
DEFAULT_WORK_DIR = os.path.join(FORTRAN_DIR, "testfiles")
DEFAULT_BASE_SETTINGS = os.path.join(REPO_ROOT, "inifiles", "params.ini")

if TESTS_DIR not in sys.path:
    sys.path.insert(0, TESTS_DIR)
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Generate high-accuracy reference outputs and compare default runs against them"
    )
    parser.add_argument("mode", choices=("generate_reference", "compare"))
    parser.add_argument("work_dir", nargs="?", default=DEFAULT_WORK_DIR, help="working directory")
    parser.add_argument("--reference_dir", "--reference-dir", required=True, help="reference root directory")
    parser.add_argument("--base_settings", "--base-settings", default=DEFAULT_BASE_SETTINGS, help="base params.ini")
    parser.add_argument("--out_files_dir", "--out-files-dir", default="test_outputs", help="output files directory")
    parser.add_argument("--clean", action="store_true", help="delete output directories before running")
    parser.add_argument("--runner", choices=("module", "command"), default="module", help="test execution backend")
    parser.add_argument("--prog", default="./camb", help="legacy executable path when --runner=command")
    parser.add_argument("--no_validate", "--no-validate", action="store_true", help="skip ini validation")
    parser.add_argument("--max_tests", "--max-tests", type=int, help="maximum number of tests to run")
    parser.add_argument("--no_sources", "--no-sources", action="store_true", help="turn off CAMB sources tests")
    parser.add_argument("--no_de", "--no-de", action="store_true", help="skip dark energy tests")
    parser.add_argument("--skip_non_flat", "--skip-non-flat", action="store_true", help="skip tests with non-zero omk")
    parser.add_argument(
        "--only_non_flat", "--only-non-flat", action="store_true", help="only run tests with non-zero omk"
    )
    parser.add_argument(
        "--override",
        "--set",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="override an ini parameter for every generated test; can be repeated",
    )
    parser.add_argument(
        "--diff_tolerance",
        "--diff-tolerance",
        type=float,
        default=1e-4,
        help="fallback numerical tolerance for compare mode",
    )
    parser.add_argument(
        "--verbose_diff_output",
        "--verbose-diff-output",
        "--verbose",
        action="store_true",
        help="print more diff details",
    )
    parser.add_argument("--num_diff", "--num-diff", action="store_true", help="use absolute diffs")

    parser.add_argument("--int_tol_boost", "--int-tol-boost", type=float, default=2.0)
    parser.add_argument("--l_sample_boost", "--l-sample-boost", type=float, default=2.0)
    parser.add_argument("--l_accuracy_boost", "--l-accuracy-boost", type=float, default=2.5)
    parser.add_argument("--lens_potential_accuracy", "--lens-potential-accuracy", type=float, default=20.0)
    parser.add_argument("--lens_margin", "--lens-margin", type=int, default=2050)
    parser.add_argument("--stability_factor", "--stability-factor", type=float, default=1.0)
    parser.add_argument("--do_late_rad_truncation", "--do-late-rad-truncation", action="store_true", default=False)
    parser.add_argument(
        "--skip_lensing_range_changes",
        "--skip-lensing-range-changes",
        action="store_true",
        help="do not modify lens_potential_accuracy or lens_margin",
    )
    parser.add_argument("--reference_ini_dir", "--reference-ini-dir", default="reference_inis")
    parser.add_argument("--reference_out_dir", "--reference-out-dir", default="test_outputs")
    return parser


def remove_tree(path):
    def onerror(func, target, exc_info):
        for _ in range(10):
            try:
                os.chmod(target, stat.S_IWRITE | stat.S_IREAD)
            except OSError:
                pass
            try:
                func(target)
                return
            except PermissionError:
                time.sleep(0.2)
        raise exc_info[1]

    shutil.rmtree(path, onerror=onerror)


def parse_override(override_text):
    if "=" not in override_text:
        raise ValueError(f"Override must use KEY=VALUE syntax: {override_text}")
    key, value = override_text.split("=", 1)
    return key.strip(), value.strip()


def apply_overrides(ini, overrides):
    for key, value in overrides:
        if key not in ini.params:
            ini.readOrder.append(key)
        ini.params[key] = value


def get_test_params(parsed_args):
    import CAMB_test_files as ctf

    ctf.args = argparse.Namespace(
        no_de=parsed_args.no_de, no_sources=parsed_args.no_sources, max_tests=parsed_args.max_tests
    )
    params = ctf.getTestParams()
    if parsed_args.skip_non_flat and parsed_args.only_non_flat:
        raise ValueError("Cannot use both --skip_non_flat and --only_non_flat")
    if parsed_args.skip_non_flat:
        return [p for p in params if not _is_non_flat(p)]
    if parsed_args.only_non_flat:
        return [p for p in params if _is_non_flat(p)]
    return params


def _is_non_flat(parameter_lines):
    for line in parameter_lines[1:]:
        if "=" not in line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        if key.lower() != "omk":
            continue
        try:
            return abs(float(value)) > 1e-12
        except ValueError:
            return True
    return False


def _apply_high_accuracy_settings(params, parsed_args):
    params.Accuracy.IntTolBoost = max(
        params.Accuracy.IntTolBoost,
        parsed_args.int_tol_boost * parsed_args.stability_factor,
    )
    params.Accuracy.lSampleBoost = max(
        params.Accuracy.lSampleBoost,
        parsed_args.l_sample_boost * parsed_args.stability_factor,
    )
    params.Accuracy.lAccuracyBoost = max(
        params.Accuracy.lAccuracyBoost,
        parsed_args.l_accuracy_boost * parsed_args.stability_factor,
    )
    params.DoLateRadTruncation = parsed_args.do_late_rad_truncation
    if not parsed_args.skip_lensing_range_changes:
        params.set_for_lmax(
            lmax=params.max_l,
            lens_potential_accuracy=parsed_args.lens_potential_accuracy * parsed_args.stability_factor,
            lens_margin=int(round(parsed_args.lens_margin * parsed_args.stability_factor)),
        )


def make_reference_inis(parsed_args, ini_dir, output_dir):
    from CAMB_test_files import write_flat_ini_file

    import camb
    from camb.inifile import IniFile

    os.makedirs(ini_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    overrides = [parse_override(item) for item in parsed_args.override]
    test_params = get_test_params(parsed_args)
    generated = []

    for parameter_lines in test_params:
        name = f"params_{parameter_lines[0]}"
        destination_ini = os.path.join(ini_dir, f"{name}.ini")
        output_root = os.path.join(output_dir, name)

        flat_lines = [
            f"output_root={output_root}",
            *parameter_lines[1:],
            f"DEFAULT({os.path.abspath(parsed_args.base_settings)})",
        ]

        with tempfile.NamedTemporaryFile("w", suffix=".ini", dir=ini_dir, delete=False, encoding="utf-8") as handle:
            temporary_ini = handle.name
        try:
            write_flat_ini_file(temporary_ini, flat_lines)
            flat_ini = IniFile(temporary_ini)
            apply_overrides(flat_ini, overrides)
            flat_ini.saveFile(temporary_ini)

            params = camb.read_ini(temporary_ini, no_validate=parsed_args.no_validate)
            _apply_high_accuracy_settings(params, parsed_args)
            camb.write_ini(params, destination_ini, validate=False)

            output_ini = IniFile(destination_ini)
            if "output_root" not in output_ini.params:
                output_ini.readOrder.insert(0, "output_root")
            output_ini.params["output_root"] = output_root
            for key in flat_ini.readOrder:
                if (
                    key == "output_root"
                    or key in PASSTHROUGH_KEYS
                    or key.startswith("output_")
                    or key.startswith("transfer_filename(")
                    or key.startswith("transfer_matterpower(")
                    or key == "transfer_interp_matterpower"
                    or key.endswith("_output_file")
                ):
                    if key not in output_ini.params:
                        output_ini.readOrder.append(key)
                    output_ini.params[key] = flat_ini.params[key]
            apply_overrides(output_ini, overrides)
            output_ini.saveFile(destination_ini)
        finally:
            os.remove(temporary_ini)

        generated.append(destination_ini)
    return generated


def run_reference_inis(reference_inis, parsed_args):
    import camb

    for ini in reference_inis:
        print(f"Running {os.path.basename(ini)}")
        if parsed_args.runner == "module":
            camb.run_ini(ini, no_validate=parsed_args.no_validate)
            continue

        import subprocess

        subprocess.run([parsed_args.prog, ini], check=True)


def run_camb_test_files(arguments):
    from CAMB_test_files import main as camb_test_files_main

    return camb_test_files_main(arguments)


def make_default_inis(parsed_args):
    from CAMB_test_files import write_flat_ini_file

    from camb.inifile import IniFile

    os.makedirs(parsed_args.work_dir, exist_ok=True)
    out_dir = os.path.join(parsed_args.work_dir, parsed_args.out_files_dir)
    os.makedirs(out_dir, exist_ok=True)

    generated = []
    overrides = [parse_override(item) for item in parsed_args.override]
    for pars in get_test_params(parsed_args):
        name = f"params_{pars[0]}"
        ini_path = os.path.join(parsed_args.work_dir, f"{name}.ini")
        write_flat_ini_file(
            ini_path,
            [
                f"output_root={os.path.join(out_dir, name)}",
                *pars[1:],
                f"DEFAULT({os.path.abspath(parsed_args.base_settings)})",
            ],
        )
        ini = IniFile(ini_path)
        apply_overrides(ini, overrides)
        ini.saveFile(ini_path)
        generated.append(ini_path)
    return generated


def run_inis(inis, parsed_args):
    import camb

    failures = 0
    for ini in inis:
        print(f"Running {os.path.basename(ini)}")
        try:
            if parsed_args.runner == "module":
                camb.run_ini(ini, no_validate=parsed_args.no_validate)
            else:
                import subprocess

                subprocess.run([parsed_args.prog, ini], check=True)
        except Exception as error:  # noqa: BLE001
            failures += 1
            print(f"  ERROR: {error}")
    return failures


def extend_common_compare_args(cmd_args, parsed_args, *, include_clean):
    cmd_args.extend(
        [
            parsed_args.work_dir,
            "--out_files_dir",
            parsed_args.out_files_dir,
            "--base_settings",
            parsed_args.base_settings,
            "--runner",
            parsed_args.runner,
            "--prog",
            parsed_args.prog,
        ]
    )
    if include_clean and parsed_args.clean:
        cmd_args.append("--clean")
    if parsed_args.no_validate:
        cmd_args.append("--no_validate")
    if parsed_args.max_tests is not None:
        cmd_args.extend(["--max_tests", str(parsed_args.max_tests)])
    if parsed_args.no_sources:
        cmd_args.append("--no_sources")
    if parsed_args.no_de:
        cmd_args.append("--no_de")
    if parsed_args.verbose_diff_output:
        cmd_args.append("--verbose_diff_output")
    if parsed_args.num_diff:
        cmd_args.append("--num_diff")
    for override in parsed_args.override:
        cmd_args.extend(["--override", override])


def _read_numeric_table(path):
    rows = []
    columns = None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                header_tokens = stripped.lstrip("#").split()
                if header_tokens:
                    columns = [token.replace("nu", "nu") for token in header_tokens]
                continue
            values = stripped.split()
            try:
                rows.append([float(item) for item in values])
            except ValueError:
                continue
    if columns is None and rows:
        columns = [str(index) for index in range(len(rows[0]))]
    return columns or [], rows


def _normalized_error(column, old_value, new_value, old_row):
    delta = abs(new_value - old_value)
    if column in {"L", "l", "ell", "k/h", "k"}:
        return None

    if "x" in column:
        first, second = column.split("x", 1)
        if first != second:
            auto_1 = old_row.get(f"{first}x{first}")
            auto_2 = old_row.get(f"{second}x{second}")
            if auto_1 is not None and auto_2 is not None:
                product = auto_1 * auto_2
                if product > 0:
                    return delta / math.sqrt(product)
    elif len(column) == 2 and column.isalpha():
        first, second = column[0], column[1]
        if first != second:
            auto_1 = old_row.get(first + first)
            auto_2 = old_row.get(second + second)
            if auto_1 is not None and auto_2 is not None:
                product = auto_1 * auto_2
                if product > 0:
                    return delta / math.sqrt(product)

    denom = abs(old_value)
    if denom > 1e-30:
        return delta / denom
    return delta


def summarize_differences(current_output_dir, reference_output_dir):
    current = Path(current_output_dir)
    reference = Path(reference_output_dir)

    current_files = {path.name: path for path in current.iterdir() if path.is_file() and path.suffix != ".ini"}
    reference_files = {path.name: path for path in reference.iterdir() if path.is_file() and path.suffix != ".ini"}
    common_names = sorted(set(current_files).intersection(reference_files))

    total_points = 0
    total_sq = 0.0
    stats = []

    for name in common_names:
        columns, current_rows = _read_numeric_table(current_files[name])
        reference_columns, reference_rows = _read_numeric_table(reference_files[name])
        if columns != reference_columns:
            min_len = min(len(columns), len(reference_columns))
            columns = columns[:min_len]
        if not current_rows or not reference_rows:
            continue

        file_points = 0
        file_sq = 0.0
        worst_rel = 0.0
        worst_abs = 0.0
        row_count = min(len(current_rows), len(reference_rows))
        for left_row, right_row in zip(current_rows[:row_count], reference_rows[:row_count]):
            col_count = min(len(left_row), len(right_row))
            if col_count == 0:
                continue
            active_columns = columns[:col_count] if columns else [str(index) for index in range(col_count)]
            old_row_map = {column: value for column, value in zip(active_columns, right_row[:col_count])}
            for column, left_value, right_value in zip(active_columns, left_row[:col_count], right_row[:col_count]):
                rel = _normalized_error(column, right_value, left_value, old_row_map)
                if rel is None:
                    continue
                file_sq += rel * rel
                file_points += 1
                worst_rel = max(worst_rel, rel)
                worst_abs = max(worst_abs, abs(left_value - right_value))

        if not file_points:
            continue

        file_rms = math.sqrt(file_sq / file_points)
        stats.append((name, file_rms, worst_rel, worst_abs, file_points))
        total_points += file_points
        total_sq += file_sq

    overall_rms = math.sqrt(total_sq / total_points) if total_points else 0.0

    print("\nDifference summary (all common numeric files):")
    print(f"  files analyzed: {len(stats)}")
    print(f"  points analyzed: {total_points}")
    print(f"  overall RMS(relative): {overall_rms:.6e}")

    if not stats:
        print("  no numeric comparable files found")
        return

    by_rms = sorted(stats, key=lambda item: item[1], reverse=True)[:10]
    by_worst = sorted(stats, key=lambda item: item[2], reverse=True)[:10]

    print("\nTop 10 files by RMS(relative):")
    for name, file_rms, worst_rel, worst_abs, points in by_rms:
        print(f"  {name}: rms={file_rms:.6e}, worst_rel={worst_rel:.6e}, worst_abs={worst_abs:.6e}, n={points}")

    print("\nTop 10 files by worst relative error:")
    for name, file_rms, worst_rel, worst_abs, points in by_worst:
        print(f"  {name}: worst_rel={worst_rel:.6e}, rms={file_rms:.6e}, worst_abs={worst_abs:.6e}, n={points}")


def generate_reference(parsed_args):
    parsed_args.reference_dir = os.path.abspath(parsed_args.reference_dir)
    reference_ini_dir = os.path.join(parsed_args.reference_dir, parsed_args.reference_ini_dir)
    reference_out_dir = os.path.join(parsed_args.reference_dir, parsed_args.reference_out_dir)

    if parsed_args.clean:
        for target in (reference_ini_dir, reference_out_dir):
            if os.path.exists(target):
                remove_tree(target)

    reference_inis = make_reference_inis(parsed_args, reference_ini_dir, reference_out_dir)
    run_reference_inis(reference_inis, parsed_args)

    print(f"High-accuracy ini directory: {reference_ini_dir}")
    print(f"High-accuracy output directory: {reference_out_dir}")
    return 0


def compare_to_reference(parsed_args):
    parsed_args.work_dir = os.path.abspath(parsed_args.work_dir)
    parsed_args.base_settings = os.path.abspath(parsed_args.base_settings)
    reference_outputs = os.path.join(os.path.abspath(parsed_args.reference_dir), parsed_args.reference_out_dir)
    use_custom_flow = parsed_args.skip_non_flat or parsed_args.only_non_flat

    if use_custom_flow:
        out_dir = os.path.join(parsed_args.work_dir, parsed_args.out_files_dir)
        if parsed_args.clean and os.path.exists(out_dir):
            remove_tree(out_dir)
        os.makedirs(out_dir, exist_ok=True)
        inis = make_default_inis(parsed_args)
        run_status = 1 if run_inis(inis, parsed_args) else 0
        diff_args = ["--diff_to", reference_outputs, "--diff_tolerance", str(parsed_args.diff_tolerance)]
        extend_common_compare_args(diff_args, parsed_args, include_clean=False)
        diff_status = run_camb_test_files(diff_args)
    else:
        run_args = ["--make_ini"]
        extend_common_compare_args(run_args, parsed_args, include_clean=True)
        run_status = run_camb_test_files(run_args)

        diff_args = ["--diff_to", reference_outputs, "--diff_tolerance", str(parsed_args.diff_tolerance)]
        extend_common_compare_args(diff_args, parsed_args, include_clean=False)
        diff_status = run_camb_test_files(diff_args)

    summarize_differences(
        os.path.join(parsed_args.work_dir, parsed_args.out_files_dir),
        reference_outputs,
    )
    return run_status or diff_status


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.work_dir = os.path.abspath(args.work_dir)
    args.base_settings = os.path.abspath(args.base_settings)

    if args.mode == "generate_reference":
        return generate_reference(args)
    return compare_to_reference(args)


if __name__ == "__main__":
    raise SystemExit(main())
