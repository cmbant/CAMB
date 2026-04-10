import argparse
import os
import shutil
import stat
import subprocess
import sys

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
FORTRAN_DIR = os.path.abspath(os.path.join(TESTS_DIR, ".."))
REPO_ROOT = os.path.abspath(os.path.join(FORTRAN_DIR, ".."))
DEFAULT_WORK_DIR = os.path.join(FORTRAN_DIR, "testfiles")
DEFAULT_BASE_SETTINGS = os.path.join(REPO_ROOT, "inifiles", "params.ini")
DEFAULT_REFERENCE_REPO = "https://github.com/cmbant/CAMB_test_outputs.git"

if TESTS_DIR not in sys.path:
    sys.path.insert(0, TESTS_DIR)


def build_parser():
    parser = argparse.ArgumentParser(description="Run the CAMB text-output regression tests from Python")
    parser.add_argument("work_dir", nargs="?", default=DEFAULT_WORK_DIR, help="working directory for generated inis")
    parser.add_argument(
        "--reference_dir",
        "--reference-dir",
        help="reference CAMB_test_outputs directory, or its test_outputs subdirectory",
    )
    parser.add_argument("--skip_diff", "--skip-diff", action="store_true", help="skip comparing against references")
    parser.add_argument("--keep_reference", "--keep-reference", action="store_true", help="keep a cloned reference")
    parser.add_argument("--clean", action="store_true", help="delete the output directory before running")
    parser.add_argument("--runner", choices=("module", "command"), default="module", help="test execution backend")
    parser.add_argument("--prog", default="./camb", help="legacy executable path when --runner=command")
    parser.add_argument("--no_validate", "--no-validate", action="store_true", help="skip ini validation")
    parser.add_argument("--base_settings", "--base-settings", default=DEFAULT_BASE_SETTINGS, help="base params.ini")
    parser.add_argument("--out_files_dir", "--out-files-dir", default="test_outputs", help="output files directory")
    parser.add_argument("--max_tests", "--max-tests", type=int, help="maximum number of tests to run")
    parser.add_argument("--no_sources", "--no-sources", action="store_true", help="turn off CAMB sources tests")
    parser.add_argument("--no_de", "--no-de", action="store_true", help="skip dark energy tests")
    parser.add_argument(
        "--diff_tolerance",
        "--diff-tolerance",
        type=float,
        default=1e-4,
        help="fallback numerical tolerance for diffing",
    )
    parser.add_argument(
        "--verbose_diff_output",
        "--verbose-diff-output",
        "--verbose",
        action="store_true",
        help="print more diff details",
    )
    parser.add_argument("--num_diff", "--num-diff", action="store_true", help="use absolute diffs")
    parser.add_argument(
        "--override",
        "--set",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="override an ini parameter for every run; can be repeated",
    )
    parser.add_argument(
        "--reference_repo",
        "--reference-repo",
        default=DEFAULT_REFERENCE_REPO,
        help="repository to clone when reference outputs are needed",
    )
    return parser


def extend_common_args(cmd_args, parsed_args):
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
    if parsed_args.clean:
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


def resolve_reference_outputs(reference_dir):
    reference_dir = os.path.abspath(reference_dir)
    test_outputs_dir = os.path.join(reference_dir, "test_outputs")
    if os.path.isdir(test_outputs_dir):
        return test_outputs_dir
    return reference_dir


def run_camb_test_files(args):
    from CAMB_test_files import main as camb_test_files_main

    return camb_test_files_main(args)


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
                import time

                time.sleep(0.2)
        raise exc_info[1]

    shutil.rmtree(path, onerror=onerror)


def clone_reference_outputs(parsed_args):
    reference_root = os.path.join(parsed_args.work_dir, "CAMB_test_outputs")
    if os.path.exists(reference_root):
        remove_tree(reference_root)
    subprocess.run(["git", "clone", "--depth=1", parsed_args.reference_repo, reference_root], check=True)
    return reference_root


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.work_dir = os.path.abspath(args.work_dir)
    args.base_settings = os.path.abspath(args.base_settings)

    run_args = ["--make_ini"]
    extend_common_args(run_args, args)
    result = run_camb_test_files(run_args)
    if result or args.skip_diff:
        return result

    cloned_reference = None
    try:
        if args.reference_dir:
            reference_outputs = resolve_reference_outputs(args.reference_dir)
        else:
            cloned_reference = clone_reference_outputs(args)
            reference_outputs = resolve_reference_outputs(cloned_reference)

        diff_args = ["--diff_to", reference_outputs, "--diff_tolerance", str(args.diff_tolerance)]
        extend_common_args(diff_args, args)
        return run_camb_test_files(diff_args)
    finally:
        if cloned_reference and not args.keep_reference and os.path.exists(cloned_reference):
            remove_tree(cloned_reference)


if __name__ == "__main__":
    raise SystemExit(main())
