"""Run findent with the repository's VS Code Fortran formatting settings."""

import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

FINDENT_ARGS = (
    "--indent=4",
    "--indent_module=0",
    "--indent_procedure=0",
    "--start_indent=4",
    "--indent_contains=0",
    "--openmp=0",
    "--indent_contains=restart",
    "--indent_select=4",
    "--indent_case=4",
    "--indent_interface=0",
    "--indent_continuation=4",
    "--indent_ampersand",
)


def hook_path(name: str) -> Path:
    path = Path(name)
    if not path.is_absolute() and not path.exists():
        candidate = SCRIPT_DIR.parent / path
        if candidate.exists():
            path = candidate
    return path


def format_path(path: Path, findent: str) -> None:
    source = path.read_bytes()
    result = subprocess.run(
        [findent, *FINDENT_ARGS],
        input=source,
        capture_output=True,
        check=False,
    )
    if result.returncode:
        error = result.stderr.decode(errors="replace").strip()
        message = f"findent failed for {path}"
        if error:
            message += f": {error}"
        raise RuntimeError(message)
    if result.stdout != source:
        path.write_bytes(result.stdout)


def main() -> int:
    findent = shutil.which("findent")
    if findent is None:
        print("findent is required to format Fortran files", file=sys.stderr)
        return 1
    try:
        for name in sys.argv[1:]:
            format_path(hook_path(name), findent)
    except (OSError, RuntimeError) as error:
        print(error, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
