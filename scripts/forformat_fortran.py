"""Run forformat with the repository's VS Code Fortran formatting settings."""

import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_DIR = SCRIPT_DIR.parent

FORFORMAT_ARGS = (
    "--full",
    "--indent=4",
    "--indent-module=0",
    "--indent-procedure=0",
    "--start-indent=4",
    "--indent-contains=0",
    "--openmp=0",
    "--indent-contains=restart",
    "--indent-select=4",
    "--indent-case=4",
    "--indent-interface=0",
    "--indent-continuation=4",
    "--indent-ampersand",
)


def hook_path(name: str) -> Path:
    path = Path(name)
    if not path.is_absolute() and not path.exists():
        candidate = REPOSITORY_DIR / path
        if candidate.exists():
            path = candidate
    return path


def format_paths(paths: list[Path], forformat: str) -> None:
    result = subprocess.run(
        [forformat, *FORFORMAT_ARGS, *(str(path) for path in paths)],
        capture_output=True,
        check=False,
        cwd=REPOSITORY_DIR,
    )
    if result.returncode:
        error = result.stderr.decode(errors="replace").strip()
        message = "forformat failed"
        if error:
            message += f": {error}"
        raise RuntimeError(message)


def main() -> int:
    paths = [hook_path(name) for name in sys.argv[1:]]
    if not paths:
        return 0

    forformat = shutil.which("forformat")
    if forformat is None:
        print("forformat is required to format Fortran files", file=sys.stderr)
        return 1
    try:
        format_paths(paths, forformat)
    except (OSError, RuntimeError) as error:
        print(error, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
