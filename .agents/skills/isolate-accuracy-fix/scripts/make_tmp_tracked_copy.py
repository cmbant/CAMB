"""Create a disposable CAMB worktree with tracked local changes replayed."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Mapping
from pathlib import Path


def run(
    args: list[str],
    *,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
    input_bytes: bytes | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(args, cwd=cwd, env=env, input=input_bytes, capture_output=True, check=check)


def git(repo: Path, args: list[str], **kwargs: object) -> subprocess.CompletedProcess[bytes]:
    return run(["git", "-C", str(repo), *args], **kwargs)


def text(output: bytes) -> str:
    return output.decode("utf-8", errors="replace").strip()


def repo_root(start: Path) -> Path:
    result = run(["git", "-C", str(start), "rev-parse", "--show-toplevel"])
    return Path(text(result.stdout))


def submodule_paths(repo: Path) -> list[str]:
    result = git(
        repo,
        [
            "submodule",
            "foreach",
            "--recursive",
            "--quiet",
            "printf '%s\\n' \"$displaypath\"",
        ],
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(text(result.stderr) or "failed to inspect submodules")
    return [line for line in text(result.stdout).splitlines() if line]


def diff_bytes(repo: Path, args: list[str]) -> bytes:
    result = git(repo, ["diff", "--binary", "--full-index", *args])
    return result.stdout


def apply_patch(repo: Path, patch: bytes, label: str) -> bool:
    if not patch:
        return False
    try:
        git(repo, ["apply", "--binary", "--whitespace=nowarn", "-"], input_bytes=patch)
    except subprocess.CalledProcessError as exc:
        message = text(exc.stderr) or f"failed to apply {label} patch"
        raise RuntimeError(message) from exc
    return True


def output_tail(output: bytes, max_lines: int = 20) -> list[str]:
    lines = text(output).splitlines()
    return lines[-max_lines:]


def build_camb(repo: Path, python: str) -> dict[str, object]:
    command = [python, "setup.py", "make"]
    forutils_path = repo / "forutils"
    env = os.environ.copy()
    env["FORUTILSPATH"] = str(forutils_path)
    result = run(command, cwd=repo, env=env, check=False)
    if result.returncode != 0:
        print(text(result.stdout), file=sys.stderr)
        print(text(result.stderr), file=sys.stderr)
        raise RuntimeError(f"build failed: {' '.join(command)}")

    library_candidates = [repo / "camb" / "camblib.so", repo / "camb" / "cambdll.dll"]
    built_library = next((path for path in library_candidates if path.exists()), None)
    return {
        "command": command,
        "forutils_path": str(forutils_path),
        "returncode": result.returncode,
        "built_library": str(built_library) if built_library else None,
        "stdout_tail": output_tail(result.stdout),
        "stderr_tail": output_tail(result.stderr),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a detached /tmp git worktree at HEAD, initialize submodules, "
            "replay staged and unstaged tracked changes from the source repo, "
            "then build CAMB in the copy. "
            "Untracked files are intentionally ignored."
        )
    )
    parser.add_argument(
        "--repo",
        type=Path,
        default=Path.cwd(),
        help="Path inside the source repository. Defaults to the current directory.",
    )
    parser.add_argument(
        "--tmp-parent",
        type=Path,
        default=Path(tempfile.gettempdir()),
        help="Parent directory for the disposable copy. Defaults to the system temp directory.",
    )
    parser.add_argument(
        "--prefix",
        default="camb-accuracy-",
        help="Name prefix for the disposable copy directory.",
    )
    parser.add_argument(
        "--ignore-dirty-submodules",
        action="store_true",
        help="Do not replay tracked staged/unstaged changes from submodules.",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python executable to use for the CAMB build. Defaults to the interpreter running this helper.",
    )
    parser.add_argument(
        "--no-build",
        action="store_true",
        help="Create the tracked source copy but do not run the CAMB build.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source = repo_root(args.repo.resolve())
    head = text(git(source, ["rev-parse", "HEAD"]).stdout)

    staged_patch = diff_bytes(source, ["--cached"])
    unstaged_patch = diff_bytes(source, [])
    submodule_patches: dict[str, tuple[bytes, bytes]] = {}
    if not args.ignore_dirty_submodules:
        for path in submodule_paths(source):
            submodule = source / path
            submodule_patches[path] = (diff_bytes(submodule, ["--cached"]), diff_bytes(submodule, []))
    destination = Path(tempfile.mkdtemp(prefix=args.prefix, dir=args.tmp_parent))

    try:
        shutil.rmtree(destination)
        run(["git", "-C", str(source), "worktree", "add", "--detach", str(destination), "HEAD"])
        git(destination, ["submodule", "update", "--init", "--recursive"])
        staged_applied = apply_patch(destination, staged_patch, "staged")
        unstaged_applied = apply_patch(destination, unstaged_patch, "unstaged")
        submodule_summary = {}
        for path, (sub_staged_patch, sub_unstaged_patch) in submodule_patches.items():
            submodule = destination / path
            submodule_summary[path] = {
                "staged_patch_applied": apply_patch(submodule, sub_staged_patch, f"{path} staged"),
                "unstaged_patch_applied": apply_patch(submodule, sub_unstaged_patch, f"{path} unstaged"),
            }
        build_summary = None if args.no_build else build_camb(destination, args.python)
    except Exception:
        shutil.rmtree(destination, ignore_errors=True)
        raise

    summary = {
        "build": build_summary,
        "copy_path": str(destination),
        "source_repo": str(source),
        "source_head": head,
        "staged_patch_applied": staged_applied,
        "submodule_patches": submodule_summary,
        "unstaged_patch_applied": unstaged_applied,
        "untracked_files": "ignored",
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
