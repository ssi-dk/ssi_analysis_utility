#!/usr/bin/env python3
"""
Extract versions of the primary executable in each Snakemake conda environment.

Assumptions
-----------
1. Snakemake stores environments under: conda/<hash>/
2. Each environment contains:
      - conda-meta/history
      - bin/ (Unix/macOS)
3. The "main tool" is approximated as:
      - the first explicitly requested package in conda-meta/history
      - excluding common infrastructure packages (python, pip, etc.)
4. Version is obtained by trying:
      - <tool> --version
      - <tool> -V
      - <tool> version

Usage
-----
    python get_env_versions.py conda/

Output
------
TSV with columns:
    environment    tool    version
"""

from pathlib import Path
import subprocess
import sys
import re

# Packages that are usually not the main tool
IGNORE = {
    "python", "pip", "setuptools", "wheel",
    "openssl", "ca-certificates", "certifi",
    "libgcc-ng", "libstdcxx-ng", "zlib",
    "xz", "bzip2", "sqlite", "tk", "ncurses",
    "readline", "tzdata", "curl", "wget",
    "mamba", "conda"
}


def parse_requested_packages(history_file: Path):
    """
    Parse explicitly requested packages from conda-meta/history.

    Looks for lines like:
        # cmd: ... create ... bwa=0.7.17 samtools=1.20
    """
    requested = []

    with history_file.open() as fh:
        for line in fh:
            if line.startswith("# cmd:"):
                # Extract package-like tokens
                tokens = line.strip().split()
                for token in tokens:
                    if token.startswith("-"):
                        continue
                    if "/" in token:
                        continue
                    if "=" in token:
                        pkg = token.split("=")[0]
                    else:
                        pkg = token
                    if re.fullmatch(r"[A-Za-z0-9_.+-]+", pkg):
                        requested.append(pkg)

    # Preserve order and remove duplicates
    seen = set()
    ordered = []
    for pkg in requested:
        if pkg not in seen:
            seen.add(pkg)
            ordered.append(pkg)

    return ordered


def choose_main_tool(packages):
    """Return the first non-ignored package."""
    for pkg in packages:
        if pkg.lower() not in IGNORE:
            return pkg
    return None


def get_version(env_path: Path, tool: str):
    """
    Run the tool inside the environment and capture version text.
    """
    exe = env_path / "bin" / tool
    if not exe.exists():
        return "executable not found"

    commands = [
        [str(exe), "--version"],
        [str(exe), "-V"],
        [str(exe), "version"],
    ]

    for cmd in commands:
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=15,
            )
            output = (result.stdout + result.stderr).strip()
            if output:
                # Return first non-empty line
                return output.splitlines()[0]
        except Exception:
            pass

    return "version not detected"


def main(conda_dir):
    conda_dir = Path(conda_dir)

    print("environment\ttool\tversion")

    for env_path in sorted(p for p in conda_dir.iterdir() if p.is_dir()):
        history = env_path / "conda-meta" / "history"
        if not history.exists():
            continue

        packages = parse_requested_packages(history)
        tool = choose_main_tool(packages)

        if tool is None:
            print(f"{env_path.name}\tNA\tno explicit tool found")
            continue

        version = get_version(env_path, tool)
        print(f"{env_path.name}\t{tool}\t{version}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python get_env_versions.py conda/")
    main(sys.argv[1])
