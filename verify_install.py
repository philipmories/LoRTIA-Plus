#!/usr/bin/env python3
"""Verify a LoRTIA Plus v1.0.0 installation.

This script intentionally depends only on the Python standard library so that
it can report missing LoRTIA Plus dependencies instead of failing on import.
"""

from __future__ import annotations

import importlib
import importlib.metadata
import os
import re
import shutil
import site
import subprocess
import sys
from dataclasses import dataclass


REFERENCE_PYTHON = (3, 10, 14)
SUPPORTED_PYTHON = (3, 10)

PYTHON_DEPENDENCIES = (
    ("NumPy", "numpy", "numpy", "1.26.4"),
    ("pandas", "pandas", "pandas", "2.2.3"),
    ("Polars", "polars", "polars", "1.8.2"),
    ("pysam", "pysam", "pysam", "0.22.1"),
    ("SciPy", "scipy", "scipy", "1.14.1"),
    ("Biopython", "biopython", "Bio", "1.88"),
)

REFERENCE_SAMTOOLS = "1.20"
REFERENCE_GNU_AWK = "5.3.1"


@dataclass
class CheckResult:
    label: str
    status: str
    detail: str


def run_command(command: list[str]) -> tuple[int, str]:
    try:
        proc = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    except OSError as exc:
        return 127, str(exc)
    return proc.returncode, proc.stdout.strip()


def check_python() -> CheckResult:
    current = sys.version_info[:3]
    if current[:2] != SUPPORTED_PYTHON:
        return CheckResult(
            "Python",
            "FAIL",
            f"{current[0]}.{current[1]}.{current[2]} installed; "
            f"LoRTIA Plus v1.0.0 supports Python 3.10.x "
            f"(reference environment: 3.10.14).",
        )
    if current != REFERENCE_PYTHON:
        return CheckResult(
            "Python",
            "WARN",
            f"{current[0]}.{current[1]}.{current[2]} installed; "
            "compatible Python 3.10.x detected, but the reference "
            "environment is pinned to 3.10.14.",
        )
    return CheckResult("Python", "PASS", "3.10.14")


def check_python_dependency(
    label: str, distribution: str, module: str, expected: str
) -> CheckResult:
    try:
        installed = importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return CheckResult(label, "FAIL", "not installed")

    try:
        importlib.import_module(module)
    except Exception as exc:
        return CheckResult(
            label,
            "FAIL",
            f"version {installed} is installed but import failed: "
            f"{type(exc).__name__}: {exc}",
        )

    if installed != expected:
        return CheckResult(
            label,
            "FAIL",
            f"{installed} installed; pinned version is {expected}",
        )
    return CheckResult(label, "PASS", installed)


def check_samtools() -> CheckResult:
    exe = shutil.which("samtools")
    if exe is None:
        return CheckResult("samtools", "FAIL", "command not found in PATH")

    rc, output = run_command([exe, "--version"])
    if rc != 0:
        return CheckResult(
            "samtools", "FAIL", f"`samtools --version` failed: {output}"
        )

    first_line = output.splitlines()[0] if output else ""
    match = re.search(r"\bsamtools\s+([0-9][^\s]*)", first_line, flags=re.I)
    installed = match.group(1) if match else first_line or "unknown"

    if installed != REFERENCE_SAMTOOLS:
        return CheckResult(
            "samtools",
            "WARN",
            f"{installed} detected; reference environment is pinned to "
            f"{REFERENCE_SAMTOOLS}",
        )
    return CheckResult("samtools", "PASS", installed)


def check_gnu_awk() -> CheckResult:
    exe = shutil.which("awk")
    if exe is None:
        return CheckResult("GNU awk", "FAIL", "`awk` command not found in PATH")

    rc, output = run_command([exe, "--version"])
    if rc != 0:
        return CheckResult(
            "GNU awk",
            "FAIL",
            "`awk --version` failed; a GNU awk installation is required "
            f"for the pinned reference environment. Output: {output}",
        )

    first_line = output.splitlines()[0] if output else ""
    if "GNU Awk" not in first_line and "GNU awk" not in first_line:
        return CheckResult(
            "GNU awk",
            "FAIL",
            f"non-GNU awk detected: {first_line or 'unknown implementation'}",
        )

    match = re.search(r"GNU [Aa]wk\s+([0-9][0-9.]*)", first_line)
    installed = match.group(1) if match else first_line

    if installed != REFERENCE_GNU_AWK:
        return CheckResult(
            "GNU awk",
            "WARN",
            f"{installed} detected; reference environment is pinned to "
            f"{REFERENCE_GNU_AWK}",
        )
    return CheckResult("GNU awk", "PASS", installed)


def check_user_site_isolation() -> CheckResult:
    env_value = os.environ.get("PYTHONNOUSERSITE")
    disabled = not bool(getattr(site, "ENABLE_USER_SITE", False))

    if env_value == "1" and disabled:
        return CheckResult(
            "Python user-site isolation",
            "PASS",
            "PYTHONNOUSERSITE=1; user-site packages disabled",
        )
    if disabled:
        return CheckResult(
            "Python user-site isolation",
            "WARN",
            "user-site packages are disabled, but PYTHONNOUSERSITE is not "
            "explicitly set to 1",
        )
    return CheckResult(
        "Python user-site isolation",
        "WARN",
        "user-site packages are enabled; the reference Conda environment "
        "sets PYTHONNOUSERSITE=1",
    )


def check_lortia_command() -> CheckResult:
    exe = shutil.which("LoRTIA")
    if exe is None:
        return CheckResult(
            "LoRTIA command",
            "FAIL",
            "command not found in PATH; run `python -m pip install .` "
            "inside the LoRTIA Plus environment",
        )

    rc, output = run_command([exe, "-h"])
    if rc != 0:
        return CheckResult(
            "LoRTIA command",
            "FAIL",
            f"`LoRTIA -h` returned exit code {rc}: {output[:500]}",
        )

    if "LoRTIA" not in output and "Long-read RNA-Seq" not in output:
        return CheckResult(
            "LoRTIA command",
            "WARN",
            "`LoRTIA -h` executed successfully, but the expected help "
            "identifier was not found",
        )

    return CheckResult("LoRTIA command", "PASS", exe)


def main() -> int:
    checks: list[CheckResult] = [check_python()]

    for args in PYTHON_DEPENDENCIES:
        checks.append(check_python_dependency(*args))

    checks.extend(
        [
            check_samtools(),
            check_gnu_awk(),
            check_user_site_isolation(),
            check_lortia_command(),
        ]
    )

    print("LoRTIA Plus v1.0.0 installation verification")
    print("=" * 48)
    for result in checks:
        print(f"[{result.status:4}] {result.label}: {result.detail}")

    failures = [r for r in checks if r.status == "FAIL"]
    warnings = [r for r in checks if r.status == "WARN"]

    print("-" * 48)
    if failures:
        print(
            f"LoRTIA Plus installation check FAILED "
            f"({len(failures)} failure(s), {len(warnings)} warning(s))."
        )
        return 1

    if warnings:
        print(
            f"LoRTIA Plus installation check passed with "
            f"{len(warnings)} warning(s)."
        )
    else:
        print("LoRTIA Plus installation check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
