"""
setup.py — builds libpapa2.so via make before packaging.

The shared library is loaded at runtime by papa2/_cdada.py (and friends) via
ctypes, so it must be present inside the installed papa2/ package directory.
We cannot use ext_modules here because the build is driven by the project
Makefile (which handles compiler flags, OpenMP, zlib, etc.).  Instead we:

  1. Run ``make libpapa2.so`` in the build_ext step so the .so is always
     freshly compiled before any wheel/sdist is assembled.
  2. Copy the resulting libpapa2.so into the papa2/ sub-package directory so
     that setuptools picks it up as package_data and installs it alongside the
     Python source files.
"""

import os
import shutil
import subprocess
from pathlib import Path

from setuptools import setup
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.build_py import build_py as _build_py

HERE = Path(__file__).parent.resolve()
SO_NAME = "libpapa2.so"
SO_SRC = HERE / SO_NAME            # built into project root by Makefile
SO_DST = HERE / "papa2" / SO_NAME  # where it must live for package_data


def _build_shared_library():
    """Run the project Makefile to (re)compile libpapa2.so."""
    print(f"[setup.py] Running: make {SO_NAME}")
    result = subprocess.run(
        ["make", SO_NAME],
        cwd=str(HERE),
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"make {SO_NAME} failed with exit code {result.returncode}. "
            "Ensure gcc/g++ and zlib-dev are installed."
        )
    if not SO_SRC.exists():
        raise FileNotFoundError(
            f"make succeeded but {SO_SRC} was not found. Check the Makefile."
        )
    # Place the .so inside the Python package so it is found at import time.
    shutil.copy2(str(SO_SRC), str(SO_DST))
    print(f"[setup.py] Copied {SO_SRC} -> {SO_DST}")


class BuildExtWithMake(_build_ext):
    """Extend build_ext to compile the shared library via make."""

    def run(self):
        _build_shared_library()
        super().run()


class BuildPyWithMake(_build_py):
    """Extend build_py so ``pip install`` (which skips build_ext) also works."""

    def run(self):
        _build_shared_library()
        super().run()


setup(
    # Core metadata lives in pyproject.toml; only packaging behaviour here.
    package_data={
        "papa2": [SO_NAME],
    },
    cmdclass={
        "build_ext": BuildExtWithMake,
        "build_py": BuildPyWithMake,
    },
)
