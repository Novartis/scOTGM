# -------------------------------------------------------
# Copyright 2023 Novartis Institutes for BioMedical Research, INC
# Author: Andac Demir
# -------------------------------------------------------
from setuptools import find_packages, setup

from src import __version__

with open("requirements.txt", "r") as f:
    required = f.read().splitlines()

setup(
    name="src",
    version=__version__,
    author="Andac Demir",
    description="sc-OTGM: Single-Cell Perturbation Modeling by Solving Optimal Mass Transport "
                "on the Manifold of Gaussian Mixtures",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    install_requires=[req.strip() for req in required],
    packages=find_packages(),
    python_requires='>=3.6',
    include_package_data=True,
)
