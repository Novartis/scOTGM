# -------------------------------------------------------
# Copyright 2023 Novartis Institutes for BioMedical Research, INC
# Author: Andac Demir
# -------------------------------------------------------
from setuptools import find_packages, setup

from src import __version__

with open("requirements.txt", "r") as f:
    required = f.readlines()


setup(
    name="src", authors="Andac Demir", version=__version__, install_requires=required, packages=find_packages(),
)
