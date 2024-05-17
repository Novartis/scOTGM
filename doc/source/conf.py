# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys

import sphinx_rtd_theme

sys.path.insert(0, os.path.abspath("../../src"))

project = 'sc-OTGM'
copyright = '2023 Novartis Institutes for BioMedical Research, INC'
author = 'Andac Demir'
release = '2023'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
        "sphinx.ext.autodoc",
        "sphinx.ext.autosummary",
        "sphinx.ext.intersphinx",
        "sphinx.ext.mathjax",
        "sphinx.ext.napoleon",
        "sphinx.ext.viewcode",
        'nbsphinx',
        ]

templates_path = ['_templates']
exclude_patterns = []

language = "English"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_static_path = ['_static']
html_logo = "logo.png"
html_css_files = [
        "custom.css"
        ]
html_show_sourcelink = False

autodoc_default_options = {
      "members": True,  # This option indicates that members should be included in the documentation.
      "undoc-members": True,  # This option indicates that members without docstrings should be included in the documentation.
      "private-members": False,  # This option indicates that private members (those starting with _) should be included in the documentation.
      "special-members": "__call__",  # This option indicates that special members (those starting and ending with __) should be included in the documentation.
      "show-inheritance": True,  # This option indicates that class inheritance should be shown in the documentation.
      "member-order": "bysource",
}
