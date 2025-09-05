# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'DEPHY-SCM'
copyright = '2025, DEPHY'
author = 'Roehrig et al.'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

import os
import sys
sys.path.insert(0, os.path.abspath(".."))  # pour trouver axis.py

extensions = [
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
#    "sphinx_copybutton",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
#    "sphinx_design",
#    "sphinxcontrib.bibtex",
    "matplotlib.sphinxext.plot_directive",
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for napoleon ----------------------------------------------------
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# -- Options for autodoc ----------------------------------------------------
autodoc_typehints = "description"
autodoc_typehints_description_target = "all"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
add_module_names = False
autosummary_generate = False

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = "sphinx_rtd_theme"
#html_theme = 'alabaster'
#html_theme = 'python_docs_theme'
html_theme = 'pydata_sphinx_theme'

html_theme_options = {
    "navbar_end": ["version-switcher", "navbar-icon-links", "theme-switcher"],
    "icon_links": [],
    "logo": {
        "text": "DEPHY-SCM",
    },
}

html_static_path = ['_static']
html_css_files = [
    "custom.css",
]


# -- Options for matplotlib plots
plot_formats = [("png", 90), "pdf"]

