# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from importlib.metadata import version as _pkg_version
sys.path.insert(0, os.path.abspath('../../'))

# -- Project information -----------------------------------------------------

project = 'htpolynet'
copyright = '2026, Cameron Abrams, Ming Huang'
author = 'Cameron Abrams, Ming Huang'

# The full version, including alpha/beta/rc tags
release = _pkg_version('htpolynet')
version = release


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.coverage', 'sphinx.ext.napoleon',
    'sphinxcontrib.bibtex', 'sphinxcontrib.mermaid', 'sphinx_copybutton', 'myst_parser',]
# extensions.append('sphinx.ext.todo')
# extensions.append('sphinx.ext.intersphinx')
# extensions.append('sphinx.ext.mathjax')
extensions.append('sphinx.ext.viewcode')

# extensions.append('sphinx.ext.graphviz')
bibtex_bibfiles = ['references.bib']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['*.xmind', 'build', 'Thumbs.db', '.DS_Store']
language = 'en'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'furo'

# The logo is line art, so it ships as two transparent variants rather than one
# opaque image: black ink for the light theme, white ink for the dark one.  An
# opaque logo would show as a white square in furo's dark mode.
html_theme_options = {
    'light_logo': 'htpolynet-logo-light.png',
    'dark_logo': 'htpolynet-logo-dark.png',
}

# The favicon, unlike the sidebar logo, keeps an opaque white disc: a browser's
# tab bar may be light or dark, and transparent line art would vanish on one of
# them.  16, 32 and 48 px in one .ico.
html_favicon = '_static/htpolynet-favicon.ico'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

numfig = True
autosummary_generate = True
todo_include_todos = True