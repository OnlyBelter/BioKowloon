# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
# import sphinx_book_theme
sys.path.insert(0, os.path.abspath(r'..'))
# sys.path.append(os.path.abspath(r'..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'BioKowloon'
copyright = '2024, OnlyBelter'
author = 'OnlyBelter'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
"myst_parser",
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    # 'sphinx.ext.doctest',
    # 'sphinx.ext.coverage',
    # 'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    # 'sphinx.ext.autosummary',
    # 'sphinx.ext.intersphinx',
    # 'matplotlib.sphinxext.plot_directive',
    # 'gallery_generator',
    # 'numpydoc',
    # 'sphinx_issues',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_static_path = ['_static']
html_theme = 'furo'
html_logo = "_static/logo.png"
html_title = "BioKowloon"
html_copy_source = True
# html_sourcelink_suffix = ""
html_favicon = "_static/logo.png"
# html_last_updated_fmt = ""
highlight_language = "python"

# ref to https://github.com/mwaskom/seaborn/blob/master/doc/conf.py
html_theme_options = {
    "source_directory": "docs/",
    "source_repository": "https://github.com/OnlyBelter/BioKowloon",
    # "source_edit_link": True,
    # "use_issues_button": True,
    # "use_repository_button": True,
    # "use_download_button": True,
}

# custom css
html_css_files = [
    '_static/custom.css',
]

# Add type of source files
source_suffix = ['.rst', '.md']

fontawesome_included = True