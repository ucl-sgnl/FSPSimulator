# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('..'))       # root directory
sys.path.insert(0, os.path.abspath('../src'))   # src directory

project = 'FSPSim'
copyright = '2023, Indigo Brownhall, Charles Constant, Santosh Bhattarai, Lee Devlin, Marek Ziebart'
author = 'Indigo Brownhall, \nCharles Constant, \nSantosh Bhattarai, \nLee Devlin, \nMarek Ziebart'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
autodoc_mock_imports = ["wget", "poliastro", "sgp4", "dotenv", "jplephem", "pyatmos"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

# -- Options for LaTeX output ------------------------------------------------
latex_engine = 'xelatex'  # xelatex is modern and supports more features
master_doc = 'index'
latex_documents = [
    (master_doc, 'FSPSim.tex', 'FSPSim Documentation',
     'Indigo Brownhall, \nCharles Constant, \nSantosh Bhattarai, \nLee Devlin, \nMarek Ziebart', 'manual'),
]