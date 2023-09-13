# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('..'))       # root directory
sys.path.insert(0, os.path.abspath('../../src'))   # src directory

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
exclude_patterns = [ 'Thumbs.db', '.DS_Store']
autodoc_mock_imports = ["wget", "poliastro", "sgp4", "dotenv", "jplephem", "pyatmos"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
import sphinx_rtd_theme

html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_theme = "sphinx_rtd_theme"

html_static_path = ['_static']
