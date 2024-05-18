import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))

# Project information
project = 'Tools Kit'
author = 'Sebastien Grunberg, Méloé Enzinger'
release = '0.1.0'

# General configuration
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]
templates_path = ['_templates']
exclude_patterns = []

# Options for HTML output
html_theme = 'alabaster'
html_static_path = ['_static']
