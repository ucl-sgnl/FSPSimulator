# FSPSimulator
This is the updated Python Version of the ODL FSP, written in python.

# Environment Installation:

''' bash
conda update -n base -c defaults conda
conda env create -f settings/fspsim_env.yml
conda activate fspsim/ source activate fspsim
'''

if you get module import errors, navigate to project root and try:  
''' bash
export PYTHONPATH="$(pwd)"    
'''