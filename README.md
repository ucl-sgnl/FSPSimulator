# FSPSimulator
This is the updated Python Version of the ODL FSP, written in python.

# Environment Installation:

``` bash
conda update -n base -c defaults conda
conda env create -f settings/fspsim_env.yml
conda activate fspsim/ source activate fspsim
```

# Installation Error Fixes:
if you get module import errors, navigate to project root and try:  
``` bash
export PYTHONPATH="$(pwd)"    
```

if you get pyatmos installation error on Mac import errors, navigate to project root and try:  
``` bash
brew install fftw  
```