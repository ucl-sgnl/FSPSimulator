# FSPSimulator
This is the updated Python Version of the ODL FSP, written in python.

# Environment Installation:

``` bash
conda update -n base -c defaults conda
conda env create -f settings/fspsim_env.yml
conda activate fspsim/ source activate fspsim
```

# Installation Error Fixes:
If you get module import errors, the first to show will be 'src' not found if you run main.py try::
Mac:
Navigate to project root:
``` bash
export PYTHONPATH="$(pwd)"    
```

Windows: 
For debugging in vscode, if you have a conda environment add the below to your .vscode/launch.json
``` json
    "module": "fspsim",
```
Add the below to main.py file
``` python
    import sys
    import os
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))) 
```

if you get pyatmos installation error on Mac import errors, navigate to project root and try:  
``` bash
brew install fftw  
```