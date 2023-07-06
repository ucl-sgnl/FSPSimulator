# FSP Simulator

Welcome to the Python version of the Orbital Debris Lab (ODL) Full Spectrum Project (FSP) Simulator. This modern adaptation provides an updated approach to running simulations.

## Getting Started

All the simulation settings should be specified within `.json` files located in `src/data/specify_simulation`.

When executing `main.py`, it will automatically run all the simulations as per the settings in the respective `.json` files.

## JSON Content Format

A typical configuration for a simulation can be represented in the following format in a JSON file:

```json
{
  "scenario_name": "simulation_name",
  "monthly_ton_capacity": "100",
  "launch_start_date": "2019-01-01",
  "startup_failure_rate": 0,
  "satellite_failure": 0,
  "sim_start_date": "2019-01-01",
  "sim_end_date": "2019-01-08",
  "output_frequency": 360,
  "integrator_step_size,": 20,
  "integrator_type": "RK45",
  "sim_object_type": "all",
  "sim_object_catalogue" : "both",
  "environment": "development",
  "repull_catalogues": false
}


# Environment Installation:

``` bash
conda update -n base -c defaults conda
conda env create -f settings/fspsim_env.yml
conda activate fspsim/ source activate fspsim
```

To run on Virtual Machine you will need miniconda to activate the environment:
``` bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda --version

# Installation Error Fixes:
If you get module import errors, the first to show will be 'src' not found if you run main.py try::
Mac:
Navigate to project root:
``` bash
export PYTHONPATH="$(pwd)" >> ~/.zshrc
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