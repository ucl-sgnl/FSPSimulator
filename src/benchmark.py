import os
import json
import time
import subprocess
import pandas as pd
import datetime
import numpy as np

def generate_simulation_specs(step_size_range, run_time_years, integrator_types, remove_operators, prediction_csv):
    # Generate step sizes from start to end with a step of 5
    step_sizes = list(range(step_size_range[0], step_size_range[1] + 1, 5))
    
    # Generate end dates for each year in the run time range
    start_date = datetime.datetime(2019, 1, 1)
    end_dates = [(start_date + datetime.timedelta(days=i*365)).strftime('%Y-%m-%d') for i in run_time_years]

    # Initialize the simulation specs list
    simulation_specs_list = []

    # Generate a simulation spec for each combination of step size, end date, integrator type, remove operator, and prediction csv
    for i, (step_size, end_date, integrator_type, remove_operator, pred_csv) in enumerate(np.array(np.meshgrid(step_sizes, end_dates, integrator_types, remove_operators, prediction_csv)).T.reshape(-1,5)):
        simulation_specs_list.append({
            "scenario_name": f"benchmark_{i+1}",
            "monthly_ton_capacity": "100",
            "launch_start_date":"2019-01-01",
            "remove_operators": remove_operator,
            "sim_start_date":"2019-01-01",
            "sim_end_date": end_date,
            "output_frequency": 360,
            "integrator_step_size": int(step_size),
            "integrator_type": integrator_type,
            "sim_object_type": "all",
            "sim_object_catalogue" : "both",
            "environment": "development",
            "repull_catalogues": bool(False),
            "satellite_predictions_csv": pred_csv
        })
        
    return simulation_specs_list

def benchmark_simulation(simulation_specs):
    """Generates a JSON file with the given simulation specs."""
    
    # Generate file name based on the specs
    file_name = f'src/data/specify_simulations/benchmarking_sims/{simulation_specs["scenario_name"]}.json'
    
    # Save the JSON file
    with open(file_name, 'w') as file:
        json.dump(simulation_specs, file)
        
def run_simulations_and_time(run_script_path):
    """Runs the main.py script and returns the time taken."""
    
    # Run the simulation and time it
    start_time = time.time()
    subprocess.run(['python', run_script_path])
    end_time = time.time()

    # Calculate and return the time taken
    return end_time - start_time

def cleanup_simulation_files(simulation_specs_list):
    """Removes the JSON files after the simulation."""
    
    for specs in simulation_specs_list:
        file_name = f'src/data/specify_simulations/benchmarking_sims/{specs["scenario_name"]}.json'
        if os.path.exists(file_name):
            os.remove(file_name)

def main():
    # Define the path to your main.py script
    run_script_path = 'src/main.py'
    step_size_range = [15, 30]
    run_time_years = [0.04, 0.03,0.02,0.01]
    integrator_types = ["RK45"] #"RK23", "DOP853", "Radau", "BDF", "LSODA"
    remove_operators = ["none"]
    prediction_csv = ["benchmark_preds_few.csv", "benchmark_preds_many.csv"]
    #TODO: will have to get the number of objects in the catalogues from the .pickle files

    simulation_specs_list = generate_simulation_specs(step_size_range, run_time_years, integrator_types, remove_operators, prediction_csv)
    print(f"Total number of benchmarks to run: {len(simulation_specs_list)}")
    # Prepare a DataFrame to store the results
    results = pd.DataFrame(columns=['scenario_name', 'step_size', 'output_frequency', 'integrator_type', 'sim_duration', 'time_taken'])

    # Run the simulations and store the results
    for specs in simulation_specs_list:
        benchmark_simulation(specs)

        time_taken = run_simulations_and_time(run_script_path)

        sim_start_date = datetime.datetime.strptime(specs["sim_start_date"], '%Y-%m-%d')
        sim_end_date = datetime.datetime.strptime(specs["sim_end_date"], '%Y-%m-%d')
        sim_duration = (sim_end_date - sim_start_date).days

        results = results.append({
            'scenario_name': specs['scenario_name'],
            'step_size': specs['integrator_step_size'],
            'output_frequency': specs['output_frequency'],
            'integrator_type': specs['integrator_type'],
            'sim_duration': sim_duration,
            'time_taken': time_taken,
            'satellite_predictions_csv': specs['satellite_predictions_csv'],
        }, ignore_index=True)

        cleanup_simulation_files([specs])
                
    # Save the results to a CSV file
    results.to_csv('src/data/results/benchmark_results.csv', index=False)

if __name__ == "__main__":
    main()