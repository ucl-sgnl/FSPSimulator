import os
import pickle
import numpy as np
def calculate_energy_from_semi_major_axis(a, m):
    """Calculate the total mechanical energy of an object in orbit around the Earth, given its semi-major axis (a) and mass (m)."""
    G = 6.67430e-11  # Gravitational constant, units: m^3 kg^-1 s^-2
    M_earth = 5.97219e24  # Mass of Earth, units: kg
    energy = -G * M_earth * m / (2 * a)
    return energy


def load_file(simName):
    # Get the absolute path of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Construct the absolute path to the pickle file
    path = os.path.join(script_dir, f'../../../src/data/catalogue/{simName}.pickle')
    
    with open(path, 'rb') as f:
        # Load the data from the pickle file
        data = pickle.load(f)
        return data
    
baseline_data = load_file('SATCAT_after_prop_baseline')
all_sats = load_file('SATCAT_after_prop_allsats')
all_sats_no_espace = load_file('SATCAT_after_prop_allsats_noespace')

data_sources = [baseline_data, all_sats, all_sats_no_espace]
data_source_names = ['Baseline Data', 'All Sats', 'All Sats No Espace']

for data, name in zip(data_sources, data_source_names):

    total_no_epochs = 0
    ephemerides = []

    # Filter satellites based on decay_date
    filtered_data = [satellite for satellite in data if pd.isna(satellite.decay_date) or satellite.decay_date > datetime.datetime(2028, 1, 1)]

    for satellite in filtered_data:
        number_of_epochs = (np.shape(satellite.ephemeris)[0])
        if number_of_epochs > total_no_epochs:
            total_no_epochs = number_of_epochs
            ephemerides = np.array(satellite.ephemeris)
    print("largest number of epochs: ", total_no_epochs)

    dates = []
    for i in range(0, total_no_epochs, 1):
        dates.append(ephemerides[i][0])
    print("dates: ", dates)

    yearly_energies = []
    for date in dates:
        kinetic_energy = 0
        potential_energy = 0
        for satellite in filtered_data:
            for epoch in satellite.ephemeris:
                if epoch[0] == date:
                    total_E += calculate_energy_from_semi_major_axis(epoch[2][0], epoch[2][1], epoch[2][2])
        yearly_energies.append(total_E)

    zipped = zip(dates, yearly_energies)
    plt.plot(dates, yearly_energies, label=name)

plt.xlabel('Date')
plt.ylabel('Total Energy (Joules)')
plt.legend()
plt.show()

plots_folder = "plots"
