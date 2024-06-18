import cantera as ct  # Cantera library for combustion simulations
import numpy as np  # Numerical Python library for mathematical operations
import matplotlib.pyplot as plt  # Library for creating plots

# Load the gas mixture from the gri30 database
gas = ct.Solution("gri30.yaml")

def calculate_adibatic_flame_temp(fuel_species, phi_min, phi_max, npoints, T, P):
    """
    Calculate adiabatic flame temperature and equilibrium composition for a given fuel.

    Parameters:
    fuel_species (str): The chemical formula of the fuel (e.g., 'C2H6' for ethane).
    phi_min (float): The minimum equivalence ratio.
    phi_max (float): The maximum equivalence ratio.
    npoints (int): The number of points to calculate between phi_min and phi_max.
    T (float): Initial temperature of the gas mixture.
    P (float): Initial pressure of the gas mixture.

    Returns:
    tuple: Arrays of equivalence ratios, adiabatic flame temperatures, and equilibrium compositions.
    """
    # Indices for the fuel, oxygen, and nitrogen species
    ifuel = gas.species_index(fuel_species)
    io2 = gas.species_index("O2")
    in2 = gas.species_index("N2")

    # Air composition: Nitrogen to Oxygen ratio
    air_N2_O2_ratio = 3.76

    # Calculate stoichiometric oxygen needed for the fuel
    if fuel_species == 'C2H6':
        stoich_O2 = 2 * gas.n_atoms(fuel_species, "C") + 0.25 * gas.n_atoms(fuel_species, "H")
    elif fuel_species == 'CH4':
        stoich_O2 = gas.n_atoms(fuel_species, "C") + 0.25 * gas.n_atoms(fuel_species, "H")

    # Create arrays to store equivalence ratios, temperatures, and compositions
    phi_values = np.linspace(phi_min, phi_max, npoints)
    tad = np.zeros(npoints)
    xeq = np.zeros((gas.n_species, npoints))

    # Loop over each equivalence ratio
    for i, phi in enumerate(phi_values):
        # Set the composition of the gas mixture
        X = np.zeros(gas.n_species)
        X[ifuel] = phi  # Fuel mole fraction
        X[io2] = stoich_O2  # Oxygen mole fraction
        X[in2] = stoich_O2 * air_N2_O2_ratio  # Nitrogen mole fraction

        # Set the gas state: temperature, pressure, and composition
        gas.TPX = T, P, X
        # Equilibrate the mixture adiabatically at constant pressure
        gas.equilibrate("HP")

        # Store the adiabatic flame temperature and equilibrium composition
        tad[i] = gas.T
        xeq[:, i] = gas.X
        print(f"At phi= {phi:.4f} Tad = {gas.T:.4f}")

    return phi_values, tad, xeq

def plot_species_mass_fractions(fuel_species, phi_values, xeq, species_names):
    """
    Plot the mass fractions of selected species against the equivalence ratio.

    Parameters:
    fuel_species (str): The chemical formula of the fuel.
    phi_values (array): Array of equivalence ratios.
    xeq (array): Array of equilibrium compositions.
    species_names (list): List of species names.
    """
    # Select species to plot
    for species in ["O2", "CO2", "CO"]:
        if species in species_names:
            index = species_names.index(species)
            plt.plot(phi_values, xeq[index, :], label=species)

    plt.xlabel("Equivalence ratio [-]")
    plt.ylabel("Mass fractions [-]")
    plt.legend(loc="best")
    plt.title(fuel_species)
    plt.savefig(f'mass_fractions_{fuel_species}.png')
    plt.show()

def plot_adibatic_flame_temp(fuel_species, phi_values, tad):
    """
    Plot the adiabatic flame temperature against the equivalence ratio.

    Parameters:
    fuel_species (str): The chemical formula of the fuel.
    phi_values (array): Array of equivalence ratios.
    tad (array): Array of adiabatic flame temperatures.
    """
    plt.plot(phi_values, tad)
    plt.xlabel("Equivalence ratio [-]")
    plt.ylabel("Adiabatic flame temperature [K]")
    plt.title(fuel_species)
    plt.savefig(f'flame_temp_{fuel_species}.png')
    plt.show()

def make_calc_for_ethane(T, P):
    """
    Perform calculations and plots for ethane (C2H6).
    """
    # Parameters for ethane
    fuel_species = "C2H6"
    phi_min = 0.3
    phi_max = 4
    npoints = 100

    # Calculate adiabatic flame temperature and equilibrium composition
    phi_values, tad, xeq = calculate_adibatic_flame_temp(fuel_species, phi_min, phi_max, npoints, T, P)
    return phi_values, tad, xeq

def make_calc_for_methane(T, P):
    """
    Perform calculations and plots for methane (CH4).
    """
    # Parameters for methane
    fuel_species = "CH4"
    phi_min = 0.3
    phi_max = 4
    npoints = 100

    # Calculate adiabatic flame temperature and equilibrium composition
    phi_values, tad, xeq = calculate_adibatic_flame_temp(fuel_species, phi_min, phi_max, npoints, T, P)
    return phi_values, tad, xeq

def make_adiabatic_plot_for_few_cases(phi_values, tad, fuel_species):
    for i in range(len(phi_values)):
        plt.plot(phi_values[i], tad[i])
    plt.xlabel("Equivalence ratio [-]")
    plt.ylabel("Adiabatic flame temperature [K]")
    plt.title(fuel_species)
    plt.legend(['T1', 'T2', "T3"])
    plt.savefig(f'flame_temp_{fuel_species}_P1.png')
    plt.show()

def make_species_plot_for_few_cases(phi_values, tad, fuel_species, species_names):
    for i in range(len(phi_values)):
        plt.figure()
        for species in ["O2", "CO2", "CO"]:
            if species in species_names:
                index = species_names.index(species)
                plt.plot(phi_values[i], xeq[i][index, :], label=species)
        plt.xlabel("Equivalence ratio [-]")
        plt.ylabel("Mass fractions [-]")
        plt.legend(loc="best")
        plt.title(f'{fuel_species}, T{i+1}')
        plt.savefig(f'mass_fractions_{fuel_species}_T{i}_P1.png')
        plt.show()

T_list = [273, 298, 323]
P_list = [101325, 202650, 303975]

phi_values = [0,0,0]
tad = [0,0,0]
xeq = [0,0,0]

# Perform calculations and plots for ethane
for i in range(len(P_list)):
    phi_values[i], tad[i], xeq[i] = make_calc_for_ethane(T_list[i],P_list[0])

make_adiabatic_plot_for_few_cases(phi_values, tad, "CH4")
make_species_plot_for_few_cases(phi_values, tad, "CH4", gas.species_names)

# Perform calculations and plots for methane
# for i in range(len(P_list)):
#     phi_values[i], tad[i], xeq[i] = make_calc_for_methane(T_list[i],P_list[1])

# make_adiabatic_plot_for_few_cases(phi_values, tad, "C2H6")
# make_species_plot_for_few_cases(phi_values, tad, "C2H6", gas.species_names)
