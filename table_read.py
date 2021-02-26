"""""""""
Import Packages and Setup Python Environment
"""""""""
import numpy as np
import argparse

b_def   = 1.5
t_def   = 0.0
k_def   = 's'
phi_def = 0.0
chi_def = 90

def read_Mach_table(tableFile, choose_k=k_def):
    # Reads a .dat file containing the maximum launch Mach numbers for a range of wind geometries and density/temperature profiles and stores in a dictionary.
    # By default returns only for spherical temperature profiles. For cylindrical call with choose_k='c'.
    table_data = np.genfromtxt(tableFile, names=True, skip_header=0, comments='#')
    Mach_data = table_data['M_b']

    # Extract density slopes b
    try:
        all_b = table_data['b']
    except ValueError:
        all_b = np.full_like(Mach_data, b_def)
    unique_b = np.unique(all_b)

    # Extract temperature slopes t
    try:
        all_t = table_data['t']
    except ValueError:
        all_t = np.full_like(Mach_data, t_def)
    unique_t = np.unique(all_t)

    # Extract temperature keys k
    try:
        all_k = table_data['T_key']
    except ValueError:
        all_k = np.array([k_def]*len(Mach_data))
    unique_k = np.unique(all_k)

    # Extract base angles phi_b
    try:
        all_phi = table_data['phi_b']
    except ValueError:
        all_phi = np.full_like(Mach_data, phi_def)
    unique_phi = np.unique(all_phi)

    # Extract launch angles chi_b
    try:
        all_chi = table_data['chi_bpi']
    except ValueError:
        all_chi = np.full_like(Mach_data, chi_def/180)
    unique_chi = np.unique(all_chi)

    # Create dictionary
    Mach_table = {}
    for potential_b in unique_b:
        Mach_table[potential_b] = {}
        for potential_t in unique_t:
            Mach_table[potential_b][potential_t] = {}
            for potential_phi in unique_phi:
                Mach_table[potential_b][potential_t][potential_phi] = {}
                for potential_chi in unique_chi:
                    Mach_table[potential_b][potential_t][potential_phi][potential_chi] = np.inf # Initialise to infinity so all true Mach numbers are less

    # Populate dictionary with lowest Mach fitting those parameters
    for b, t, k, phi_b, chi_b, Mach_b in zip(all_b, all_t, all_k, all_phi, all_chi, Mach_data):
        if Mach_b < Mach_table[b][t][phi_b][chi_b] and k==choose_k:
            # If less than existing record, replace it
            Mach_table[b][t][phi_b][chi_b] = Mach_b

    return Mach_table

def search_Mach_table(Mach_table, b, t, phi_b, chi_b):
    try:
        return Mach_table[b][t][phi_b][chi_b]
    except:
        return np.nan

def main():
    """""""""
    Argument Setup
    """""""""
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", "-f", type=str, default='launch_Mach.dat', help='File containing launch Mach numbers')
    parser.add_argument("--b", "-b", type=float, nargs='+', default=[b_def], help='Density Profile Index'    )
    parser.add_argument("--t", "-t", type=float, nargs='+', default=[t_def], help='Temperature Profile Index')
    parser.add_argument("--phi_b", "-p", type=float, nargs='+', default=[phi_def],  help='Angle of the wind base (degrees) [Default: 0]' )
    parser.add_argument("--chi_b", "-c", type=float, nargs='+', default=[chi_def], help='Angle of the wind base (degrees) [Default: 90]')
    args = parser.parse_args()

    lengths = [len(args.b), len(args.t), len(args.phi_b), len(args.chi_b)]
    if not all(l==lengths[0] for l in lengths):
        raise AssertionError("The lists of parameters provided do not all have the same lengths.")

    Mach_table = read_Mach_table(args.file)

    for b, t, phi_b, chi_b in zip(args.b, args.t, args.phi_b, args.chi_b):
        Mach_b = search_Mach_table(Mach_table, b, t, phi_b, chi_b/180)
        if not np.isnan(Mach_b):
            print("b = {:3.2f}, t = {:3.2f}, phi_b = {:2.0f}, chi_b = {:3.3f}: \t M_b = {:3.3f}".format(b, t, phi_b, chi_b, Mach_b))
        else:
            print("Mach number could not be found for b = {:3.2f}, t = {:3.2f}, phi_b = {:2.0f}, chi_b = {:3.3f}".format(b, t, phi_b, chi_b))

    return 1

if __name__=="__main__":
    main()

