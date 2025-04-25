"""

The probability of each bond breaking is p1 and can change during the simulation.
The bonds are in a column with the force perpendicular to the column.
All bonds breaking would be a full rapture.

"""


#Importing import libraries
import numpy as np
import pandas as pd
from Bonds import bonds # Class bonds in Bonds.py
from Bonds import add_bonds_to_dic


# Parameters
model = 'p_increase' # 'p_const', 'p_increase', 'gradiant' , 'zipper'
#The types of model:
#    1) 'p_const '- The force on each bond is constant throughout the simulation even after a bond breaks.
#    2) 'p_increase' - The force on a bond increases after a bond breaks F/not_broken_bonds.
#    3) 'gradiant' - The force is higher at the top bond and decreases further down (like bending).
#    4) 'zipper'- The bonds brake one by one with all the force on one bond (like a zipper).
rebind_on = False  # If True the bonds can rebind after breaking.
Rebind_P1 = 0.00001 #The probability of a bond to rebind.
number_of_bonds = 10
number_of_s_changes = 1 #Number of stress changes. At each stress the probability of a bond breaking will change
number_of_repetitions = 100 # Number of repetitions for each stress (or probability)
starting_p0 = 0.0001 #The probability of a bond breaking for the first stress simulation
p_change = 0.0001 # The change in the probability of a bond breaking (if modeling for several stresses)
f_const = 1. # The Force/KbT. In the 'p_increase' and the 'gradiant' model f_const is used to calculate
# the change in the probability of a bond breaking after a bond is broken or rebinds

#File Names
"""
To save the data from the simulations there are 3 file:
1) file_name_bonds: Will save the info of when each bond breaks or rebinds. 
This could be a lot of data. Therefore, it will only be used if write_bonds is True
2) file_name_ns: Will save the data of the number of cycles when all the bonds break for each stress.
This could also be a lot of data and will only be used if write_ns is true.
3) file_name_stress: Will save the data with the mean, median, and standard deviation at each stress. 
"""

directory = 'results/'
date = '2025-04-25'
write_bonds = False # If True will write file #1
file_name_b = 'all_breaks' # A unique str to add to file_name_bonds

write_ns = False # If True will write file #2
file_name_alln = 'broke_at' # A unique str to add to file_name_ns

file_name_s = 'S_N' # A unique str to add to file_name_stress

# Dictionary for the mean, median, and standard deviation at each stress (or P).
all_stress = {'P': [],
              'Mean': [],
              'Median': [],
              'STD': []
              }

if write_ns:
    # If true will write to file_name_ns
    file_name_ns = directory + date + '_' + model + '_reb' + str(rebind_on) + file_name_alln + '.csv'
    # Writes the details of the simulation
    with open(file_name_ns, 'w') as file_ns:
        file_ns.write(
            f"Total Bonds = {number_of_bonds}, Model = {model}, Rebind = {rebind_on} , Rebind_P1 = {Rebind_P1}\n")

for j in range (number_of_s_changes):
    """
    This loop is designed to measure at increasing stress by changing f_const and starting_p (for S-N curve).
    """
    f_const = f_const  # Can be used to calculate starting_p or the change in p for each stress (will adjust at a future time)
    starting_p = starting_p0 + p_change*j # Increasing the breaking probability of a single bond
    all_stress['P'].append(starting_p) # Appending the breaking probability of a single bond to the dictionary
    all_n = np.array([]) #array to collect the N (number of cycle) where all bonds break for each repetitions

    if write_bonds:
        # If true will create file with file_name_bonds
        file_name_bonds = directory + date + '_' + model + '_reb' + str(rebind_on) + '_P-' + str(starting_p).replace('.', '_') + file_name_b + '.csv'
        with open(file_name_bonds, 'w') as file_b:
            file_b.write("")

    for i in range (number_of_repetitions):
        """
        Loop to run number_of_repetitions for each stress (to get average, SD, ETC.)
        """
        # Initializing  the class bonds for each step
        material = bonds(number_of_bonds, starting_p, f_const, model = model, rebind_on = rebind_on, Rebind_P1 = Rebind_P1)

        n = 0 #Counting the number of cycles in each repetition
        new_broken = np.array([]) #

        if write_bonds:
            # Dictionary to collect info of when each bond breaks or rebinds
            dic_bonds = {'Cycle': [],
                          'Broken': [],
                          'Broke/Rebind this step': []
                          }

        while material.not_broken > 0:
            """
            Loop to run cycles till all bonds are broken or n > 10^8
            """
            n += 1
            if n > 10**8: # Stops the while loop if run more than 10^8 times.
                print("n = {n}")
                break

            if material.rebind_on:
                # If rebind option is on will test to see if a broken bond can rebind by calling bonds.rebind() function
                did_rebind, new_bind= material.rebind()
                if did_rebind and write_bonds:
                    # If a bond rebinds and writ_bonds is True it will add the detail to the dictionary
                    dic_bonds=add_bonds_to_dic(dic_bonds, n, material.broken, new_bind)

            did_break, new_broken = material.break_bond() # One cycle to check if any bond breaks.
            if did_break and write_bonds:
                # If a bond breaks and write_bonds is true will add the detail to the dictionary
                dic_bonds = add_bonds_to_dic(dic_bonds, n, material.broken, new_broken)

        all_n = np.append(all_n, n) # Once all bonds are broken and the while loop ends will add n to the list

        if write_bonds:
            # If true will write to file_name_bonds
            # Writes the details of the simulation to a file with name 'ile_name_bonds'
            with open(file_name_bonds, 'a') as file_b:
                file_b.write(f"Repetition = {i+1}, Total Bonds = {material.NB}, P1 = {material.p1}, F_const = {material.F_const}, Model = {material.MODEL}, Rebind = {rebind_on} , Rebind_P1 = {Rebind_P1}\n")
            # Writes the details of dic_bonds
            df_bonds = pd.DataFrame(dic_bonds) # Converts the dictionary to a Data Frame for easier writing
            df_bonds.to_csv(file_name_bonds, sep = ',', mode = 'a', index = False)

    if write_ns:
        #If true will write to file_name_ns
        df_all_n = pd.DataFrame(all_n) # Converts the list to a Data Frame for easier writing
        df_all_n.columns = [f'Breaking Cycle [P1 = {material.p1}, F_const = {material.F_const}]']
        # Writes the number of cycles till all bonds break for all 'number_of_repetitions' to file 'file_name_ns'
        df_all_n.T.to_csv(file_name_ns, sep=',', mode='a')

    # Calculates the mean, median and std for each stress from the number of cycles of all the repetitions
    # and adds it to the dictionary 'all_stress'
    all_stress['Mean'].append(round(all_n.mean(), 2))
    all_stress['Median'].append(np.median(all_n))
    all_stress['STD'].append(round(all_n.std(), 2))

# Writing file_name_stress
df_all_stress = pd.DataFrame(all_stress) # Converts the dictionary to a Data Frame for easier writing
file_name_stress = directory + date + '_' + model + '_reb' + str(rebind_on) +  file_name_s + '.csv'
# Writes the details of the simulation
# and 'df_all_stress' to the file with name 'file_name_stress'
with open(file_name_stress, 'w') as file_s:
    file_s.write(f"Total Bonds = {material.NB}, F_const = {material.F_const}, Model = {material.MODEL}, Rebind = {rebind_on} , Rebind_P1 = {Rebind_P1}\n")
df_all_stress.to_csv(file_name_stress, sep=',', mode='a', index = False)








