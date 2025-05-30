"""
Simulating mechanical fatigue of NB bonds.
The probability of each bond breaking is p1 and can change during the simulation.
The bonds are in a column with the force perpendicular to the column.
All bonds breaking would be a full rapture.

"""


#Importing libraries
import numpy as np
import pandas as pd
import math
from Bonds import bonds # Class bonds in Bonds.py
from Bonds import add_bonds_to_dic
from Bonds import plot_hist
from Bonds import plot_p_n
from tqdm import tqdm #Time bar
import os


# Parameters
model = 'p_increase' # 'p_const', 'p_increase',  'zipper'
#The types of model:
#    1) 'p_const '- The force on each bond is constant throughout the simulation even after a bond breaks.
#    2) 'p_increase' - The force on a bond increases after a bond breaks F/not_broken_bonds.
#    3) 'zipper'- The bonds brake one by one with all the force on one bond (like a zipper).
rebind_on = True  # If True the bonds can rebind after breaking.
Rebind_P1 = 0.0001 #The probability of a bond to rebind.
number_of_bonds = 13
number_of_s_changes =  10 # Number of stress changes. At each stress the probability of a bond breaking will change
number_of_repetitions = 100 # Number of repetitions for each stress (or probability)
starting_p0 = 2.67*10**-8 #The probability of a bond breaking with zero external force
f_change = 3 # The change in the force of a bond breaking (if modeling for several stresses)
f_const0 = 40 # The Total Force*r in kbT unit for the first stress.
f_const = f_const0 # The Total Force*r that changes for each stress change.
# In the 'p_increase'
#File Names
"""
To save the data from the simulations there are 3 file:
1) file_name_bonds: Will save the info of when each bond breaks or rebinds. 
This could be a lot of data. Therefore, it will only be used if write_bonds is True
2) file_name_ns: Will save the data of the number of cycles when all the bonds break for each stress.
This could also be a lot of data and will only be used if write_ns is true.
3) file_name_stress: Will save the data with the mean, median, and standard deviation at each stress. 
"""

directory = 'results/p_increase/'
date = '2025-05-29'
write_bonds = False # If True will write file #1
file_name_b = 'all_breaks_00' # A unique str to add to file_name_bonds

write_ns = True # If True will write file #2
file_name_alln = 'broke_at_N_00' # A unique str to add to file_name_ns

file_name_s = 'S_N_00' # A unique str to add to file_name_stress

if not os.path.exists(directory):
    # if directory does not exist it will create it
    os.makedirs(directory)

#Option to show plots of the results
to_plot_hist = False # If true will plot a histogram of the number of cycles for a full break for each probability
to_plot_p_n = False  # If true will plot a probability vs median number of cycles for a full break



# Dictionary for the mean, median, and standard deviation at each stress (or P).
all_stress = {'P': [],
              'Fr': [],
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

for j in tqdm (range (number_of_s_changes), desc = f"P range: {starting_p0 * math.exp(f_const0/number_of_bonds)} to {starting_p0 * math.exp((f_const0 + f_change*(number_of_s_changes-1))/number_of_bonds)}"):
    """
    This loop is designed to measure at increasing stress by changing f_const and starting_p (for S-N curve).
    tqdm is used to show a progress bar.
    """
    f_const = f_const0 + f_change*j  # Is used to calculate starting_p or the change in p for each stress
    if model == 'zipper':
        starting_p = starting_p0 * math.exp(f_const)  #If the model is 'zipper' all the force is on one bond
    else:
        starting_p = starting_p0 * math.exp(f_const/number_of_bonds)  # Increasing the breaking probability of a single bond

    all_stress['P'].append(starting_p) # Appending the breaking probability of a single bond to the dictionary
    all_stress['Fr'].append(f_const)  # Appending the force to the dictionary
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
                #print(f"{i}) n: {n}")  #print loop number and n = 10^8
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
            # Writes the details of the simulation to a file with name 'file_name_bonds'
            with open(file_name_bonds, 'a') as file_b:
                file_b.write(f"Repetition = {i+1}, Total Bonds = {material.NB}, P1 = {material.p1}, F_const = {material.F_const}, Model = {material.MODEL}, Rebind = {rebind_on} , Rebind_P1 = {Rebind_P1}\n")
            # Writes the details of dic_bonds
            df_bonds = pd.DataFrame(dic_bonds) # Converts the dictionary to a Data Frame for easier writing
            df_bonds.to_csv(file_name_bonds, sep = ',', mode = 'a', index = False)

    if to_plot_hist:
        # If true will plot a histogram
        plot_hist(all_n, number_of_repetitions//250, starting_p)


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
    file_s.write(f"Total Bonds = {material.NB}, Model = {material.MODEL}, Rebind = {rebind_on} , Rebind_P1 = {Rebind_P1}\n")
df_all_stress.to_csv(file_name_stress, sep=',', mode='a', index = False)

if to_plot_p_n:
    # If true will plot a p_n curve
    plot_p_n(df_all_stress)


