import numpy as np
import random
import math
import matplotlib.pyplot as plt


class bonds:
# This class sets up an object of NB bonds and details the different operation that can be done.

    def __init__(self, NB, p1, F_const = 1.0, model = 'p_const', rebind_on = False, Rebind_P1 = 0.0):
        """
        This function initiates the class and generates a vector with NB bonds and a vector with
        the corresponding probabilities.

        Args:
        NB (int): Number of bonds.
        p1 (float): Initial probability of a bond breaking.
        F_const (float): A constant to calculate the probability change when using 'p_increase'.
        model (str): The type of model to run:
        1) 'p_const '- The force on each bond is constant throughout the simulation even after a bond breaks.
        2) 'p_increase' - The force on a bond increases after a bond breaks F/not_broken_bonds.
        3) 'zipper'- The bonds break one by one with all the force on one bond (like a zipper).
        rebind_on (boolean): If True the bonds can rebind after breaking.
        Rebind_P1 (float): The probability of a bond to rebind.

        """
        #self parameters
        self.NB = NB
        self.p1 = p1
        self.F_const = F_const
        self.MODEL = model
        self.rebind_on = rebind_on
        self.Rebind_P1 = Rebind_P1

        self.broken = 0 #Number of boroken bonds
        self.not_broken = self.NB  # Number of not broken bonds

        self.allbonds = np.ones(shape=(self.NB)) #Bond vector. 1 for bond 0 for broken bond
        self.rebind_ps = np.zeros(shape=(self.NB))  # Vector with the probability of rebinding for each bond

        #Seting up the initial probability vector for each model type
        if self.MODEL == 'zipper':
            self.p_start_zipper()
        else:
            self.p_start_same()



    def break_bond (self):
        """
                This function tests if a bond should break and if it does it will adjust the
                bond and probability vectors.

                Returns:
                did_break(boolean): If a bond is broken it will return true
                (np_array (float)) If bond is broken return the number of bond that broke.
        """

        did_break = False

        # If a random number (0-1) is larger the the probability of the bond breaking it will break.
        rand_array = np.random.rand(self.NB) - self.ps
        broke = (rand_array<0).sum()

        if broke > 0:
            # If a bond broke it will reset allbonds with 0 at the broken location and reset the probabilities self.ps
            self.allbonds = np.where(rand_array < 0, 0, self.allbonds)
            self.not_broken -= broke
            self.broken += broke
            self.new_p1() #Set the new probabilities after a bond broke.
            if self.rebind_on:
                self.rebind_ps = np.where(rand_array < 0, self.Rebind_P1, self.rebind_ps)
            did_break = True
            broken_bond_location = np.transpose(np.argwhere(rand_array < 0)+1)
            return did_break, broken_bond_location

        return did_break, []


    def new_p1(self):
        #Adjstes the probability vector ps with the new probabilities for each bond after a bond breaks.
        if self.not_broken == 0:
            self.ps[:] = 0
        elif self.MODEL == 'p_increase':
            self.cal_p_increase()
        elif self.MODEL == 'zipper':
            self.ps[self.broken] = self.p1
            self.ps[self.broken-1] = 0
        elif self.MODEL == 'p_const':
            self.ps = np.where(self.allbonds == 0, 0, self.p1)


    def p_start_same (self):
        #sets the vector ps with all probabilities the same
        self.ps = np.repeat(self.p1, self.NB)


    def p_start_zipper(self):
        # sets the vector ps with probabilities of breaking only for the first bond
        self.ps = np.zeros(shape=(self.NB))
        self.ps[0] = self.p1



    def rebind(self):
        # Checks if a bond rebinds for each model type and calculates the
        did_rebind = False
        rebind_bond_location = []
        if self.MODEL == 'zipper':
            if random.random() < self.Rebind_P1:
                self.bonds_rebinds(1)
                self.ps[self.broken] = self.p1
                self.ps[self.broken + 1] = 0
                did_rebind = True
                return did_rebind, [self.broken+1]

        else:
            rand_rebind_array = np.random.rand(self.NB) - self.rebind_ps
            rebinds = (rand_rebind_array < 0).sum()
            if rebinds > 0:
                self.bonds_rebinds(rebinds)
                #Update the vectors
                self.allbonds = np.where(rand_rebind_array < 0, 1, self.allbonds)
                self.rebind_ps = np.where(rand_rebind_array < 0, 0, self.rebind_ps)
                rebind_bond_location = np.transpose(np.argwhere(rand_rebind_array < 0) + 1)
                did_rebind = True

                if self.MODEL == 'p_increase':
                    #Update self.p1
                    self.cal_p_increase()
                    self.ps = np.where(rand_rebind_array < 0, self.p1, self.ps)
                else:
                    self.ps = np.where(rand_rebind_array < 0, self.p1, self.ps)

        return did_rebind, rebind_bond_location



    def bonds_rebinds (self, n_rebinds):
        # Adjusts the bond counter if a bond rebinds.
        self.not_broken += n_rebinds
        self.broken -= n_rebinds

    def cal_p_increase(self):
        # Calculates the probability of each bond breaking after a bond breaks/rebinds
        # and the force increases/decreases on all other bonds.
        newp = self.p1 * math.exp(self.F_const * (1 / float(self.not_broken) - 1 / float(self.NB)))
        self.ps = np.where(self.allbonds == 0, 0, newp)



def add_bonds_to_dic(dic_bonds, n, broken, new_changes):
    """
    This function will append the cycle number, the number of broken bonds, broken bonds,
    and the location of the broken/rebind bonds to the dictionary dic_bonds and return the dictionary

   Args:
       dic_bonds (dictionary)
        n (int): Cycle number
        broken (int): The number of broken bonds
        new_changes (list of ints): The location of the broken/rebind bonds

    Returns:
        dic_bonds (dictionary)

    """

    dic_bonds['Cycle'].append(n)
    dic_bonds['Broken'].append(broken)
    dic_bonds['Broke/Rebind this step'].append(new_changes)

    return (dic_bonds)


def plot_hist (n_s, bins, p1_start):
    """
        This function will plot a histogram

       Args:
            n_s (int): Cycle number
            bins (int): The number bins
            p1_start (float): The probability of the bond breaking

        """
    plt.hist(n_s, bins = bins)
    plt.title(f'Histogram for p = {p1_start}')
    plt.xlabel('n')
    plt.ylabel('Counts')
    plt.show()



def plot_p_n (df_p_n):
    """
        This function will plot a probability (p)
        vs the median of the number of cycles for all bonds to break (n)
        on a log log scale

       Args:
            df_p_n (pandas df)

        """
    plt.loglog(df_p_n['Median'], df_p_n['P'], 'o')
    plt.title('P vs N')
    plt.xlabel('n')
    plt.ylabel('p')
    plt.show()
