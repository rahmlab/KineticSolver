from datareader import *
from constants import R, R_kcal, kB, h
import numpy as np
import sys
import decimal

def ResetAll(normaldist=False,scale=0.0):
    ResetKineticConstants(normaldist,scale)
    ResetConcs()
    ResetDerivatives()


def ResetKineticConstants(normaldist=False,scale=0.0):

    if normaldist:
        r_shifts=np.random.normal(loc=0.0,scale=scale,size=N_reactions)
        TS_shifts=np.random.normal(loc=0.0,scale=scale,size=N_reactions)
        print(f'TS shifts:\t {TS_shifts}')
        print(f'reac shifts:\t {r_shifts}')

    for n in range(N_reactions):
        if normaldist:
            DeltaGTS_curr = reactions[n].DeltaG_TS + TS_shifts[n]
            DeltaGr_curr = reactions[n].DeltaG_reaction + r_shifts[n]
            if reactions[n].ref[0] != -1:
                for ref_n in reactions[n].ref:
                    DeltaGTS_curr = DeltaGTS_curr - r_shifts[ref_n]
                    DeltaGr_curr = DeltaGr_curr - r_shifts[ref_n]

            if DeltaGr_curr < 0.0 and DeltaGTS_curr < 0:
                DeltaGTS_curr = 0.0
            elif DeltaGr_curr > 0.0 and DeltaGr_curr > DeltaGTS_curr:
                DeltaGTS_curr = DeltaGr_curr + 0.0
        else:
            DeltaGTS_curr = reactions[n].DeltaG_TS
            DeltaGr_curr = reactions[n].DeltaG_reaction
        reactions[n].kinetic_constant_fw = Eyring(DeltaGTS_curr,T)
        reactions[n].kinetic_constant_bw = Eyring(DeltaGTS_curr - DeltaGr_curr,T)
    
        print(f'Reaction {n}: {DeltaGTS_curr}, {DeltaGr_curr}')
"""
def ResetConcs():
    for sp in species:
        sp.conc = 0.0
        if sp.name in initial_species:
            sp.conc = initial_species[sp.name] 

"""
def ResetConcs():
    concs = np.zeros(len(species))
    for i, sp in enumerate(species):
        if sp.name in initial_species:
            concs[i] = initial_species[sp.name]
        sp.conc = concs[i]



"""
def UpdateConcProd():
    for rxn in reactions:
            rxn.conc_product_fw = 1
            rxn.conc_product_bw = 1

            for rn in range(len(rxn.reactants)):
                rxn.conc_product_fw *= rxn.reactants[rn].conc ** rxn.r_stoic[rn]

            for pn in range(len(rxn.products)):
                rxn.conc_product_bw *= rxn.products[pn].conc ** rxn.p_stoic[pn]

"""
"""
def UpdateDerivatives(c):
    
    #Vectorized version of UpdateDerivatives.
    #Updates the derivatives of species concentrations based on reaction kinetics.
   
    #Parameters:
     #   c (array): Current concentrations of all species.
    
    # Update species concentrations
    for i, sp in enumerate(species):
        sp.conc = c[i]
        sp.derivative = 0.0  # Reset derivatives

    # Prepare arrays for vectorized calculations
    num_reactions = len(reactions)
    num_species = len(species)

    reactant_matrix = np.zeros((num_reactions, num_species))  # Stoichiometries of reactants
    product_matrix = np.zeros((num_reactions, num_species))   # Stoichiometries of products
    forward_rates = np.zeros(num_reactions)                   # Forward reaction rates
    backward_rates = np.zeros(num_reactions)                  # Backward reaction rates

    # Populate stoichiometry matrices and rate constants
    for i, rxn in enumerate(reactions):
        for sp, stoic in zip(rxn.reactants, rxn.r_stoic):
            reactant_matrix[i, sp.index] = stoic
        for sp, stoic in zip(rxn.products, rxn.p_stoic):
            product_matrix[i, sp.index] = stoic

        # Calculate forward and backward rates
        reactant_concs = c ** reactant_matrix[i]
        product_concs = c ** product_matrix[i]

        rxn.conc_product_fw = np.prod(reactant_concs)
        rxn.conc_product_bw = np.prod(product_concs)

        forward_rates[i] = rxn.kinetic_constant_fw * rxn.conc_product_fw
        backward_rates[i] = rxn.kinetic_constant_bw * rxn.conc_product_bw

    # Net reaction rates
    net_rates = forward_rates - backward_rates

    # Update derivatives
    for i in range(num_species):
        derivatives = -np.sum(reactant_matrix[:, i] * net_rates) + np.sum(product_matrix[:, i] * net_rates)
        species[i].derivative = derivatives
"""


def UpdateDerivatives(c):
    for i in range(len(c)):
        species[i].conc = c[i]
        species[i].derivative = 0.0

    for rxn in reactions:
        reactant_concs = np.array([sp.conc for sp in rxn.reactants])
        product_concs = np.array([sp.conc for sp in rxn.products])
        r_stoic = np.array(rxn.r_stoic)
        p_stoic = np.array(rxn.p_stoic)

        rxn.conc_product_fw = np.prod(reactant_concs ** r_stoic)
        rxn.conc_product_bw = np.prod(product_concs ** p_stoic)

        fw_term = rxn.kinetic_constant_fw * rxn.conc_product_fw
        bw_term = rxn.kinetic_constant_bw * rxn.conc_product_bw

        for sp, stoic in zip(rxn.reactants, r_stoic):
            sp.derivative += bw_term * stoic - fw_term * stoic

        for sp, stoic in zip(rxn.products, p_stoic):
            sp.derivative += fw_term * stoic - bw_term * stoic

"""
def UpdateDerivatives():
    for sp in species:
        sp.derivative = 0.0

    for rxn in reactions:
        for n_r in range(len(rxn.reactants)):
            #e.g.2A -> B. r = -1/2 d[A]/dt = k[A]^2 => d[A]/dt = -2k[A]^2
            #    2A <- B  r = d[B]/dt = 1/2 d[A]/dt = k[B] => d[A]/dt = 2k[B]
            rxn.reactants[n_r].derivative -= rxn.kinetic_constant_fw * rxn.conc_product_fw * rxn.r_stoic[n_r]
            rxn.reactants[n_r].derivative += rxn.kinetic_constant_bw * rxn.conc_product_bw * rxn.r_stoic[n_r]

        for n_p in range(len(rxn.products)):
            rxn.products[n_p].derivative += rxn.kinetic_constant_fw * rxn.conc_product_fw * rxn.p_stoic[n_p]
            rxn.products[n_p].derivative -= rxn.kinetic_constant_bw * rxn.conc_product_bw * rxn.p_stoic[n_p]

"""

def ResetDerivatives():
    for sp in species:
        sp.derivative = 0.0

    for rxn in reactions:
        reactant_concs = np.array([sp.conc for sp in rxn.reactants])
        product_concs = np.array([sp.conc for sp in rxn.products])
        r_stoic = np.array(rxn.r_stoic)
        p_stoic = np.array(rxn.p_stoic)

        fw_term = rxn.kinetic_constant_fw * rxn.conc_product_fw
        bw_term = rxn.kinetic_constant_bw * rxn.conc_product_bw

        for sp, stoic in zip(rxn.reactants, r_stoic):
            sp.derivative += bw_term * stoic - fw_term * stoic

        for sp, stoic in zip(rxn.products, p_stoic):
            sp.derivative += fw_term * stoic - bw_term * stoic


def CalcRate(c,rxn_list):
    N_times = len(c[:,0])
    rate = np.zeros(N_times)
    for rxn_n in rxn_list:
        r_concprod = np.ones(N_times)
        for r_n in range(len(reactions[rxn_n].reactants)):
            r_concprod *= c[:,reactions[rxn_n].reactants[r_n].index] ** reactions[rxn_n].r_stoic[r_n]
        p_concprod = np.ones(N_times)
        for p_n in range(len(reactions[rxn_n].products)):
            p_concprod *= c[:,reactions[rxn_n].products[p_n].index] ** reactions[rxn_n].p_stoic[p_n]

    
        rate += reactions[rxn_n].kinetic_constant_fw * r_concprod 
        rate -= reactions[rxn_n].kinetic_constant_bw * p_concprod
    return rate


def Get_AdeFormationRate(c):
    #rxn 5: p11(5) <==> Ade(22)
    #rxn 13: p16(11) <==> Ade(22) + HCN(0)
    rate_rxn5 =     reactions[5].kinetic_constant_fw * c[5]     - reactions[5].kinetic_constant_bw * c[22]
    rate_rxn13 =    reactions[13].kinetic_constant_fw * c[11]   - reactions[13].kinetic_constant_bw * c[22] * c[0]
    return rate_rxn5 + rate_rxn13


def CalcDRC(rate_perturb,rate_unp,pct):
    # Degree of rate-control: X = ki/r * dr/dki
    # Note that ki/dki = 1/pct
    #return (rate_perturb - rate_unp)/rate_unp/(pct*0.01)
    #return (rate_perturb/rate_unp - 1.00)/(pct*0.01)
    return (np.log(rate_perturb)-np.log(rate_unp)) / np.log(1+pct*0.01)

def CalcDRCfromConc(c,c_unp,pct):
    return (np.log(c)-np.log(c_unp)) / np.log(1+pct*0.01)

def Eyring(DeltaG,T):
    return kB*T/h*np.exp(-DeltaG/(R_kcal*T))

def CalcBarrier(k):
    return -R_kcal*T*np.log(k*h/(kB*T))

def AlmostLogSpacedIndexes(N_points):       #Why not log? Because it don't want repeated indexes. You can get: 1 1 1 2 2 4 8 ...
    base = 11**(1/N_points)

    steps = []

    for i in range(N_points):
        steps.append(base**i)
    indexes = []

    for st in steps:
        indexes.append(int((st-1)*N_timesteps/10))

    return indexes

def powers_of_tens_array(num,step=1.0):

    i_max = int(decimal.Decimal(num).log10())
    #i_max = int(decimal.Decimal(num).ln())
    
    if i_max <= 0:
        i = i_max - 1
    else:
        i = 0.0

    powers = []

    while i <= i_max:
        powers.append(10**i)
        #powers.append(np.e**i)
        i+=step

    return powers


