#datareader.py: reads reactions.txt and initialize all variables

from constants import *
import numpy as np
import sys
import argparse
import os

class Reaction(object):     #Reaction object
    def __init__(self,num,DeltaG_TS,DeltaG_reaction,kinetic_constant_fw,kinetic_constant_bw,reactants,r_stoic,products,p_stoic,conc_product_fw,conc_product_bw,ref):
        self.num = num                                  # int, reaction number (from 0)
        self.DeltaG_TS = DeltaG_TS              # Gibbs energy of the barrier (TS)
        self.DeltaG_reaction = DeltaG_reaction            # Gibbs energy of the reaction
        self.kinetic_constant_fw = kinetic_constant_fw  # float, forward kinetic constant
        self.kinetic_constant_bw = kinetic_constant_bw  # float, backward kinetic constant
        self.reactants = reactants                      # list of Species objects of reactants
        self.products = products                        # list of Species objects of products
        self.r_stoic = r_stoic                          # array of int, stoichiometries of reactants
        self.p_stoic = p_stoic                          # array of int, stoichiometries of products
        self.conc_product_fw = conc_product_fw          # float, conc. product of reactants as in kinetic law, e.g. 2A + B -> C,  conc_product_fw = [A]^2 * [B]
        self.conc_product_bw = conc_product_bw          # float, conc. product of products as in kinetic law,  e.g. 2A + B -> 3C, conc_product_bw = [C]^3
        self.ref = ref


class Species(object):  #Species object
    def __init__(self,index,name,conc=0.0):
        self.index = index
        self.name = name        #string, species name
        self.conc = conc        #float, concentration
        self.derivative = 0.0   #float, concentration derivative wrt time
        self.Gibbs_energy = 0.0


def Eyring(DeltaG,T):   #Returns kinetic constant with kappa = 1
    return kB*T/h*np.exp(-DeltaG/(R_kcal*T))

def parse_line(line):   #Read and interpret "reactions.txt" file
    r_stoic = []        #stoic of reactants
    p_stoic = []        #stoic of products
    rxn=line.split('|') 
    species=rxn[0].split('>')
    reactants_names = [s.strip() for s in species[0].split('+')]
    for i,reac in enumerate(reactants_names):   #parse reactants
        if reac[0].isdigit():               #if it starts with a number, the number is interpreted as a stoic coefficient
            r_stoic.append(int(reac[0]))
            reactants_names[i] = reac[1:]   #takes the species name only
        else:
            r_stoic.append(1)


    products_names = [s.strip() for s in species[1].split('+')] 
    for i,prod in enumerate(products_names): #parse products
        if prod[0].isdigit():
            p_stoic.append(int(prod[0]))
            products_names[i] = prod[1:]
        else:
            p_stoic.append(1)
  
    energies = [s.strip() for s in rxn[1].split(',')]   #extract energies
    G_TS = float(energies[0])
    G_r = float(energies[1])
    ref = np.array(energies[2:],dtype=int)
    return reactants_names,r_stoic,products_names,p_stoic,G_TS,G_r,ref

def PrintReadReactions():
    print('')
    print('Reactions:')
    for rxn in reactions:
        print(f'{rxn.num}.',end='\t')
        n_reac = len(rxn.reactants)
        n_prod = len(rxn.products)
        for n_r in range(n_reac):
            print(f'{rxn.r_stoic[n_r]}{rxn.reactants[n_r].name}',end='\t')
            if n_r != len(rxn.reactants)-1:
                print(f'+',end='\t')
            elif n_reac == 2:
                print(f'<=>',end='\t')
            else:
                print(f'\t\t<=>',end='\t')
        for n_p in range(len(rxn.products)):
            print(f'{rxn.p_stoic[n_p]}{rxn.products[n_p].name}',end='\t')
            if n_p != len(rxn.products)-1:
                print(f'+',end='\t')
            elif n_prod == 1:
                print(f'', end='\t\t')
        print(f'|\t{rxn.DeltaG_TS},\t{rxn.DeltaG_reaction},\tRef is {rxn.ref}')
    print('')
    print('Initial concentrations:')
    for sp in species:
        if sp.name in initial_species:
            print(f'[{sp.name}] =\t{initial_species[sp.name]}')
    print('')
    print(f'Temperature is {T}')


def CalcSpeciesEnergies():

    reac_side = np.zeros(N_species) 
    prod_side = np.zeros(N_species)
    
    for rxn in reactions:
        for reac in rxn.reactants:
            if reac_side[reac.index] == 0:
                reac_side[reac.index] = 1
        for prod in rxn.products:
            if prod_side[prod.index] == 0:
                prod_side[prod.index] = 1


# ------------------- Program starts ---------------------------------


parser = argparse.ArgumentParser(prog='Kinetic simulator by Marco Cappelletti',description='Perform kinetic simulation of a set of reactions',epilog='You can change the default values in constants.py file. Read the README file!')
parser.add_argument('-f','--reactions-file', type=str, help=f'Name of the input (reactions) file. Default = {reactionsfile}.', default=reactionsfile)
parser.add_argument('-R','--run_type', type=str, help='Type of run.  Default = kinetics.', choices=['kinetics','DRC','distribution'], default='kinetics')
parser.add_argument('-o','--output-file', type=str, help=(f"Name of the csv output file. Default =\n" 
                                                        f"\t{output_kinetics} if run_type = kinetics;\n"
                                                        f"\t{output_DRC} if run_type = DRC;\n"
                                                        f"\t{output_distribution} if run_type = distribution."), default=None)
parser.add_argument('-T', type=float, help=f'Temperature of the simulation. Beware! Gibbs energies depend on temperature. Default = {T}.', default=T)
parser.add_argument('-t','--t_max', type=float, help=f'Maximum time of the simulation. Default = {t_max}.', default=t_max)
parser.add_argument('-N','--N-data', type=int, help=f'Number of data points printed in the output file. Default = {N_data}.', default=N_data)
parser.add_argument('-D','--DRCprod', type=str, help=f'Product for which to calculate the Degree of Rate Control (DRC). Default = {DRCprod}.', default=DRCprod)
parser.add_argument('-c','--percent-change', type=float, help=f'Change of kinetics constant for DRC, in percent. Default = {pct}.', default=pct)
parser.add_argument('-S','--time-steps', type=int, help=f'Number of data (time steps) to solve the ODEs. Even numbers are more efficient. Default = {N_timesteps}.', default=N_timesteps)
parser.add_argument('-a','--atol',type=float,help=f'odeint relative error tolerance. Default = {atol}', default=atol)
parser.add_argument('-r','--rtol',type=float,help=f'odeint absolute error tolerance. Default = {rtol}', default=rtol)
args = parser.parse_args()


reactionsfile = args.reactions_file
run_type = args.run_type
output_file = args.output_file

T = args.T
t_max = args.t_max
N_data = args.N_data
DRCprod = args.DRCprod
pct = args.percent_change
N_timesteps = args.time_steps
atol = args.atol
rtol = args.rtol

if output_file is None:
    if run_type == 'kinetics':
        output_file = output_kinetics
    elif run_type == 'DRC':
        output_file = f'DRCs/DRCs_{DRCprod}.csv'
    elif run_type == 'distribution':
        output_file = output_distribution
    else:
        try:
            print("Something went wrong with the output file names")
        except:
            sys.exit(1)


if not os.path.isfile(reactionsfile):
    try:
        print(f'Error: {reactionsfile} does not exist!')
    except:
        sys.exit(1)
else:
    print(f'Reactions are read from {reactionsfile}')

init_species = []
init_concs = []

with open(reactionsfile) as fi:       #this is to extract intial concentrations
    for l in fi.readlines():
        l = l.strip()
        if l and l[0] == "[":       #detect inital concentration lines  
            init_species.append(l.split('[', 1)[1].split(']')[0].strip())
            init_concs.append(float(l.split('=', 1)[1].strip()))
    
    if not init_concs:
        try:
            print("Error: no intial conditions!")
        except:
            sys.exit(1)

initial_species = dict(zip(init_species,init_concs))

all_species_names=[]    #list of string with names, to check if a species has already shown up
species=[]              #list of Species objects
reactions = []          #list of Reaction objects
DeltaGr = []
DeltaGTS = []
i=0                     #reaction index
sp_idx = 0              #species index

with open(reactionsfile) as fi:
    for line in fi.readlines():
        line = line.strip()
        if not line or line[0] in '#%![':   
            # skip blank and comment lines, and initial conditions lines
            continue

        reactants,r_stoic,products,p_stoic,G_TS,G_r,ref = parse_line(line)

        DeltaGr.append(G_r)
        DeltaGTS.append(G_TS)

        kfw = Eyring(G_TS,T)        #compute forward kinetic constant
        kbw = Eyring(G_TS-G_r,T)    #compute backward kinetic constant

        #append Reaction object
        reactions.append(Reaction(num=i,DeltaG_TS=G_TS,DeltaG_reaction=G_r,kinetic_constant_fw=kfw,kinetic_constant_bw=kbw,\
                                  reactants=[],r_stoic=[],products=[],p_stoic=[],conc_product_fw=0.0,conc_product_bw=0.0,ref=ref))

        for reac in reactants:
            if reac not in all_species_names:   #check if it is new
                if reac in initial_species:     #check if it has an initial concentration
                    new_species = Species(sp_idx,reac,conc=initial_species[reac])  
                else:
                    new_species = Species(sp_idx,reac) #conc = 0.0 by default
                species.append(new_species)     #append Species to species array
                all_species_names.append(reac)
                sp_idx +=1

            for sp in species:
                if sp.name == reac:
                    reactions[i].reactants.append(sp)  

        for prod in products:
            if prod not in all_species_names:   #check if it is new
                if prod in initial_species:     #check if it has an initial concentration
                    new_species = Species(sp_idx,prod,conc=initial_species[prod])
                else:
                    new_species = Species(sp_idx,prod) #conc = 0.0 by default
                species.append(new_species)     #append Species to species array
                all_species_names.append(prod)
                sp_idx +=1
        
            for sp in species:
                if sp.name == prod:
                    reactions[i].products.append(sp)
        
        reactions[i].r_stoic = r_stoic
        reactions[i].p_stoic = p_stoic

        i+=1


for rxn in reactions:   #Compute concentration products
    conc_product_fw = 1     #see above. Initialized to 1
    conc_product_bw = 1

    for rn in range(len(rxn.reactants)):
        if rxn.reactants[rn].conc != 0: #if either one reactant is zero, then 
            conc_product_fw *= rxn.reactants[rn].conc ** rxn.r_stoic[rn]    #e.g. 2A -> B, r =  k * [A]**2
        else:
            conc_product_fw = 0

    for pn in range(len(rxn.products)):
        if rxn.products[pn].conc != 0:
            conc_product_bw *= rxn.products[pn].conc ** rxn.p_stoic[pn]     #e.g. 2A -> B, r =  k * [A]**2 
        else:
            conc_product_bw = 0

    rxn.conc_product_fw = conc_product_fw
    rxn.conc_product_bw = conc_product_bw


N_reactions = len(reactions)
N_species = len(all_species_names)

for rxn in reactions:   #Compute derivatives
    for n_r in range(len(rxn.reactants)):   
        rxn.reactants[n_r].derivative -= rxn.kinetic_constant_fw * rxn.conc_product_fw * rxn.r_stoic[n_r]   #e.g.  2A -> B. r = -2d[A]/dt = d[B]/dt
        rxn.reactants[n_r].derivative += rxn.kinetic_constant_bw * rxn.conc_product_bw * rxn.r_stoic[n_r]

    for n_p in range(len(rxn.products)):
        rxn.products[n_p].derivative += rxn.kinetic_constant_fw * rxn.conc_product_fw * rxn.p_stoic[n_p]
        rxn.products[n_p].derivative -= rxn.kinetic_constant_bw * rxn.conc_product_bw * rxn.p_stoic[n_p]


PrintReadReactions()

#CalcSpeciesEnergies()

if __name__ == "__main__":
    print(N_reactions)

    for sp in species:
        print(sp.name,sp.conc)
        
    for rxn in reactions:
        print("rxn", rxn.num,rxn.conc_product_fw,rxn.conc_product_bw)
        for j in range(len(rxn.reactants)):
            print(rxn.r_stoic[j],rxn.reactants[j].name,rxn.reactants[j].conc,rxn.reactants[j].derivative)

        print(">")

        for k in range(len(rxn.products)):
            print(rxn.p_stoic[k],rxn.products[k].name,rxn.products[k].conc,rxn.products[k].derivative)
        print("------")

