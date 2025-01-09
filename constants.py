
#Physics constants -----
Na = 6.02214076E23  # Avogadro's number mol^-1
kB = 1.380649E-23           # Boltzmann constant in J/K
R =  6.02214076 * 1.380649    # Gas constant in J/(mol*K)
R_kcal = R / 4184    # Gas constant in kcal/(mol*K)
h = 6.62607015E-34  # Planck's Constant in J*s


#Default values -------

#Temperature
T = 278.0

#Number of timesteps to solve ODE
N_timesteps = 1000000 

#Max time of the simulation: must be a power of 10 for now!
t_max = 1e12

#Number of data points that will be printed (<= N_timesteps)
#This is not precise, and it may vary by +/- log10(t_max)
N_data = 100

#Product for which compute DRCs
DRCprod = 'Ade'

#Path and name of the outputfile
output_kinetics='results/concs-ode_S-1e5.csv'
output_DRC=f'DRCs/DRCs_{DRCprod}.csv'
output_distribution='distributions/concs_distribution.csv'

#Path and name of the input file, where to read the reactions
reactionsfile = 'reactions.txt'

#Change of kinetics constant for DRC, in percent
pct = 1.0

#Tolerance to solve (default for odeint is ~1.5e-8, but use at least 1e-10 if you have a lot of reactions)
atol = 1e-9
rtol = 1e-8

