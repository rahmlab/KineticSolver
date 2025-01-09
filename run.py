from kinetics_solver import *
from utils import *
import numpy as np
import sys
from datetime import datetime


t_max_arr = powers_of_tens_array(t_max,step=1.0)
t0 = 0

N_data_per_timestep = int(N_data/len(t_max_arr))
N_data_eff = N_data_per_timestep * len(t_max_arr) + 1

indexes = AlmostLogSpacedIndexes(N_data_per_timestep)

print(indexes)

ResetAll()
#ResetAll(normaldist=True,scale=1.0)

init_conc = np.empty(N_species)
for sp_idx in range(N_species):
    init_conc[sp_idx] = species[sp_idx].conc

match run_type:
    case 'kinetics':    #-------------------------------------------------------

        print(f'Performing kinetics simulation up to t = {t_max}')    
        print(f'Number of timesteps to solve ode is {N_timesteps}')
        print(f'A number of {N_data_eff} data points will be printed in \'{output_file}\'')
        print(f'Tolerances are {atol}, {rtol}')

        times = np.empty(N_data_eff)
        concs = np.empty(shape=(N_data_eff,N_species))

        k=0

        header = "t," + ",".join(all_species_names)
        
        for curr_t_max in t_max_arr:
            print(f't is {curr_t_max}, t0 is {t0}')
            t, c = SolveODE(init_conc,curr_t_max,return_all_t=True)

            if c is None:
                print("FATAL ERROR: ODE solver failed")
                print("No output has been produced")
                sys.exit(1)

            for idx in indexes:
                times[k] = t[idx] + t0
                concs[k] = c[idx]
                k += 1

            t0 = t[-1]+t0
            init_conc = c[-1]


        times[-1] = t0
        concs[-1] = init_conc


        np.savetxt(output_file, np.column_stack([times,concs]), delimiter=',', header = header, comments='')


    case 'DRC':    #-------------------------------------------------------

        print(f'Performing DRC calculations of all {N_reactions} reactions with pct = {pct}')
        print(f'Kinetics simulation are run up to t = {t_max}')
        print(f'Number of timesteps to solve ode is {N_timesteps-1}')
        print(f'A number of {N_data_eff} data points will be printed in \'{output_file}\'')

        header = "t," + ",".join(f'X_rc{str(x)}' for x in range(N_reactions))

        DRCs = np.empty((N_data_eff,N_reactions))
        rate_unp = np.empty(N_data_eff)
        rates = np.empty((N_data_eff,N_reactions))
        DRCprodconc_unp = np.empty(N_data_eff)
        DRCprodconcs = np.empty((N_data_eff,N_reactions))

        rxn_list = []

        for rxn in reactions:
            for prod in rxn.products:
                if prod.name == DRCprod:
                    rxn_list.append(rxn.num)

        if not rxn_list:
            try:
                print("ERROR: no reactions to form desired product!")
            except:
                sys.exit(1)
        else:
            print(f'product is formed by reactions {rxn_list}')

        for sp in species:
            if sp.name == DRCprod:
                    DRCprod_idx = sp.index


        for i in range(-1,N_reactions):
            ResetAll()
            ResetKineticConstants()

            t0 = 0

            init_conc = np.empty(N_species)
            for sp_idx in range(N_species):
                init_conc[sp_idx] = species[sp_idx].conc

            times = np.empty(N_data_eff)
            concs = np.empty(shape=(N_data_eff,N_species))

            k = 0

            if i >= 0:     #Nraise increases the barrier of reaction i. It increases k_fw and k_bw by pct%
                reactions[i].kinetic_constant_fw *= 1. + pct/100
                reactions[i].kinetic_constant_bw *= 1. + pct/100
          

            for curr_t_max in t_max_arr:

                t, c = SolveODE(init_conc,curr_t_max,return_all_t=True)

                for idx in indexes:
                    times[k] = t[idx] + t0
                    concs[k] = c[idx]
                    k += 1
                
                t0 = t[-1]+t0
                init_conc = c[-1]
            
            times[-1] = t0
            concs[-1] = init_conc

            if  i == -1:
                rate_unp = CalcRate(concs,rxn_list)
                #np.savetxt(f'DRCs/rate_unp.csv', rate_unp, delimiter=',', header = f't,rate_unp', comments='')
                #c_unp = concs[:,DRCprod_idx]
                DRCprodconc_unp = concs[:,DRCprod_idx]
            else:
                rates[:,i] = CalcRate(concs,rxn_list)
                DRCs[:,i] = CalcDRC(rates[:,i],rate_unp,pct) 
                DRCprodconcs[:,i] = concs[:,DRCprod_idx]
                #DRCs[:,i] = CalcDRCfromConc(concs[:,DRCprod_idx],c_unp,pct) 
                #np.savetxt(f'DRCs/rate_r{i}.csv', rate_curr, delimiter=',', header = f't,rate_{i}', comments='')
            
            #header_species = "t," + ",".join(all_species_names)
            #with open("DRCs/progress.csv", "ab") as pf:
            #    np.savetxt(pf, [i])

        
        #np.savetxt(f'{output_file}_rates', np.column_stack([times,rate_unp,rates]), delimiter=',', header = header, comments='')
        np.savetxt(f'{output_file}_concs', np.column_stack([times,DRCprodconc_unp,DRCprodconcs]), delimiter=',', header = header, comments='')
        np.savetxt(output_file, np.column_stack([times,DRCs]), delimiter=',', header = header, comments='')



    case 'distribution':    #-------------------------------------------------------

        filename = os.path.splitext(output_file)[0]

        t_max_arr = powers_of_tens_array(t_max,step=1.0)

        print(f'Performing concentration distribution calculations with standard deviation of 1')
        print(f'Number of timesteps to solve ode is {N_timesteps-1}')
        print(f'A number of {N_data_eff} data points will be printed in \'{filename}_t-N.csv\' files, with N = log10(t_max)')
        print(f'Automated calculation of means and std.dev is work in progress.')

        header = ",".join(all_species_names) + ",k3_fw,k3_bw,k5_fw,k5_bw,k10_fw,k10_bw,k13_fw,k13_bw,k22_fw"

        for idx in range(N_data_eff):
            curr_filename=f'{filename}-wks_t-{idx}.csv'
            if not os.path.isfile(curr_filename):
                with open(curr_filename, 'ab') as fd:
                    np.savetxt(fd, [], header=header, comments='')

        for i in range(100):
            ResetAll(normaldist=True,scale=1.0)

            k5_fw = reactions[5].kinetic_constant_fw
            k5_bw = reactions[5].kinetic_constant_bw
            k13_fw = reactions[13].kinetic_constant_fw
            k13_bw = reactions[13].kinetic_constant_bw
            k22_fw = reactions[22].kinetic_constant_fw

            k3_fw = reactions[3].kinetic_constant_fw
            k3_bw = reactions[3].kinetic_constant_bw
            k10_fw = reactions[10].kinetic_constant_fw
            k10_bw = reactions[10].kinetic_constant_bw

            init_conc = [sp.conc for sp in species]
            t0=0
            times = np.empty(N_data_eff)
            concs = np.empty(shape=(N_data_eff,N_species))
            k=0

            for curr_t_max in t_max_arr:
                t, c = SolveODE(init_conc,curr_t_max,return_all_t=True)

                if c is None:
                    print(f"FATAL ERROR: ODE solver failed at time {t0}, iteration {i}")
                    print("No output has been produced")
                    break

                for idx in indexes:
                    times[k] = t[idx] + t0
                    concs[k] = c[idx]
                    k += 1

                t0 = t[-1]+t0
                init_conc = c[-1]

                print(t0)
                c = np.empty(N_species)

            times[-1] = t0
            concs[-1] = init_conc

            for idx in range(N_data_eff):
                curr_filename=f'{filename}-wks_t-{idx}.csv'
                with open(curr_filename, 'ab') as fd:
                    np.savetxt(fd, np.reshape(np.append(concs[idx,:],[k3_fw,k3_bw,k5_fw,k5_bw,k10_fw,k10_bw,k13_fw,k13_bw,k22_fw]),(1,N_species+9)),delimiter=',')


            #np.savetxt(curr_filename, concs, delimiter=',', header = header, comments='')



        """
        for curr_t_max in t_max_arr:
            curr_filename=f'{filename}_t-{int(100*np.log10(curr_t_max))}.csv'
            if not os.path.isfile(curr_filename):
                with open(curr_filename, 'ab') as fd:
                    np.savetxt(fd, [], header=header, comments='')

        for i in range(100):
            ResetAll(normaldist=True,scale=np.sqrt(3))

            init_conc = [sp.conc for sp in species]
            t0=0

            for curr_t_max in t_max_arr:

                t, c = SolveODE(init_conc,curr_t_max,return_all_t=False)
                
                if c is None:
                    print(f"FATAL ERROR: ODE solver failed at time {t0}, iteration {i}")
                    print(f"Output has been interrupted.")
                    break

                curr_filename=f'{filename}_t-{int(100*np.log10(curr_t_max))}.csv'

                with open(curr_filename, 'ab') as fd:
                    np.savetxt(fd, c.reshape(1,c.shape[0]), delimiter=',')

                t0 = t+t0
                init_conc = c[:]   

                print(t0)

                c = np.empty(N_species)

                """
