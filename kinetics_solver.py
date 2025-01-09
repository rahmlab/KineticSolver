import numpy as np
from scipy.integrate import ode
from utils import *


def GetDerivatives(t,c):
    UpdateDerivatives(c)
    return np.array([sp.derivative for sp in species])


def SolveODE(init_conc,Delta_t,return_all_t=True):
    dt = Delta_t/N_timesteps
    times=np.empty(N_timesteps+1)
    Concs=np.empty(shape=(N_timesteps+1,N_species))
    o = ode(GetDerivatives)
    
    o.set_integrator('lsoda',method='bdf',rtol=rtol, atol=atol)
    #o.set_integrator('dopri5')
    #o.set_integrator('rk8pd')
    o.set_initial_value(init_conc)

    Concs[0] = init_conc
    times[0] = 0.0
    
    i=0

    while o.successful() and o.t < Delta_t-dt/100:
        i+=1
        Concs[i] = o.integrate(o.t+dt)
        times[i] = o.t

    if o.successful():
        UpdateDerivatives(Concs[-1])
        if return_all_t:
            return times, Concs
        else:
            return times[-1], Concs[-1]
    else:
        return None, None



