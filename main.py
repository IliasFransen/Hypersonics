import numpy as np
import matplotlib.pyplot as plt
from GeneticAlgorithm import GeneticAlgorithmOptimization
#from State_System[TEMP] import StateUpdate
from scipy.integrate import solve_ivp
from math import atan2

# Atmosphere parameters
rho0 = 1.225 # density of air at sea level [kg/m^3]

atm_params = [rho0]

# Spacecraft parameters (temp values)
R0 = 0.25 # nose radius (m)
S = 0.8 # reference area [m^2]
m = 2000 # mass [kg]

sc_params = [R0, S]

# ICs (temps values)
h0 = 80000 # initial altitude [m]
Vx0 = 0 # initial x velocity [m/s]
Vy0 = -1800 # initial y velocity [m/s]

ICs = [Vx0, Vy0, h0] # initial conditions vector 

# control var

alpha = 0 # AoA [rad], to be optimized

g = 9.81

def EOMs(t, x, *args):
    
    alpha, N = args
    
    [Vx, Vy, h] = x
    
    Beta = atan2(Vy, Vx)
    
    eta = Beta - alpha

    dVx_dt = N * np.cos(eta - np.pi) / m
    dVy_dt = -g + N * np.sin(eta - np.pi) / m
    dh_dt = Vy
    
    return [dVx_dt, dVy_dt, dh_dt]

if __name__ == "__main__":
    
    GA = GeneticAlgorithmOptimization()
    
    t = np.arange(0, 10, 1)
    
    x = []
    x.append(ICs)

    for i in range(0, len(t) - 1):
        
        t1 = t[i]
        t2 = t[i+1]
        tspan = [t1, t2]

        opt = GA.getSolution(x, alpha, atm_params, sc_params) # run genetic algorithm to get optimum
        
        alpha = opt[0]
        N = opt[1]
        
        # ODE stuff with optimum

        x0 = x[-1] # initialize with previous solution
        
        args = [alpha, N]

        sol = solve_ivp(EOMs, tspan, x0, args=args)
        
        x.append( [sol.y[0][-1], sol.y[1][-1], sol.y[2][-1]] ) 
            
        
    x = list(zip(*x))
    
    Vx = x[0]
    Vx2 = tuple(x**2 for x in Vx)
    Vy = x[1]
    Vy2 = tuple(x**2 for x in Vy)
    V = tuple(np.sqrt(x + y) for x,y in zip(Vx2, Vy2))
    h = x[2]

    plt.figure(1)
    plt.plot(t, h)
    
    plt.figure(2)
    plt.plot(V, h)
    plt.show()