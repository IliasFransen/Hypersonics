import numpy as np
from GeneticAlgorithm import GeneticAlgorithmOptimization

# Atmosphere parameters
rho0 = 1.225 # density of air at sea level [kg/m^3]

atm_params = [rho0]

# Spacecraft parameters (temp values)
R0 = 0.25 # nose radius (m)
S = 0.8 # reference area [m^2]

sc_params = [R0, S]

# ICs (temps values)
h0 = 100000 # initial altitude [m]
Vx0 = 100000 # initial x velocity [m/s]
Vy0 = 100000 # initial y velocity [m/s]

ics = [Vx0, Vy0, h0] 

# control var

alpha = 0 # AoA [rad], to be optimized

if __name__ == "__main__":
    
    GA = GeneticAlgorithmOptimization()
    
    t = np.arange(0, 10, 1)
    
    x = []
    x.append(ics)

    for i in range(0, len(t)):
        
        x = ics[-1]

        opt = GA.getSolution() # run genetic algorithm to get optimum
                
        alpha = opt[0]
        N = opt[1]
        
        # ODE stuff with optimum

        sol = []
        
        x.append(sol)
        