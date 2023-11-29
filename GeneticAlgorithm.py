from inspect import _void
import numpy as np
import random
from deap import base, creator, tools, algorithms

class GeneticAlgorithmOptimization:
           
    # GA setup

    ngen = 10 # number of generations
    nind = 10 # number of individuals
    eta = 5.0 
    MUTPB = 0.5 # probability of mutation
    CXPB = 0.5 # probability of crossover
    hof = 5 # hall of fame size
    
    # control variable bounds

    alpha_bounds = np.array([0, 2 * np.pi / 180]) # AoA bounds [rad] 
    alpha_min = alpha_bounds[0]
    alpha_max = alpha_bounds[1]
    
    # constraints

    q_max = 45 # max. stag. heat [W/cm^3]

    def __init__(self):

        alpha_min = self.alpha_min
        alpha_max = self.alpha_max
        eta = self.eta
        
        creator.create("FitnessMin", base.Fitness, weights=(-1.0))
        creator.create("Individual", list, fitness=creator.FitnessMin)
        
        self.toolbox = base.Toolbox()
        toolbox = self.toolbox
        evaluate = self.evalulate        
        check_feasibility = self.check_feasibility
        penalty = self.penalty

        toolbox.register("alpha", random.uniform, alpha_min, alpha_max)
        toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.alpha, n = 1)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        toolbox.register("evaluate", evaluate)
        toolbox.register("mate", tools.cxSimulatedBinaryBounded, low=alpha_min, up=alpha_max, eta=eta)
        toolbox.mutate("mutate", tools.mutPolynomialBounded, low=alpha_min, up=alpha_max, eta=eta, indpb=0.08)
        toolbox.register("select", tools.selNSGA2)    
        toolbox.decorate("evaulate", tools.DeltaPenalty(check_feasibility, 1000, penalty))
        
    def solve(self, alpha):
      
        rho0 = self.atm_params[0]
        beta = 6700 # scale height [m]
        z = self.x[0] # height [m]
        rho = rho0 * np.exp(-beta * z) # density at height [kg/m^3]

        V = self.x[1]

        R0 = self.sc_params[0]
        S = self.sc_params[1]

        """
        Cd = 2 * pow(np.sin(alpha), 3) # Newtonian theory coeff. of drag

        D = 0.5 * rho * pow(V, 2) * Cd * S

        Cl = 2 * pow(np.sin(alpha), 2) * np.cos(alpha) # Newtonian theory coeff. of lift

        L = 0.5 * rho * pow(V, 2) * Cl * S
        """
        
        Cp = 2 * pow(np.sin(alpha), 2) # coeff. of pressure Newtonian theory
        q_dyn = 0.5 * rho * pow(V, 2) # dynamic pressure [Pa]
        
        N = Cp * q_dyn * S # normal force [N]

        q = self.k * np.sqrt((self.rho / self.R0)) * pow(self.V, 3) # Sutton-Graves stagnation point heat flux approximation [W/m^2]
        
        # add other constraints here if needed
        
        return [alpha, N, q]

    def getSolution(self, x, atm_params, sc_params):
         
        # Inputs x = [h, ]

        [h, v] = x[-1]
        alpha = 1
        [R0, S] = sc_params

        hof = self.solve() 

        alpha = hof[0][0]
        
        sol = self.solve(alpha)

        return sol
    
    def optimize(self):
        
        toolbox = self.toolbox
        pop = toolbox.population(n = self.nind)
        NGEN = self.ngen
        CXPB = self.CXPB
        MUTPB = self.MUTPB    
        hof = self.hof

        algorithms.eaSimple(pop, toolbox, cxpb=CXPB, mutpb=MUTPB, ngen=NGEN, halloffame=hof, verbose=True)
        
        return hof
    
    def check_feasibility(self, params):
        
        alpha = params[0]

        q = self.getSolution(alpha)[3]

        if q > self.q_max:
            return False
        else:
            return True
        
    def evaluate(self):
        
        weights = [0.5, 0.5]
        
    def penalty(self):
        
        return True