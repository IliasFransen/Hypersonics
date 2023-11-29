import numpy as np
import random
from deap import base, creator, tools, algorithms
from Atmosphereic_Conditions import Get_Density
from Modified_Newton import Get_Normal

class GeneticAlgorithmOptimization:
           
    # GA setup

    ngen = 5 # number of generations
    nind = 5 # number of individuals
    eta = 5.0 
    mutpb = 0.5 # probability of mutation
    cxpb = 0.5 # probability of crossover
    nhof = 5 # hall of fame size
    
    # control variable bounds

    alpha_bounds = [0, 2 * np.pi / 180] # AoA bounds [rad] 
    alpha_min = alpha_bounds[0]
    alpha_max = alpha_bounds[1]
    
    # constraints

    q_max = 45 * 10**6 # max. stag. heat [W/m^2]
    
    # Sutton-Graves stagnation point heat transfer coefficient for earth
    k = 1.7415*10**(-4)

    def __init__(self):

        alpha_min = self.alpha_min
        alpha_max = self.alpha_max
        eta = self.eta
        
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMin)
        
        self.toolbox = base.Toolbox()
        toolbox = self.toolbox
        evaluate = self.evaluate        
        check_feasibility = self.check_feasibility

        toolbox.register("alpha", random.uniform, alpha_min, alpha_max)
        toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.alpha, n = 1)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        toolbox.register("evaluate", evaluate)
        toolbox.register("mate", tools.cxSimulatedBinaryBounded, low=alpha_min, up=alpha_max, eta=eta)
        toolbox.register("mutate", tools.mutPolynomialBounded, low=alpha_min, up=alpha_max, eta=eta, indpb=0.08)
        toolbox.register("select", tools.selNSGA2)    
        toolbox.decorate("evaluate", tools.DeltaPenalty(check_feasibility, 1000))
        
    def solve(self, alpha):
      
        """
        rho0 = self.atm_params[0]
        beta = 6700 # scale height [m]
        z = self.x[0] # height [m]
        rho = rho0 * np.exp(-beta * z) # density at height [kg/m^3]
        """
        
        h = self.x[2]

        rho = Get_Density(h).tolist()[0]
        
        Vx = self.x[0]
        Vy = self.x[1]

        R0 = self.sc_params[0]
        S = self.sc_params[1]

        """
        Cd = 2 * pow(np.sin(alpha), 3) # Newtonian theory coeff. of drag

        D = 0.5 * rho * pow(V, 2) * Cd * S

        Cl = 2 * pow(np.sin(alpha), 2) * np.cos(alpha) # Newtonian theory coeff. of lift

        L = 0.5 * rho * pow(V, 2) * Cl * S
        """
        
        # Some BS calculations i know they arent right was just testing GA
        Cp = 2 * pow(np.sin(alpha), 2) # coeff. of pressure Newtonian theory
        V = np.sqrt(Vx**2 + Vy**2) # velocity magnitude
        q_dyn = 0.5 * rho * pow(V, 2) # dynamic pressure [Pa]
        
        N = Cp * q_dyn * S # normal force [N]

        
        #N = Get_Normal(Vx, Vy, h, gamma, x_lst, y_lst, alpha)
        
        q = self.k * np.sqrt((rho / R0)) * pow(V, 3) # Sutton-Graves stagnation point heat flux approximation [W/m^2]
        
        # add other constraints here if needed
        
        return [alpha, N, q]

    def getSolution(self, x, alpha, atm_params, sc_params):
         
        # state x = [h, Vx, Vy]

        [Vx, Vy, h] = x[-1]
        
        [R0, S] = sc_params
        
        [rho0] = atm_params
        
        self.x = x[-1]
        self.sc_params = sc_params
        self.atm_params = atm_params

        hof = self.optimize() 

        alpha = hof[0][0]
        
        sol = self.solve(alpha)

        return sol
    
    # Evolutionary algorithm

    def optimize(self):
        
        toolbox = self.toolbox
        pop = toolbox.population(n = self.nind)
        NGEN = self.ngen
        CXPB = self.cxpb
        MUTPB = self.mutpb    
        hof = tools.HallOfFame(self.nhof)

        algorithms.eaSimple(pop, toolbox, cxpb=CXPB, mutpb=MUTPB, ngen=NGEN, halloffame=hof, verbose=True)
        
        return hof
    
    # Check feasibility of solution

    def check_feasibility(self, params):
        
        alpha = params[0]

        q = self.solve(alpha)[2]

        if q > self.q_max:
            return False
        else:
            return True
        
    # Evaluate solution with objective function

    def evaluate(self, params):
        
        alpha = params[0]

        weights = [0.5]
        
        sol = self.solve(alpha)
        alpha = sol[0]
        q = sol[2]

        objectives = [q]
        
        objfcn = sum(x * y for x, y in zip(weights, objectives)) / sum(weights) # weighted sum average objective function
        
        return objfcn,
        