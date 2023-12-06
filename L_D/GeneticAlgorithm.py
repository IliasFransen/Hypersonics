from math import isnan
import numpy as np
import random
from deap import base, creator, tools, algorithms
from Atmosphereic_Conditions import Get_Density
from Modified_Newton import Get_Lift, Get_Drag, Get_Tangential, Get_Normal
from State_System import StateUpdate
from scipy.integrate import odeint


class GeneticAlgorithmOptimization:
    # GA setup

    ngen = 3  # number of generations
    nind = 12  # number of individuals
    eta = 5.0  # SBX crossover operator
    mutpb = 0.5  # probability of mutation
    cxpb = 0.6  # probability of crossover
    nhof = 2  # hall of fame size

    # control variable bounds

    alpha_bounds = np.array([-10, 10]) * np.pi / 180  # AoA bounds [rad]
    alpha_min = alpha_bounds[0]
    alpha_max = alpha_bounds[1]

    # constraints

    q_max = 2.5 * 10 ** 6  # max. stag. heat [W/m^2]
    ng_max = 6.0  # max. g-load [g]

    # Sutton-Graves stagnation point heat transfer coefficient for earth
    k = 1.7415 * 10 ** (-4)

    def __init__(self):

        random.seed(64)

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
        toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.alpha, n=1)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        toolbox.register("evaluate", evaluate)
        toolbox.register("mate", tools.cxSimulatedBinaryBounded, low=alpha_min, up=alpha_max, eta=eta)
        toolbox.register("mutate", tools.mutPolynomialBounded, low=alpha_min, up=alpha_max, eta=eta, indpb=0.2)
        toolbox.decorate("evaluate", tools.DeltaPenalty(check_feasibility, 1E10))
        toolbox.register("select", tools.selNSGA2)

    def solve(self, alpha, V, h, fpa):

        g = self.atm_params[0]
        gamma = self.atm_params[1]

        rho = Get_Density(h)

        R0 = self.sc_params[0]
        m = self.sc_params[1]
        x_lst = self.sc_params[2]
        y_lst = self.sc_params[3]

        N = Get_Normal(V, h, gamma, x_lst, y_lst, alpha)
        T = Get_Tangential(V, h, gamma, x_lst, y_lst, alpha)

        L = Get_Lift(N, T, V, h, alpha, x_lst, y_lst)
        D = Get_Drag(N, T, V, h, alpha, x_lst, y_lst)

        # constraints
        ng = np.sqrt(L** 2 + D** 2) / (m * g)  # load deceleration [g's]

        q = self.k * np.sqrt((rho / R0)) * pow(V, 3)  # Sutton-Graves stagnation point heat flux approximation [W/m^2]

        return [alpha, L, D, q, ng]

    def getSolution(self, x, atm_params, sc_params, tspan):

        self.x = x[-1]
        self.sc_params = sc_params
        self.atm_params = atm_params
        self.tspan = tspan

        hof = self.optimize()

        alpha = hof[0][0]

        g = self.atm_params[0]
        gamma = self.atm_params[1]

        m = self.sc_params[1]
        x_lst = self.sc_params[2]
        y_lst = self.sc_params[3]

        x0 = self.x

        args = (g, m, alpha, gamma, x_lst, y_lst)
        sol = odeint(StateUpdate, x0, tspan, args=args)

        V = sol[-1, 0]
        fpa = sol[-1, 1]
        h = sol[-1, 2]

        opt = self.solve(alpha, V, h, fpa)

        return opt, sol

    # Evolutionary algorithm

    def optimize(self):
        
        toolbox = self.toolbox
        pop = toolbox.population(n=self.nind)
        NGEN = self.ngen
        CXPB = self.cxpb
        MUTPB = self.mutpb
        
        """
        hof = tools.HallOfFame(self.nhof)
        
        algorithms.eaSimple(pop, toolbox, cxpb=CXPB, mutpb=MUTPB, ngen=NGEN, halloffame=hof, verbose=False)
        
        return hof
        """
        
        pareto = tools.ParetoFront()       
        pareto.clear()

        invalid_ind = [ind for ind in pop if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
            
        pop = toolbox.select(pop, len(pop))
        
        for gen in range(1, NGEN):

            offspring = tools.selTournamentDCD(pop, len(pop))
            offspring = [toolbox.clone(ind) for ind in offspring]
            
            for ind1, ind2 in zip(offspring[::2], offspring[1::2]):
                
                if random.random() <= CXPB:
                    toolbox.mate(ind1, ind2)
                    del ind1.fitness.values, ind2.fitness.values
            
                for mutant in offspring:
                    if random.random() < MUTPB:
                        toolbox.mutate(mutant)
                        del mutant.fitness.values
                    
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit
            
            pop = toolbox.select(pop + offspring, self.nind)
            
        pareto.update(pop)
        
        return pareto
        
        
    # Check feasibility of solution

    def check_feasibility(self, params):

        alpha = params[0]

        g = self.atm_params[0]
        gamma = self.atm_params[1]

        m = self.sc_params[1]
        x_lst = self.sc_params[2]
        y_lst = self.sc_params[3]

        tspan = self.tspan

        x0 = self.x

        args = (g, m, alpha, gamma, x_lst, y_lst)

        sol = odeint(StateUpdate, x0, tspan, args=args)

        V = sol[-1, 0]
        fpa = sol[-1, 1]
        h = sol[-1, 2]

        opt = self.solve(alpha, V, h, fpa)

        q = opt[3]
        ng = opt[4]

        if isnan(q):
            return False
        if ng > self.ng_max:
            print('ng failed: ng = %.5f' % ng)
            return False
        if q > self.q_max:
            print('q failed: q = %.5f' % q)
            return False
        else:
            print('pass')
            return True

    # Evaluate solution with objective function

    def evaluate(self, params):

        alpha = params[0]

        g = self.atm_params[0]
        gamma = self.atm_params[1]

        m = self.sc_params[1]
        x_lst = self.sc_params[2]
        y_lst = self.sc_params[3]

        x0 = self.x
        tspan = self.tspan

        args = (g, m, alpha, gamma, x_lst, y_lst)

        sol = odeint(StateUpdate, x0, tspan, args=args)

        V = sol[-1, 0]
        fpa = sol[-1, 1]
        h = sol[-1, 2]

        weights = [0.4, 0.6]

        opt = self.solve(alpha, V, h, fpa)

        q = opt[3]
        ng = opt[4]

        objectives = [q, ng]

        objfcn = sum(x * y for x, y in zip(weights, objectives)) / sum(weights)  # weighted sum average objective function
        
        return objfcn,
