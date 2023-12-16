from asyncio.windows_events import NULL
import numpy as np
import random
from deap import base, creator, tools
from Atmosphereic_Conditions import Get_Density, Get_SoundSpeed
from Modified_Newton import Get_Lift, Get_Drag
from State_System import StateUpdate, getGForce, getStagHeatFlux, getStagHeatLoad
from scipy.integrate import odeint


class GlobalGeneticAlgorithmOptimization:
    # GA setup

    ngen = 3 # number of generations
    npop = 12  # number of populations
    eta = 10.0  # SBX crossover operator
    mutpb = 0.1  # probability of mutation
    cxpb = 0.9 # probability of crossover

    # constraints

    q_max = 2.5E6 # max. stagnation point heat flux [W/m^2]
    Q_max = 5.1E8 # max. stagnation point heat load [J/m^2]
    ng_max = 8.0  # max. g-load [g]

    def __init__(self):

        random.seed(64)

        creator.create("FitnessMulti", base.Fitness, weights=(-1.0, -1.0, -1.0))
        creator.create("Individual", list, fitness=creator.FitnessMulti)

        self.toolbox = base.Toolbox()
        toolbox = self.toolbox
        evaluate = self.evaluate
        check_feasibility = self.check_feasibility

        toolbox.register("evaluate", evaluate)
        toolbox.decorate("evaluate", tools.DeltaPenalty(check_feasibility, 1E10))
        toolbox.register("select", tools.selNSGA2)

    def solve(self, alphas, t):
        
        x0 = self.x
        
        g = self.atm_params[0]
        gamma = self.atm_params[1]
        Pr = self.atm_params[2]

        R0 = self.sc_params[0]
        m = self.sc_params[1]
        x_lst = self.sc_params[2]
        y_lst = self.sc_params[3]
        S = self.sc_params[4]
        
        L = np.zeros(len(t))
        D = np.zeros(len(t))
        ng = np.zeros(len(t))
        q = np.zeros(len(t))
        Q = np.zeros(len(t))
        V = np.zeros(len(t))
        fpa = np.zeros(len(t))
        h = np.zeros(len(t))
        AoA = np.zeros(len(t))
        M = np.zeros(len(t))
        s = np.zeros(len(t))

        for i in range(1, len(t) - 1):

            x_sol = x0[-1]

            args = (g, m, gamma, x_lst, y_lst, S, alphas[i])
            sol = odeint(StateUpdate, x_sol, [t[i], t[i + 1]], args=args)

            V[i] = sol[-1, 0]
            fpa[i] = sol[-1, 1]
            h[i] = sol[-1, 2]
            s[i] = sol[-1,3]
            
            rho = Get_Density(h[i])

            AoA[i] = alphas[i]
        
            L[i] = Get_Lift(V[i], h[i], gamma, x_lst, y_lst, S, alphas[i])
            D[i] = Get_Drag(V[i], h[i], gamma, x_lst, y_lst, S, alphas[i])

            dVdt = StateUpdate([V[i], fpa[i], h[i], s[i]], NULL, g, m, gamma, x_lst, y_lst, S, alphas[i])[0]

            a = Get_SoundSpeed(sol[-1, 2])
            M[i] = sol[-1, 0] / a   # Mach number

            # constraints
            ng[i] = getGForce(dVdt, g)  # load deceleration [g's]

            q[i] = getStagHeatFlux(h[i], M[i], gamma, Pr, R0)  # Sutton-Graves stagnation point heat flux approximation [W/m^2]
            
            Q[i] = getStagHeatLoad(q, t, i)  # heat load
            
            x0.append([V[i], fpa[i], h[i], s[i]])

            # Break solver loop if Mach < 2
            if M[i] < 3:
                t = t[0:i + 1]
                q = q[0:i + 1]
                L = L[0:i + 1]
                D = D[0:i + 1]
                Q = Q[0:i + 1]
                ng = ng[0:i + 1]
                V = V[0:i + 1]
                h = h[0:i + 1]
                AoA = AoA[0:i + 1]
                M = M[0:i + 1]
                s = s[0:i + 1]
                self.x = [x0[0]]
                break
            
        self.x = [x0[0]]
        return [alphas, L, D, q, ng, Q, V, h, M]

    def getSolution(self, x, atm_params, sc_params, t, alpha_min, alpha_max, dalpha):

        self.x = x
        self.sc_params = sc_params
        self.atm_params = atm_params
        self.t = t
        self.alpha_min = alpha_min
        self.alpha_max = alpha_max
        self.dalpha = dalpha
            
        alphas = self.optimize(alpha_min, alpha_max, dalpha, t)

        opt = self.solve(alphas, t)

        return opt

    # Evolutionary algorithm
    """
    def alpha_init(self, alpha_min, alpha_max, t):
        
        alphas = np.linspace(alpha_min, alpha_max, len(t))

        alphas = alphas.tolist()

        return alphas
    """
    def optimize(self, alpha_min, alpha_max, dalpha, t):
        
        eta = self.eta

        toolbox = self.toolbox
        
        alpha_seq = np.linspace(alpha_min, alpha_max, len(t))
        alpha_seq = alpha_seq.tolist()

        toolbox.register("alpha", random.uniform, alpha_min, alpha_max)
        #toolbox.register("alpha", self.alpha_init, alpha_min, alpha_max, t)
        toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.alpha, n=len(t))
        #toolbox.register("individual", tools.initIterate, list, self.alpha_init(alpha_min, alpha_max, t))
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        toolbox.register("mate", tools.cxSimulatedBinaryBounded, low=alpha_min, up=alpha_max, eta=eta)
        toolbox.register("mutate", tools.mutPolynomialBounded, low=alpha_min, up=alpha_max, eta=eta, indpb=1/len(t))

        pop = toolbox.population(n=self.npop)
        NGEN = self.ngen
        CXPB = self.cxpb
        #MUTPB = self.mutpb
        MUTPB = 1 / len(t)
        pareto = tools.ParetoFront()

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

            pop = toolbox.select(pop + offspring, self.npop)

        pareto.update(pop)
        
        final_pareto = pareto.items
        
        best = min(final_pareto, key=lambda ind: ind.fitness.values[0])

        return best

    # Check feasibility of solution

    def check_feasibility(self, params):

        alphas = params

        t = self.t

        opt = self.solve(alphas, t)
        
        q = opt[3]
        ng = opt[4]
        Q = opt[5]
        h = opt[7]
        
        if any(np.isnan(x) for x in q) or any(np.isnan(x) for x in ng) or any(np.isnan(x) for x in Q):
            print('Failed: NaNs in solution')
            return False
        if any(x > 1.5E5 for x in h):
            print('wrong way')
            return False
        #if max(dalphas) > self.dalpha:
        #    print("dalpha failed: dalpha = %.5f" % max(dalphas))
        #    return False
        if max(ng) > self.ng_max:
            print('ng failed: ng = %.5f' % max(ng))
            return False
        if max(q) > self.q_max:
            print('q failed: q = %.5f' % max(q))
            return False
        if max(Q) > self.Q_max:
            print('Q failed: Q = %.5f' % max(q))
            return False
        else:
            print('pass')
            return True

    # Evaluate solution with objective function

    def evaluate(self, params):

        alphas = params

        opt = self.solve(alphas, self.t)

        q_max = max(opt[3])
        ng_max = max(opt[4])
        Q_max = max(opt[5])

        return q_max, ng_max, Q_max
