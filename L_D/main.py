import numpy as np
import matplotlib.pyplot as plt
from GeneticAlgorithm import GeneticAlgorithmOptimization
from Atmosphereic_Conditions import Get_SoundSpeed
from State_System import Get_VInit
from scipy.integrate import simpson
from heatshieldpoints import generate_heatshield_points

# Spacecraft parameters (temp values)
dhs = 3.9116  # heatshield diameter (m)
# R0hs = 4.69392  # heatshield radius of curvature
hhs = 0.635  # heatshield height (m)
R0hs = (dhs/2)**2/hhs  # heatshield radius of curvature

m = 5357 / np.pi  # mass [kg]
# x_lst = np.arange(-2, 2.1, 0.1) # x coords
# y_lst = [(x/5)**2 for x in x_lst] # y coords
x_lst, y_lst = generate_heatshield_points(dhs, hhs)

sc_params = [R0hs, m, x_lst, y_lst]

# ICs (temps values)
h0 = 120000  # initial altitude [m]
fpa0 = 6.5 * np.pi/180 + np.pi  # reentry angle [rad]
V0 = 11200  # initial velocity magnitude [m/s]
ICs = [V0, fpa0, h0]  # initial conditions vector

# control var

alpha = 0  # AoA [rad], to be optimized

g = 9.80665  # acceleration due to gravity [m/s^2]
gamma = 1.4  # ratio of specific heats

atm_params = [g, gamma]

if __name__ == "__main__":
    
    GA = GeneticAlgorithmOptimization()
    
    t = np.linspace(0, 500, 101)
    q = np.zeros(len(t))
    Q = np.zeros(len(t))
    L = np.zeros(len(t))
    D = np.zeros(len(t))
    ng = np.zeros(len(t))
    AOA = np.zeros(len(t))
    M = np.zeros(len(t))

    x = [ICs]

    for i in range(0, len(t) - 1):
        
        print(t[i])
        t1 = t[i]
        t2 = t[i+1]
        tspan = [t1, t2]

        opt, sol = GA.getSolution(x, atm_params, sc_params, tspan)  # run genetic algorithm to get optimum
        a = Get_SoundSpeed(sol[-1,2])

        AOA[i] = opt[0]
        alpha = opt[0]
        L[i] = opt[1]
        D[i] = opt[2]
        q[i] = opt[3]
        ng[i] = opt[4]
        M[i] = sol[-1, 0]/a

        Q[i] = simpson(q[0:i + 1], t[0:i + 1])

        x.append([sol[-1, 0], sol[-1, 1], sol[-1, 2]])

        # At Mach 2 Newtonian breaks so stop simulation
        if M[i] < 2:
            t = t[0:i + 1]
            q = q[0:i + 1]
            L = L[0:i + 1]
            D = D[0:i + 1]
            Q = Q[0:i + 1]
            ng = ng[0:i + 1]
            AOA = AOA[0:i + 1]
            M = M[0:i+1]
            break
        # # Parachute height = 8 km, exit loop if h < 8 km
        # if sol[-1, 2] < 8000:
        #     t = t[0:i + 1]
        #     q = q[0:i + 1]
        #     L = L[0:i + 1]
        #     D = D[0:i + 1]
        #     Q = Q[0:i + 1]
        #     ng = ng[0:i + 1]
        #     AOA = AOA[0:i + 1]
        #
        #     break
        
    x = list(zip(*x))

    V = x[0]
    fpa = x[1]
    h = x[2]

    if len(h) == len(t)+1:
        V = V[0:-1]
        fpa = fpa[0:-1]
        h = h[0:-1]

    Vx = np.multiply(V, np.cos(fpa))
    Vy = np.multiply(V, np.sin(fpa))

    plt.figure(1)
    plt.plot(V, h)
    plt.xlabel('V [m/s]')
    plt.ylabel('Altitude [m]')
    plt.tight_layout()

    plt.figure(2)
    plt.plot(t, h)
    plt.xlabel('t [s]')
    plt.ylabel('Altitude [m]')
    plt.tight_layout()

    plt.figure(3)
    plt.plot(t, V)
    plt.xlabel('t [s]')
    plt.ylabel('Velocity [m/s]')
    plt.tight_layout()

    plt.figure(4)
    plt.plot(ng, h)
    plt.xlabel('Deceleration Load [g]')
    plt.ylabel('Altitude [m]')

    plt.figure(5)
    plt.plot(t, ng)
    plt.xlabel('t [s]')
    plt.ylabel('Deceleration Load [g]')
    plt.tight_layout()

    plt.figure(6)
    plt.plot(t, M)
    plt.xlabel('t [s]')
    plt.ylabel('Mach number [-]')
    plt.tight_layout()

    plt.figure(7)
    plt.plot(q, h)
    plt.xlabel(r'Stagnation point heat flux [W/m$^2$]')
    plt.ylabel('Altitude [m]')

    plt.figure(8)
    plt.plot(Q, h)
    plt.xlabel(r'Stagnation point heat load [J/m$^2$]')
    plt.ylabel('Altitude [m]')


    plt.show()
