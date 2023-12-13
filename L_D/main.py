import numpy as np
import matplotlib.pyplot as plt
from GlobalGeneticAlgorithm import GlobalGeneticAlgorithmOptimization
from Modified_Newton import Get_LD
from State_System import getGForce, getStagHeatFlux, getStagHeatLoad
from heatshieldpoints import generate_heatshield_points
from Atmosphereic_Conditions import Get_Density, Get_SoundSpeed

# Spacecraft parameters (temp values)
dhs = 3.9116  # heatshield diameter (m)
# R0hs = 4.69392  # heatshield radius of curvature
hhs = 0.635  # heatshield height (m)
R0hs = (dhs/2)**2/hhs
S = np.pi*(dhs/2)**2

m = 5357  # mass [kg]
# x_lst = np.arange(-2, 2.1, 0.1) # x coords
# y_lst = [(x/5)**2 for x in x_lst] # y coords
x_lst, y_lst = generate_heatshield_points(dhs, hhs)


sc_params = [R0hs, m, x_lst, y_lst, S]

# ICs (temps values)
V0 = 11000  # initial velocity magnitude [m/s]
fpa0 = -6.5 * np.pi/180 # reentry angle [rad]
h0 = 120000  # initial altitude [m]
x0 = 0

ICs = [V0, fpa0, h0, x0]  # initial conditions vector

alpha0 = 0 * np.pi/180 # initial AoA [rad]

alpha_min = -5 * np.pi/180
alpha_max = 5 * np.pi/180

g = 9.80665  # acceleration due to gravity [m/s^2]
gamma = 1.4  # ratio of specific heats

Pr = 0.74 # Prandtl number

atm_params = [g, gamma, Pr]

if __name__ == "__main__":
    
    GA = GlobalGeneticAlgorithmOptimization()

    dt = 1
    t = np.arange(0, 1500, dt)
    
    dalpha = 0.1 * np.pi / 180 * dt

    # solve t = 0 stuff

    AoA0 = alpha0
    L0, D0 = Get_LD(V0, h0, gamma, x_lst, y_lst, alpha0)
    M0 = V0 / Get_SoundSpeed(h0)
    q0 = getStagHeatFlux(h0, M0, gamma, Pr, R0hs)
    Q0 = getStagHeatLoad([q0], t, 0)
    ng0 = getGForce(0, g)
    
    x = [ICs]

    opt = GA.getSolution(x, atm_params, sc_params, t, alpha_min, alpha_max, dalpha)

    AoA = opt[0]; AoA.insert(0, AoA0); AoA = np.array(AoA)
    L = opt[1]; L[0] = L0
    D = opt[2]; D[0] = D0
    q = opt[3]; q[0] = q0
    ng = opt[4]; ng[0] = ng0
    Q = opt[5]; Q[0] = Q0
    V = opt[6]; V[0] = V0
    h = opt[7]; h[0] = h0
    M = opt[8]; M[0] = M0

    t = t[0:len(ng)]
    AoA = AoA[0:len(t)]

    plt.figure(1)
    plt.plot(t, h)
    plt.xlabel('t [s]')
    plt.ylabel('Altitude [m]')

    plt.figure(2)
    plt.plot(t, V)
    plt.xlabel('t [s]')
    plt.ylabel('Velocity [m/s]')

    plt.figure(3)
    plt.plot(t, np.gradient(V, t))
    plt.xlabel('t [s]')
    plt.ylabel('Grad(Velocity)=a [m/s^2]')

    plt.figure(4)
    plt.plot(V, h)
    plt.xlabel('V [m/s]')
    plt.ylabel('Altitude [m]')
    
    plt.figure(5)
    plt.plot(q, h)
    plt.xlabel(r'Stagnation point heat flux [W/m$^2$]')
    plt.ylabel('Altitude [m]')
    
    plt.figure(6)
    plt.plot(Q, h)
    plt.xlabel(r'Stagnation point heat load [J/m$^2$]')
    plt.ylabel('Altitude [m]')
    
    plt.figure(7)
    plt.plot(ng, h)
    plt.xlabel('Deceleration Load [g]')
    plt.ylabel('Altitude [m]')
    
    plt.figure(8)
    plt.plot(t, ng)
    plt.xlabel('t [s]')
    plt.ylabel('Deceleration Load [g]')
    
    plt.figure(9)
    plt.plot(M, L/D)
    plt.xlabel('Mach Number')
    plt.ylabel(r'$\frac{L}{D}$')
    
    plt.figure(10)
    plt.plot(t, AoA * 180/np.pi)
    plt.ylim([-10, 10])
    plt.xlabel('t [s]')
    plt.ylabel(r'$\alpha$ [deg]')

    plt.figure(11)
    plt.plot(t, np.gradient(AoA * 180/np.pi))
    plt.xlabel('t [s]')
    plt.ylabel(r'$\frac{d\alpha}{dt}$ [deg/s]')

    plt.show()
