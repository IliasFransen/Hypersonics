import numpy as np
import matplotlib.pyplot as plt
<<<<<<< Updated upstream
from GeneticAlgorithm import GeneticAlgorithmOptimization
from State_System import Get_VInit
from scipy.integrate import simpson
from heatshieldpoints import generate_heatshield_points

# Spacecraft parameters (temp values)
dhs = 3.9116  # heatshield diameter (m)
R0hs = 4.69392  # heatshield radius of curvature
hhs = 0.635  # heatshield height (m)

m = 5357 / np.pi  # mass [kg]
# x_lst = np.arange(-2, 2.1, 0.1) # x coords
# y_lst = [(x/5)**2 for x in x_lst] # y coords
x_lst, y_lst = generate_heatshield_points(R0hs, dhs, hhs)

sc_params = [R0hs, m, x_lst, y_lst]

# ICs (temps values)
h0 = 120000  # initial altitude [m]
fpa0 = 6.5 * np.pi/180 + np.pi  # reentry angle [rad]
V0 = 11200  # initial velocity magnitude [m/s]
=======
from scipy.integrate import odeint
from Modified_Newton import Get_Lift, Get_Drag
from State_System import StateUpdate, getGForce, getStagHeatFlux, getStagHeatLoad
from heatshieldpoints import generate_heatshield_points
from Atmosphereic_Conditions import Get_Density, Get_SoundSpeed
from GlobalGeneticAlgorithm import GlobalGeneticAlgorithmOptimization

# Spacecraft parameters (temp values)
dhs = 3.9116  # heatshield diameter (m)
hhs = 0.635  # heatshield height (m)
R0hs = (dhs / 2) ** 2 / hhs
S = np.pi * (dhs / 2) ** 2

m = 5357  # mass [kg]
x_lst, y_lst = generate_heatshield_points(dhs, hhs)

sc_params = [R0hs, m, x_lst, y_lst, S]

# ICs 
h0 = 120000  # initial altitude [m]
fpa0 = 6.5 * np.pi / 180 # reentry angle [rad]
V0 = 11000  # initial velocity magnitude [m/s]
x0 = 0
>>>>>>> Stashed changes

ICs = [V0, fpa0, h0]  # initial conditions vector

<<<<<<< Updated upstream
# control var

alpha = 0  # AoA [rad], to be optimized

g = 9.80665  # acceleration due to gravity [m/s^2]
gamma = 1.4  # ratio of specific heats

atm_params = [g, gamma]
=======
g = 9.80665  # acceleration due to gravity [m/s^2]
gamma = 1.4  # ratio of specific heats

Pr = 0.74  # Prandtl number

atm_params = [g, gamma, Pr]
>>>>>>> Stashed changes

dt = 1
t = np.arange(0, 1500, dt)

alpha_min = -30 * np.pi/180
alpha_max = -25 * np.pi/180 
dalpha = 0
alpha0 = -25 * np.pi/180

AoA0 = alpha0
L0 = Get_Lift(V0, h0, gamma, x_lst, y_lst, S, alpha0)
D0 = Get_Drag(V0, h0, gamma, x_lst, y_lst, S, alpha0)
M0 = V0 / Get_SoundSpeed(h0)
q0 = getStagHeatFlux(h0, M0, gamma, Pr, R0hs)
Q0 = getStagHeatLoad([q0], t, 0)
ng0 = getGForce(0, g)
    
<<<<<<< Updated upstream
    GA = GeneticAlgorithmOptimization()
    
    t = np.linspace(0, 500, 101)
    q = np.zeros(len(t))
    Q = np.zeros(len(t))
    L = np.zeros(len(t))
    D = np.zeros(len(t))
    ng = np.zeros(len(t))
    AOA = np.zeros(len(t))
    
    x = [ICs]

    for i in range(0, len(t) - 1):
        
        print(t[i])
        t1 = t[i]
        t2 = t[i+1]
        tspan = [t1, t2]

        opt, sol = GA.getSolution(x, atm_params, sc_params, tspan)  # run genetic algorithm to get optimum
        
        AOA[i] = opt[0]
        alpha = opt[0]
        L[i] = opt[1]
        D[i] = opt[2]
        q[i] = opt[3]
        ng[i] = opt[4]
        
        Q[i] = simpson(q[0:i + 1], t[0:i + 1])

        x.append([sol[-1, 0], sol[-1, 1], sol[-1, 2]])

        # Parachute height = 8 km, exit loop if h < 8 km
        if sol[-1, 2] < 8000:
            t = t[0:i + 1]
            q = q[0:i + 1]
            L = L[0:i + 1]
            D = D[0:i + 1]
            Q = Q[0:i + 1]
            ng = ng[0:i + 1]
            AOA = AOA[0:i + 1]

            break
        
    x = list(zip(*x))

    V = x[0]
    fpa = x[1]
    h = x[2]

    Vx = np.multiply(V, np.cos(fpa))
    Vy = np.multiply(V, np.sin(fpa))
    # if len(h) == len(t)+1:
    #     V = V[0:-1]
    #     fpa = fpa[0:-1]
    #     h = h[0:-1]

    plt.figure(1)
    plt.plot(t, h)
    plt.xlabel('t [s]')
    plt.ylabel('Altitude [m]')
    plt.tight_layout()

    plt.figure(2)
    plt.plot(t, V)
    plt.xlabel('t [s]')
    plt.ylabel('Velocity [m/s]')
    plt.tight_layout()

    plt.figure(3)
    plt.plot(t, Vy)
    plt.xlabel('t [s]')
    plt.ylabel('Velocity y [m/s]')
    plt.tight_layout()

    plt.figure(4)
    plt.plot(t, Vx)
    plt.xlabel('t [s]')
    plt.ylabel('Velocity x [m/s]')
    plt.tight_layout()

    plt.figure(5)
    plt.plot(t, np.gradient(V, t))
    plt.xlabel('t [s]')
    plt.ylabel('Grad(Velocity)=a [m/s^2]')
    plt.tight_layout()

    plt.figure(6)
    plt.plot(V, h)
    plt.xlabel('V [m/s]')
    plt.ylabel('Altitude [m]')
    plt.tight_layout()

    # plt.figure(5)
    # plt.plot(q, h)
    # plt.xlabel(r'Stagnation point heat flux [W/m$^2$]')
    # plt.ylabel('Altitude [m]')
    #
    # plt.figure(6)
    # plt.plot(Q, h)
    # plt.xlabel(r'Stagnation point heat load [J/m$^2$]')
    # plt.ylabel('Altitude [m]')
    #
    # plt.figure(7)
    # plt.plot(ng, h)
    # plt.xlabel('Deceleration Load [g]')
    # plt.ylabel('Altitude [m]')
=======
x = [ICs]
    
GA = GlobalGeneticAlgorithmOptimization()

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

# solve t = 0 stuff
"""
X0 = [V0, fpa0, h0, x0]  # initial conditions vector
arg = (g, m, gamma, x_lst, y_lst, S)
X = odeint(StateUpdate, X0, t, arg)
V = X[:, 0]
fpa = X[:, 1]
h = X[:, 2]
x = X[:, 3]
>>>>>>> Stashed changes

dVdt = np.zeros(len(V))
dfpadt = np.zeros(len(V))
dhdt = np.zeros(len(V))
dxdt = np.zeros(len(V))
M = np.zeros(len(V))
q = np.zeros(len(V))
Q = np.zeros(len(V))
L = np.zeros(len(V))
D = np.zeros(len(V))

for i in range(len(V)):
    [dVdt[i], dfpadt[i], dhdt[i], dxdt[i]] = StateUpdate([V[i], fpa[i], h[i], x[i]], t, g, m, gamma, x_lst, y_lst, S)
    M[i] = V[i] / Get_SoundSpeed(h[i])
    q[i] = getStagHeatFlux(h[i], M[i], gamma, Pr, R0hs)
    L[i] = Get_Lift(V[i], h[i], gamma, x_lst, y_lst, S)
    D[i] = Get_Drag(V[i], h[i], gamma, x_lst, y_lst, S)

    if M[i] < 3:
        break

M = M[M>=3]
for i in range(len(M)):
    Q[i] = getStagHeatLoad(q, t, i)

ng = -dVdt/g

V = V[0:len(M)]
fpa = fpa[0:len(M)]
h = h[0:len(M)]
x = x[0:len(M)]

dVdt = dVdt[0:len(M)]
dfpadt = dfpadt[0:len(M)]
dhdt = dhdt[0:len(M)]
dxdt = dxdt[0:len(M)]

D = D[0:len(M)]
L = L[0:len(M)]

q = q[0:len(M)]
Q = q[0:len(M)]
ng = ng[0:len(M)]
t = t[0:len(M)]
"""

plt.figure(1)
plt.plot(t, h)
plt.xlabel('t [s]')
plt.ylabel('Altitude [m]')

plt.figure(2)
plt.plot(t, V)
plt.xlabel('t [s]')
plt.ylabel('Velocity [m/s]')

plt.figure(3)
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
plt.plot(t, AoA * 180/np.pi)
plt.xlabel('t [s]')
plt.ylabel(r'$\alpha$ [deg]')

plt.figure(10)
plt.plot(t, np.gradient(AoA) * 180/np.pi)
plt.xlabel('t [s]')
plt.ylabel(r'$\frac{d\alpha}{dt}$ [deg/s]')              

plt.show()
