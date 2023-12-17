import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from Modified_Newton import Get_Lift, Get_Drag, Get_Tangential, Get_Normal
from State_System import StateUpdate, getGForce, getStagHeatFlux, getStagHeatLoad
from heatshieldpoints import generate_heatshield_points
from Atmosphereic_Conditions import Get_Density, Get_SoundSpeed

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
fpa0 = 6.5 * np.pi / 180  # reentry angle [rad]
V0 = 11200  # initial velocity magnitude [m/s]
x0 = 0

ICs = [V0, fpa0, h0, x0]  # initial conditions vector

alpha0 = 1 * np.pi/180  # initial AoA [rad]

g = 9.80665  # acceleration due to gravity [m/s^2]
gamma = 1.4  # ratio of specific heats

Pr = 0.74  # Prandtl number

atm_params = [g, gamma, Pr]

dt = 1
t = np.arange(0, 270, dt)

# solve t = 0 stuff

X0 = [V0, fpa0, h0, x0]  # initial conditions vector
arg = (g, m, gamma, x_lst, y_lst, alpha0, S)
X = odeint(StateUpdate, X0, t, arg)
V = X[:, 0]
fpa = X[:, 1]
h = X[:, 2]
x = X[:, 3]

dVdt = np.zeros(len(V))
dfpadt = np.zeros(len(V))
dhdt = np.zeros(len(V))
dxdt = np.zeros(len(V))
M = np.zeros(len(V))
q = np.zeros(len(V))
Q = np.zeros(len(V))
L = np.zeros(len(V))
D = np.zeros(len(V))
Ct = np.zeros(len(V))
Cn = np.zeros(len(V))


for i in range(len(V)):
    [dVdt[i], dfpadt[i], dhdt[i], dxdt[i]] = StateUpdate([V[i], fpa[i], h[i], x[i]], t, g, m, gamma, x_lst, y_lst, alpha0, S)
    M[i] = V[i] / Get_SoundSpeed(h[i])
    q[i] = getStagHeatFlux(h[i], M[i], gamma, Pr, R0hs)
    L[i] = Get_Lift(V[i], h[i], gamma, x_lst, y_lst, alpha0, S)
    D[i] = Get_Drag(V[i], h[i], gamma, x_lst, y_lst, alpha0, S)
    Ct[i] = Get_Tangential(V[i], h[i], gamma, x_lst, y_lst, alpha0)
    Cn[i] = Get_Normal(V[i], h[i], gamma, x_lst, y_lst, alpha0)

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

Ct = Ct[0:len(M)]
Cn = Cn[0:len(M)]


q = q[0:len(M)]

t = t[0:len(M)]


# plt.figure(1)
# plt.plot(t, h)
# plt.xlabel('t [s]')
# plt.ylabel('Altitude [m]')

# plt.figure(2)
# plt.plot(t, V)
# plt.xlabel('t [s]')
# plt.ylabel('Velocity [m/s]')
#
plt.figure(3)
plt.plot(V, h)
plt.xlabel('V [m/s]')
plt.ylabel('Altitude [m]')

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
#
# plt.figure(8)
# plt.plot(t, ng)
# plt.xlabel('t [s]')
# plt.ylabel('Deceleration Load [g]')

plt.show()
