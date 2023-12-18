import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from Modified_Newton import Get_Lift, Get_Drag, Get_Ca, Get_Cn
from State_System import StateUpdate, getGForce, getStagHeatFlux, getStagHeatLoad
from heatshieldpoints import generate_heatshield_points
from Atmosphereic_Conditions import Get_Density, Get_SoundSpeed

# Spacecraft parameters (temp values)
Rhs = 1.9558  # heatshield radius (m)
R0hs = 4.6939  # heatshield radius of curvature (m)
delta = 33 * np.pi / 180
S = np.pi*Rhs**2
Aref = np.pi*R0hs**2

m = 5357  # mass [kg]
# x_lst, y_lst = generate_heatshield_points(dhs, hhs)

sc_params = [R0hs, m, delta, S]

# ICs 
h0 = 120000  # initial altitude [m]
fpa0 = 6.5 * np.pi / 180  # reentry angle [rad]
V0 = 11200  # initial velocity magnitude [m/s]
x0 = 0

alpha0 = 20 * np.pi / 180  # initial AoA [rad]

g = 9.80665  # acceleration due to gravity [m/s^2]
gamma = 1.4  # ratio of specific heats

Pr = 0.74  # Prandtl number

atm_params = [g, gamma, Pr]

dt = 1
t = np.arange(0, 600, dt)

# solve t = 0 stuff

X0 = [V0, fpa0, h0, x0]  # initial conditions vector
arg = (g, m, gamma, R0hs, delta, alpha0, S)
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
Cn = np.zeros(len(V))
Ca = np.zeros(len(V))

for i in range(len(V)):
    [dVdt[i], dfpadt[i], dhdt[i], dxdt[i]] = StateUpdate([V[i], fpa[i], h[i], x[i]], t, g, m, gamma, R0hs, delta,
                                                         alpha0, S)
    M[i] = V[i] / Get_SoundSpeed(h[i])
    q[i] = getStagHeatFlux(h[i], M[i], gamma, Pr, R0hs)
    L[i] = Get_Lift(V[i], h[i], gamma, R0hs, delta, alpha0, S)[0]
    D[i] = Get_Drag(V[i], h[i], gamma, R0hs, delta, alpha0, S)[0]
    Cn[i] = Get_Cn(V[i], h[i], gamma, R0hs, delta, alpha0, S)
    Ca[i] = Get_Ca(V[i], h[i], gamma, R0hs, delta, alpha0, S)

    if M[i] < 3:
        break

M = M[M >= 3]
for i in range(len(M)):
    Q[i] = getStagHeatLoad(q, t, i)

V = V[0:len(M)]
fpa = fpa[0:len(M)]
h = h[0:len(M)]
x = x[0:len(M)]

dVdt = dVdt[0:len(M)]
dfpadt = dfpadt[0:len(M)]
dhdt = dhdt[0:len(M)]
dxdt = dxdt[0:len(M)]

ng = -dVdt / g

D = D[0:len(M)]
L = L[0:len(M)]

Cn = Cn[0:len(M)]
Ca = Ca[0:len(M)]

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
# plt.figure(3)
# plt.plot(V, h)
# plt.xlabel('V [m/s]')
# plt.ylabel('Altitude [m]')

# plt.figure(5)
# plt.plot(q, h)
# plt.xlabel(r'Stagnation point heat flux [W/m$^2$]')
# plt.ylabel('Altitude [m]')

# plt.figure(6)
# plt.plot(Q, h)
# plt.xlabel(r'Stagnation point heat load [J/m$^2$]')
# plt.ylabel('Altitude [m]')

# plt.figure(7)
# plt.plot(ng, h)
# plt.xlabel('Deceleration Load [g]')
# plt.ylabel('Altitude [m]')
#
plt.figure(8)
plt.plot(t, ng)
plt.xlabel('t [s]')
plt.ylabel('Deceleration Load [g]')

plt.show()
