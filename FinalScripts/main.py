import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from Modified_Newton import Get_Lift, Get_Drag, Get_Ca, Get_Cn
from State_System import StateUpdate, getGForce, getStagHeatFlux, getStagHeatLoad
from heatshieldpoints import generate_heatshield_points
from Atmosphereic_Conditions import Get_Density, Get_SoundSpeed
from GlobalGeneticAlgorithm import GlobalGeneticAlgorithmOptimization

# Spacecraft parameters (temp values)
Rhs = 1.9558  # heatshield radius (m)
R0hs = 4.6939  # heatshield radius of curvature (m)
delta = 33 * np.pi / 180
S = np.pi*Rhs**2
Aref = np.pi * R0hs**2
m = 5470  # mass [kg]
# x_lst, y_lst = generate_heatshield_points(dhs, hhs)

sc_params = [R0hs, m, delta, S, Aref]

# ICs
h0 = 121920  # initial altitude [m]
fpa0 = 6.49 * np.pi / 180  # reentry angle [rad]
V0 = 11032  # initial velocity magnitude [m/s]
x0 = 0

ICs = [V0, fpa0, h0, x0]

g = 9.80665  # acceleration due to gravity [m/s^2]
gamma = 1.4  # ratio of specific heats

Pr = 0.74  # Prandtl number

atm_params = [g, gamma, Pr]

dt = 1
t = np.arange(0, 900, dt)

GA = GlobalGeneticAlgorithmOptimization()

alpha_min = 18 * np.pi/180
alpha_max = 24 * np.pi/180
alpha_0 = 18 * np.pi/180

dalpha = 0.1*np.pi/180

AoA0 = alpha_0
L0, CL0 = Get_Lift(V0, h0, gamma, R0hs, delta, alpha_0, S)
D0, CD0 = Get_Drag(V0, h0, gamma, R0hs, delta, alpha_0, S)
CN0 = Get_Cn(V0, h0, gamma, R0hs, delta, alpha_0, S)
CA0 = Get_Ca(V0, h0, gamma, R0hs, delta, alpha_0, S)
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
s = opt[9]; s[0] = x0
CL = opt[10]; CL[0] = CL0
CD = opt[11]; CD[0] = CD0
CA = opt[12]; CA[0] = CA0
CN = opt[13]; CN[0] = CN0

t = t[0:len(ng)]
AoA = AoA[0:len(t)]

# Plotting

image_format = 'svg'

plt.figure(1)
plt.plot(t, h/1000)

plt.grid()

plt.xlabel('t [s]')
plt.ylabel('Altitude [km]')

plt.tight_layout()
image_name = 'h_vs_t.svg'
plt.savefig(image_name, format=image_format,dpi=1200)
#plt.savefig('h_vs_t.png')

plt.figure(2)
plt.plot(t, V/1000)

plt.grid()

plt.xlabel('t [s]')
plt.ylabel('Velocity [km/s]')

plt.tight_layout()
image_name = 'V_vs_t.svg'
plt.savefig(image_name, format=image_format, dpi=1200)
#plt.savefig('V_vs_t.png')

plt.figure(3)
plt.plot(V/1000, h/1000)

plt.grid()

plt.xlabel('V [km/s]')
plt.ylabel('Altitude [km]')

plt.tight_layout()
image_name = 'h_vs_V.svg'
plt.savefig(image_name, format=image_format, dpi=1200)
#plt.savefig('V_vs_h.png')

plt.figure(5)
plt.plot(q/1E6, h/1000)

plt.grid()

plt.xlabel(r'Stagnation point heat flux [MW/m$^2$]')
plt.ylabel('Altitude [km]')

plt.tight_layout()
image_name = 'h_vs_qw.svg'
plt.savefig(image_name, format=image_format, dpi=1200)
#plt.savefig('qw_vs_h.png')

plt.figure(6)
plt.plot(Q/1E6, h/1000)

plt.grid()

plt.xlabel(r'Stagnation point heat load [MJ/m$^2$]')
plt.ylabel('Altitude [km]')

plt.tight_layout()
image_name = 'h_vs_Q.svg'
plt.savefig(image_name, format=image_format, dpi=1200)
#plt.savefig('Q_vs_h.png')

plt.figure(7)
plt.plot(ng, h/1000)

plt.grid()

plt.xlabel('Deceleration Load [g]')
plt.ylabel('Altitude [km]')

plt.tight_layout()
image_name = 'h_vs_ng.svg'
plt.savefig(image_name, format=image_format, dpi=1200)
#plt.savefig('ng_vs_h.png')

plt.figure(8)
plt.plot(t, ng)

plt.grid()

plt.xlabel('t [s]')
plt.ylabel('Deceleration Load [g]')

plt.tight_layout()
image_name = 'ng_vs_t.svg'
plt.savefig(image_name, format=image_format, dpi=1200)
#plt.savefig('ng_vs_t.png')

plt.figure(9)
plt.plot(M, h/1000)

plt.grid()

plt.xlabel('Mach Number [-]')
plt.ylabel('Altitude [km]')

plt.tight_layout()
image_name = 'h_vs_M.svg'
plt.savefig(image_name, format=image_format, dpi=1200)
#plt.savefig('h_vs_M.png')

plt.figure(10)
plt.plot(t, q/1E6)

plt.grid()

plt.xlabel('t [s]')
plt.ylabel(r'Stagnation point heat flux [MW/m$^2$]')

plt.tight_layout()
image_name = 'qw_vs_t.svg'
plt.savefig(image_name, format=image_format, dpi=1200)
#plt.savefig('qw_vs_t.png')

plt.figure(11)
plt.plot(AoA*180/np.pi, h/1000)

plt.grid()

plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel('Altitude [km]')

plt.tight_layout()
image_name = 'h_vs_AoA.svg'
plt.savefig(image_name, format=image_format, dpi=1200)
#plt.savefig('h_vs_AoA.png')

plt.figure(12)
plt.plot(s/1000, h/1000)

plt.grid()

plt.xlabel('Range [km]')
plt.ylabel('Altitude [km]')

plt.tight_layout()
image_name = 'h_vs_x.svg'
plt.savefig(image_name, format=image_format, dpi=1200)
#plt.savefig('h_vs_x.png')

plt.figure(13)
plt.plot(CA, h/1000)

plt.grid()

plt.xlabel(r'$C_A$ [-]')
plt.ylabel('Altitude [km]')

plt.tight_layout()
image_name = 'h_vs_CA.svg'
plt.savefig(image_name, format=image_format, dpi=1200)
#plt.savefig('h_vs_CA.png')

plt.figure(14)
plt.plot(CN, h/1000)

plt.grid()

plt.xlabel(r'$C_N$ [-]')
plt.ylabel('Altitude [km]')

plt.tight_layout()
image_name = 'h_vs_CN.svg'
plt.savefig(image_name, format=image_format, dpi=1200)
#plt.savefig('h_vs_CN.png')

plt.figure(15)
plt.plot(CL, h/1000)

plt.grid()

plt.xlabel(r'$C_L$ [-]')
plt.ylabel('Altitude [km]')

plt.tight_layout()
image_name = 'h_vs_CL.svg'
plt.savefig(image_name, format=image_format, dpi=1200)
#plt.savefig('h_vs_CL.png')

plt.figure(16)
plt.plot(CD, h/1000)

plt.grid()

plt.xlabel(r'$C_D$ [-]')
plt.ylabel('Altitude [km]')

plt.tight_layout()
image_name = 'h_vs_CD.svg'
plt.savefig(image_name, format=image_format, dpi=1200)
#plt.savefig('h_vs_CD.png')

plt.figure(17)
plt.plot(t, AoA*180/np.pi)

plt.grid()

plt.xlabel('t [s]')
plt.ylabel(r'$\alpha$ [deg]')

plt.tight_layout()
image_name = 'AoA_vs_t.svg'
plt.savefig(image_name, format=image_format, dpi=1200)

plt.figure(18)
plt.plot(t, CA)

plt.grid()

plt.xlabel('t [s]')
plt.ylabel(r'$C_A$')

plt.tight_layout()
image_name = 'CA_vs_t.svg'
plt.savefig(image_name, format=image_format, dpi=1200)

plt.figure(19)
plt.plot(t, CN)

plt.grid()

plt.xlabel('t [s]')
plt.ylabel(r'$C_N$')

plt.tight_layout()
image_name = 'CN_vs_t.svg'
plt.savefig(image_name, format=image_format, dpi=1200)

plt.figure(20)
plt.plot(t, CD)

plt.grid()

plt.xlabel('t [s]')
plt.ylabel(r'$C_D$')

plt.tight_layout()
image_name = 'CD_vs_t.svg'
plt.savefig(image_name, format=image_format, dpi=1200)

plt.figure(21)
plt.plot(t, CL)

plt.grid()

plt.xlabel('t [s]')
plt.ylabel(r'$C_L$')

plt.tight_layout()
image_name = 'CL_vs_t.svg'
plt.savefig(image_name, format=image_format, dpi=1200)

#plt.show()
