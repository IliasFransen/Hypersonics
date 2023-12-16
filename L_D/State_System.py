import numpy as np
<<<<<<< Updated upstream
from Modified_Newton import Get_Drag, Get_Lift, Get_Tangential, Get_Normal
=======
from Modified_Newton import Get_Drag, Get_Lift
from Atmosphereic_Conditions import Get_Density, Get_SoundSpeed, Get_Temperature
from scipy.integrate import simpson
>>>>>>> Stashed changes


# StateUpdate is gives the system of differential equations that needs to be S

<<<<<<< Updated upstream
def StateUpdate(states: list, t: list, g: float, m: float, alpha: float, gamma: float, x_lst: list, y_lst: list):
    V, fpa, h = states

    r = 6378 * 1000 + h

    N = Get_Normal(V, h, gamma, x_lst, y_lst, alpha)
    T = Get_Tangential(V, h, gamma, x_lst, y_lst, alpha)

    L = Get_Lift(N, T, alpha)
    D = Get_Drag(N, T, alpha)

    dVdt = -D / m - g * np.sin(fpa)
    dfpadt = L / (V * m) - (g/V - V / r) * np.cos(fpa)
    dhdt = V*np.sin(fpa)

    return [dVdt, dfpadt, dhdt]
=======
def StateUpdate(states: list, t: list, g: float, m: float, gamma: float, x_lst: list, y_lst: list, S : float, alpha : float):
    V, fpa, h, x = states

    Re = 6378 * 1000
    r = Re + h
    g = g * (Re/r)**2

    D = Get_Drag(V, h, gamma, x_lst, y_lst, S, alpha)
    L = Get_Lift(V, h, gamma, x_lst, y_lst, S, alpha)

    dVdt = -D / m + g * np.sin(fpa)
    dfpadt = -L / (V * m) + (g / V - V / r) * np.cos(fpa)
    dhdt = -V*np.sin(fpa)
    dxdt = V*np.cos(fpa)
>>>>>>> Stashed changes


<<<<<<< Updated upstream
# Get_VInit gives intial velocities, fpa must be measured from POSITIVE X AXIS

def Get_VInit(V0: float, fpa0: float):
    V0 = abs(V0)

    Vx0 = V0 * np.cos(fpa0)
    Vy0 = V0 * np.sin(fpa0)

    return Vx0, Vy0

=======
def getViscosity(T : float):
>>>>>>> Stashed changes

"""
def Solver(V0 : float, beta0 : float, h0 : float, t: list, g : float, m : float, alpha : list, gamma : float, x_lst : list, y_lst : list):
    Vx0 , Vy0 = Get_VInit(V0, beta0)
    State0 = [Vx0, Vy0, h0]

<<<<<<< Updated upstream
    param = (g, m, alpha, gamma, x_lst, y_lst) #doubt I can set array as a parameter
=======
    T_ref = 288
    mu_ref = 1.789 * 10**(-5)
    S = 110
    
    mu = mu_ref * (T / T_ref)**1.5 * (T_ref + S) / (T + S)
    
    return mu
>>>>>>> Stashed changes

    sol = odeint(StateUpdate,State0 , t, args=(param) )

    return sol


#TEST
print(Solver(10000, 3.25, 100000, [0,1,2,3,4,5,6,7,8,9], 9.81, 1000, [0,0,0,0,0,0,0,0,0,0], 1.4, [-1,0,1], [0,0,0]))
"""
