import numpy as np
from scipy.integrate import odeint
from FBD_Stuff import Get_Angles
from Modified_Newton import Get_Normal, Get_Tangential

#StateUpdate is gives the system of differential equations that needs to be S

def StateUpdate( states: list, t : list,g : float, m : float, alpha : float, gamma : float, x_lst : list, y_lst : list):
    Vx, Vy, h = states
    R = 6378*1000
    eta = Get_Angles(Vx, Vy, alpha)[1]
    N = Get_Normal(Vx, Vy, h, gamma, x_lst, y_lst, alpha)
    T = Get_Tangential(Vx, Vy, h, gamma, x_lst, y_lst, alpha)

    dVxdt = N*np.cos(eta-np.pi)/m + T*np.cos(eta-np.pi*0.5)/m
    dVydt = -g - Vx**2/(R+h) + N*np.sin(eta-np.pi)/m + T*np.sin(eta-np.pi*0.5)/m
    dhdt = Vy

    return [dVxdt, dVydt, dhdt]


#Get_VInit gives intial velocities eta must be measured from POSITIVE X AXIS

def Get_VInit (V0 : float, beta0 : float):
    V0 = abs(V0)

    Vx0 = V0 * np.cos(beta0)
    Vy0 = V0 * np.sin(beta0)

    return Vx0, Vy0

"""
def Solver(V0 : float, beta0 : float, h0 : float, t: list, g : float, m : float, alpha : list, gamma : float, x_lst : list, y_lst : list):
    Vx0 , Vy0 = Get_VInit(V0, beta0)
    State0 = [Vx0, Vy0, h0]

    param = (g, m, alpha, gamma, x_lst, y_lst) #doubt I can set array as a parameter

    sol = odeint(StateUpdate,State0 , t, args=(param) )

    return sol


#TEST
print(Solver(10000, 3.25, 100000, [0,1,2,3,4,5,6,7,8,9], 9.81, 1000, [0,0,0,0,0,0,0,0,0,0], 1.4, [-1,0,1], [0,0,0]))
"""
