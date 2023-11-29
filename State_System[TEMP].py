import numpy as np
from scipy.integrate import solve_ivp
from FBD_Stuff import Get_Angles
from Modified_Newton import Get_Normal

#StateUpdate is gives the system of differential equations that needs to be S

def StateUpdate(t : list, states: list, g : float, m : float, alpha : list, gamma : float, x_lst : list, y_lst : list):
    Vx, Vy, h = states

    eta = Get_Angles(Vx, Vy, alpha(t))
    N = Get_Normal(Vx, Vy, h, gamma, x_lst, y_lst, alpha(t))

    dVxdt = N*np.cos(eta-np.pi)/m
    dVydt = -g + N*np.sin(eta-np.pi)/m
    dhdt = Vy

    return [dVxdt, dVydt, dhdt]


#Get_VInit gives intial velocities eta must be measured from POSITIVE X AXIS

def Get_VInit (V0 : float, eta0 : float):
    V0 = abs(V0)

    Vx0 = V0 * np.cos(eta0)
    Vy0 = V0 * np.sin(eta0)

    return Vx0, Vy0


def Solver(V0 : float, eta0 : float, h0 : float, t0 : float, t_end : float, g : float, m : float, alpha : list, gamma : float, x_lst : list, y_lst : list):
    Vx0 , Vy0 = Get_VInit(V0, eta0)
    State0 = [Vx0, Vy0, h0]

    aoa = lambda t: alpha

    param = (g, m, aoa, gamma, x_lst, y_lst) #doubt I can set array as a parameter

    sol = solve_ivp(StateUpdate,(t0,t_end) ,State0, param )

    return sol

print(Solver(3000, 3.5, 100000, 0, 100, 9.81, 100 , np.zeros(1000), 1.4, (-1,0,1), (0,0,0)))

