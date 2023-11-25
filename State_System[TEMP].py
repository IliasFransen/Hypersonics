import numpy as np
from FBD_Stuff import Get_Angles
from Modified_Newton import Get_Normal

#StateUpdate is gives the system of differential equations that needs to be S

def StateUpdate(N : list, g : float, m : float, alpha : float, states: list, gamma : float, x_lst : list, y_lst : list):
    Vx, Vy, h = states

    eta = Get_Angles(Vx, Vy, alpha)
    N = Get_Normal(Vx, Vy, h, gamma, x_lst, y_lst, alpha)

    dVxdt = N*np.cos(eta-np.pi)/m
    dVydt = -g + N*np.sin(eta-np.pi)/m
    dhdt = Vy

    return [dVxdt, dVydt, dhdt]