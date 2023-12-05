import numpy as np
from FBD_Stuff import Get_Angles
from Modified_Newton import Get_Normal, Get_Tangential


# StateUpdate is gives the system of differential equations that needs to be S

def StateUpdate(states: list, t: list, g: float, m: float, alpha: float, gamma: float, x_lst: list, y_lst: list):
    Vx, Vy, h = states

    r = 6378 * 1000 + h
    eta = Get_Angles(Vx, Vy, alpha)[1]
    N = Get_Normal(Vx, Vy, h, gamma, x_lst, y_lst, alpha)
    T = Get_Tangential(Vx, Vy, h, gamma, x_lst, y_lst, alpha)

    dVxdt = N * np.cos(eta - np.pi) / m + T * np.cos(eta - np.pi * 0.5) / m
    dVydt = -g - Vx ** 2 / r + N * np.sin(eta - np.pi) / m + T * np.sin(eta - np.pi * 0.5) / m
    dhdt = Vy

    return [dVxdt, dVydt, dhdt]


# Get_VInit gives intial velocities eta must be measured from POSITIVE X AXIS
def Get_VInit(V0: float, beta0: float):
    V0 = abs(V0)

    Vx0 = V0 * np.cos(beta0)
    Vy0 = V0 * np.sin(beta0)

    return Vx0, Vy0
