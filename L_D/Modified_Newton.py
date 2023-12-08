from tkinter import N
import numpy as np
from scipy.integrate import simpson, trapezoid
import csv

from scipy.integrate._ivp.radau import C
from Atmosphereic_Conditions import Get_SoundSpeed, Get_Density


# x and y lists must be points specified from left to right on sideview
# nr of points must be odd

# Get_Theta gives the local incidence angle at each panel on the surface

def Get_Theta(x_lst: list, y_lst: list, alpha: float):
    theta = np.array([])

    for i in range(len(x_lst) - 1):

        slope = (y_lst[i + 1] - y_lst[i]) / (x_lst[i + 1] - x_lst[i])

        if i + 1 > len(x_lst) / 2:
            theta_temp = np.arctan(abs(slope)) - alpha
        else:
            theta_temp = np.arctan(abs(slope)) + alpha
        
        theta = np.append(theta, theta_temp)
        
    return theta


def Get_MidPoints(x_lst: list, y_lst: list):  # Get_MidPoints computes the x and y coodinate of the midpoints of each panel

    x_mid = np.array([])
    y_mid = np.array([])

    for i in range(len(x_lst) - 1):
        x_mid = np.append(x_mid, (x_lst[i] + x_lst[i + 1]) / 2)
        y_mid = np.append(y_mid, (y_lst[i] + y_lst[i + 1]) / 2)

    return x_mid, y_mid


# Get_Sin2Int gives the integral of sin^2 of these angles

"for now 2D, might just be multiplying with pi to get 3D, not sure"


def Get_Sin2Int(x_lst: list, y_lst: list, alpha: float):
    theta = Get_Theta(x_lst, y_lst, alpha)

    x = Get_MidPoints(x_lst, y_lst)[0]

    # have to use midpoints of panels, lose one point
    integral = simpson(np.sin(theta) ** 2, x)

    return integral


# Get_CpMax gives the cpmax from modified newtonian theory

def Get_CpMax(V: float, h: float, gamma: float):
    V_inf = V
    a = Get_SoundSpeed(h)
    M_inf = V_inf / a

    part1 = (1 - gamma + 2 * gamma * M_inf ** 2) / (gamma + 1)

    part2 = (((gamma + 1) ** 2 * M_inf ** 2) / (4 * M_inf ** 2 * gamma - 2 * (gamma - 1))) ** (gamma / (gamma - 1))

    Cp_max = 2 / (gamma * M_inf ** 2) * (part1 * part2 - 1)

    return Cp_max


# Get_Normal gives the normal force working on the SC

def Get_Normal(V: float, h: float, gamma: float, x_lst: list, y_lst: list, alpha: float):

    Cp_max = Get_CpMax(V, h, gamma)
    sin2th = Get_Sin2Int(x_lst, y_lst, alpha)

    Cn = Cp_max * sin2th

    N = 0.5 * Cn * Get_Density(h) * V**2 * np.pi * 3.9116/2

    return Cn


def Get_Tangential(V: float, h: float, gamma: float, x_lst: list, y_lst: list, alpha: float):

    CP_max = Get_CpMax(V, h, gamma)
    theta = Get_Theta(x_lst, y_lst, alpha)

    y = Get_MidPoints(x_lst, y_lst)[1]

    Cp_local = CP_max * np.sin(theta) ** 2

    integral_left = simpson(np.flip(Cp_local[:len(Cp_local) // 2]), np.flip(y[:len(Cp_local) // 2]))

    integral_right = simpson(Cp_local[len(Cp_local) // 2:], y[len(Cp_local) // 2:])

    T = 0.5 * (integral_right - integral_left) * Get_Density(h) * V**2 * np.pi * 3.9116/2

    return integral_right - integral_left


def Get_Lift(N: float, T: float, V : float, h : float, alpha: float, x_lst : list, y_lst : list):
    L = -N*np.sin(alpha) + T*np.cos(alpha)
    #L = -T * np.cos(np.pi - alpha) - N * np.sin(np.pi - alpha)
    #Cpmax = Get_CpMax(V, h, 1.4)
    #sin2int = Get_Sin2Int(x_lst, y_lst, alpha)
    #CL = Cpmax * sin2int * np.cos(alpha)
    #L = 0.5 * CL * Get_Density(h) * V**2
    #Cl = 2 * np.sin(np.pi - alpha)**2 * np.cos(np.pi - alpha)
    #L = 0.5 * Get_Density(h) * V**2 * Cl * 12
    return L


def Get_Drag(N: float, T: float, V : float, h : float, alpha: float, x_lst : list, y_lst : list):
    D = N * np.cos(alpha) + T * np.sin(alpha)
    #D = T * np.sin(np.pi - alpha) - N * np.cos(np.pi - alpha)
    #Cpmax = Get_CpMax(V, h, 1.4)
    #sin2int = Get_Sin2Int(x_lst, y_lst, alpha)
    #CD = Cpmax * sin2int * np.sin(alpha)
    #D = 0.5 * CD * Get_Density(h) * V**2
    #Cd = 2 * np.sin(np.pi - alpha)**3
    #D = 0.5 * Get_Density(h) * V**2 * Cd * 12
    return D


def Get_LD(V : float, h : float, alpha : float):

    Cpmax = Get_CpMax(V, h, 1.4)
    q = 0.5 * Get_Density(h) * V**2
    
    CD = Cpmax * np.sin(alpha)**2 * np.cos(alpha)
    CL = Cpmax * np.sin(alpha)**3
    
    L = q * CL * 12
    D = q * CD * 12
    
    return L, D