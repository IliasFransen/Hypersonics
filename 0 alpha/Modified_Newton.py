import numpy as np
from scipy.integrate import simpson, trapezoid
from Atmosphereic_Conditions import Get_SoundSpeed, Get_Density


# x and y lists must be points specified from left to right on sideview
# nr of points must be odd

# Get_Theta gives the local incidence angle at each panel on the surface

def Get_Theta(x_lst: list, y_lst: list):
    theta = np.array([])

    for i in range(len(x_lst) - 1):

        slope = (y_lst[i + 1] - y_lst[i]) / (x_lst[i + 1] - x_lst[i])
        theta_temp = np.arctan(abs(slope))
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


def Get_Sin2Int(x_lst: list, y_lst: list):
    theta = Get_Theta(x_lst, y_lst)

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

def Get_Drag(V: float, h: float, gamma: float, x_lst: list, y_lst: list, S):

    Cp_max = Get_CpMax(V, h, gamma)
    sin2th = Get_Sin2Int(x_lst, y_lst)

    b = x_lst[-1] - x_lst[0]
    Cn = 1/b * Cp_max * sin2th
    D = 1/2 * Cn * Get_Density(h) * V**2 * S
    return D


def Get_Lift(V: float, h: float, gamma: float, x_lst: list, y_lst: list, S):

    CP_max = Get_CpMax(V, h, gamma)
    theta = Get_Theta(x_lst, y_lst)

    y = Get_MidPoints(x_lst, y_lst)[1]

    Cp_local = CP_max * np.sin(theta) ** 2

    integral_left = simpson(np.flip(Cp_local[:len(Cp_local) // 2]), np.flip(y[:len(Cp_local) // 2]))

    integral_right = simpson(Cp_local[len(Cp_local) // 2:], y[len(Cp_local) // 2:])

    c = y_lst[-1] - y_lst[0]
    Ct = (integral_right - integral_left) / c

    L = 1/2 * Ct * Get_Density(h) * V**2 * S
    return Ct

def Get_length(x_lst: list, y_lst: list):
    x,y = Get_MidPoints(x_lst, y_lst)

    lengths = np.array([])
    for i in range(len(x)-1):
        leng = ((x[i]-x[i+1])**2 + (y[i]-y[i+1])**2)**0.5
        lengths = np.append(lengths,leng)

    middle_ind = int(np.floor(len(lengths)/2))
    lengths = np.delete(lengths,middle_ind)

    return np.sum(lengths)
