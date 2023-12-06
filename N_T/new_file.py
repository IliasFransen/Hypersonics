import numpy as np
from scipy.integrate import simpson

from Atmosphereic_Conditions import Get_SoundSpeed, Get_Density


# x and y lists must be points specified from left to right on sideview
# nr of points must be odd

# Get_Theta gives the local incidence angle at each panel on the surface

def Get_Theta(x_lst: list, y_lst: list, alpha: float):
    theta = np.array([])

    for i in range(len(x_lst) - 1):

        slope = (y_lst[i + 1] - y_lst[i]) / (x_lst[i + 1] - x_lst[i])
        
        if i+1 > len(x_lst) / 2:
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

def Get_Sin2Int(x_lst: list, y_lst: list, alpha: float):
    theta = Get_Theta(x_lst, y_lst, alpha)

    x = Get_MidPoints(x_lst, y_lst)[0]

    # have to use midpoints of panels, lose one point
    integral = simpson(np.sin(theta) ** 2, x)

    return integral


# Get_CpMax gives the cpmax from modified newtonian theory

def Get_CpMax(Vx: float, Vy: float, h: float, gamma: float):
    V_inf = (Vx ** 2 + Vy ** 2) ** 0.5
    a = Get_SoundSpeed(h)
    M_inf = V_inf / a

    part1 = (1 - gamma + 2 * gamma * M_inf ** 2) / (gamma + 1)

    part2 = (((gamma + 1) ** 2 * M_inf ** 2) / (4 * M_inf ** 2 * gamma - 2 * (gamma - 1))) ** (gamma / (gamma - 1))

    Cp_max = 2 / (gamma * M_inf ** 2) * (part1 * part2 - 1)

    return Cp_max


# Get_Normal gives the normal force working on the SC

def Get_Normal(Vx: float, Vy: float, h: float, gamma: float, x_lst: list, y_lst: list, alpha: float):
    Cp_max = Get_CpMax(Vx, Vy, h, gamma)
    sin2th = Get_Sin2Int(x_lst, y_lst, alpha)

    Cn = Cp_max * sin2th

    N = Cn * Get_Density(h) * (Vx ** 2 + Vy ** 2)*0.5

    return N


def Get_Tangential(Vx: float, Vy: float, h: float, gamma: float, x_lst: list, y_lst: list, alpha: float):

    CP_max = Get_CpMax(Vx, Vy, h, gamma)
    theta = Get_Theta(x_lst, y_lst, alpha)

    y = Get_MidPoints(x_lst, y_lst)[1]

    Cp_local = CP_max * np.sin(theta) ** 2

    integral_left = simpson(Cp_local[:len(Cp_local) // 2], y[:len(Cp_local) // 2])

    integral_right = simpson(Cp_local[len(Cp_local) // 2:], y[len(Cp_local) // 2 :])

    T = (integral_right - integral_left) * Get_Density(h) * (Vx ** 2 + Vy ** 2)*0.5
    return T

def Get_length(x_lst: list, y_lst: list):
    x,y = Get_MidPoints(x_lst, y_lst)

    lengths = np.array([])
    for i in range(len(x)-1):
        len = ((x[i]-x[i+1])**2 + (y[i]-y[i+1])**2)**0.5
        lengths = np.append(lengths,len)

    middle_ind = np.floor(len(lengths))
    lengths = np.delete(lengths,middle_ind)

    return np.sum(lengths)

def Get_CLCD(Vx: float, Vy: float, h: float, gamma: float, x_lst: list, y_lst: list, alpha: float):
    T = Get_Tangential(Vx, Vy, h, gamma, x_lst, y_lst, alpha)
    N = Get_Normal(Vx, Vy, h, gamma, x_lst, y_lst, alpha)

    L = N*np.cos(alpha) - T*np.sin(alpha)
    D = T*np.cos(alpha) + N * np.sin(alpha)

    S = Get_length(x_lst,y_lst)

    CL = 2*L/(Get_Density(h)*(Vx**2+Vy**2)*S)
    CD = 2*D/(Get_Density(h)*(Vx**2+Vy**2)*S)

    return CL, CD
