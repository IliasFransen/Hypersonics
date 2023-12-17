import numpy as np
from scipy.integrate import simpson, trapezoid
from Atmosphereic_Conditions import Get_SoundSpeed, Get_Density

"""
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


def Get_MidPoints(x_lst: list,
                  y_lst: list):  # Get_MidPoints computes the x and y coodinate of the midpoints of each panel

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
"""

# Get_CpMax gives the cpmax from modified newtonian theory

def Get_CpMax(V: float, h: float, gamma: float):
    V_inf = V
    a = Get_SoundSpeed(h)
    M_inf = V_inf / a

    part1 = (1 - gamma + 2 * gamma * M_inf ** 2) / (gamma + 1)

    part2 = (((gamma + 1) ** 2 * M_inf ** 2) / (4 * M_inf ** 2 * gamma - 2 * (gamma - 1))) ** (gamma / (gamma - 1))

    Cp_max = 2 / (gamma * M_inf ** 2) * (part1 * part2 - 1)

    return Cp_max


"""
def Get_Normal(V: float, h: float, gamma: float, x_lst: list, y_lst: list, alpha: float):

    Cp_max = Get_CpMax(V, h, gamma)
    sin2th = Get_Sin2Int(x_lst, y_lst, alpha)

    b = abs(x_lst[-1] - x_lst[0])
    Cn = 1/b * Cp_max * sin2th

    return Cn

def Get_Tangential(V: float, h: float, gamma: float, x_lst: list, y_lst: list, alpha: float):
    CP_max = Get_CpMax(V, h, gamma)
    theta = Get_Theta(x_lst, y_lst, alpha)

    y = Get_MidPoints(x_lst, y_lst)[1]

    Cp_local = CP_max * np.sin(theta) ** 2

    integral_left = simpson(np.flip(Cp_local[:len(Cp_local) // 2]), np.flip(y[:len(Cp_local) // 2]))

    integral_right = simpson(Cp_local[len(Cp_local) // 2:], y[len(Cp_local) // 2:])

    c = max(y_lst) - min(y_lst)
    Ct = (integral_right - integral_left) / c

    return Ct

def Get_Drag(V: float, h: float, gamma: float, x_lst: list, y_lst: list, alpha: float, S):
    Cn = Get_Normal(V, h, gamma, x_lst, y_lst, alpha)
    Ct = Get_Tangential(V, h, gamma, x_lst, y_lst, alpha)

    Cd = Ct*np.sin(alpha)+Cn*np.cos(alpha)

    D = 1/2 * Cd * Get_Density(h) * V**2 * S
    return D


def Get_Lift(V: float, h: float, gamma: float, x_lst: list, y_lst: list, alpha: float, S):
    Cn = Get_Normal(V, h, gamma, x_lst, y_lst, alpha)
    Ct = Get_Tangential(V, h, gamma, x_lst, y_lst, alpha)

    Cl = Ct * np.cos(alpha) - Cn * np.sin(alpha)

    L = 1/2 * Cl * Get_Density(h) * V**2 * S
    return L
"""


def Get_Ca(V: float, h: float, gamma: float, R0: float, delta: float, alpha: float, S: float):
    CP_max = Get_CpMax(V, h, gamma)
    d = delta
    if 0 <= alpha <= d:
        C = np.pi / 2 * (
                    np.sin(alpha) ** 2 * np.cos(d) ** 4 / 2 - np.cos(alpha) ** 2 * np.sin(d) ** 4 + np.cos(alpha) ** 2)
    elif d < alpha <= np.pi - d:
        C1 = np.cos(alpha) * np.arccos(np.sin(d)/np.sin(alpha))
        C2 = (np.sin(alpha) ** 2 * np.cos(d) ** 4 / 2 - np.cos(alpha) ** 2 * np.sin(d) ** 4 + np.cos(alpha) ** 2) * (np.pi / 2 + np.arcsin(np.tan(d)/np.tan(alpha)))
        C3 = np.cos(alpha) * np.sin(d)/2 * (1-3*np.sin(d)**2) * np.sqrt(np.sin(alpha)**2 - np.sin(d)**2)
        C = 1 / 2 * (C1 + C2 + C3)
    elif np.pi - d < alpha <= np.pi:
        C = 0
    else:
        print("error")
        exit()
    Ca = CP_max / S * R0 ** 2 * C
    return Ca


def Get_Cn(V: float, h: float, gamma: float, R0: float, delta: float, alpha: float, S: float):
    CP_max = Get_CpMax(V, h, gamma)
    d = delta
    if 0 <= alpha <= d:
        C = np.pi / 2 * np.cos(alpha) * np.sin(alpha) * (np.cos(d)) ** 4
    elif d < alpha <= np.pi - d:
        C1 = np.arccos(np.sin(d) / np.sin(alpha))
        C2 = np.cos(alpha) * np.cos(d) ** 4 * (np.pi / 2 + np.arcsin(np.tan(d) / np.tan(alpha)))
        C3 = np.sin(d) / 3 * np.sqrt(np.sin(alpha) ** 2 - np.sin(d) ** 2) * (
                np.sin(d) ** 2 * (3 - 1 / np.sin(alpha) ** 2) - 5)
        C = np.sin(alpha) / 2 * (C1 + C2 + C3)
    elif np.pi - d < alpha <= np.pi:
        C = 0
    else:
        print("error")
        exit()
    Cn = CP_max / S * R0 ** 2 * C
    return Cn


def Get_Drag(V: float, h: float, gamma: float, R0: float, delta: float, alpha: float, S):
    Cn = Get_Cn(V, h, gamma, R0, delta, alpha, S)
    Ca = Get_Ca(V, h, gamma, R0, delta, alpha, S)

    Cd = Cn*np.sin(alpha)+Ca*np.cos(alpha)

    D = 1/2 * Cd * Get_Density(h) * V**2 * S
    return D, Cd


def Get_Lift(V: float, h: float, gamma: float, R0: float, delta: float, alpha: float, S):
    Cn = Get_Cn(V, h, gamma, R0, delta, alpha, S)
    Ca = Get_Ca(V, h, gamma, R0, delta, alpha, S)

    Cl = Cn * np.cos(alpha) - Ca * np.sin(alpha)

    L = 1/2 * Cl * Get_Density(h) * V**2 * S
    return L, Cl
