import numpy as np
from scipy.integrate import simpson
<<<<<<< Updated upstream

=======
>>>>>>> Stashed changes
from Atmosphereic_Conditions import Get_SoundSpeed, Get_Density


# x and y lists must be points specified from left to right on sideview
# nr of points must be odd

# Get_Theta gives the local incidence angle at each panel on the surface

def Get_Theta(x_lst: list, y_lst: list, alpha : float):
    theta = np.array([])

    for i in range(len(x_lst) - 1):

        slope = (y_lst[i + 1] - y_lst[i]) / (x_lst[i + 1] - x_lst[i])

<<<<<<< Updated upstream
        if alpha > 0:
            if i > len(x_lst) / 2:
                theta_temp = np.arctan(abs(slope)) - alpha
            else:
                theta_temp = np.arctan(abs(slope)) + alpha
        else:
            if i > len(x_lst) / 2:
                theta_temp = np.arctan(abs(slope)) + alpha
            else:
                theta_temp = np.arctan(abs(slope)) - alpha

=======
        if i + 1 < len(x_lst)/2:
            theta_temp = np.arctan(abs(slope)) - alpha
        else:
            theta_temp = np.arctan(abs(slope)) + alpha
>>>>>>> Stashed changes
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


def Get_Sin2Int(x_lst: list, y_lst: list, alpha : float):
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

def Get_Drag(V: float, h: float, gamma: float, x_lst: list, y_lst: list, S : float, alpha : float):
    """
    Cp_max = Get_CpMax(V, h, gamma)
    sin2th = Get_Sin2Int(x_lst, y_lst)

<<<<<<< Updated upstream
    Cn = Cp_max * sin2th

    N = Cn * Get_Density(h) * V**2

    return N
=======
    b = x_lst[-1] - x_lst[0]
    Cn = 1/b * Cp_max * sin2th
    """
    Ct, Cn = Get_CtCn(V, h, gamma, x_lst, y_lst, alpha)
    Cd = Ct * np.sin(alpha) + Cn * np.cos(alpha)
    D = 1/2 * Cd * Get_Density(h) * V**2 * S
    return D
>>>>>>> Stashed changes


def Get_Lift(V: float, h: float, gamma: float, x_lst: list, y_lst: list, S : float, alpha : float):
    """
    CP_max = Get_CpMax(V, h, gamma)
    theta = Get_Theta(x_lst, y_lst)

    y = Get_MidPoints(x_lst, y_lst)[1]

    Cp_local = CP_max * np.sin(theta) ** 2

    integral_left = simpson(Cp_local[:len(Cp_local) // 2 + 1], y[:len(Cp_local) // 2 + 1])

    integral_right = simpson(Cp_local[len(Cp_local) // 2 + 1:], y[len(Cp_local) // 2 + 1:])

    T = (integral_right - integral_left) * Get_Density(h) * V**2

    return T

    c = max(y_lst) - min(y_lst)
    Ct = (integral_right - integral_left) / c
    """
    Ct, Cn = Get_CtCn(V, h, gamma, x_lst, y_lst, alpha)
    Cl = Ct * np.cos(alpha) - Cn * np.sin(alpha)
    L = 1/2 * Cl * Get_Density(h) * V**2 * S
    return L

<<<<<<< Updated upstream
def Get_Lift(N: float, T: float, alpha: float):
    L = -N*np.sin(alpha) + T*np.cos(alpha)
    return L


def Get_Drag(N: float, T: float, alpha: float):
    D = N * np.cos(alpha) + T * np.sin(alpha)
    return D
=======
def Get_CtCn(V : float, h : float, gamma : float, x_lst : list, y_lst : list, alpha : float):
    
    Cp_max = Get_CpMax(V, h, gamma)
    #theta = Get_Theta(x_lst, y_lst, alpha)
    """
    y = Get_MidPoints(x_lst, y_lst)[1]

    Cp_local = Cp_max * np.sin(theta) ** 2

    integral_left = simpson(np.flip(Cp_local[:len(Cp_local) // 2]), np.flip(y[:len(Cp_local) // 2]))

    integral_right = simpson(Cp_local[len(Cp_local) // 2:], y[len(Cp_local) // 2:])

    c = max(y_lst) - min(y_lst)
    Cn = (integral_right - integral_left) / c

    sin2th = Get_Sin2Int(x_lst, y_lst, alpha)
    
    b = x_lst[-1] - x_lst[0]
    Ct = 1/b * Cp_max * sin2th
    """
    S = 12
    R0 = 4.6
    
    d = 33*np.pi/180
    
    Cn = np.pi/2 * np.cos(alpha) * np.sin(alpha) * np.cos(d)**4
    Ct = np.pi/2 * (np.sin(alpha)**2 * np.cos(d)**4/2 - np.cos(alpha) * np.sin(d)**4 + np.cos(alpha)**2)

    return Cn, Ct

# from heatshieldpoints import generate_heatshield_points
#
# # Spacecraft parameters (temp values)
# dhs = 3.9116  # heatshield diameter (m)
# hhs = 0.635  # heatshield height (m)
# S = np.pi * (dhs / 2) ** 2
#
#
# x_lst, y_lst = generate_heatshield_points(dhs, hhs)
# print(Get_Lift(10000, 100000, 1.4, x_lst, y_lst, S))


def Get_length(x_lst: list, y_lst: list):
    x,y = Get_MidPoints(x_lst, y_lst)

    lengths = np.array([])
    for i in range(len(x)-1):
        leng = ((x[i]-x[i+1])**2 + (y[i]-y[i+1])**2)**0.5
        lengths = np.append(lengths,leng)

    middle_ind = int(np.floor(len(lengths)/2))
    lengths = np.delete(lengths,middle_ind)

    return np.sum(lengths)
>>>>>>> Stashed changes
