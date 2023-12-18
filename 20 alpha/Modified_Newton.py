import numpy as np
from scipy.integrate import simpson, trapezoid
from Atmosphereic_Conditions import Get_SoundSpeed, Get_Density


# Get_CpMax gives the cpmax from modified newtonian theory

def Get_CpMax(V: float, h: float, gamma: float):
    V_inf = V
    a = Get_SoundSpeed(h)
    M_inf = V_inf / a

    part1 = (1 - gamma + 2 * gamma * M_inf ** 2) / (gamma + 1)

    part2 = (((gamma + 1) ** 2 * M_inf ** 2) / (4 * M_inf ** 2 * gamma - 2 * (gamma - 1))) ** (gamma / (gamma - 1))

    Cp_max = 2 / (gamma * M_inf ** 2) * (part1 * part2 - 1)

    return Cp_max


def Get_Ca(V: float, h: float, gamma: float, R0: float, delta: float, alpha: float, S: float):
    CP_max = Get_CpMax(V, h, gamma)
    d = delta
    if 0 <= alpha <= d:
        C = np.pi / 2 * (
                np.sin(alpha) ** 2 * np.cos(d) ** 4 / 2 - np.cos(alpha) ** 2 * np.sin(d) ** 4 + np.cos(alpha) ** 2)
    elif d < alpha <= np.pi - d:
        C1 = np.cos(alpha) * np.arccos(np.sin(d) / np.sin(alpha))
        C2 = (np.sin(alpha) ** 2 * np.cos(d) ** 4 / 2 - np.cos(alpha) ** 2 * np.sin(d) ** 4 + np.cos(alpha) ** 2) * (
                np.pi / 2 + np.arcsin(np.tan(d) / np.tan(alpha)))
        C3 = np.cos(alpha) * np.sin(d) / 2 * (1 - 3 * np.sin(d) ** 2) * np.sqrt(np.sin(alpha) ** 2 - np.sin(d) ** 2)
        C = 1 / 2 * (C1 + C2 + C3)
    elif np.pi - d < alpha <= np.pi:
        C = 0
    else:
        print("error")
        exit()
    Ca = CP_max / np.pi * C
    return Ca


def Get_Cn(V: float, h: float, gamma: float, R0: float, delta: float, alpha: float, S: float):
    CP_max = Get_CpMax(V, h, gamma)
    d = delta
    if 0 <= alpha <= d:
        C = np.pi / 2 * np.cos(alpha) * np.sin(alpha) * np.cos(d) ** 4
    elif d < alpha <= np.pi - d:
        C1 = np.arccos(np.sin(d) / np.sin(alpha))
        C2 = np.cos(alpha) * np.cos(d) ** 4 * (np.pi / 2 + np.arcsin(np.tan(d) / np.tan(alpha)))
        C3 = np.sin(d) / 3 * (np.sin(d) ** 2 * (3 - 1 / np.sin(alpha) ** 2) - 5) * np.sqrt(
            np.sin(alpha) ** 2 - np.sin(d) ** 2)
        C = np.sin(alpha) / 2 * (C1 + C2 + C3)
    elif np.pi - d < alpha <= np.pi:
        C = 0
    else:
        print("error")
        exit()
    Cn = CP_max / np.pi * C
    return Cn


def Get_Drag(V: float, h: float, gamma: float, R0: float, delta: float, alpha: float, S):
    Cn = Get_Cn(V, h, gamma, R0, delta, alpha, S)
    Ca = Get_Ca(V, h, gamma, R0, delta, alpha, S)

    Cd = Cn * np.sin(alpha) + Ca * np.cos(alpha)

    D = 1 / 2 * Cd * Get_Density(h) * V ** 2 * S
    return D, Cd


def Get_Lift(V: float, h: float, gamma: float, R0: float, delta: float, alpha: float, S):
    Cn = Get_Cn(V, h, gamma, R0, delta, alpha, S)
    Ca = Get_Ca(V, h, gamma, R0, delta, alpha, S)

    Cl = Cn * np.cos(alpha) - Ca * np.sin(alpha)

    L = 1 / 2 * Cl * Get_Density(h) * V ** 2 * S
    return L, Cl
