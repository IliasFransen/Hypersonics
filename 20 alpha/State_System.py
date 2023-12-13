import numpy as np
from Modified_Newton import Get_Drag, Get_Lift
from Atmosphereic_Conditions import Get_Density, Get_SoundSpeed, Get_Temperature
from scipy.integrate import simpson


# StateUpdate is gives the system of differential equations that needs to be S

def StateUpdate(states: list, t: list, g: float, m: float, gamma: float, x_lst: list, y_lst: list, alpha: float, S):
    V, fpa, h, x = states

    Re = 6378 * 1000
    r = Re + h
    g = g * (Re / r) ** 2

    L = Get_Lift(V, h, gamma, x_lst, y_lst, alpha, S)[0]
    D = Get_Drag(V, h, gamma, x_lst, y_lst, alpha, S)[0]

    dVdt = -D / m + g * np.sin(fpa)
    dfpadt = -L / (V * m) + (g / V - V / r) * np.cos(fpa)
    dhdt = -V * np.sin(fpa)
    dxdt = V * np.cos(fpa)

    return [dVdt, dfpadt, dhdt, dxdt]


def getViscosity(T: float):
    # Sutherland's Law

    T_ref = 288
    mu_ref = 1.789 * 10 ** (-5)
    S = 110

    mu = mu_ref * (T / T_ref) ** 1.5 * (T_ref + S) / (T + S)

    return mu


def normalShockProperties(h: float, M: float, gamma: float):
    rho1 = Get_Density(h)
    rho2 = rho1 * ((gamma + 1) * M ** 2) / (2 + (gamma - 1) * M ** 2)

    T1 = Get_Temperature(h)
    T2 = T1 * (1 + 2 * gamma / (gamma + 1) * (M ** 2 - 1)) * (2 + (gamma - 1) * M ** 2) / ((gamma + 1) * M ** 2)

    mu2 = getViscosity(T2)

    return rho2, mu2


def getStagHeatFlux(h: float, M: float, gamma: float, Pr: float, R0: float):
    rho_inf = Get_Density(h)

    rho_e, mu_e = normalShockProperties(h, M, gamma)

    a = Get_SoundSpeed(h)

    U_inf = M * a

    q = 2 ** (-3 / 4) * 0.76 * Pr ** (-0.6) / np.sqrt(R0) * np.sqrt(rho_e ** 0.5 * rho_inf ** 0.5 * mu_e) * U_inf ** 2.5

    return q


def getStagHeatLoad(q: list, t: list, i: int):
    Q = simpson(q[0:i + 1], t[0:i + 1])

    return Q


def getGForce(dVdt: float, g: float):
    ng = abs(dVdt) / g

    return ng
