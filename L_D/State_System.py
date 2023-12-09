import numpy as np
from Modified_Newton import Get_Drag, Get_Lift, Get_Tangential, Get_Normal
from Atmosphereic_Conditions import Get_Density
from scipy.integrate import simpson


# StateUpdate is gives the system of differential equations that needs to be S

def StateUpdate(states: list, t: list, g: float, m: float, alpha: float, gamma: float, x_lst: list, y_lst: list):
    V, fpa, h = states
    
    r = 6378 * 1000 + h

    N = Get_Normal(V, h, gamma, x_lst, y_lst, alpha)
    T = Get_Tangential(V, h, gamma, x_lst, y_lst, alpha)

    L = Get_Lift(N, T, V, h, alpha, x_lst, y_lst)
    D = Get_Drag(N, T, V, h, alpha, x_lst, y_lst)
    
    """
    dVdt = -D / m - g * np.sin(fpa)
    dfpadt = L / (V * m) - (g / V - V / r) * np.cos(fpa)
    dhdt = V*np.sin(fpa)
    """
    
    B = 372 # Ballistic coeff.
    rho = Get_Density(h)

    dVdt = -rho * V**2 / (2 * B) + g * np.sin(fpa)
    dfpadt = -rho * V * L/D / (2 * B) + g * np.cos(fpa) / V - V * np.cos(fpa) / r
    dhdt = -V * np.sin(fpa)
    
    return [dVdt, dfpadt, dhdt]

# Get_VInit gives intial velocities, fpa must be measured from POSITIVE X AXIS

def Get_VInit(V0: float, fpa0: float):
    V0 = abs(V0)

    Vx0 = V0 * np.cos(fpa0)
    Vy0 = V0 * np.sin(fpa0)

    return Vx0, Vy0

def getStagHeatFlux(k : float, R0 : float, rho : float, V : float):
    
    q = k * np.sqrt(rho/R0) * V**3
    
    return q

def getStagHeatLoad(q : list, t : list, i : int):
    
    Q = simpson(q[0:i + 1], t[0:i + 1])
    
    return Q

def getGForce(dVdt : float, g : float):
    
    ng = abs(dVdt) / g
    
    return ng