import numpy as np
from math import atan2

#Get_Angles takes the velocity components and the angle of attack in RADIANS and returns the angle between velocity vectors and coordianate system and the angle between SC and coodinate system

def Get_Angles(Vx: float, Vy: float, Beta : float, alpha: float):
    #Beta = atan2(Vy, Vx)
    eta = Beta - alpha
    return eta
