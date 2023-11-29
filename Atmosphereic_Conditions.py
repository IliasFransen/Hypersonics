import numpy as np
from ambiance import Atmosphere as atmos

#Get_SoundSpeed takes the altitude in meters as a float and returns the speed of sound at said altitude in meters per second

def Get_SoundSpeed (h: float):
    a = atmos(h).speed_of_sound
    return float(a)


#Get_Density gives the density at a certain altitiude

def Get_Density(h: float):
    rho = atmos(h).density
    return float(rho)

