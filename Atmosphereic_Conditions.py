import numpy as np
import isacalc as isa

#Get_SoundSpeed takes the altitude in meters as a float and returns the speed of sound at said altitude in meters per second

def Get_SoundSpeed (h: float):
    atmosphere = isa.get_atmosphere()
    T, P, d, a, mu = isa.calculate_at_h(h, atmosphere)
    return float(a)


#Get_Density gives the density at a certain altitiude

def Get_Density(h: float):
    atmosphere = isa.get_atmosphere()
    T, P, rho, a, mu = isa.calculate_at_h(h, atmosphere)
    return float(rho)

