import numpy as np
import fluids.atmosphere as isa
import matplotlib.pyplot as plt

#Get_SoundSpeed takes the altitude in meters as a float and returns the speed of sound at said altitude in meters per second

def Get_SoundSpeed (h: float):
    
    properies = isa.ATMOSPHERE_1976(h, dT=0.0)
    a = (287*1.4*properies.T)**0.5
    return float(a)


#Get_Density gives the density at a certain altitiude

def Get_Density(h: float):
    properies = isa.ATMOSPHERE_1976(h, dT=0.0)
    rho = properies.rho
    return float(rho)
