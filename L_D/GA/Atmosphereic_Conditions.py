import numpy as np
import fluids.atmosphere as isa
import matplotlib.pyplot as plt

#Get_SoundSpeed takes the altitude in meters as a float and returns the speed of sound at said altitude in meters per second

def Get_SoundSpeed (h: float):
    properties = isa.ATMOSPHERE_1976(h)
    a = (287*1.4*properties.T)**0.5
    return float(a)

#Get_Density gives the density at a certain altitiude

def Get_Density(h: float):
    properties = isa.ATMOSPHERE_1976(h)
    rho = properties.rho
    return float(rho)

def Get_Temperature(h : float):
    
    properties = isa.ATMOSPHERE_1976(h)
    T = properties.T
    
    return T