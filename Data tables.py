import numpy as np
from scipy.optimize import fsolve

#Tables for compressor chics
def Coldchic(p = None, mdot = None): #cm3/s to Pa
    if p == None: 
        poly = np.poly1d([-2.03241506e-07,  2.18666233e-04, -1.17200431e-01, -4.34105885e+01,
  7.07758580e+04])
        return poly(mdot)
    
    elif mdot == None:
        poly = np.poly1d([-1.12242099e-17,  1.58646220e-12, -1.69351833e-07, -2.28742118e-03,
  7.29640030e+02])
        return poly(p)
    
    else:
        raise("Invalid input, only input one of pressure or mass flow")
          

def Hotchic(p = None, mdot = None): #cm3/s to Pa
    if p == None: 
        poly = np.poly1d([1.30073533e-06, -1.00398975e-03,  9.91494147e-02, -8.16106902e+01,
  6.23799629e+04])
        return poly(mdot)
    
    elif mdot == None:
        poly = np.poly1d([-8.70199042e-18, -4.64584399e-13,  6.03677597e-08, -7.97101169e-03,
  5.06006301e+02])
        return poly(p)
    
    else:
        raise("Invalid input, only input one of pressure or mass flow")

#Tables for friction factor
def Friction(Re = None, epsilon = 0.002, D = 8 ): #https://www.engineeringtoolbox.com/surface-roughness-ventilation-ducts-d_209.html
    #Using Haarland equation
    return (-1.8 * np.log10((epsilon / D / 3.7)**1.11 + 6.9 / Re))**-2


#