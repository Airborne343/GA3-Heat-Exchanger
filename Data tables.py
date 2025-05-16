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

def Kc(sigma, Re):
    #These fitting is determined via Figure 8 in the handout
    poly_m = np.poly1d([-3.4112395e9,  1.3174895e6,   -97.6365546,  -0.4])
    poly_c = np.poly1d([1.15879727e10, -7.49109769e6,  1633.23004,  0.4])
    
    x = 1.0 / Re
    return poly_m(x) * sigma + poly_c(x)

def Ke(sigma, Re):
    #These fitting is determined via Figure 8 in the handout
    # Cubic-in-x coefficients for each sigma^k term (descending powers of x)
    # d3(x) multiplies sigma^3
    _coef_d3 = np.array([-1.11383048e+11,   4.92744792e+07,  -4.37578675e+03,  -5.95238095e-02])
    # d2(x) multiplies sigma^2
    _coef_d2 = np.array([ 1.93887205e+11,  -8.58331318e+07,   7.75224537e+03,   1.05357143e+00])
    # d1(x) multiplies sigma^1
    _coef_d1 = np.array([-9.66470141e+10,   4.53015098e+07,  -5.10931576e+03,  -1.99404762e+00])
    # d0(x) multiplies sigma^0 (constant term)
    _coef_d0 = np.array([-5.18889235e-04,   2.65109503e-07,  -3.10673330e-11,   1.00000000e+00])

    σ, R = np.broadcast_arrays(sigma, Re)
    # x = 1/Re, but x=0 when Re==inf to hit the asymptote curves
    x = np.where(np.isinf(R), 0.0, 1.0 / R)

    # Evaluate each cubic d_i(x) via polyval
    d3 = np.polyval(_coef_d3, x)
    d2 = np.polyval(_coef_d2, x)
    d1 = np.polyval(_coef_d1, x)
    d0 = np.polyval(_coef_d0, x)

    # Return the combined cubic‐in‐sigma polynomial
    return d3*σ**3 + d2*σ**2 + d1*σ + d0
