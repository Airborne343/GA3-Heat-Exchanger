import numpy as np
from Datatables import Friction, Kc, Ke, Hotchic, Coldchic
from HXobj import HeatExchanger 
from sympy import symbols, Eq, solve
import matplotlib.pyplot as plt

x = HeatExchanger(tube_count = 10, baffle_count = 9, type = "tri")


y = np.linspace(0, 1/1000, 100)
HXchic = []
Compressorchic = []
   
def Pressuredropcalc(Volflow):
    V_pipe = Volflow/x.area_pipe
    q_pipe = 1/2 * x.density * V_pipe ** 2
    Re_pipe = x.density * V_pipe * x.tube_ID / x.dynamic_viscosity

    P_pipe = q_pipe * x.length / x.tube_ID * Friction(Re_pipe)
    P_ends = q_pipe * (Kc(x.sigma, Re_pipe) + Ke(x.sigma, Re_pipe))
    P_nozzle = 0.5 * x.density * (Volflow/x.area_nozzle)**2 * 2

    P_total = P_pipe + P_ends + P_nozzle

    return P_total - Hotchic(Volflow)


