import numpy as np
from Datatables import Friction, Kc, Ke
from HXobj import HeatExchanger 

x = HeatExchanger(tube_count = 13, baffle_count = 9, type = "tri")


Qguess = 0.45/1000

V_pipe = Qguess/x.pipearea
q_pipe = 1/2 * 1000 * V_pipe ** 2
Re_pipe = x.density * V_pipe * x.tube_ID / x.dynamic_viscosity

P_pipe = q_pipe * x.length / x.tube_ID * Friction(Re_pipe)
P_ends = q_pipe * (Kc(x.sigma, Re_pipe) + Ke(x.sigma, Re_pipe))





print(P_ends , "hi")