import numpy as np
from Datatables import Friction, Kc, Ke, Hotchic
from HXobj import HeatExchanger 
from sympy import symbols, Eq, solve

x = HeatExchanger(tube_count = 13, baffle_count = 9, type = "tri")


Volflow= 0.45/1000

V_pipe = Volflow/x.area_pipe
q_pipe = 1/2 * x.density * V_pipe ** 2
Re_pipe = x.density * V_pipe * x.tube_ID / x.dynamic_viscosity

P_pipe = q_pipe * x.length / x.tube_ID * Friction(Re_pipe)
P_ends = q_pipe * (Kc(x.sigma, Re_pipe) + Ke(x.sigma, Re_pipe))
P_nozzle = 0.5 * x.density * (Volflow/x.area_nozzle)**2 * 2

Ploss = P_pipe + P_ends + P_nozzle

x = symbols('x')
poly_expr = sum(c * x**i for i, c in enumerate(reversed(Hotchic(Ploss))))
roots = solve(Eq(poly_expr, Ploss/(Volflow) * x), x)

print(roots)