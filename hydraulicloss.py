import numpy as np
from Datatables import Friction
from HXobj import HeatExchanger 

x = HeatExchanger(tube_count = 13, baffle_count = 3, type = "tri")

Qguess = 0.3/1000
Pipevelocity = Qguess/x.pipearea
Pipedynpressure = 1/2 * 1000 * Pipevelocity ** 2
Re_pipe = 

P_pipe = Pipedynpressure * x.length / x.tube_ID * Friction(Re_pipe)

print(P_pipe)