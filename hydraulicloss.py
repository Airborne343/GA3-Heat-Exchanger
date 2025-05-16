import numpy as np
from HXobj import HeatExchanger 

x = HeatExchanger(tube_count = 13, baffle_count = 3, type = "tri")

Qguess = 100 #Cm3/s
Pipevelocity = Qguess/x.