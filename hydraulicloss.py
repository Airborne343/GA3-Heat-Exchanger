import numpy as np
from Datatables import Friction, Kc, Ke, Hotchic, Coldchic
from HXobj import HeatExchanger 
import matplotlib.pyplot as plt

Hx = HeatExchanger(length = 0.35, tube_count = 1, baffle_count = 9, type = "triangle", passes = 2, N_shell = 1)


def P_drop_hot(mhot, Hx): #Function to calculate hot side pressure drop

    #Working out velocity, dynamic pressure and reynolds number inside tube
    Volflow = mhot/Hx.density
    V_pipe = Volflow/(Hx.area_tubes/Hx.passes)
    q_pipe = 1/2 * Hx.density * V_pipe ** 2
    Re_pipe = Hx.density * V_pipe * Hx.tube_ID / Hx.dynamic_viscosity

    #Working out three sources of pressure drop, tube friction, end losses for tubes and inlet nozzle drops
    P_pipe = q_pipe * Hx.length / Hx.tube_ID * Friction(Re_pipe)
    P_ends = q_pipe * (Kc(Hx.sigma, Re_pipe) + Ke(Hx.sigma, Re_pipe)) * Hx.passes
    P_nozzle = 0.5 * Hx.density * (Volflow/Hx.area_nozzle)**2 * 2

    #Summing pressure losses and combining
    P_total = P_pipe + P_ends + P_nozzle
    return [P_total - Hotchic(qdot = Volflow), P_total, Hotchic(qdot = Volflow)]

def P_drop_cold(mcold, Hx):
    #Calculating velocity, characteristic diameter and shell reynolds number
    Volflow = mcold/Hx.density
    V_shell = Volflow/Hx.area_shell

    if Hx.type == "triangle":
        D_e_const1 = 1.10
        D_e_const2 = 0.917
    if Hx.type == "square":
        D_e_const1 = 1.27
        D_e_const2 = 0.785
    else:
        raise("invalid type given, should be triangle or square")

    charc_D_shell = D_e_const1 * (Hx.pitch ** 2 - D_e_const2 * Hx.tube_OD **2) / Hx.tube_OD
    Re_shell = (Hx.density * V_shell * charc_D_shell)/Hx.dynamic_viscosity

    #Working out pressure losses due to shell and nozzle exit
    ploss_shell = 4 * Hx.a * ((Re_shell)**(-0.15)) * Hx.tube_count * Hx.density * ((V_shell)**2)
    V_nozzle = mcold/(Hx.density * Hx.area_nozzle)
    ploss_nozzle = 2 * 0.5 * Hx.density * (V_nozzle**2)


    ploss_cold_tot = ploss_nozzle + ploss_shell

    return [ploss_cold_tot - Coldchic(qdot = Volflow), ploss_cold_tot, Coldchic(qdot = Volflow)]

def iteration(pressurefunction, Hx, initialmass = 0.45, tol = 0.005, maxiter = 15):
    massflow = initialmass

    for i in range(maxiter):
        #Working out 
        p = pressurefunction(massflow, Hx)[0]
        dp_dq = (pressurefunction(massflow * 1.05, Hx)[0] - pressurefunction(massflow * 0.95, Hx)[0]) / (2 * 0.05 * massflow)

        #Newton raphson method to calculate pressurefunction drop
        massflow_new = massflow - p / dp_dq
        
        if abs(massflow - massflow_new) < massflow * tol:
            return massflow

        massflow = massflow_new

    raise("Maximum iterations reached without convergence")


mhot = iteration(P_drop_hot, Hx)
mcold = iteration(P_drop_cold, Hx)
print(mhot, mcold)
print(Hotchic(mhot/Hx.density), Coldchic(mcold/Hx.density))