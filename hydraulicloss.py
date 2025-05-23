import numpy as np
from Datatables import Friction, Kc, Ke, Hotchic, Coldchic
from HXobj import HeatExchanger 
import matplotlib.pyplot as plt

Hx = HeatExchanger(length = 0.35, tube_count = 10, baffle_count = 9, type = "60", passes = 2, N_shell = 1, pitch = 12/1000, bundle_height = 0.4)


def P_drop_hot(mhot, Hx): #Function to calculate hot side pressure drop

    #Working out velocity, dynamic pressure and reynolds number inside tube
    Volflow = mhot/Hx.density
    V_tube = Volflow/(Hx.area_tubes/Hx.passes)
    q_tube = 1/2 * Hx.density * V_tube ** 2
    Re_pipe = Hx.density * V_tube * Hx.tube_ID / Hx.dynamic_viscosity

    #Working out three sources of pressure drop, tube friction, end losses for tubes and inlet nozzle drops
    P_tube = q_tube * Hx.length / Hx.tube_ID * Friction(Re_pipe)
    P_ends = q_tube * (Kc(Hx.sigma, Re_pipe) + Ke(Hx.sigma, Re_pipe)) * Hx.passes
    P_nozzle = 0.5 * Hx.density * (Volflow/Hx.area_nozzle)**2 * 2
    P_pipe = 9320 * (1000*Volflow/0.436)**2 #This is taken from the max flow point of the hot chic

    #Summing pressure losses and combining
    P_total = P_tube + P_ends + P_nozzle + P_pipe
    return [P_total - Hotchic(qdot = Volflow), P_total, Hotchic(qdot = Volflow)]

def P_drop_cold_old(mcold, Hx):
    #Calculating velocity, characteristic diameter and shell reynolds number
    Volflow = mcold/Hx.density
    V_shell = Volflow/Hx.area_shell

    if Hx.type == "triangle":
        D_e_const1 = 1.10
        D_e_const2 = 0.917
    elif Hx.type == "square":
        D_e_const1 = 1.27
        D_e_const2 = 0.785

    charc_D_shell = D_e_const1 * (Hx.pitch ** 2 - D_e_const2 * Hx.tube_OD **2) / Hx.tube_OD
    Re_shell = (Hx.density * V_shell * charc_D_shell)/Hx.dynamic_viscosity

    #Working out pressure losses due to shell and nozzle exit
    ploss_shell = 4 * Hx.a * ((Re_shell)**(-0.15)) * Hx.tube_count * Hx.density * ((V_shell)**2)
    V_nozzle = mcold/(Hx.density * Hx.area_nozzle)
    ploss_nozzle = 2 * 0.5 * Hx.density * (V_nozzle**2)
    P_pipe = 15840 * (1000*Volflow/0.6580)**2    #This is taken from the max flow point of the cold chic
    ploss_cold_tot = ploss_nozzle + ploss_shell + P_pipe

    return [ploss_cold_tot - Coldchic(qdot = Volflow), ploss_cold_tot, Coldchic(qdot = Volflow)]

def P_drop_cold(mcold, Hx):
    #Calculating velocity, characteristic diameter and shell reynolds number
    massvelocity = mcold / Hx.area_shell

    if Hx.type == "60":
        D_e_const1 = 1.10
        D_e_const2 = 0.917
    elif Hx.type == "0" or Hx.type == "45":
        D_e_const1 = 1.27
        D_e_const2 = 0.785

    charc_D_shell = D_e_const1 * (Hx.pitch ** 2 - D_e_const2 * Hx.tube_OD **2) / Hx.tube_OD
    Re_shell = (massvelocity * charc_D_shell)/Hx.dynamic_viscosity

    euler_coeff_dict = {"0": [0.267, 0.249e4, -0.927e7, 0.10e11, 0], "45" :[0.245, 0.339e4, -0.938e7, 0.132e11, -0.599e13], "60" :[0.245, 0.339e4, -0.938e7, 0.132e11, -0.599e13]}
    eulercoeff = euler_coeff_dict[Hx.type][0] + euler_coeff_dict[Hx.type][1]/ (Re_shell) + euler_coeff_dict[Hx.type][2]/ (Re_shell**2) + euler_coeff_dict[Hx.type][3]/ (Re_shell**3) + euler_coeff_dict[Hx.type][4]/ (Re_shell**4) 

    maxdynamicpressure = 0.5*massvelocity**2 / Hx.density

    P_drop_crossflow = maxdynamicpressure*eulercoeff*(Hx.baffle_count + 1) * Hx.N_shell * Hx.rows

    windowarea = (Hx.D_shell**2/8)*(2*np.arccos(1-2*Hx.baffle_height) - np.sin(2*np.arccos(1-2*Hx.baffle_height)))
    P_drop_window = 2*maxdynamicpressure*(Hx.area_shell/windowarea) * (1 + 0.3*(1-Hx.crossflow_prop)) * Hx.baffle_count * Hx.N_shell

    

    V_nozzle = mcold/(Hx.density * Hx.area_nozzle)
    ploss_nozzle = 2 * 0.5 * Hx.density * (V_nozzle**2)
    P_pipe = 15840 * (1000*mcold/(0.6580*Hx.density))**2    #This is taken from the max flow point of the cold chic

    ploss_cold_tot = P_drop_crossflow +P_drop_window + ploss_nozzle +  P_pipe


    return [ploss_cold_tot - Coldchic(qdot = mcold/Hx.density), ploss_cold_tot, Coldchic(qdot = mcold/Hx.density)]

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

P_drop_cold(0.5, Hx)


mhot = iteration(P_drop_hot, Hx)
mcold = iteration(P_drop_cold, Hx)
print(mhot, mcold)
print(Hotchic(mhot/Hx.density), Coldchic(mcold/Hx.density))