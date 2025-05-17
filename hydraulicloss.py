import numpy as np
from Datatables import Friction, Kc, Ke, Hotchic, Coldchic
from HXobj import HeatExchanger 
import matplotlib.pyplot as plt

Hx = HeatExchanger(tube_count = 6, baffle_count = 9, type = "triangle")


def P_drop_hot(mhot, Hx): #Function to calculate hot side pressure drop

    #Working out velocity, dynamic pressure and reynolds number inside tube
    Volflow = mhot/Hx.density
    V_pipe = Volflow/Hx.area_tubes
    q_pipe = 1/2 * Hx.density * V_pipe ** 2
    Re_pipe = Hx.density * V_pipe * Hx.tube_ID / Hx.dynamic_viscosity

    #Working out three sources of pressure drop, tube friction, end losses for tubes and inlet nozzle drops
    P_pipe = q_pipe * Hx.length / Hx.tube_ID * Friction(Re_pipe)
    P_ends = q_pipe * (Kc(Hx.sigma, Re_pipe) + Ke(Hx.sigma, Re_pipe))
    P_nozzle = 0.5 * Hx.density * (Volflow/Hx.area_nozzle)**2 * 2

    #Summing pressure losses and combining
    P_total = P_pipe + P_ends + P_nozzle
    return [P_total - Hotchic(qdot = Volflow), P_total, Hotchic(qdot = Volflow)]


def iteration(pressure, Hx, initialmass = 0.45, tol = 0.005, maxiter = 15):
    massflow = initialmass

    for i in range(maxiter):
        #Working out 
        p = pressure(massflow, Hx)[0]
        dp_dq = (pressure(massflow * 1.05, Hx)[0] - pressure(massflow * 0.95, Hx)[0]) / (2 * 0.05 * massflow)

        #Newton raphson method to 
        massflow_new = massflow - p / dp_dq
        
        print(massflow)
        if abs(massflow - massflow_new) < massflow * tol:
            return massflow

        massflow = massflow_new

    raise("Maximum iterations reached without convergence")


print(iteration(P_drop_hot, Hx))

flows =  np.linspace(0.01,0.6, 50)
HXpressures = []
Pumppressures = []
for k in flows:
    HXpressures.append(P_drop_hot(k, Hx)[1])
    Pumppressures.append(P_drop_hot(k, Hx)[2])



plt.plot(flows, HXpressures)
plt.plot(flows, Pumppressures)

plt.show()


import numpy as np
from HXobj import HeatExchanger
from scipy.optimize import fsolve
from Datatables import Friction, Kc, Ke, Coldchic

hex = HeatExchanger(tube_count = 13, baffle_count = 9, type = "square")


def calculate_pressure_loss(mcold):
    V_shell = mcold/(hex.density * hex.area_shell)
    charc_D_shell = hex.D_shell * (hex.area_shell/hex.area_pipe)
    Re_shell = (hex.density * V_shell * charc_D_shell)/hex.dynamic_viscosity
    ploss_shell = 4 * hex.a * ((Re_shell)**(-0.15)) * hex.tube_count * hex.density * ((V_shell)**2)

    V_nozzle = mcold/(hex.density * hex.area_nozzle)
    ploss_nozzle = 2 * 0.5 * hex.density * (V_nozzle**2)

    ploss_cold_tot = ploss_nozzle + ploss_shell
    return ploss_cold_tot

def compressor_pressure(mcold):
    qcold = mcold/hex.density
    coeffs = Coldchic(qdot=None)
    pcomp = np.polyval(coeffs, qcold)
    return pcomp

# Iterative Optimization Parameters
mcold = 0.5
convergence_limit = 1e-3
step_factor = 0.01
max_iterations = 100
iteration = 0

while iteration < max_iterations:
    ploss = calculate_pressure_loss(mcold)
    pcomp = compressor_pressure(mcold)
    difference = ploss - pcomp
    
    print(f"Iteration {iteration+1}: Mass Flow Rate = {mcold:.4f} kg/s, "
          f"Pressure Loss = {ploss:.2f} Pa, Compressor Pressure = {pcomp:.2f} Pa, "
          f"Difference = {difference:.4f}")
    
    if abs(difference) <= convergence_limit:
        print("\nConvergence achieved!")
        break

    mcold -= step_factor * difference
    
    iteration += 1

if iteration >= max_iterations:
    print("\nMaximum iterations reached without convergence.")

print(f"\nOptimal Mass Flow Rate: {mcold:.4f} kg/s")
print(f"Final Pressure Loss: {ploss:.2f} Pa")
print(f"Final Compressor Pressure: {pcomp:.2f} Pa")

