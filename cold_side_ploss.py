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