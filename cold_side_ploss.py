import numpy as np
from HXobj import HeatExchanger
from Datatables import Friction, Kc, Ke, Coldchic

hex = HeatExchanger(tube_count = 13, baffle_count = 9, type = "square")


def calculate_pressure_loss(mcold):
    try:
        V_shell = mcold / (hex.density * hex.area_shell)
        charc_D_shell = hex.D_shell * (hex.area_shell / hex.area_pipe)
        Re_shell = (hex.density * V_shell * charc_D_shell) / hex.dynamic_viscosity

        if Re_shell <= 0:
            return float('inf')

        f_shell = 4 * hex.a * (Re_shell ** -0.15)
        ploss_shell = f_shell * hex.tube_count * hex.density * (V_shell ** 2)
        
        V_nozzle = mcold / (hex.density * hex.area_nozzle)
        ploss_nozzle = 2 * 0.5 * hex.density * (V_nozzle ** 2)
        ploss_cold_tot = ploss_nozzle + ploss_shell
        return ploss_cold_tot
    
    except Exception as e:
        print(f"Error in pressure loss calculation: {e}")
        return float('inf')

def compressor_pressure(mcold):
    qcold = mcold/hex.density
    coeffs = Coldchic(qdot=None)
    pcomp = np.polyval(coeffs, qcold)
    return pcomp

def compressor_massflow(pressure):
    coeffs = Coldchic(p=None)
    mcold_comp = np.polyval(coeffs, pressure)
    return mcold_comp

mcold_guess = 0.5
limit_factor = 0.3
n_iter = 1
convergence_thresh = 100

#to prevent instability
upper_bound = 100
lower_bound = 0.01

for i in range(n_iter):
    ploss_cold = calculate_pressure_loss(mcold_guess)
    pcomp = compressor_pressure(mcold_guess)
    error = abs(pcomp - ploss_cold)

    if error < convergence_thresh:
        print(f"Convergence achieved in {i+1} iterations:")
        print(f"Compressor Pressure Rise: {pcomp:.2f} Pa")
        print(f"Cold Side Pressure Loss: {ploss_cold:.2f} Pa")
        print(f"Cold Side Mass Flow Rate: {mcold_guess:.4f} kg/s")
        break

    new_mcold = compressor_massflow((pcomp + ploss_cold) / 2)
    mcold_guess = mcold_guess + limit_factor * (new_mcold - mcold_guess)
    mcold_guess = max(min(mcold_guess, upper_bound), lower_bound)

else:
    print("No convergence after maximum iterations.")
    print(f"Final estimated mass flow rate: {mcold_guess:.4f} kg/s")