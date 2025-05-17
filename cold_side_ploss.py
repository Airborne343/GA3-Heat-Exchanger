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
    return qcold, pcomp

mcold_guess = 0.5
qcold_guess, comp_ploss = compressor_pressure(mcold_guess)
ploss_total = calculate_pressure_loss(mcold_guess)

print(qcold_guess, comp_ploss, ploss_total)