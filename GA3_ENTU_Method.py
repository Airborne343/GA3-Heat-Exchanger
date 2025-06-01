import numpy as np
import math
from HXobj import HeatExchanger
from Datatables import find_eff

def heat_transfer_coefficient(mhot, mcold, Hx):
    
    #inner heat convection coefficient
    m_tube = mhot/Hx.tube_count
    V_tube = m_tube/(Hx.density * Hx.area_tube)
    Re_tube = (Hx.density * V_tube * Hx.tube_ID)/Hx.dynamic_viscosity
    Nu_inner = 0.023 * (Re_tube ** 0.8) * (Hx.Prandtl_no ** 0.3)
    conv_coeff_inner = (Nu_inner * Hx.water_heat_conductivity)/Hx.tube_ID

    #outer heat convection coefficient
    windowarea = (Hx.D_shell**2/8)*(2*np.arccos(1-2*(1-Hx.baffle_height)) - np.sin(2*np.arccos(1-2*(1-Hx.baffle_height))))
    V_window = mcold/(windowarea * Hx.density)
    V_shell = mcold/(Hx.density * Hx.area_shell)
    Re_shell = (Hx.density * (0.9*V_shell + 0.1*V_window) * Hx.tube_OD)/Hx.dynamic_viscosity
    Nu_outer = Hx.c * (Re_shell ** 0.6) * (Hx.Prandtl_no ** 0.3)
    #conv_coeff_outer_factor
    #J_i = np.exp(Hx.A + Hx.B*np.log(Re_shell) + Hx.C*(np.log(Re_shell)**2) + Hx.D*(Re_shell**4) + Hx.E*(np.log(Re_shell)**5)) #tube_arrangement_correction
    J_c = 0.55 + 0.72*Hx.crossflow_prop
    conv_coeff_outer = ((Nu_outer * Hx.water_heat_conductivity)/Hx.tube_OD) * J_c

    #heat conduction through walls
    conv_walls = (Hx.inner_surface_area * np.log(Hx.tube_OD/Hx.tube_ID))/(2*np.pi*Hx.tube_heat_conductivity *Hx.length)

    H = 1/(conv_walls + (1/conv_coeff_inner) + (1/conv_coeff_outer)*(Hx.inner_surface_area/Hx.outer_surface_area))
    return H

#initialise variables
# Hx = HeatExchanger(length = 0.35, pitch = 0.012, tube_count = 13, baffle_count = 9, type = "60", passes = 1, N_shell = 1)
# m_hot = 0.47
# m_cold = 0.5
# A_ht = Hx.tube_count * np.pi * Hx.length * Hx.tube_ID
# H = heat_transfer_coefficient(m_hot, m_cold, Hx)

#ENTU Method
def effective_NTU(HX, H, m_hot, m_cold, Thot_in = 60, Tcold_in = 20):
    C_hot = HX.heat_cap * m_hot
    C_cold = HX.heat_cap * m_cold
    C_min = min(C_hot, C_cold)
    C_max = max(C_hot, C_cold)
    C_r = C_min/C_max

    
    A_ht = HX.tube_count * np.pi * HX.length * HX.tube_ID
    NTU = (H*A_ht)/C_min

    eff = find_eff(NTU, C_r, HX.N_shell)

    q_max = C_min * (Thot_in - Tcold_in)
    q_abs = eff * q_max 

    Thot_out = Thot_in - q_abs/C_hot
    Tcold_out = Tcold_in + q_abs/C_cold

    return [NTU, eff, Thot_out, Tcold_out, q_abs]

# result = effective_NTU(Hx, H, m_hot, m_cold, Hx.temp_hot, Hx.temp_cold)

# print(f"ENTU NTU: {result[0]:.3f}")
# print(f"ENTU Effectiveness: {result[1]:.3f}")
# print(f"ENTU Hot Side Outlet Temperature: {result[2]:.2f} C")
# print(f"ENTU Cold Side Outlet Temperature: {result[3]:.2f} C")
# print(f"ENTU Heat Transfer: {result[4]:.2f} W")