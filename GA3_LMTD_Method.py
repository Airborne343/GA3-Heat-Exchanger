import numpy as np
import math
from HXobj import HeatExchanger

def heat_transfer_coefficient(mhot, mcold, Hx):
    
    #inner heat convection coefficient
    m_tube = mhot/Hx.tube_count
    V_tube = m_tube/(Hx.density * Hx.area_tube)
    Re_tube = (Hx.density * V_tube * Hx.tube_ID)/Hx.dynamic_viscosity
    Nu_inner = 0.023 * (Re_tube ** 0.8) * (Hx.Prandtl_no ** 0.3)
    conv_coeff_inner = (Nu_inner * Hx.water_heat_conductivity)/Hx.tube_ID

    #outer heat convection coefficient
    V_shell = mcold/(Hx.density * Hx.area_shell)
    Re_shell = (Hx.density * V_shell * Hx.charc_D_shell)/Hx.dynamic_viscosity
    c = 0.15 #constant
    Nu_outer = c * (Re_shell ** 0.6) * (Hx.Prandtl_no ** 0.3)
    conv_coeff_outer = (Nu_outer * Hx.water_heat_conductivity)/Hx.tube_OD

    #heat conduction through walls
    conv_walls = (Hx.inner_surface_area * np.log(Hx.tube_OD/Hx.tube_ID))/(2*np.pi*Hx.tube_heat_conductivity *Hx.length)

    H = 1/(conv_walls + (1/conv_coeff_inner) + (1/conv_coeff_outer)*(Hx.inner_surface_area/Hx.outer_surface_area))
    return H

#initialise variables
Hx = HeatExchanger(tube_count = 13, baffle_count = 9, type = "triangle")
m_hot = 0.47
m_cold = 0.5
A_ht = Hx.tube_count * np.pi * Hx.length * Hx.tube_ID
F = 1 
H = heat_transfer_coefficient(m_hot, m_cold, Hx)

#LMTD Method
n_iter = 100
convergence_thresh = 1e-4

C_hot = m_hot * Hx.heat_cap
C_cold = m_cold * Hx.heat_cap

Thot_out_init = 59.99

for i in range(n_iter):
    Q_val = C_hot * (Hx.temp_hot - Thot_out_init)
    Tcold_out = Hx.temp_cold + Q_val/C_cold

    Delta_T1 = Hx.temp_hot - Tcold_out
    Delta_T2 = Thot_out_init - Hx.temp_cold

    if Delta_T1 <=0 or Delta_T2 <=0:
        raise ValueError("Unrealistic Guess for Delta T")
    
    LMTD = (Delta_T1 - Delta_T2)/np.log(Delta_T1/Delta_T2)
    Q_LMTD = H * A_ht * F * LMTD
    Thot_out_LMTD = Hx.temp_hot - Q_LMTD/C_hot

    if abs(Thot_out_LMTD - Thot_out_init) < convergence_thresh:
        break

    Thot_out_init = 0.1 * Thot_out_LMTD + 0.9 * Thot_out_init

Q_fn = (Hx.temp_hot - Thot_out_LMTD) * C_hot
Tcold_out_fn = Hx.temp_cold + Q_fn/C_cold

print(f"Convergence in {i+1} iterations:")
print(f"Hot Side Outlet Temperature: {Thot_out_LMTD:.2f} C")
print(f"Cold Side Outlet Temperature: {Tcold_out_fn:.2f} C")
print(f"Heat Transfer Coefficient {H}")
print(f"LMTD {LMTD:.2f} C")
print(f"Heat Transfer {Q_fn:.2f} W")