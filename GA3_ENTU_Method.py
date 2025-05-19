import numpy as np
import math

#initialise variables
cp_w = 4179
rho_w = 990.1
k_w = 0.632
k_tube = 386
L = 0.35
d_inner = 0.006
Pr = 4.31
mu_w = 6.51*(10^-4)

Thot_in = 60
Tcold_in = 20
m_hot = 0.47
m_cold = 0.5
H = 3927
N_tubes = 13
A_ht = N_tubes * np.pi * L * d_inner
F = 1 

#ENTU Method
def effective_NTU(HX, H, m_hot, m_cold, Thot_in = 60 ,Tcold_in = 20):
    C_hot = HX.heat_cap * m_hot
    C_cold = HX.heat_cap * m_cold
    C_min = min(C_hot, C_cold)
    C_max = max(C_hot, C_cold)
    C_r = C_min/C_max

    
    A_ht = HX.tube_count * np.pi * HX.length * HX.tube_ID
    NTU = (H*A_ht)/C_min

    eff_1 = 2/(1 + C_r + np.sqrt(1 + C_r**2)*((1 + np.exp(-NTU * np.sqrt(1 + C_r**2)))/(1 - np.exp(-NTU * np.sqrt(1 + C_r**2)))))

    if HX.N_shell == 1:
        eff = eff_1
    
    else:
        eff = (((1-eff_1*C_r)/(1-eff_1))*HX.N_shell - 1)/(((1-eff_1*C_r)/(1-eff_1))**HX.N_shell - C_r)

    q_max = C_min * (Thot_in - Tcold_in)
    q_abs = eff * q_max 

    Thot_out = Thot_in - q_abs/C_hot
    Tcold_out = Tcold_in + q_abs/C_cold

    return [NTU, eff, Thot_out, Tcold_out, q_abs]

result = effective_NTU(H, N_tubes, m_hot, Thot_in, m_cold, Tcold_in, cp_w, 1)

print(f"ENTU NTU: {result[0]:.3f}")
print(f"ENTU Effectiveness: {result[1]:.3f}")
print(f"ENTU Hot Side Outlet Temperature: {result[2]:.2f} C")
print(f"ENTU Cold Side Outlet Temperature: {result[3]:.2f} C")
print(f"ENTU Heat Transfer: {result[4]:.2f} W")