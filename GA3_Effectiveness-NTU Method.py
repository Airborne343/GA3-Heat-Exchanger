import numpy as np
import math
import matplotlib as plt
from scipy.optimize import differential_evolution

#initialise variables
Thot_in = 60
Tcold_in = 20
cp_w = 4179
rho_w = 990.1
k_w = 0.632
k_tube = 386
mu_w = 6.51*(10^-4)
Pr = 4.31

def effective_NTU(H, A_ht, m_hot, Thot_in, m_cold, Tcold_in, cp_w, N_shell):
    C_hot = cp_w * m_hot
    C_cold = cp_w * m_cold
    C_min = min(C_hot, C_cold)
    C_max = max(C_hot, C_cold)
    C_r = C_min/C_max
    
    NTU = (H*A_ht)/C_min

    eff_1 = 2/(1 + C_r + np.sqrt(1 + C_r**2)*((1 + np.exp(-NTU * np.sqrt(1 + C_r**2)))/(1 - np.exp(-NTU * np.sqrt(1 + C_r**2)))))

    if N_shell == 1:
        eff = eff_1
    
    else:
        eff = (((1-eff_1*C_r)/(1-eff_1))**N_shell - 1)/(((1-eff_1*C_r)/(1-eff_1))**N_shell - C_r)

    q_max = C_min * (Thot_in - Tcold_in)
    q_abs = eff * q_max 

    Thot_out = Thot_in - q_abs/C_hot
    Tcold_out = Tcold_in - q_abs/C_cold

    return {NTU, eff, q_abs, Thot_out, Tcold_out}

def maximum_effectiveness(x):
    H, N_tubes, m_hot, m_cold, N_shell = x
    N_tubes = int(round(N_tubes))
    N_shell = int(round(N_shell))

    #safety
    if H <=0 or N_tubes <=1 or m_cold <=0 or m_hot <=0 or N_shell <=1:
        return 1e6
    
    C_hot = cp_w * m_hot
    C_cold = cp_w * m_cold
    C_min = min(C_hot, C_cold)
    C_max = max(C_hot, C_cold)
    C_r = C_min/C_max

    A_ht = N_tubes*np.pi*0.35*0.006
    
    NTU = (H*A_ht)/C_min

    eff_1 = 2/(1 + C_r + np.sqrt(1 + C_r**2)*((1 + np.exp(-NTU * np.sqrt(1 + C_r**2)))/(1 - np.exp(-NTU * np.sqrt(1 + C_r**2)))))

    if N_shell == 1:
        eff = eff_1
    
    else:
        eff = (((1-eff_1*C_r)/(1-eff_1))**N_shell - 1)/(((1-eff_1*C_r)/(1-eff_1))**N_shell - C_r)

    return -eff

bounds = [
    (100, 5000),  # H (W/m².K)
    (1, 20),      # N_tubes (integer)
    (0.3, 1.5),   # m_cold (kg/s)
    (0.3, 1.5),   # m_hot (kg/s)
    (1, 3)       # N_shell (integer)
]

result = differential_evolution(maximum_effectiveness, bounds, strategy='best1bin', maxiter=1000, tol=1e-7)

opt_H, opt_N_tubes, opt_m_cold, opt_m_hot, opt_N_shell = result.x
opt_N_shell = int(round(opt_N_shell))
max_eff = -result.fun

print(f"Optimal H: {opt_H:.2f} W/m²·K")
print(f"Optimal N_tubes: {opt_N_tubes}")
print(f"Optimal m_cold: {opt_m_cold:.2f} kg/s")
print(f"Optimal m_hot: {opt_m_hot:.2f} kg/s")
print(f"Optimal N_shell: {opt_N_shell}")
print(f"Maximum Effectiveness: {max_eff:.4f}")