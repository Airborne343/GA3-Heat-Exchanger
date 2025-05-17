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

#LMTD Method
n_iter = 100
convergence_thresh = 1e-4

C_hot = m_hot * cp_w
C_cold = m_cold * cp_w

Thot_out_init = 59.99

for i in range(n_iter):
    Q_val = C_hot * (Thot_in - Thot_out_init)
    Tcold_out = Tcold_in + Q_val/C_cold

    Delta_T1 = Thot_in - Tcold_out
    Delta_T2 = Thot_out_init - Tcold_in

    if Delta_T1 <=0 or Delta_T2 <=0:
        raise ValueError("Unrealistic Guess for Delta T")
    
    LMTD = (Delta_T1 - Delta_T2)/np.log(Delta_T1/Delta_T2)
    Q_LMTD = H * A_ht * F * LMTD
    Thot_out_LMTD = Thot_in - Q_LMTD/C_hot

    if abs(Thot_out_LMTD - Thot_out_init) < convergence_thresh:
        break

    Thot_out_init = 0.1 * Thot_out_LMTD + 0.9 * Thot_out_init

Q_fn = (Thot_in - Thot_out_LMTD) * C_hot
Tcold_out_fn = Tcold_in + Q_fn/C_cold

print(f"Convergence in {i+1} iterations:")
print(f"Hot Side Outlet Temperature: {Thot_out_LMTD:.2f} C")
print(f"Cold Side Outlet Temperature: {Tcold_out_fn:.2f} C")
print(f"LMTD {LMTD:.2f} C")
print(f"Heat Transfer {Q_fn:.2f} W")