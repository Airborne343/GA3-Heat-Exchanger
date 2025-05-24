import numpy as np
from HXobj import HeatExchanger

Hx = HeatExchanger(length = 0.35, tube_count = 10, baffle_count = 9, type = "60", passes = 2, N_shell = 1, pitch = 12/1000, baffle_height= 0.7, bundle_height=0.7, rows=4)

mass_limit = 1.20 #kg

#Design Variables for Heat Exchanger
# Length of Heat Exchanger -> max 0.35m
#we need to relax the constraint on length as a design variable instead of keeping it at 0.35m constant

#non-negotiable masses
mass_shell_pul = 0.650 #kg/m (pul = per unit length)
mass_tube_pul = 0.20 #kg/m
mass_ABS_pua = 2.39 #kg/m^2
rho_resin = 1150 #kg/m^3 for tubesheets and endplates
o_rings011 = 0.8/1000 #kg for each ring
o_rings036 = 5.3/1000 #kg for each ring

# area of ABS sheet per baffle = pi/4 * (64 * 10^-3)^2 - (N * pi/4 * (8*10^-3)^2)
windowarea = (Hx.D_shell**2/8)*(2*np.arccos(1-2*(1-Hx.baffle_height)) - np.sin(2*np.arccos(1-2*(1-Hx.baffle_height))))
baffle_area = (Hx.baffle_count) * ((np.pi/4 * (Hx.D_shell)**2) - Hx.tube_count *(np.pi/4 * (Hx.tube_OD)**2) - windowarea)
#assuming the 2 parts of 2 baffles can be combined into 1 (no you cant and i need to change this)
tube_length = Hx.length - (4 * 0.025) - (2 * 0.005) #0.025 due to constraint on nozzle and 0.005 as each end plate adds 0.005 to overall length

#resin volume (idk what thickness is yet)
tubesheet_thickness = 0.009
endplate_thickness = 0.007
resin_vol =  (2 * ((np.pi/4 * (Hx.D_shell)**2) - Hx.tube_count *(np.pi/4 * (Hx.tube_OD)**2)) * tubesheet_thickness) + (2 * (np.pi/4 * (Hx.D_shell)**2) * endplate_thickness)

mass_resin = rho_resin * resin_vol
mass_baffle = baffle_area * mass_ABS_pua
mass_nozzles = 0.025 * 4
mass_shell = mass_shell_pul * Hx.length
mass_tube = Hx.tube_count * mass_tube_pul * tube_length
mass_rings = (2 * o_rings011) + (2 * o_rings036) #no clue which one to use where (ask demonstrator)

total_mass = mass_resin + mass_baffle + mass_nozzles + mass_shell + mass_tube + mass_rings
#note:leave a margin of error so the max we can go should be like 1.17, not 1.20 on the dot

print(total_mass)