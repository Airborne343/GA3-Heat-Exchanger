import numpy as np
from HXobj import HeatExchanger

#config list: N_tubes, shape (0/45 - square, 60 - triangle), pitch, N_shell, passes
iterations = [[10, 'square', 12, 2, 2],
[16, 'square', 12, 2, 4],
[12, 'square', 12, 2, 2],
[12, 'square', 12, 2, 4],
[10, 'square', 14, 2, 2],
[16, 'triangle', 12, 2, 2],
[10, 'triangle', 14, 2, 2],
[18, 'triangle', 12, 2, 2],
[16, 'triangle', 12, 2, 4],
[10, 'triangle', 12, 2, 2]]

#constants
mass_shell_pul = 0.650
mass_tube_pul = 0.20
mass_ABS_pua = 2.39
rho_resin = 1150
o_rings011 = 0.8 / 1000
o_rings036 = 5.3 / 1000
mass_nozzles = 0.025 * 4
thickness = 0.005
mass_limit = 1.17

valid_designs = []

for config in iterations:
    N_tubes, shape, pitch_mm, max_shells, max_passes = config
    pitch_m = pitch_mm/1000 #in meters
    tube_length = 3.5/N_tubes #for maximum heat transfer
    tube_length = min(3.5 / N_tubes, 0.35 - 0.11)
    Hx_length = tube_length + 0.11
    
    for baffles in range(6, 11):
        for N_shell in [1, 2]:
            for passes in [2, 4]:
                if N_shell > max_shells or passes > max_passes:
                    continue

                try:
                    Hx = HeatExchanger(length = Hx_length, pitch = pitch_m, tube_count = N_tubes,  baffle_count = baffles, type = shape, passes = passes, N_shell = N_shell)
                    baffle_area = baffles * ((np.pi/4 * (Hx.D_shell)**2) - N_tubes *(np.pi/4 * (Hx.tube_OD)**2)) * 0.7
                    resin_vol = 2 * ((np.pi / 4 * Hx.D_shell**2) -
                                    N_tubes * (np.pi / 4 * Hx.tube_OD**2)) * thickness \
                                + 2 * (np.pi / 4 * Hx.D_shell**2) * thickness
                    
                    mass_resin = rho_resin * resin_vol
                    mass_baffle = baffle_area * mass_ABS_pua
                    mass_shell = mass_shell_pul * Hx.length
                    mass_tube = N_tubes * mass_tube_pul * tube_length
                    mass_rings = 2 * o_rings011 + 2 * o_rings036

                    
                    total_mass = mass_resin + mass_baffle + mass_nozzles + mass_shell + mass_tube + mass_rings

                    if total_mass <= mass_limit:
                        valid_designs.append({
                            'tubes': N_tubes,
                            'shape': shape,
                            'pitch': pitch_mm,
                            'baffles': baffles,
                            'passes': passes,
                            'shells': N_shell,
                            'tube_length': round(tube_length, 3),
                            'hx_length': round(Hx_length, 3),
                            'total_mass': round(total_mass, 4)
                        })
                
                except Exception as e:
                    print(f"Error with config {config}: {e}")

valid_designs.sort(key=lambda x: x['tube_length'], reverse=True)

for design in valid_designs:
    print(design)