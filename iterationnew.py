import numpy as np
import pandas as pd
import os
from HXobj import HeatExchanger
from hydraulicloss import P_drop_cold, P_drop_hot, iteration
from GA3_ENTU_Method import heat_transfer_coefficient, effective_NTU

#config list: N_tubes, shape (0/45 - square, 60 - triangle), pitch, N_shell, passes
iterations = [[10, '45', 12, 2, 2, 5, 0.66, "Rect-C-12-45"],
[16, '0', 12, 2, 4, 4, 0.69, "Squ-NC-12"],
[12, '45', 12, 2, 2, 5, 0.66, "Dia-NC-12"],
[12, '0', 12, 2, 4, 4, 0.781, "Squ-NC-14"],
[10, '45', 14, 2, 2, 5, 0.75, "Dia-NC-14"],
[16, '60', 12, 2, 2, 5, 0.781, "Tri-C-12"],
[10, '60', 14, 2, 2, 3, 0.5, "Tri-C-12"],
[18, '60', 12, 2, 2, 5, 0.875, "Tri-NC-12"],
[16, '60', 12, 2, 4, 5, 0.875, "Tri-NC-12"],
[10, '60', 12, 2, 2, 4, 0.812, "Tri-NC-15"]]

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
hydraulic_results = []
ENTU_results = [] #final

def massfunction(t_length, Hx):
    Hx.length = t_length + 0.11

    windowarea = (Hx.D_shell**2/8)*(2*np.arccos(1-2*(1-Hx.baffle_height)) - np.sin(2*np.arccos(1-2*(1-Hx.baffle_height))))
    baffle_area = (Hx.baffle_count) * ((np.pi/4 * (Hx.D_shell)**2) - Hx.tube_count *(np.pi/4 * (Hx.tube_OD)**2) - windowarea)
    resin_vol = (2 * ((np.pi/4 * (Hx.D_shell)**2) - Hx.tube_count *(np.pi/4 * (Hx.tube_OD)**2)) * 0.009) + (2 * (np.pi/4 * (Hx.D_shell)**2) * 0.007)
    
    mass_resin = rho_resin * resin_vol
    mass_baffle = baffle_area * mass_ABS_pua
    mass_shell = mass_shell_pul * Hx.length
    mass_tube = N_tubes * mass_tube_pul * t_length
    mass_rings = 2 * o_rings011 + 2 * o_rings036
    mass_splitter = t_length * 0.15

    total_mass = mass_resin + mass_baffle + mass_nozzles + mass_shell + mass_tube + mass_rings + mass_splitter

    return [total_mass - 1.1, total_mass]



for config in iterations:
    N_tubes, shape, pitch_mm, max_shells, max_passes, rows, bundle_height, config_name = config
    pitch_m = pitch_mm/1000 #in meter

    baffleheights = np.arange(0.5, 0.85, 0.025)

    for baffles in range(4, 15):
        for baffle_height in baffleheights:
            for N_shell in [1, 2]:
                for passes in [2, 4]:
                    if N_shell > max_shells or passes > max_passes:
                        continue
                    
                
                    try:

                        Hx = HeatExchanger(length = 0.35, pitch = pitch_m, tube_count = N_tubes,  baffle_count = baffles, type = shape, passes = passes, N_shell = N_shell, rows = rows, bundle_height = bundle_height, baffle_height = baffle_height)
                        
                        t_length = iteration(massfunction, Hx, initialmass=0.25)
                        tube_length = min(3.5 / N_tubes, t_length, 0.25)

                        Hx_length = 0.11 + tube_length  # Adjust length with calculated tube_length
                        Hx.length = Hx_length  # Update Hx with the new length

                        total_mass = massfunction(tube_length, Hx)[1]
                    
                        if total_mass <= mass_limit:
                            valid_designs.append({
                                'config_name' : config_name,
                                'tubes': N_tubes,
                                'shape': shape,
                                'pitch': pitch_m,
                                'baffles': baffles,
                                'passes': passes,
                                'shells': N_shell,
                                'tube_length': round(tube_length, 5),
                                'Hx_length': round(Hx_length, 5),
                                'total_mass': round(total_mass, 5),
                                'baffle_height' : baffle_height,
                                'bundle_height' : bundle_height,
                                'rows': rows
                            })
                    
                    except Exception as e:
                        print(f"Error with config {config}: {e}")

valid_designs.sort(key=lambda x: x['tube_length'], reverse=True)

# for design in valid_designs:
#     print(design)

for design in valid_designs:
    try:
        Hx = HeatExchanger(
            length = design['Hx_length'],
            pitch = design['pitch'],
            tube_count = design['tubes'],
            baffle_count = design['baffles'],
            type = design['shape'],
            passes = design['passes'],
            N_shell = design['shells'],
            baffle_height = design['baffle_height'],
            bundle_height = design['bundle_height'],
            rows = design['rows']
        )

        mhot = iteration(P_drop_hot, Hx)
        mcold = iteration(P_drop_cold, Hx)
        Q_hot = mhot * Hx.heat_cap * (Hx.temp_hot - Hx.temp_cold)
        Q_cold = mcold * Hx.heat_cap * (Hx.temp_hot - Hx.temp_cold)
        Q_min = min(Q_hot, Q_cold)

        hydraulic_results.append({
            **design,
            'mhot': round(mhot, 4),
            'mcold': round(mcold, 4),
            'Q_hot': round(Q_hot, 2),
            'Q_cold': round(Q_cold, 2),
            'Q_min': round(Q_min, 2),
        })

    except Exception as e:
        print(f'Failed for design {design}: {e}')
    
hydraulic_results.sort(key=lambda x: x['Q_min'], reverse=True)

# for result in hydraulic_results:
#     print(result)

for design in hydraulic_results:
    try:
        Hx = HeatExchanger(
            length = design['Hx_length'],
            pitch = design['pitch'],
            tube_count = design['tubes'],
            baffle_count = design['baffles'],
            type = design['shape'],
            passes = design['passes'],
            N_shell = design['shells'],
            baffle_height = design['baffle_height'],
            bundle_height = design['bundle_height'],
            rows = design['rows']
        )

        m_hot = design['mhot']
        m_cold = design['mcold']

        H = heat_transfer_coefficient(m_hot, m_cold, Hx)
        NTU, eff, Thot_out, Tcold_out, q_abs = effective_NTU(Hx, H, m_hot, m_cold, Hx.temp_hot, Hx.temp_cold)

        ENTU_results.append({
            **design,
            'NTU': NTU,
            'Effectiveness': eff,
            'Thot_out': Thot_out,
            'Tcold_out': Tcold_out,
            'Q_abs': q_abs
        })

    except Exception as e:
        print(f"ENTU calculation failed for design {design}: {e}")

ENTU_results.sort(key=lambda x: x['Q_abs'], reverse=True)

for result in ENTU_results:
    print(result)

df = pd.DataFrame(ENTU_results)
df.to_excel('GA3_HeatExchanger_Optimisation (with length optimiser (capped length)).xlsx', index = False)
print("Results Exported! :D")
print(os.getcwd())