import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from HXobj import HeatExchanger
from hydraulicloss import P_drop_cold, P_drop_hot, iteration
from GA3_ENTU_Method import heat_transfer_coefficient, effective_NTU

final_result = []

#constants
mass_shell_pul = 0.650
mass_tube_pul = 0.20
mass_ABS_pua = 2.39
rho_resin = 1150
o_rings011 = 0.8 / 1000
o_rings036 = 5.3 / 1000
mass_nozzles = 0.025 * 4

def total_mass(Hx):
    windowarea = (Hx.D_shell**2/8)*(2*np.arccos(1-2*(1-Hx.baffle_height)) - np.sin(2*np.arccos(1-2*(1-Hx.baffle_height))))
    baffle_area = (Hx.baffle_count) * ((np.pi/4 * (Hx.D_shell)**2) - Hx.tube_count *(np.pi/4 * (Hx.tube_OD)**2) - windowarea)
    resin_vol = (2 * ((np.pi/4 * (Hx.D_shell)**2) - Hx.tube_count *(np.pi/4 * (Hx.tube_OD)**2)) * 0.009) + (2 * (np.pi/4 * (Hx.D_shell)**2) * 0.007)

    mass_resin = rho_resin * resin_vol
    mass_baffle = baffle_area * mass_ABS_pua
    mass_shell = mass_shell_pul * Hx.length
    mass_tube = Hx.tube_count * mass_tube_pul * Hx.tube_length
    mass_rings = 2 * o_rings011 + 2 * o_rings036
    mass_splitter = Hx.tube_length * 0.15

    totalmass = mass_resin + mass_baffle + mass_nozzles + mass_shell + mass_tube + mass_rings + mass_splitter
    return totalmass


result_repository = {
    "Group A - 2025": {
        "year": 2025,
        "group_name": "Group A",
        "HX_Length": 0.325,
        "HX_tube_length": 0.250,
        "HX_pitch": 14.5/1000,
        "HX_tube_count": 14,
        "HX_baffle_count": 5,
        "HX_shape": "60",
        "HX_passes": 2,
        "HX_shells": 2,
        "HX_baffle_height": 0.793,
        "HX_bundle_height": 0.744,
        "rows": 7
    },

    "Group B - 2025": {
        "year": 2025,
        "group_name": "Group B",
        "HX_Length": 0.325,
        "HX_tube_length": 0.260,
        "HX_pitch": 12/1000,
        "HX_tube_count": 12,
        "HX_baffle_count": 6,
        "HX_shape": "60",
        "HX_passes": 4,
        "HX_shells": 2,
        "HX_baffle_height": 0.657,
        "HX_bundle_height": 0.659,
        "rows": 6
    },

    "Group C - 2025": {
        "year": 2025,
        "group_name": "Group C",
        "HX_Length": 0.267,
        "HX_tube_length": 0.195,
        "HX_pitch": 12/1000,
        "HX_tube_count": 16,
        "HX_baffle_count": 7,
        "HX_shape": "60",
        "HX_passes": 2,
        "HX_shells": 2,
        "HX_baffle_height": 0.748,
        "HX_bundle_height": 0.775,
        "rows": 5
    },

    "Group D - 2025": {
        "year": 2025,
        "group_name": "Group D",
        "HX_Length": 0.343,
        "HX_tube_length": 0.266,
        "HX_pitch": 15/1000,
        "HX_tube_count": 12,
        "HX_baffle_count": 8,
        "HX_shape": "60",
        "HX_passes": 2,
        "HX_shells": 2,
        "HX_baffle_height": 0.700,
        "HX_bundle_height": 0.805,
        "rows": 4
    },

    "Group E - 2025": {
        "year": 2025,
        "group_name": "Group E",
        "HX_Length": 0.339,
        "HX_tube_length": 0.262,
        "HX_pitch": 13/1000,
        "HX_tube_count": 12,
        "HX_baffle_count": 6,
        "HX_shape": "60",
        "HX_passes": 2,
        "HX_shells": 2,
        "HX_baffle_height": 0.750,
        "HX_bundle_height": 0.672,
        "rows": 5
    }
}



for result in result_repository.values():
    try:
        Hx = HeatExchanger(
            length = result['HX_Length'],
            tube_length = result['HX_tube_length'],
            pitch = result['HX_pitch'],
            tube_count = result['HX_tube_count'],
            baffle_count = result['HX_baffle_count'],
            type = result['HX_shape'],
            passes = result['HX_passes'],
            N_shell = result['HX_shells'],
            baffle_height = result['HX_baffle_height'],
            bundle_height = result['HX_bundle_height'],
            rows = result['rows']
        )

        totalmass = total_mass(Hx)
        mhot = iteration(P_drop_hot, Hx)
        mcold = iteration(P_drop_cold, Hx)
        [pdrop_hot_diff, p_hotside_total, p_drop_hotchic] = P_drop_hot(mhot, Hx)
        [pdrop_cold_diff, p_coldside_total, p_drop_xflow, p_drop_window] = P_drop_cold(mcold, Hx)
        
        Q_hot = mhot * Hx.heat_cap * (Hx.temp_hot - Hx.temp_cold)
        Q_cold = mcold * Hx.heat_cap * (Hx.temp_hot - Hx.temp_cold)
        Q_min = min(Q_hot, Q_cold)
        H = heat_transfer_coefficient(mhot, mcold, Hx)
        NTU, eff, Thot_out, Tcold_out, q_abs = effective_NTU(Hx, H, mhot, mcold, Hx.temp_hot, Hx.temp_cold)

        filtered_result = {k: v for k, v in result.items() if k not in ('hotpoly', 'coldpoly')}

        filtered_result.update({
            'mass': round(totalmass, 4),
            'mhot': round(mhot, 4),
            'mcold': round(mcold, 4),
            'Pdrop_hot': round(p_hotside_total, 4),
            'Pdrop_cold': round(p_coldside_total, 4),
            'Q_hot': round(Q_hot, 2),
            'Q_cold': round(Q_cold, 2),
            'Q_min': round(Q_min, 2),
            'NTU': NTU,
            'Effectiveness': eff,
            'Thot_out': Thot_out,
            'Tcold_out': Tcold_out,
            'Q_abs': q_abs
        })

        final_result.append(filtered_result)
        print("Done!")

    except Exception as e:
        print(f"Calculations failed for design {result}: {e}")

# print(final_result)
df = pd.DataFrame(final_result)
df.to_excel('2025_heat_exchanger_results.xlsx', index=False)