import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from HXobj import HeatExchanger
from hydraulicloss import iteration, Friction, Kc, Ke
from GA3_ENTU_Method import heat_transfer_coefficient, effective_NTU


final_result = []

def P_drop_hot(mhot, Hx, hotpoly): #Function to calculate hot side pressure drop

    #Working out velocity, dynamic pressure and reynolds number inside tube
    Volflow = mhot/Hx.density
    V_pipe = Volflow/(Hx.area_tubes/Hx.passes)
    q_pipe = 1/2 * Hx.density * V_pipe ** 2
    Re_pipe = Hx.density * V_pipe * Hx.tube_ID / Hx.dynamic_viscosity

    #Working out three sources of pressure drop, tube friction, end losses for tubes and inlet nozzle drops
    P_pipe = q_pipe * Hx.length / Hx.tube_ID * Friction(Re_pipe)
    P_ends = q_pipe * (Kc(Hx.sigma, Re_pipe) + Ke(Hx.sigma, Re_pipe)) * Hx.passes
    P_nozzle = 0.5 * Hx.density * (Volflow/Hx.area_nozzle)**2 * 2

    #Summing pressure losses and combining
    P_total = P_pipe + P_ends + P_nozzle
    return [P_total - hotpoly(Volflow), P_total, hotpoly(Volflow)]

def P_drop_cold(mcold, Hx, coldpoly):
    #Calculating velocity, characteristic diameter and shell reynolds number
    Volflow = mcold/Hx.density
    V_shell = Volflow/Hx.area_shell

    if Hx.type == "60":
        D_e_const1 = 1.10
        D_e_const2 = 0.917
    elif Hx.type == "0" or Hx.type == "45":
        D_e_const1 = 1.27
        D_e_const2 = 0.785
    else:
        raise Exception("invalid type given, should be '60', '0', or '45'")

    charc_D_shell = D_e_const1 * (Hx.pitch ** 2 - D_e_const2 * Hx.tube_OD **2) / Hx.tube_OD
    Re_shell = (Hx.density * V_shell * charc_D_shell)/Hx.dynamic_viscosity

    #Working out pressure losses due to shell and nozzle exit
    ploss_shell = 4 * Hx.a * ((Re_shell)**(-0.15)) * Hx.tube_count * Hx.density * ((V_shell)**2)
    V_nozzle = mcold/(Hx.density * Hx.area_nozzle)
    ploss_nozzle = 2 * 0.5 * Hx.density * (V_nozzle**2)


    ploss_cold_tot = ploss_nozzle + ploss_shell

    return [ploss_cold_tot - coldpoly(Volflow), ploss_cold_tot, coldpoly(Volflow)]

def iteration(pressurefunction, Hx, initialmass = 0.45, tol = 0.005, maxiter = 15):
    massflow = initialmass

    for i in range(maxiter):
        #Working out 
        p = pressurefunction(massflow, Hx)[0]
        dp_dq = (pressurefunction(massflow * 1.05, Hx)[0] - pressurefunction(massflow * 0.95, Hx)[0]) / (2 * 0.05 * massflow)

        #Newton raphson method to calculate pressurefunction drop
        massflow_new = massflow - p / dp_dq
        
        if abs(massflow - massflow_new) < massflow * tol:
            return massflow

        massflow = massflow_new

    raise("Maximum iterations reached without convergence")

##2024
# COLD SIDE data
cold_flow_2024 = np.array([0.6333, 0.6083, 0.5750, 0.5083, 0.4250, 0.3583, 0.3083, 0.2417, 0.1917, 0.1583]) / 1000
cold_pressure_2024 = np.array([0.1024, 0.1444, 0.1870, 0.2717, 0.3568, 0.4203, 0.4626, 0.5152, 0.5597, 0.5776]) * 100000

# HOT SIDE data
hot_flow_2024 = np.array([0.4826, 0.4340, 0.3924, 0.3507, 0.3021, 0.2535, 0.1979, 0.1493, 0.1111, 0.0694]) / 1000
hot_pressure_2024 = np.array([0.0944, 0.1662, 0.2297, 0.2820, 0.3294, 0.3856, 0.4447, 0.5006, 0.5311, 0.5615]) * 100000

# Fit 4th-order polynomials
cold_poly_coeffs_2024 = np.polyfit(cold_flow_2024, cold_pressure_2024, 3)
hot_poly_coeffs_2024 = np.polyfit(hot_flow_2024, hot_pressure_2024, 3)

# Create polynomial functions
cold_poly_2024 = np.poly1d(cold_poly_coeffs_2024)
hot_poly_2024 = np.poly1d(hot_poly_coeffs_2024)

##2023
# COLD SIDE data
cold_flow_2023 = np.array([0.7083, 0.6417, 0.5750, 0.5083, 0.4250, 0.3583, 0.3083, 0.2417, 0.1917, 0.1583]) / 1000
cold_pressure_2023 = np.array([0.1310, 0.2017, 0.2750, 0.3417, 0.4038, 0.4503, 0.4856, 0.5352, 0.5717, 0.5876]) * 100000

# HOT SIDE data
hot_flow_2023 = np.array([0.4722, 0.4340, 0.3924, 0.3507, 0.3021, 0.2535, 0.1979, 0.1493, 0.1111, 0.0694]) / 1000
hot_pressure_2023 = np.array([0.0538, 0.1192, 0.1727, 0.2270, 0.2814, 0.3366, 0.3907, 0.4456, 0.4791, 0.5115]) * 100000

# Fit 4th-order polynomials
cold_poly_coeffs_2023 = np.polyfit(cold_flow_2023, cold_pressure_2023, 3)
hot_poly_coeffs_2023 = np.polyfit(hot_flow_2023, hot_pressure_2023, 3)

# Create polynomial functions
cold_poly_2023 = np.poly1d(cold_poly_coeffs_2023)
hot_poly_2023 = np.poly1d(hot_poly_coeffs_2023)

##2022
# COLD SIDE data
cold_flow_2022 = np.array([0.5833, 0.5083, 0.4750, 0.4250, 0.3792, 0.3417, 0.2958, 0.2583, 0.2125, 0.1708]) / 1000
cold_pressure_2022 = np.array([0.1113, 0.2157, 0.2538, 0.3168, 0.3613, 0.4031, 0.4511, 0.4846, 0.5181, 0.5573]) * 100000

# HOT SIDE data
hot_flow_2022 = np.array([0.4583, 0.4236, 0.4010, 0.3611, 0.3125, 0.2639, 0.2222, 0.1597, 0.1181, 0.0694]) / 1000
hot_pressure_2022 = np.array([0.1333, 0.1756, 0.2024, 0.2577, 0.3171, 0.3633, 0.4233, 0.4784, 0.5330, 0.5715]) * 100000

# Fit 4th-order polynomials
cold_poly_coeffs_2022 = np.polyfit(cold_flow_2022, cold_pressure_2022, 3)
hot_poly_coeffs_2022 = np.polyfit(hot_flow_2022, hot_pressure_2022, 3)

# Create polynomial functions
cold_poly_2022 = np.poly1d(cold_poly_coeffs_2022)
hot_poly_2022 = np.poly1d(hot_poly_coeffs_2022)

##2019
# COLD SIDE data
cold_flow_2019 = np.array([0.6917, 0.6750, 0.6292, 0.5917, 0.5458, 0.5083, 0.4625, 0.4250, 0.3792, 0.3417, 0.2958, 0.2542, 0.2125, 0.1708]) / 1000
cold_pressure_2019 = np.array([0.1475, 0.1619, 0.2178, 0.2607, 0.3041, 0.3417, 0.3756, 0.4118, 0.4423, 0.4711, 0.5031, 0.5297, 0.5561, 0.5823]) * 100000

# HOT SIDE data
hot_flow_2019 = np.array([0.5382, 0.5278, 0.4931, 0.4549, 0.4201, 0.3854, 0.3507, 0.3160, 0.2813, 0.2465, 0.2118, 0.1771, 0.1424, 0.1076, 0.0694]) / 1000
hot_pressure_2019 = np.array([0.1101, 0.1315, 0.1800, 0.2185, 0.2537, 0.2999, 0.3440, 0.3780, 0.4149, 0.4547, 0.5005, 0.5271, 0.5677, 0.5971, 0.6045]) * 100000

# Fit 4th-order polynomials
cold_poly_coeffs_2019 = np.polyfit(cold_flow_2019, cold_pressure_2019, 3)
hot_poly_coeffs_2019 = np.polyfit(hot_flow_2019, hot_pressure_2019, 3)

# Create polynomial functions
cold_poly_2019 = np.poly1d(cold_poly_coeffs_2019)
hot_poly_2019 = np.poly1d(hot_poly_coeffs_2019)

##2018
# COLD SIDE data
cold_flow_2018 = np.array([0.4426, 0.4255, 0.4055, 0.3913, 0.3799, 0.3628, 0.3485, 0.3286, 0.3058, 0.2801, 0.2573, 0.2317, 0.2060, 0.1861, 0.1576, 0.1319, 0.1034, 0.0806, 0.0664, 0.0521]) / 1000
cold_pressure_2018 = np.array([0.1068, 0.1418, 0.1779, 0.2056, 0.2382, 0.2601, 0.2858, 0.3187, 0.3627, 0.4037, 0.4426, 0.4845, 0.5213, 0.5569, 0.6036, 0.6412, 0.6838, 0.7121, 0.7343, 0.7744]) * 100000

# HOT SIDE data
hot_flow_2018 = np.array([0.4954, 0.4805, 0.4640, 0.4475, 0.4310, 0.4145, 0.3980, 0.3815, 0.3650, 0.3485, 0.3320, 0.3155, 0.2990, 0.2825, 0.2660, 0.2495, 0.2330, 0.2165, 0.2000, 0.1819, 0.1670, 0.1472, 0.1307, 0.1142, 0.1010, 0.0845, 0.0680, 0.0515]) / 1000
hot_pressure_2018 = np.array([0.0989, 0.1245, 0.1541, 0.1827, 0.2083, 0.2339, 0.2625, 0.2880, 0.3115, 0.3330, 0.3575, 0.3800, 0.4014, 0.4249, 0.4503, 0.4647, 0.4900, 0.5134, 0.5337, 0.5470, 0.5703, 0.5966, 0.6068, 0.6150, 0.6242, 0.6304, 0.6375, 0.6457]) * 100000

# Fit 4th-order polynomials
cold_poly_coeffs_2018 = np.polyfit(cold_flow_2018, cold_pressure_2018, 3)
hot_poly_coeffs_2018 = np.polyfit(hot_flow_2018, hot_pressure_2018, 3)

# Create polynomial functions
cold_poly_2018 = np.poly1d(cold_poly_coeffs_2018)
hot_poly_2018 = np.poly1d(hot_poly_coeffs_2018)

##2017
# COLD SIDE data
cold_flow_2017 = np.array([0.4967, 0.4739, 0.4511, 0.4312, 0.4055, 0.3856, 0.3628, 0.3371, 0.3115, 0.2887, 0.2659, 0.2431, 0.2231, 0.2003, 0.1775, 0.1547, 0.1376, 0.1120, 0.0664]) / 1000
cold_pressure_2017 = np.array([0.0674, 0.1309, 0.2043, 0.2654, 0.3279, 0.3829, 0.4391, 0.4903, 0.5415, 0.5824, 0.6203, 0.6541, 0.6848, 0.7105, 0.7361, 0.7577, 0.7721, 0.7916, 0.8183]) * 100000

# HOT SIDE data
hot_flow_2017 = np.array([0.4937, 0.4789, 0.4640, 0.4475, 0.4294, 0.4129, 0.3980, 0.3782, 0.3650, 0.3452, 0.3320, 0.3188, 0.3007, 0.2792, 0.2644, 0.2479, 0.2330, 0.2165, 0.2000, 0.1835, 0.1670, 0.1489, 0.1340, 0.1175, 0.0994, 0.0845, 0.0664]) / 1000
hot_pressure_2017 = np.array([0.0579, 0.0845, 0.1091, 0.1317, 0.1644, 0.1890, 0.2125, 0.2401, 0.2575, 0.2891, 0.3055, 0.3239, 0.3494, 0.3750, 0.3953, 0.4137, 0.4360, 0.4564, 0.4727, 0.4970, 0.5073, 0.5276, 0.5478, 0.5570, 0.5632, 0.5734, 0.5805]) * 100000

# Fit 4th-order polynomials
cold_poly_coeffs_2017 = np.polyfit(cold_flow_2017, cold_pressure_2017, 3)
hot_poly_coeffs_2017 = np.polyfit(hot_flow_2017, hot_pressure_2017, 3)

# Create polynomial functions
cold_poly_2017 = np.poly1d(cold_poly_coeffs_2017)
hot_poly_2017 = np.poly1d(hot_poly_coeffs_2017)

result_repository = {
    "Group A - 2024": {
        "year": 2024,
        "group_name": "Group A",
        "HX_Length": 0.340,
        "HX_pitch": 12/1000,
        "HX_tube_count": 12,
        "HX_baffle_count": 8,
        "HX_shape": "60",
        "HX_passes": 2,
        "HX_shells": 1,
        "HX_baffle_height": 0.785,
        "HX_bundle_height": 0.489,
        "rows": 3,
        "hotpoly": hot_poly_2024,
        "coldpoly": cold_poly_2024
    },

    "Group B - 2024": {
        "year": 2024,
        "group_name": "Group B",
        "HX_Length": 0.329,
        "HX_pitch": 16/1000,
        "HX_tube_count": 12,
        "HX_baffle_count": 8,
        "HX_shape": "60",
        "HX_passes": 2,
        "HX_shells": 1,
        "HX_baffle_height": 0.750,
        "HX_bundle_height": 0.879,
        "rows": 5,
        "hotpoly": hot_poly_2024,
        "coldpoly": cold_poly_2024
    },

    "Group C - 2024": {
        "year": 2024,
        "group_name": "Group C",
        "HX_Length": 0.360,
        "HX_pitch": 12/1000,
        "HX_tube_count": 12,
        "HX_baffle_count": 8,
        "HX_shape": "60",
        "HX_passes": 2,
        "HX_shells": 2,
        "HX_baffle_height": 0.754,
        "HX_bundle_height": 0.597,
        "rows": 6,
        "hotpoly": hot_poly_2024,
        "coldpoly": cold_poly_2024
    },

    "Group D - 2024": {
        "year": 2024,
        "group_name": "Group D",
        "HX_Length": 0.320,
        "HX_pitch": 12/1000,
        "HX_tube_count": 14,
        "HX_baffle_count": 6,
        "HX_shape": "60",
        "HX_passes": 2,
        "HX_shells": 1,
        "HX_baffle_height": 0.750,
        "HX_bundle_height": 0.649,
        "rows": 5,
        "hotpoly": hot_poly_2024,
        "coldpoly": cold_poly_2024
    },

    "Group E - 2024": {
        "year": 2024,
        "group_name": "Group E",
        "HX_Length": 0.320,
        "HX_pitch": 13/1000,
        "HX_tube_count": 15,
        "HX_baffle_count": 7,
        "HX_shape": "60",
        "HX_passes": 3,
        "HX_shells": 1,
        "HX_baffle_height": 0.702,
        "HX_bundle_height": 0.932,
        "rows": 8,
        "hotpoly": hot_poly_2024,
        "coldpoly": cold_poly_2024
    }

    # "Group A - 2023": {
    #     "year": 2023,
    #     "group_name": "Group A",
    #     "HX_Length": 0.334,
    #     "HX_pitch": 15.2/1000,
    #     "HX_tube_count": 14,
    #     "HX_baffle_count": 10,
    #     "HX_shape": "triangle",
    #     "HX_passes": 2,
    #     "HX_shells": 1,
    #     "hotpoly": hot_poly_2023,
    #     "coldpoly": cold_poly_2023
    # },

    # "Group B - 2023": {
    #     "year": 2023,
    #     "group_name": "Group B",
    #     "HX_Length": 0.271,
    #     "HX_pitch": 16/1000,
    #     "HX_tube_count": 12,
    #     "HX_baffle_count": 8,
    #     "HX_shape": "triangle",
    #     "HX_passes": 2,
    #     "HX_shells": 1,
    #     "hotpoly": hot_poly_2024,
    #     "coldpoly": cold_poly_2024
    # },

    # "Group C - 2024": {
    #     "year": 2024,
    #     "group_name": "Group C",
    #     "HX_Length": 0.360,
    #     "HX_pitch": 12/1000,
    #     "HX_tube_count": 12,
    #     "HX_baffle_count": 8,
    #     "HX_shape": "triangle",
    #     "HX_passes": 2,
    #     "HX_shells": 2,
    #     "hotpoly": hot_poly_2024,
    #     "coldpoly": cold_poly_2024
    # },

    # "Group D - 2024": {
    #     "year": 2024,
    #     "group_name": "Group D",
    #     "HX_Length": 0.320,
    #     "HX_pitch": 12/1000,
    #     "HX_tube_count": 14,
    #     "HX_baffle_count": 6,
    #     "HX_shape": "triangle",
    #     "HX_passes": 2,
    #     "HX_shells": 1,
    #     "hotpoly": hot_poly_2024,
    #     "coldpoly": cold_poly_2024
    # },

    # "Group E - 2024": {
    #     "year": 2024,
    #     "group_name": "Group E",
    #     "HX_Length": 0.320,
    #     "HX_pitch": 13/1000,
    #     "HX_tube_count": 15,
    #     "HX_baffle_count": 7,
    #     "HX_shape": "triangle",
    #     "HX_passes": 3,
    #     "HX_shells": 1,
    #     "hotpoly": hot_poly_2024,
    #     "coldpoly": cold_poly_2024
    # },

}

for result in result_repository.values():
    try:
        Hx = HeatExchanger(
            length = result['HX_Length'],
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

        mhot = iteration(lambda m, Hx: P_drop_hot(m, Hx, result['hotpoly']), Hx)
        mcold = iteration(lambda m, Hx: P_drop_cold(m, Hx, result['coldpoly']), Hx)
        Q_hot = mhot * Hx.heat_cap * (Hx.temp_hot - Hx.temp_cold)
        Q_cold = mcold * Hx.heat_cap * (Hx.temp_hot - Hx.temp_cold)
        Q_min = min(Q_hot, Q_cold)
        H = heat_transfer_coefficient(mhot, mcold, Hx)
        NTU, eff, Thot_out, Tcold_out, q_abs = effective_NTU(Hx, H, mhot, mcold, Hx.temp_hot, Hx.temp_cold)

        filtered_result = {k: v for k, v in result.items() if k not in ('hotpoly', 'coldpoly')}

        filtered_result.update({
            'mhot': round(mhot, 4),
            'mcold': round(mcold, 4),
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
df.to_excel('past_year_heat_exchanger_results_edited.xlsx', index=False)