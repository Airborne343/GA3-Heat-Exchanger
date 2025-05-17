import numpy as np
import matplotlib.pyplot as plt

# COLD SIDE data
cold_flow = np.array([0.6580, 0.6290, 0.5830, 0.5380, 0.4670, 0.3920, 0.3210, 0.2790, 0.2210, 0]) / 1000
cold_pressure = np.array([0.1584, 0.1958, 0.2493, 0.3127, 0.3723, 0.4436, 0.4950, 0.5318, 0.5739, 0.7077]) * 100000

# HOT SIDE data
hot_flow = np.array([0.4360, 0.3870, 0.3520, 0.3110, 0.2600, 0.2290, 0.1670, 0.1180, 0.0690, 0.0010]) / 1000
hot_pressure = np.array([0.0932, 0.1688, 0.2209, 0.2871, 0.3554, 0.4041, 0.4853, 0.5260, 0.5665, 0.6239]) * 100000

# Fit 4th-order polynomials
cold_poly_coeffs = np.polyfit(cold_pressure, cold_flow, 3)
hot_poly_coeffs = np.polyfit(hot_pressure, hot_flow, 3)

# Create polynomial functions
cold_poly = np.poly1d(cold_poly_coeffs)
hot_poly = np.poly1d(hot_poly_coeffs)

# Print coefficients
print("Cold side 4th-order polynomial coefficients:")
print(cold_poly_coeffs)

print("\nHot side 4th-order polynomial coefficients:")
print(hot_poly_coeffs)

# Optional: plot fits
x_cold = np.linspace(0, max(cold_pressure)*1.2, 200)
x_hot = np.linspace(0, max(hot_pressure)*1.2, 200)

plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.plot(cold_pressure, cold_flow, 'bo', label='Cold data')
plt.plot(x_cold, cold_poly(x_cold), 'b-', label='Cold fit')
plt.title('Cold Side Fit')
plt.xlabel('Flowrate (litres/s)')
plt.ylabel('Pressure rise (bar)')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(hot_pressure, hot_flow, 'ro', label='Hot data')
plt.plot(x_hot, hot_poly(x_hot), 'r-', label='Hot fit')
plt.title('Hot Side Fit')
plt.xlabel('Flowrate (litres/s)')
plt.ylabel('Pressure rise (bar)')
plt.legend()

plt.tight_layout()
plt.show()