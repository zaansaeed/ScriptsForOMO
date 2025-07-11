import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define your decay function
def decay_func(x, a, b):
    return 1 / (1 + a * x**b)

# Your target points
x_data = np.array([10, 20, 30])
y_data = np.array([0.9, 0.7, 0.1])

# Fit the curve
popt, _ = curve_fit(decay_func, x_data, y_data, p0=[1e-8, 7])
a_fit, b_fit = popt

# Show the result
print(f"Fitted: a = {a_fit:.3e}, b = {b_fit:.3f}")

# Plot
x = np.linspace(0, 100, 500)
y = decay_func(x, a_fit, b_fit)

plt.plot(x, y, label=fr'$f(x) = \frac{{1}}{{1 + {a_fit:.2e} \cdot x^{{{b_fit:.2f}}}}}$', linewidth=2)
plt.scatter(x_data, y_data, color='red', zorder=10, label="Target Points")
plt.axvline(10, linestyle='--', color='gray', linewidth=1)
plt.axvline(20, linestyle='--', color='gray', linewidth=1)
plt.axvline(30, linestyle='--', color='gray', linewidth=1)
plt.legend()
plt.grid(True)
plt.title("Fit Decay Function to 3 Specific Points")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.show()
