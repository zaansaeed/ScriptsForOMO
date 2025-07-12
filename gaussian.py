import numpy as np
import matplotlib.pyplot as plt

def dynamic_exponent(x):
    base = 2.7112
    max_exp = 10
    ramp = 1 / (1 + np.exp(-(x - 30)))  # smooth ramp after x=30
    return base + (max_exp - base) * ramp

def smooth_sharpened_function(x):
    A = 0.9776682
    B = -540981.1
    C = 4103.594
    exp = dynamic_exponent(x)
    return B + (A - B) / (1 + (x / C)**exp)

# Plot
x = np.linspace(0, 60, 1000)
y = smooth_sharpened_function(x)

plt.plot(x, y, label="Sharpened Function")
plt.axvline(30, color='gray', linestyle='--', label='x = 30')
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Smoothed Steep Drop After x = 30")
plt.grid(True)
plt.legend()
plt.show()