import numpy as np
import matplotlib.pyplot as plt

def lorenz63(t, state, sigma, rho, beta):
    x, y, z = state
    dxdt = sigma * (y - x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta * z
    return np.array([dxdt, dydt, dzdt])

def runge_kutta4(func, y0, t0, t1, dt, *args):
    t = np.arange(t0, t1, dt)
    n = len(t)
    y = np.zeros((n, len(y0)))
    y[0] = y0
    for i in range(n - 1):
        k1 = dt * func(t[i], y[i], *args)
        k2 = dt * func(t[i] + 0.5 * dt, y[i] + 0.5 * k1, *args)
        k3 = dt * func(t[i] + 0.5 * dt, y[i] + 0.5 * k2, *args)
        k4 = dt * func(t[i] + dt, y[i] + k3, *args)
        y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return t, y

# Parameters
sigma = 10.0
rho = 28.0
beta = 8.0 / 3.0

# Initial conditions
y0 = np.array([1.0, 1.0, 1.0])

# Time settings
t0 = 0
t1 = 50
dt = 0.01

# Solve the Lorenz-63 model
t, y = runge_kutta4(lorenz63, y0, t0, t1, dt, sigma, rho, beta)

# Plot time series
plt.figure(figsize=(10, 5))
plt.plot(t, y[:, 0], label='x')
plt.plot(t, y[:, 1], label='y')
plt.plot(t, y[:, 2], label='z')
plt.xlabel('Time')
plt.ylabel('Values')
plt.title('Lorenz-63 Time Series')
plt.legend()
plt.grid(True)
plt.show()

# Plot phase space
plt.figure(figsize=(8, 6))
plt.plot(y[:, 0], y[:, 1], label='xy')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Lorenz-63 Phase Space')
plt.grid(True)
plt.show()

