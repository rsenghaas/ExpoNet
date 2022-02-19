import numpy as np


def fwd_euler(f, y, dt, theta):
    return dt*f(y, theta)


def rk4(f, y, dt, theta, expo):
    k1 = f(y, theta, expo)
    if k1[0] == 0:
        return np.array([0,0,0])
    k2 = f(y + dt * k1 / 2, theta, expo)
    if k2[0] == 0:
        return np.array([0,0,0])
    k3 = f(y + dt * k2 / 2, theta, expo)
    if k3[0] == 0:
        return np.array([0,0,0])
    k4 = f(y + dt * k3, theta, expo)
    if k4[0] == 0:
        return np.array([0,0,0])
    else:
        return 1 / 6 * dt * (k1 + 2 * k2 + 2 * k3 + k4)
