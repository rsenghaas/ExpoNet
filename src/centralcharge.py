import numpy as np
from scipy.integrate import quad

def a0(q):
    return 1

def a1(q):
    return -6*q + 45*q**2 - 560*q**3 + 17325 / 2 * q**4

def a2(q):
    return -18*q + 423 / 2 * q**2 - 2972* q**3 + 389415 / 8 * q**4

def Pi_0(q):
    return 1

def Pi_1(q):
    return 1 / (2 * np.pi * 1j) * (a0(q) * np.log(q) + a1(q))

def Pi_2(q):
    return 1 / (2 * np.pi * 1j)**2 *(a0(q) * np.log(q)**2 + 2 * a1(q) *  np.log(q) + a2(q))


def Z(k, z):
    if k == 0:
        return Pi_2(z) + Pi_1(z) + 1/2
    elif k == - 1:
        return Pi_2(z) - Pi_1(z) + 1/2
    else:
        return Pi_2(z) + (2 * k + 1) + (k**2 + k + 1/2)