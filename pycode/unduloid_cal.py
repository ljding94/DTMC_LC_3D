import numpy as np
from scipy.special import ellipeinc
from scipy import integrate

# small tools

## unduloid related:
def unduloid_r(a, e, theta, phi):
    rsqrtcomm = np.sqrt((1 - e * np.cos(theta)) / (1 + e * np.cos(theta)))
    x = a * np.sqrt(1 - e * e) * rsqrtcomm * np.cos(phi)
    y = a * np.sqrt(1 - e * e) * rsqrtcomm * np.sin(phi)
    z = a * e * np.sin(theta) * rsqrtcomm - a * ellipeinc(theta - np.pi / 2, e * e)
    r = np.array([x, y, z])
    return r


def unduloid_sqrtg(a, e, theta, phi=0):
    sqrtg = a * a * (1 - e * e)
    sqrtg *= np.sqrt((1 - e * np.cos(theta)) / np.power((1 + e * np.cos(theta)), 3))
    return sqrtg


def unduloid_int_theta_sqrtg(theta, a, e):
    return unduloid_sqrtg(a, e, theta, 0)


def unduloid_A(a, e, theta1):
    # A
    # area of unduloid from theta\in(-theta1,theta1)
    A = integrate.quad(unduloid_int_theta_sqrtg, a=-theta1, b=theta1, args=(a, e))[0]
    A *= 2 * np.pi
    return A


def unduloid_l(a, e, theta1):
    # z(-theta1)-z(theta1)
    zcomm = e * np.sin(theta1) * np.sqrt((1 - e * np.cos(theta1)) / (1 + e * np.cos(theta1)))
    z0 = -zcomm
    z0 -= ellipeinc(-theta1 - np.pi / 2, e * e)
    z1 = zcomm
    z1 -= ellipeinc(theta1 - np.pi / 2, e * e)
    l = a * (z0 - z1)
    return l


def unduloid_R(a, e, theta1):
    # R/a
    R = a * np.sqrt(1 - e * e) * np.sqrt((1 - e * np.cos(theta)) / (1 + e * np.cos(theta)))
    return R


def unduloid_E(a, e, theta1, lam, kar, C0):
    E = 4 * np.pi * lam * unduloid_R_unit_a(a, e, theta1)
    E += 0.5 * kar * unduloid_A(a, e, theta1) * np.power(1 / a - C0, 2)
    return E

def obj_unduloid_E(xs,lam,kar,C0):
    a,e,theta1 = xs
    return unduloid_E(a, e, theta1, lam, kar, C0)

def nlc_unduloid(A,lf):
    pass
    #TODO: finish implementation for unduloid


