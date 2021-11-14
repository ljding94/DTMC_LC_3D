import numpy as np
from scipy.special import ellipeinc
from scipy import integrate
from scipy.optimize import minimize, NonlinearConstraint, Bounds

# small tools

## unduloid related:
def unduloid_rs(a, e, theta, phi):
    rsqrtcomm = np.sqrt((1 - e * np.cos(theta)) / (1 + e * np.cos(theta)))
    x = a * np.sqrt(1 - e * e) * rsqrtcomm * np.cos(phi)
    y = a * np.sqrt(1 - e * e) * rsqrtcomm * np.sin(phi)
    z = a * e * np.sin(theta) * rsqrtcomm - a * ellipeinc(theta - np.pi / 2, e * e)
    r = np.array([x, y, z])
    return r


def unduloid_z(a, e, theta):
    rsqrtcomm = np.sqrt((1 - e * np.cos(theta)) / (1 + e * np.cos(theta)))
    z = a * e * np.sin(theta) * rsqrtcomm - a * ellipeinc(theta - np.pi / 2, e * e)
    return z


def unduloid_r(a, e, theta):
    r = a * np.sqrt(1 - e * e) * np.sqrt((1 - e * np.cos(theta)) / (1 + e * np.cos(theta)))
    return r


def unduloid_sqrtg(a, e, theta, phi=0):
    sqrtg = a * a * (1 - e * e)
    sqrtg *= np.sqrt((1 - e * np.cos(theta)) / np.power((1 + e * np.cos(theta)), 3))
    return sqrtg


def unduloid_int_theta_sqrtg(theta, a, e):
    return unduloid_sqrtg(a, e, theta, 0)


def unduloid_A(a, e, theta1):
    # A with bump in the middle at theta=pi
    # area of unduloid from theta\in(pi-theta1,pi+theta1)
    A = integrate.quad(unduloid_int_theta_sqrtg, a=-theta1, b=theta1, args=(a, e))[0]
    A *= 2 * np.pi
    return A


def unduloid_l(a, e, theta1):
    l = np.abs(unduloid_z(a, e, -theta1) - unduloid_z(a, e, theta1))
    return l


def unduloid_E(a, e, theta1, lam, kar, C0, mu, A):
    E = 4 * np.pi * lam * unduloid_r(a, e, theta1)
    E += 0.5 * kar * unduloid_A(a, e, theta1) * np.power(1 / a - C0, 2)
    E += mu * np.power(unduloid_A(a, e, theta1) - A, 2)
    return E


def obj_unduloid_E(xs, lam, kar, C0, mu, A):
    a, e, theta1 = xs
    return unduloid_E(a, e, theta1, lam, kar, C0, mu, A)


def nlc_unduloid(A, lf):
    con1 = lambda x: unduloid_A(x[0], x[1], x[2])
    nlc1 = NonlinearConstraint(con1, 0.5*A, np.inf)
    con2 = lambda x: unduloid_l(x[0], x[1], x[2])
    nlc2 = NonlinearConstraint(con2, lf, np.inf)
    #con3 = lambda x: x[1]*x[1]
    #nlc3 = NonlinearConstraint(con3, 0.0, 0.8)
    return (nlc1, nlc2)#, nlc3)


def opt_unduloid_E(lam, kar, C0, mu, A, lf):
    paras = (lam, kar, C0, mu, A)
    xs0 = (np.sqrt(A)/7, 0.0, 3*np.pi)  # a,e,theta1
    print("unduloid_A(x[0], x[1], x[2])",unduloid_A(xs0[0], xs0[1], xs0[2]))
    print("unduloid_l(x[0], x[1], x[2])",unduloid_l(xs0[0], xs0[1], xs0[2]))
    boundsx=[(0,np.sqrt(A)),(-0.9,0.9),(0,10*np.pi)]
    # boundsx=[(),(-1,1),()]
    nlcs = nlc_unduloid(A, lf)
    res = minimize(obj_unduloid_E, xs0, args=paras, method="trust-constr", constraints=nlcs,bounds=boundsx)
    return res
    #method="trust-constr"

