import numpy as np
from scipy.optimize import minimize, NonlinearConstraint, Bounds


def catelinder_A(b, h, l):
    # b: neck radius
    # h: catenoid total height
    # l: catelinder seperation can be larger than pulling length lf
    A = 2 * (l - h)
    A += h + b * np.sinh(h / b)
    A *= np.pi * b
    return A


def catelinder_E(lam, kar, mu, A, l, h, b):
    E = lam * 4 * np.pi * b * np.cosh(h / (2 * b)) + kar * np.pi * (l - h) / b + mu * np.power(catelinder_A(b, h, l) - A, 2)
    return E


# elements for optimization
def obj_catelinder_E(xs, lam, kar, mu, A):
    b, h, l = xs
    # optimize all b,h, and l
    return catelinder_E(lam, kar, mu, A, l, h, b)


def nlc_catelinder(A, lf):
    con1 = lambda x: catelinder_A(x[0], x[1], x[2])
    # nlc1 = NonlinearConstraint(con1, 0.5 * A, 2 * A)
    nlc1 = NonlinearConstraint(con1, 0.5*A, 1.5*A)
    con2 = lambda x: (x[1] - x[2])
    nlc2 = NonlinearConstraint(con2, -np.inf, 0)
    con3 = lambda x: x[2]
    nlc3 = NonlinearConstraint(con3, lf, np.inf)
    return (nlc1, nlc2, nlc3)
    # (nlc1, nlc2)


def opt_catelinder_E(lam, kar, mu, lf, A):
    paras = (lam, kar, mu, A)
    bhl0 = (A / (2 * np.pi * lf), 0, lf)
    # boundsbh = [(0.1, 1.2*A / (2 * np.pi * lf)), (0, lf)]
    nlcbh = nlc_catelinder(A, lf)
    res = minimize(obj_catelinder_E, bhl0, args=paras, method="trust-constr", constraints=nlcbh)  # , bounds=boundsbh)
    return res
