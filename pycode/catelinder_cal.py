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

def catelinder_E(lam,kar,l,h,b):
    E = lam*4*np.pi*b*np.cosh(h/(2*b))+kar*np.pi*(l-h)/b
    return E

# elements for optimization
def obj_catelinder_E(xs,lam,kar):
    b,h,l=xs
    # optimize all b,h, and l
    return catelinder_E(lam,kar,l,h,b)

def nlc_catelinder(A):
    con1 = lambda x: catelinder_A(x[0], x[1], x[2])
    nlc1 = NonlinearConstraint(con1, 0.5*A,  2*A)
    con2 = lambda x: (x[1]-x[2])
    nlc2 = NonlinearConstraint(con2,-np.inf,0)
    return (nlc1,nlc2)


def opt_catelinder_E(lam,kar,lf,A):
    paras=(lam,kar)
    bhl0=(A/(2*np.pi*lf),lf,lf)
    boundsbh=[(0.1,A/(2*np.pi*lf)),(0,lf),(lf,2*lf)]
    nlcbh = nlc_catelinder(A)
    res = minimize(obj_catelinder_E, bhl0,args=paras,constraints=nlcbh,bounds=boundsbh)
    return res



