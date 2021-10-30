import numpy as np
from scipy.special import ellipeinc
from scipy import integrate
# small tools

## unduloid related:
def unduloid_r(a,e,theta,phi):
    rsqrtcomm=np.sqrt((1-e*np.cos(theta))/(1+e*np.cos(theta)))
    x=a*np.sqrt(1-e*e)*rsqrtcomm*np.cos(phi)
    y=a*np.sqrt(1-e*e)*rsqrtcomm*np.sin(phi)
    z=a*e*np.sin(theta)*rsqrtcomm-a*ellipeinc(theta-np.pi/2,e*e)
    r = np.array([x,y,z])
    return r

def unduloid_sqrtg(a,e,theta,phi=0):
    sqrtg = a*a*(1-e*e)
    sqrtg *= np.sqrt((1-e*np.cos(theta))/np.power((1+e*np.cos(theta)),3))
    return sqrtg
def unduloid_int_theta_sqrtg(theta,a,e):
    return unduloid_sqrtg(a,e,theta,0)

def unduloid_A_unit_a(e,theta1):
    #A/a^2
    # area of unduloid from theta\in(-theta1,theta1)
    # as unit of a^2
    #print("e,theta1=",e,theta1)
    A = integrate.quad(unduloid_int_theta_sqrtg,a=0,b=theta1,args=(1,e))[0]
    A *= 4*np.pi
    return A

def unduloid_lf_unit_a(e,theta1):
    # z1/a
    zcomm=e*np.sin(theta1)*np.sqrt((1-e*np.cos(theta1))/(1+e*np.cos(theta1)))
    z0 = -zcomm
    z0 -= ellipeinc(-theta1-np.pi/2,e*e)
    z1 = zcomm
    z1 -= ellipeinc(theta1-np.pi/2,e*e)
    lf=z0-z1
    return lf

#FIXME: z is not symmetric about theta=0, thus need to find lf another way

def unduloid_R_unit_a(e,theta1):
    #R/a
    R = np.sqrt(1-e*e)*np.sqrt((1-e*np.cos(theta))/(1+e*np.cos(theta)))
    return R

def unduloid_E(a,A,z1,lam,kar,C0):
    # firstly find e and theta1 from A/a^2 and z1/a
    # TODO: get e theta1 from (A,lf) first

    E = 4*np.pi*lam*a*unduloid_R_unit_a(e, theta1)

