import numpy as np
from scipy import optimize, integrate
import matplotlib.pyplot as plt

# cholesteric on cylinder of radius 1 and length l


# beta = m/2*(phi-tan(alpha)*z/R)
# uc = cos(beta)*n+sin(beta)*(cos(alpha)*z+sin(alpha)*phi)
# un = cos(beta)*n-sin(beta)phi
# interpolate between these two modes
# uitp = cos(beta)*n+sin(beta)*((1-gamma)cos(alpha)*\hat{z}+((1-gamma)sin(alpha)+gamma)*\hat{phi})
# u  = uitp/(sqrt(uitp.uitp))


def n(m, alpha, gamma, phi, z):
    nx = np.cos(phi)
    ny = np.sin(phi)
    nz = np.zeros(len(phi))
    return np.array([nx, ny, nz])

# previous variable under tan(alpha)z , no z/R definition
'''
def u(m, alpha, gamma, phi, z):
    ux = (np.cos(phi) * np.cos((m * (phi - z * np.tan(alpha))) / 2.0) - (-1 + gamma + gamma * np.sin(alpha)) * np.sin(phi) * np.sin((m * (phi - z * np.tan(alpha))) / 2.0)) / np.sqrt(1 + (-1 + gamma) * gamma - (-1 + gamma) * gamma * (-np.sin(alpha) + np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha))))
    uy = (np.cos((m * (phi - z * np.tan(alpha))) / 2.0) * np.sin(phi) + np.cos(phi) * (-1 + gamma + gamma * np.sin(alpha)) * np.sin((m * (phi - z * np.tan(alpha))) / 2.0)) / np.sqrt(1 + (-1 + gamma) * gamma - (-1 + gamma) * gamma * (-np.sin(alpha) + np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha))))
    uz = (gamma * np.cos(alpha) * np.sin((m * (phi - z * np.tan(alpha))) / 2.0)) / np.sqrt(1 + (-1 + gamma) * gamma - (-1 + gamma) * gamma * (-np.sin(alpha) + np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha))))
    return np.array([ux, uy, uz])


def S_R1(m, alpha, gamma, phi, z):
    S = -0.5 * (np.cos((m * (phi - z * np.tan(alpha))) / 2.0) * (-2 + m - m * gamma - 2 * (-1 + gamma) * gamma + (2 - 2 * gamma) * gamma * np.sin(alpha) + 2 * (-1 + gamma) * gamma * np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha)))) / (1 + (-1 + gamma) * gamma - (-1 + gamma) * gamma * (-np.sin(alpha) + np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha)))) ** 1.5
    return S


def T_R1(m, alpha, gamma, phi, z):
    T = -0.25 * (2 * m * gamma / np.cos(alpha) + 2 * gamma * np.cos(alpha) * (-1 + gamma - np.cos(m * (phi - z * np.tan(alpha))) * (-1 + gamma + gamma * np.sin(alpha))) + (2 * m * (-1 + gamma) + gamma ** 2 + gamma ** 2 * np.cos(2 * alpha)) * np.tan(alpha)) / (-1 + gamma - gamma ** 2 + (-1 + gamma) * gamma * (-np.sin(alpha) + np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha))))
    return T


def B_R1(m, alpha, gamma, phi, z):
    Bx = (
        np.cos(phi)
        * (4 - 2 * m * (-1 + gamma) * (-1 + (2 - 2 * gamma) * gamma) + 2 * gamma * (-6 + 5 * gamma) * (1 + (-1 + gamma) * gamma) + gamma * (-2 * (6 - 2 * m) * np.sin(alpha) + 2 * gamma * (14 + 2 * m * (-2 + gamma) + gamma * (-15 + 7 * gamma)) * np.sin(alpha) - 4 * (-1 + gamma) * np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha)) * (-1 + gamma + gamma * np.sin(alpha)) ** 2 - 2 * gamma * np.cos(2 * alpha) * (3 + gamma * (-5 + 3 * gamma) + (-1 + gamma) * gamma * np.sin(alpha))))
        * np.sin((m * (phi - z * np.tan(alpha))) / 2.0) ** 2
        + (-1 + gamma + gamma * np.sin(alpha)) * (2 + m * (-1 + gamma) + 2 * (-1 + gamma) * gamma + gamma * (-2 + 2 * gamma) * np.sin(alpha) - 2 * (-1 + gamma) * gamma * np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha))) * np.sin(phi) * np.sin(m * (phi - z * np.tan(alpha)))
    ) / (4.0 * (1 + (-1 + gamma) * gamma - (-1 + gamma) * gamma * (-np.sin(alpha) + np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha)))) ** 2)

    By = (
        (4 - 2 * m * (-1 + gamma) * (-1 + (2 - 2 * gamma) * gamma) + 2 * gamma * (-6 + 5 * gamma) * (1 + (-1 + gamma) * gamma) + gamma * (-2 * (6 - 2 * m) * np.sin(alpha) + 2 * gamma * (14 + 2 * m * (-2 + gamma) + gamma * (-15 + 7 * gamma)) * np.sin(alpha) - 4 * (-1 + gamma) * np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha)) * (-1 + gamma + gamma * np.sin(alpha)) ** 2 - 2 * gamma * np.cos(2 * alpha) * (3 + gamma * (-5 + 3 * gamma) + (-1 + gamma) * gamma * np.sin(alpha))))
        * np.sin(phi)
        * np.sin((m * (phi - z * np.tan(alpha))) / 2.0) ** 2
        + np.cos(phi) * (-1 + gamma + gamma * np.sin(alpha)) * (-2 + m - m * gamma - 2 * (-1 + gamma) * gamma + (2 - 2 * gamma) * gamma * np.sin(alpha) + 2 * (-1 + gamma) * gamma * np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha))) * np.sin(m * (phi - z * np.tan(alpha)))
    ) / (4.0 * (1 + (-1 + gamma) * gamma - (-1 + gamma) * gamma * (-np.sin(alpha) + np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha)))) ** 2)

    Bz = (m * (1 - gamma) * gamma * np.cos(alpha) * np.sin(m * (phi - z * np.tan(alpha)))) / (4.0 * (1 + (-1 + gamma) * gamma - (-1 + gamma) * gamma * (-np.sin(alpha) + np.cos(m * (phi - z * np.tan(alpha))) * (1 + np.sin(alpha)))) ** 2)

    B = np.array([Bx, By, Bz])
    return B


def S(m, alpha, gamma, phi, z, R):
    S = -0.5*(np.cos((m*(phi - z*np.tan(alpha)))/2.)*(-2 + m - m*gamma - 2*(-1 + gamma)*gamma + (2 + m*(-1 + R) - 2*gamma)*gamma*np.sin(alpha) + 2*(-1 + gamma)*gamma*np.cos(m*(phi - z*np.tan(alpha)))*(1 + np.sin(alpha))))/ (R*(1 + (-1 + gamma)*gamma - (-1 + gamma)*gamma*(-np.sin(alpha) + np.cos(m*(phi - z*np.tan(alpha)))*(1 + np.sin(alpha))))**1.5)
    return S

def T(m, alpha, gamma, phi, z, R):
    T = -0.25*(2*m*R*gamma/np.cos(alpha) + 2*gamma*np.cos(alpha)*(-1 + m - m*R + gamma - np.cos(m*(phi - z*np.tan(alpha)))*(-1 + gamma + gamma*np.sin(alpha))) +(2*m*R*(-1 + gamma) + gamma**2 + gamma**2*np.cos(2*alpha))*np.tan(alpha))/(R*(-1 + gamma - gamma**2 + (-1 + gamma)*gamma*(-np.sin(alpha) + np.cos(m*(phi - z*np.tan(alpha)))*(1 + np.sin(alpha)))))
    return T

def B(m, alpha, gamma, phi, z, R):
    Bx = (np.cos(phi)*(4 + 2*gamma*(-6 + 5*gamma)*(1 + (-1 + gamma)*gamma) - 2*m*(-1 + gamma)*(-1 + gamma*(2 + (-3 + R)*gamma)) +gamma*(-2*(6 + m*(-3 + R))*np.sin(alpha) + 2*gamma*(14 + gamma*(-15 + 7*gamma) + 2*m*(-3 + R + 2*gamma - R*gamma))*np.sin(alpha) - 4*(-1 + gamma)*np.cos(m*(phi - z*np.tan(alpha)))*(1 + np.sin(alpha))*(-1 + gamma + gamma*np.sin(alpha))**2 - 2*gamma*np.cos(2*alpha)*(3 + gamma*(-5 + 3*gamma) + m*(-1 + R + gamma - R*gamma) + (-1 + gamma)*gamma*np.sin(alpha))))*np.sin((m*(phi - z*np.tan(alpha)))/2.)**2 + (-1 + gamma + gamma*np.sin(alpha))*(2 + m*(-1 + gamma) + 2*(-1 + gamma)*gamma + gamma*(-2 + m - m*R + 2*gamma)*np.sin(alpha) - 2*(-1 + gamma)*gamma*np.cos(m*(phi - z*np.tan(alpha)))*(1 + np.sin(alpha)))*np.sin(phi)*np.sin(m*(phi - z*np.tan(alpha))))/(4.*R*(1 + (-1 + gamma)*gamma - (-1 + gamma)*gamma*(-np.sin(alpha) + np.cos(m*(phi - z*np.tan(alpha)))*(1 + np.sin(alpha))))**2)

    By = ((4 + 2*gamma*(-6 + 5*gamma)*(1 + (-1 + gamma)*gamma) - 2*m*(-1 + gamma)*(-1 + gamma*(2 + (-3 + R)*gamma)) + gamma*(-2*(6 + m*(-3 + R))*np.sin(alpha) + 2*gamma*(14 + gamma*(-15 + 7*gamma) + 2*m*(-3 + R + 2*gamma - R*gamma))*np.sin(alpha) - 4*(-1 + gamma)*np.cos(m*(phi - z*np.tan(alpha)))*(1 + np.sin(alpha))*(-1 + gamma + gamma*np.sin(alpha))**2 - 2*gamma*np.cos(2*alpha)*(3 + gamma*(-5 + 3*gamma) + m*(-1 + R + gamma - R*gamma) + (-1 + gamma)*gamma*np.sin(alpha))))*np.sin(phi)*np.sin((m*(phi - z*np.tan(alpha)))/2.)**2 + np.cos(phi)*(-1 + gamma + gamma*np.sin(alpha))*(-2 + m - m*gamma - 2*(-1 + gamma)*gamma + (2 + m*(-1 + R) - 2*gamma)*gamma*np.sin(alpha) + 2*(-1 + gamma)*gamma*np.cos(m*(phi - z*np.tan(alpha)))*(1 + np.sin(alpha)))*np.sin(m*(phi - z*np.tan(alpha))))/(4.*R*(1 + (-1 + gamma)*gamma - (-1 + gamma)*gamma*(-np.sin(alpha) + np.cos(m*(phi - z*np.tan(alpha)))*(1 + np.sin(alpha))))**2)

    Bz = (m*gamma*np.cos(alpha)*(1 - gamma + (-1 + R)*gamma*np.sin(alpha))*np.sin(m*(phi - z*np.tan(alpha))))/(4.*R*(1 + (-1 + gamma)*gamma - (-1 + gamma)*gamma*(-np.sin(alpha) + np.cos(m*(phi - z*np.tan(alpha)))*(1 + np.sin(alpha))))**2)

    B = np.array([Bx, By, Bz])
    return B

'''

def beta_cal(m,alpha,R,phi,z):
    return m/2*(phi-np.tan(alpha)*z/R)

def u_cal(m, alpha, gamma, R, phi, z):
    beta = beta_cal(m,alpha,R,phi,z)

    ux = (np.cos(phi)*np.cos(beta) - (-1 + gamma + gamma*np.sin(alpha))*np.sin(phi)*np.sin(beta))/np.sqrt(np.cos(beta)**2 + (1 + 2*(-1 + gamma)*gamma + 2*(-1 + gamma)*gamma*np.sin(alpha))*np.sin(beta)**2)

    uy =  (np.cos(beta)*np.sin(phi) + np.cos(phi)*(-1 + gamma + gamma*np.sin(alpha))*np.sin(beta))/np.sqrt(np.cos(beta)**2 + (1 + 2*(-1 + gamma)*gamma + 2*(-1 + gamma)*gamma*np.sin(alpha))*np.sin(beta)**2),

    uz =  (gamma*np.cos(alpha)*np.sin(beta))/np.sqrt(np.cos(beta)**2 + (1 + 2*(-1 + gamma)*gamma + 2*(-1 + gamma)*gamma*np.sin(alpha))*np.sin(beta)**2)

    return np.array([uz,uy,uz])

def S_cal(m, alpha, gamma, R, phi, z):
    beta = beta_cal(m,alpha,R,phi,z)

    S = (np.cos(beta)*(2 + m*(-1 + gamma) + 2*(-1 + gamma)*gamma - 2*(-1 + gamma)*gamma*(np.cos(2*beta) - 2*np.sin(alpha)*np.sin(beta)**2)))/(2.*R*(np.cos(beta)**2 + (1 + 2*(-1 + gamma)*gamma + 2*(-1 + gamma)*gamma*np.sin(alpha))*np.sin(beta)**2)**1.5)

    return S

def T_cal(m, alpha, gamma, R, phi, z):
    beta = beta_cal(m,alpha,R,phi,z)

    T = -(((m*gamma*np.cos(alpha))/2. + (-1 + gamma + gamma*np.sin(alpha))*(gamma*np.cos(alpha)*np.sin(beta)**2 + (m*np.tan(alpha))/2.))/(R*(-1 + gamma - gamma**2 + (-1 + gamma)*gamma*(np.cos(2*beta) - 2*np.sin(alpha)*np.sin(beta)**2))))

    return T

def B_cal(m, alpha, gamma, R, phi, z):
    beta = beta_cal(m,alpha,R,phi,z)

    Bx = (np.sin(beta)*((m*(-1 + gamma)*(np.cos(beta)*(-1 + gamma + gamma*np.sin(alpha))*np.sin(phi) + np.cos(phi)*np.sin(beta) + 2*(-1 + gamma)*gamma*np.cos(phi)*(1 + np.sin(alpha))*np.sin(beta)))/2. + (-1 + gamma + gamma*np.sin(alpha))*(np.cos(beta)*np.sin(phi) + np.cos(phi)*(-1 + gamma + gamma*np.sin(alpha))*np.sin(beta))*(1 + (-1 + gamma)*gamma - (-1 + gamma)*gamma*(np.cos(2*beta) - 2*np.sin(alpha)*np.sin(beta)**2))))/(R*(np.cos(beta)**2 + (1 + 2*(-1 + gamma)*gamma + 2*(-1 + gamma)*gamma*np.sin(alpha))*np.sin(beta)**2)**2)

    By = (np.sin(beta)*((m*(-1 + gamma)*(-(np.cos(phi)*np.cos(beta)*(-1 + gamma + gamma*np.sin(alpha))) + np.sin(phi)*np.sin(beta) + 2*(-1 + gamma)*gamma*(1 + np.sin(alpha))*np.sin(phi)*np.sin(beta)))/2. + (-1 + gamma + gamma*np.sin(alpha))*(-(np.cos(phi)*np.cos(beta)) + (-1 + gamma + gamma*np.sin(alpha))*np.sin(phi)*np.sin(beta))*(1 + (-1 + gamma)*gamma - (-1 + gamma)*gamma*(np.cos(2*beta) - 2*np.sin(alpha)*np.sin(beta)**2))))/(R*(np.cos(beta)**2 + (1 + 2*(-1 + gamma)*gamma + 2*(-1 + gamma)*gamma*np.sin(alpha))*np.sin(beta)**2)**2)

    Bz = -0.25*(m*(-1 + gamma)*gamma*np.cos(alpha)*np.sin(2*beta))/(R*(np.cos(beta)**2 + (1 + 2*(-1 + gamma)*gamma + 2*(-1 + gamma)*gamma*np.sin(alpha))*np.sin(beta)**2)**2)

    return np.array([Bx, By, Bz])


# integral

def intSS_unit_length(m, alpha, gamma, bn_phi, bn_z, R):
    phi, z = np.meshgrid(2 * np.pi * np.linspace(0.5 / bn_phi, 1 - 0.5 / bn_phi, bn_phi), np.linspace(0.5 / bn_z, 1 - 0.5 / bn_z, bn_z))
    phi, z = phi.flatten(), z.flatten()
    dphi = 2 * np.pi / bn_phi
    dz = 1 / bn_z
    #s = S(m, alpha, gamma, phi, z, R)
    s = S_cal(m, alpha, gamma, R, phi, z)
    res = s * s * R * dphi * dz
    return res.sum()


def intTT_unit_length(m, alpha, gamma, bn_phi, bn_z, R):
    phi, z = np.meshgrid(2 * np.pi * np.linspace(0.5 / bn_phi, 1 - 0.5 / bn_phi, bn_phi), np.linspace(0.5 / bn_z, 1 - 0.5 / bn_z, bn_z))
    phi, z = phi.flatten(), z.flatten()
    dphi = 2 * np.pi / bn_phi
    dz = 1 / bn_z
    #t = T(m, alpha, gamma, phi, z, R)
    t = T_cal(m, alpha, gamma, R, phi, z)
    res = t * t * R* dphi * dz
    return res.sum()


def intT_unit_length(m, alpha, gamma, bn_phi, bn_z, R):
    phi, z = np.meshgrid(2 * np.pi * np.linspace(0.5 / bn_phi, 1 - 0.5 / bn_phi, bn_phi), np.linspace(0.5 / bn_z, 1 - 0.5 / bn_z, bn_z))
    phi, z = phi.flatten(), z.flatten()
    dphi = 2 * np.pi / bn_phi
    dz = 1 / bn_z
    #t = T(m, alpha, gamma, phi, z, R)
    t = T_cal(m, alpha, gamma, R, phi, z)
    res = t * R*dphi * dz
    return res.sum()


def intBB_unit_length(m, alpha, gamma, bn_phi, bn_z, R):
    phi, z = np.meshgrid(2 * np.pi * np.linspace(0.5 / bn_phi, 1 - 0.5 / bn_phi, bn_phi), np.linspace(0.5 / bn_z, 1 - 0.5 / bn_z, bn_z))
    phi, z = phi.flatten(), z.flatten()
    dphi = 2 * np.pi / bn_phi
    dz = 1 / bn_z
    #b = B(m, alpha, gamma, phi, z, R)
    b = B_cal(m, alpha, gamma, R, phi, z)
    bb = b[0] * b[0] + b[1] * b[1] + b[2] * b[2]
    res = bb * R*dphi * dz
    return res.sum()


def intC_unit_length(m, alpha, gamma, bn_phi, bn_z, R):
    phi, z = np.meshgrid(2 * np.pi * np.linspace(0.5 / bn_phi, 1 - 0.5 / bn_phi, bn_phi), np.linspace(0.5 / bn_z, 1 - 0.5 / bn_z, bn_z))
    phi, z = phi.flatten(), z.flatten()
    dphi = 2 * np.pi / bn_phi
    dz = 1 / bn_z
    ui = u_cal(m, alpha, gamma, R, phi, z)
    ni = n(m, alpha, gamma, phi, z)
    un = ui[0] * ni[0] + ui[1] * ni[1] + ui[2] * ni[2]
    res = (1 - un * un) * R* dphi * dz
    return res.sum()


def Ftot_unit_length(K, C, q, m, alpha, gamma, bn_phi, bn_z, R):
    intSS = intSS_unit_length(m, alpha, gamma, bn_phi, bn_z, R)
    #intTT = intTT_unit_length(m, alpha, gamma, bn_phi, bn_z, R)
    intT = intT_unit_length(m, alpha, gamma, bn_phi, bn_z, R)
    intBB = intBB_unit_length(m, alpha, gamma, bn_phi, bn_z, R)
    intC = intC_unit_length(m, alpha, gamma, bn_phi, bn_z, R)
    Ftot = 0.5 * K * (intSS + (intT - q)**2  + intBB) + 0.5 * C * intC
    return Ftot


def obj_Ftot(xs, K, C, q, m, bn_phi, bn_z, R=1):
    alpha, gamma = xs
    return Ftot_unit_length(K, C, q, m, alpha, gamma, bn_phi, bn_z, R)


def opt_alpha_gamma_fun(fun, K, C, q, m, bn_phi, bn_z, R, method):
    paras = (K, C, q, m, bn_phi, bn_z, R)
    alpha_gamma_0 = (0.25*np.pi, 0.5)
    if method:
        opt = optimize.minimize(fun, alpha_gamma_0, args=paras, bounds=((0, 0.5 * np.pi), (0, 1)), method=method)
    else:
        opt = optimize.minimize(fun, alpha_gamma_0, args=paras)
    return opt
