import numpy as np
from scipy import optimize, integrate
import matplotlib.pyplot as plt

# cholesteric on cylinder of radius 1 and length l

def inte_SS(m, alpha):
    # ans = 2*m*np.pi+np.sin(m*z*np.tan(alpha))+np.sin(2*m*np.pi-m*z*np.tan(alpha))
    # ans /=2*m
    return np.pi


def inte_Tq(q, m, alpha):
    # ansT = np.sin(2*alpha)*(-2*m*np.pi+np.sin(m*z*np.tan(alpha))+np.sin(2*m*np.pi-m*z*np.tan(alpha)))/(4*m)
    # ansT += m*np.pi/np.cos(alpha)
    ansT = m * np.pi / np.cos(alpha) + np.pi * np.cos(alpha) * np.sin(alpha)
    # ansTT = 4*m*np.pi*(3-3*np.cos(4*alpha)+16*m*m/np.cos(alpha)**2)/(128*m)+m*np.pi*np.sin(alpha)
    ansTT = 3 - 3 * np.cos(4 * alpha) + 16 * m * m / (np.cos(alpha)) ** 2 + 32 * m * np.sin(alpha)
    ansTT *= np.pi / 32
    ans = ansTT - 2 * q * ansT + 2 * np.pi * q * q
    return ans


def inte_BB(m, alpha):
    # ans = np.cos(alpha)**2*(4*m*np.pi-np.sin(2*m*z*np.tan(alpha))-np.sin(4*m*np.pi-2*m*z*np.tan(alpha)))
    # ans -= 8*np.sin(alpha)**2*(-2*m*np.pi+np.sin(m*z*np.tan(alpha))+np.sin(2*m*np.pi-m*z*np.tan(alpha)))
    # ans *= np.sin(alpha)**2/(16*m)
    ans = np.pi / 8 * (5 - 3 * np.cos(2 * alpha)) * np.sin(alpha) ** 2
    return ans


def inte_C(m, alpha):
    return np.pi


def inte_ftot(alpha, K, q, C, m):
    ans = 0.5 * K * (inte_SS(m, alpha) + inte_Tq(q, m, alpha) + inte_BB(m, alpha)) + 0.5 * C * inte_C(m, alpha)
    return ans



def Dalpha_inte_ftot(alpha, K, q, C, m):
    ans = 1 / 4 * K * np.pi * np.cos(alpha) ** 2
    ans *= -4 * q + 2 * m / np.cos(alpha) + 3 * np.cos(2 * alpha) * np.tan(alpha) + 2 * m / np.cos(alpha) ** 3 * (-2 * q + m / np.cos(alpha)) * np.tan(alpha) + 4 * q * np.tan(alpha) ** 2
    return ans


def opt_alpha_inte_ftot(K,q,C,m):
    pass


def find_alpha():
    K = 1
    qs = np.arange(0, 1.0, 0.1)
    C = 1
    m = 2
    alphas = []
    for i in range(len(qs)):
        root = optimize.fsolve(Dalpha_inte_ftot, 0.1, (K, qs[i], C, m))
        alphas.append(root[0])
    print("q:",qs)
    print("alpha",alphas)
    plt.figure()
    for m in [1,2,3]:
        alphas = []
        for i in range(len(qs)):
            root = optimize.fsolve(Dalpha_inte_ftot, 0.1, (K, qs[i], C, m))
            alphas.append(root[0])
        plt.plot(qs,np.cos(alphas),label="m=%.0f"%m)
    plt.xlabel("q",fontsize=12)
    plt.ylabel(r"$\alpha$",fontsize=12)
    plt.legend()
    plt.show()
    plt.close()
    return 0

def find_alpha_opt():
    print("find_alpha_opt()\n")
    K = 1
    qs = np.arange(0.0, 2.0, 0.1)
    C = 1
    plt.figure()
    axalpha = plt.subplot2grid((2, 2), (0, 0))
    axE = plt.subplot2grid((2, 2), (1, 0),sharex=axalpha)
    axalpha1 = plt.subplot2grid((2, 2), (0, 1))
    axE1 = plt.subplot2grid((2, 2), (1, 1),sharex=axalpha1)
    for m in [0,1,2,3,4]:
        alphas = []
        Es = []
        for i in range(len(qs)):
            res = optimize.minimize(inte_ftot,0.1, args=(K,qs[i],C,m),method="Powell")
            Es.append(res.fun)
            alphas.append(res.x)
            print("q:",qs)
            print("alpha",alphas)
        axalpha.plot(qs,4/(2*np.pi)*np.cos(alphas),label="m=%.0f"%m)
        axE.plot(qs,Es,"--",label="m=%.0f"%m)
        # compare with kc ~ 3/2 tan(4pi q*del_l) where del_l is the bead-bread distance for the triangle mech model
        del_l=0.01
        kc = 3/2*np.tan(4*np.pi*del_l*qs)
        axalpha1.plot(kc,4/(2*np.pi)*np.cos(alphas),label="m=%.0f"%m)
        axE1.plot(kc,Es,"--",label="m=%.0f"%m)
        #   $axalpha1.plot(3*np.tan(2*qs),np.cos(alphas),label="m=%.0f"%m)
        #   axE1.plot(3*np.tan(2*qs),Es,"--",label="m=%.0f"%m)
    axE.set_xlabel("q",fontsize=12)
    axE1.set_xlabel(r"$k_c$",fontsize=12)
    axalpha.set_ylabel(r"$\frac{4}{2\pi}\cos(\alpha)$",fontsize=12)
    axE.set_ylabel(r"$E$",fontsize=12)
    plt.legend()
    plt.show()
    plt.close()
    return 0