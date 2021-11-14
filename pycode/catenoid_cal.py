import numpy as np
import matplotlib.pyplot as plt

def catenoid_A(b,z1,alp):
    A = 2*np.pi*np.sqrt(2)*b*b
    A *=np.sinh(z1/b)+np.sinh(z1*alp/b)
    return A
def catenoid_L(b,z1,alp):
    L = 2*np.pi*b*(np.cosh(z1/b)+np.cosh(z1*alp/b))
    return L

def catenoid_IK(b,z1,alp):
    # untegral of gaussian curvature
    pass


def catenoid_A_L_plot():
    alps = np.arange(0.1,0.91,0.2)
    z1s = np.arange(0.1,2,0.1)
    As = []
    Ls = []
    b=1
    for i in range(len(alps)):
        print("alp=",alps[i])
        As.append([])
        Ls.append([])
        for j in range(len(z1s)):
            print("z1=",z1s[j])
            As[i].append(catenoid_A(b, z1s[j], alps[i]))
            Ls[i].append(catenoid_L(b, z1s[j], alps[i]))
    plt.figure()
    for i in range(len(alps)):
        plt.plot(Ls[i],As[i],"-",label=r"$z_1/z_2=%.1f$"%alps[i])
    plt.plot(2*np.sqrt(np.pi*np.array(As[-1])),As[-1],"--",label="disk")
    plt.xlabel(r"$L/b$")
    plt.ylabel(r"$A/b^2$")
    plt.legend()
    plt.savefig("catenoid_A_L_plot.pdf",format="pdf")
    plt.close()