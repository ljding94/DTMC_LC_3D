import numpy as np
import matplotlib.pyplot as plt
from unduloid_cal import *
from catelinder_plot import F_smooth

def unduloid_F_plot():
    print("plotting unduloid pulling force")

    lam = 2.0
    kar = 40
    C0s = [0.0, 0.05,0.1]
    A = 300
    lfs = np.arange(10, 20, 0.2)

    Es = []
    Fs = []
    As = []
    optas = []
    optes = []
    opttht1s = []

    for i in range(len(C0s)):
        print("calculating: C0=", C0s[i])
        Es.append([])
        As.append([])
        optas.append([])
        optes.append([])
        opttht1s.append([])
        for j in range(len(lfs)):
            res = opt_unduloid_E(lam, kar, C0s[i], lfs[j], A)
            Es[i].append(res.fun)
            optas[i].append(res.x[0])
            optes[i].append(res.x[1])
            opttht1s[i].append(res.x[2])
            As[i].append(unduloid_A(res.x[0], res.x[1], res.x[2]))
        Fs.append(np.gradient(Es[i])/ (lfs[1] - lfs[0]))
    optas = np.array(optas)
    optes = np.array(optes)
    opttht1s = np.array(opttht1s)
    ppi = 72
    plt.figure()
    na=6
    fig, axs = plt.subplots(na, 2, figsize=(246 / ppi * 2, 246 / ppi * na/2), sharex=True)
    for i in range(len(C0s)):
        axs[0,0].plot(lfs,Es[i],"x",label=r"$C_0=%.2f$" % C0s[i])
        axs[1,0].plot(lfs,Fs[i],"x",label=r"$C_0=%.2f$" % C0s[i])
        axs[1,1].plot(lfs,F_smooth(Fs[i],2),"x",label=r"$C_0=%.2f$" % C0s[i])
        axs[2,1].plot(lfs,F_smooth(Fs[i],4),"x",label=r"$C_0=%.2f$" % C0s[i])
        axs[2,0].plot(lfs,As[i],"x",label=r"$C_0=%.2f$" % C0s[i])
        axs[3,0].plot(lfs,optas[i],"x",label=r"$C_0=%.2f$" % C0s[i])
        axs[4,0].plot(lfs,optes[i],"x",label=r"$C_0=%.1f$" % C0s[i])
        axs[5,0].plot(lfs,opttht1s[i],"x",label=r"$C_0=%.1f$" % C0s[i])
    axs[0,0].set_ylabel("E")
    axs[1,0].set_ylabel("F")
    axs[1,1].set_ylabel("F sm=2")
    axs[2,1].set_ylabel("F sm=4")
    axs[2,0].set_ylabel("A")
    axs[3,0].set_ylabel("a")
    axs[4,0].set_ylabel("e")
    axs[5,0].set_ylabel("theta1")
    axs[5,0].set_xlabel("lf")
    axs[0,0].legend()
    plt.tight_layout(pad=0.1)
    plt.savefig("unduloid_F.pdf", format="pdf")
    plt.close()

def unduloid_config_plot(a,e,theta1):
    thetas = np.linspace(-theta1,theta1,50)
    zs = unduloid_z(a, e, thetas)
    rs = unduloid_r(a, e, thetas)

    plt.figure()
    ax=plt.subplot2grid((1,1), (0,0))
    ax.axis("equal")
    ax.plot(zs,rs,"b-")
    ax.plot(zs,-rs,"b-")
    plt.savefig("unduloid_a%.1f_e%.1f_theta1%.1f.pdf"%(a,e,theta1),format="pdf")
    plt.close()