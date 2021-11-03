import matplotlib.pyplot as plt
from catelinder_cal import *


def F_smooth(n):
    # smooth force curve with near points
    pass

def catelinder_F_plot():
    print("plotting catelinder extension force")

    kars = [20, 40, 60, 80]
    lam = 2.0
    A = 300
    lfs = np.arange(10, 20, 0.2)
    Es = []
    Fs = []
    As = []
    optbs = []
    opths = []
    optls = []

    for i in range(len(kars)):
        print("calculating: kar=", kars[i])
        Es.append([])
        As.append([])
        optbs.append([])
        opths.append([])
        optls.append([])
        for j in range(len(lfs)):
            res = opt_bh_E(lam, kars[i], lfs[j], A)
            Es[i].append(res.fun)
            optbs[i].append(res.x[0])
            opths[i].append(res.x[1])
            optls[i].append(res.x[2])
            As[i].append(catelinder_A(res.x[0], res.x[1], res.x[2]))
        Fs.append(np.gradient(Es[i]) / 1)
    optbs=np.array(optbs)
    opths=np.array(opths)
    optls=np.array(optls)

    plt.figure()
    fig, axs = plt.subplots(8, 2, figsize=(246 / ppi * 2, 246 / ppi * 4), sharex=True)
    for i in range(len(kars)):
        axs[0,0].plot(lfs, Es[i], "+", label=r"$\kappa=%.0f$" % kars[i])
        axs[1,0].plot(lfs, Fs[i], "+", label=r"$\kappa=%.0f$" % kars[i])
        axs[2,0].plot(lfs, As[i], "+", label=r"$\kappa=%.0f$" % kars[i])
        axs[3,0].plot(lfs, optbs[i], "+", label=r"$\kappa=%.0f$" % kars[i])
        axs[4,0].plot(lfs, opths[i], "+", label=r"$\kappa=%.0f$" % kars[i])
        axs[5,0].plot(lfs, optls[i], "+", label=r"$\kappa=%.0f$" % kars[i])
        axs[6,0].plot(lfs, optbs[i]*4*np.pi*np.cosh(0.5*opths[i]/optbs[i]), "+", label=r"$\kappa=%.0f$" % kars[i])
        axs[7,0].plot(lfs, 2*np.pi*(optls[i]-opths[i])/optbs[i], "+", label=r"$\kappa=%.0f$" % kars[i])
    axs[4,0].plot(lfs, lfs, "--", label=r"$lf=lf$")
    axs[5,0].plot(lfs, lfs, "--", label=r"$lf=lf$")

    axs[0,0].set_ylabel("E")
    axs[1,0].set_ylabel("F")
    axs[2,0].set_ylabel("A")
    #axs[2,0].set_ylim(0, 800)

    axs[3,0].set_ylabel("b")
    axs[4,0].set_ylabel("h")
    axs[5,0].set_ylabel("l")
    axs[6,0].set_ylabel(r"$\int ds$")
    axs[7,0].set_ylabel(r"$\int (2H)^2 dA$")
    axs[7,0].set_xlabel("lf")
    axs[0,0].legend()
    plt.tight_layout(pad=0.1)

    plt.savefig("caetelinder_F.pdf", format="pdf")
    plt.close()


def catelinder_config_plot(b, h, l):
    # elongating along x direction
    x1 = np.linspace(-l / 2, -l / 2 + h / 2, 50)
    x2 = np.linspace(-l / 2 + h / 2, l / 2 - h / 2, 50)
    x3 = np.linspace(l / 2 - h / 2, l / 2, 50)
    r1 = b * np.cosh((x1+l/2-h/2) / b)
    r2 = b *np.ones(len(x2))
    r3 = b * np.cosh((x3-l/2+h/2) / b)

    x=np.concatenate((x1,x2,x3))
    r=np.concatenate((r1,r2,r3))
    #phi = np.linspace(0, 2 * np.pi, 20)
    #暂时没有必要

    plt.figure()
    ax=plt.subplot2grid((1,1), (0,0))
    ax.axis("equal")
    ax.plot(x,r,"b-")
    ax.plot(x,-r,"b-")
    #ax.plot([x[0],x[0]],[-r[0],r[0]],"b-")
    #ax.plot([x[-1],x[0-1]],[-r[-1],r[-1]],"b-")
    ax.plot([x2[0],x2[0]],[-b,b],"g--")
    ax.plot([x2[-1],x2[-1]],[-b,b],"g--")
    plt.savefig("catelinder_b%.1f_h%.1f_l%.1f.pdf"%(b,h,l),format="pdf")
    plt.close()
