from unduloid_cal import *
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np
from mpl_toolkits import mplot3d
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle
from matplotlib import cm
from matplotlib.colors import Normalize
from nematic_cholesteric_2mod_cylinder_cal import *

def cylinder_2mod_plot(foldername, pars):
    all_data = []
    for i in range(len(pars)):
        K, C, m, R = pars[i]
        filename = foldername + "/optFtot_K%.1f_C%.1f_m%.0f_R%.1f_qs.csv" % (K, C, m, R)
        data = np.loadtxt(filename, skiprows=1, delimiter=",", unpack=True)
        all_data.append(data)

    ppi = 72
    plt.figure()
    #plt.rc("text", usetex=True)
    na = 3
    fig, axs = plt.subplots(na, 1, figsize=(246 / ppi * 1, 246 / ppi * 0.8 * na), sharex=True)  # , sharex=True
    # cpar_aj = cpar-np.outer([2.8, 2.0, 1.5, 0.8, 0], np.ones(len(cpar[0])))
    for i in range(len(all_data)):
        q, optFtot, optalpha, optgamma = all_data[i]
        K, C, m , R = pars[i]
        axs[0].plot(q, K*np.pi*np.ones(len(q)),"k:",label="m=0")

        axs[0].plot(q, optFtot, "-o", mfc="None",label="K%.1f_C%.1f_m%.0f_R%.1f" % (K, C, m, R))
        axs[1].plot(q, np.tan(optalpha), "-o", mfc="None",label="K%.1f_C%.1f_m%.0f_R%.1f" % (K, C, m, R))
        axs[2].plot(q, optgamma, "-o", mfc="None",label="K%.1f_C%.1f_m%.0f" % (K, C, m))

    axs[2].set_xlabel("q")
    axs[0].set_ylabel(r"$F_{tot}$")
    axs[0].set_yscale('symlog', linthreshy=1)
    axs[1].set_ylabel(r"$\tan(\alpha)$")
    axs[2].set_ylabel(r"$\gamma$")
    axs[1].legend()
    plt.tight_layout(pad=0.5)
    plt.savefig(filename[:-4] + ".pdf")
    plt.close()

def cylinder_2mod_plot_normlq(foldername, pars, normqlim=0.5):
    all_data = []
    for i in range(len(pars)):
        K, C, m, R = pars[i]
        filename = foldername + "/optFtot_K%.1f_C%.1f_m%.0f_R%.1f_qs.csv" % (K, C, m, R)
        data = np.loadtxt(filename, skiprows=1, delimiter=",", unpack=True)
        all_data.append(data)

    ppi = 72
    plt.figure()
    na = 3
    fig, axs = plt.subplots(na, 1, figsize=(246 / ppi * 1, 246 / ppi * 0.8 * na), sharex=True)  # , sharex=True
    # cpar_aj = cpar-np.outer([2.8, 2.0, 1.5, 0.8, 0], np.ones(len(cpar[0])))
    for i in range(len(all_data)):
        q, optFtot, optalpha, optgamma = all_data[i]
        K, C, m , R = pars[i]
        axs[0].plot(q/C, optFtot, "-o", mfc="None",label="K%.1f_C%.1f_m%.0f_R%.1f" % (K, C, m, R))
        axs[1].plot(q/C, np.tan(optalpha), "-o", mfc="None",label="K%.1f_C%.1f_m%.0f_R%.1f" % (K, C, m, R))
        axs[2].plot(q/C, optgamma, "-o", mfc="None",label="K%.1f_C%.1f_m%.0f" % (K, C, m))
    for ax in axs:
        ax.set_xlim(0,normqlim)
    print("??\n")
    axs[2].set_xlabel("q/C")

    axs[0].set_ylabel(r"$F_{tot}$")
    axs[1].set_ylabel(r"$\tan(\alpha)$")
    axs[2].set_ylabel(r"$\gamma$")
    axs[0].legend()
    plt.tight_layout(pad=0.5)
    plt.savefig(filename[:-4] + "_normlq.pdf")
    plt.close()

def find_critical_q(q,y,kind):
    trshd_dict = {"optalpha": 0.5, "optgamma": 0.9}
    trshd = trshd_dict[kind]
    print("trshd=",trshd)
    crtc_q, q_err = None, None
    for i in range(len(q)-1):
        if(y[i]<trshd and y[i+1]>trshd):
            crtc_q = 0.5*(q[i]+q[i+1])
            q_err = q[i+1]-q[i]
            break
    return [crtc_q,q_err]


def cylinder_critical_q_C(foldername, pars):
    all_data = []
    crtc_qs, q_errs = [], []
    Cs = []
    for i in range(len(pars)):
        K, C, m, R = pars[i]
        filename = foldername + "/optFtot_K%.1f_C%.1f_m%.0f_R%.1f_qs.csv" % (K, C, m, R)
        data = np.loadtxt(filename, skiprows=1, delimiter=",", unpack=True)
        q, optFtot, optalpha, optgamma = data
        print("q",q)
        print("optalpha",optalpha)
        crtc_q,q_err = find_critical_q(q,optalpha,"optalpha")
        if(crtc_q):
            Cs.append(C)
            crtc_qs.append(crtc_q)
            q_errs.append(q_err)

    ppi = 72
    plt.figure()
    fig, ax = plt.subplots(1, 1, figsize=(246 / ppi * 1, 246 / ppi * 0.8))  # , sharex=True
    ax.errorbar(Cs,crtc_qs,yerr=q_errs)
    plt.show()


def cylinder_c_q_hist(foldername, pars):
    all_data = []
    optalphas, optgammas = [], []
    for i in range(len(pars)):
        K, C, m, R = pars[i]
        filename = foldername + "/optFtot_K%.1f_C%.1f_m%.0f_R%.1f_qs.csv" % (K, C, m, R)
        data = np.loadtxt(filename, skiprows=1, delimiter=",", unpack=True)
        all_data.append(data)
        optalphas.append(data[2])
        optgammas.append(data[3])
    #q, optFtot, optalpha, optgamma = all_data[i]
    #K, C, m , R = pars[i]


    ppi = 72
    plt.figure()
    fig, axs = plt.subplots(1, 2, figsize=(246 / ppi * 1, 246 / ppi * 0.5), sharey=True, sharex=True)  # , sharex=True

    axs[0].pcolor(np.array(optalphas))
    axs[1].pcolor(np.array(optgammas))

    plt.show()




def A_z1_e_theta1_demo():
    e, theta1 = np.meshgrid(np.linspace(0.1, 0.9, 100), np.linspace(0.1 * np.pi, 0.9 * np.pi, 100))
    A = []
    for i in range(len(e)):
        A.append([])
        for j in range(len(e[0])):
            A[i].append(unduloid_A_unit_a(e[i][j], theta1[i][j]))
    A = np.array(A)
    lf = unduloid_lf_unit_a(e, theta1)
    print("np.shape(A),np.shape(lf)", np.shape(A), np.shape(lf))
    Af = A.flatten()
    lff = lf.flatten()

    fig = plt.figure()
    ax1 = fig.add_subplot(131, projection="3d")
    ax2 = fig.add_subplot(132, projection="3d")
    ax3 = fig.add_subplot(133, projection="3d")

    ax1.plot_surface(e, lf, A)
    ax1.set_zlabel("A")
    ax2.plot_surface(e, lf, theta1)
    ax2.set_zlabel("theta1")

    f = interpolate.interp2d(e, lf, theta1, kind="cubic")

    plt.show()
    plt.close()


def director_field_plot(m, alpha, gamma, bn_phi, bn_z):
    #u(m, alpha, gamma, phi, z)
    #phi, z = np.meshgrid(2 * np.pi * np.linspace(0.5 / bn_phi, 1 - 0.5 / bn_phi, bn_phi), np.linspace(0.5 / bn_z, 1 - 0.5 / bn_z, bn_z))
    phi, z = np.meshgrid(2 * np.pi * np.linspace(0.5 / bn_phi, 1 + 0.5 / bn_phi, bn_phi), np.linspace(0.5 / bn_z, 1 + 0.5 / bn_z, bn_z))
    X,Y,Z=np.cos(phi),np.sin(phi),z
    phi, z = phi.flatten(), z.flatten()


    plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot_surface(X,Y,Z,linewidth=0.5,shade=0,color="gray",edgecolor="black",alpha=0.5,rstride=5, cstride=5)
    #ax.plot_surface(X,Y,Z,linewidth=0.5,shade=0,color="gray",edgecolor="black",alpha=0.5)

    x,y = np.cos(phi),np.sin(phi)
    ux,uy,uz=u(m, alpha, gamma, phi, z)
    d = 0.5
    deg=np.arccos(np.abs(ux*x+uy*y))
    nu = np.abs(ux*x+uy*y)
    print("np.abs(ux*x+uy*y)",nu)
    print("mean=",np.average(nu))
    print("std=",np.std(nu))
    #deg=np.arccos(np.abs(uz))

    cmap = cm.get_cmap("jet_r")
    norm=Normalize(vmin=0,vmax=0.5*np.pi)

    for i in range(len(x))[::5]:
        #ax.plot3D([x[i]-0.5*d*ux[i],x[i]+0.5*d*ux[i]],[y[i]-0.5*d*uy[i],y[i]+0.5*d*uy[i]],[z[i]-0.5*d*uz[i],z[i]+0.5*d*uz[i]],"-",linewidth=1,color=cmap(norm(abs_un[i])),label=r"$u$")
        ax.plot3D([x[i]-0.5*d*ux[i],x[i]+0.5*d*ux[i]],[y[i]-0.5*d*uy[i],y[i]+0.5*d*uy[i]],[z[i]-0.5*d*uz[i],z[i]+0.5*d*uz[i]],"-",linewidth=1,color=cmap(norm(deg[i])),label=r"$u$")
    ax.set_frame_on(False)
    ax.set_axis_off()
    plt.show()
    plt.close()

