import numpy as np
import matplotlib.pyplot as plt
from Oplot import *


def E_compare(foldername, pars1, pars2, par_nm, par_dg, mode):
    colors, alphas = None, None
    data1, data2, O_label1, O_label2 = [], [], [], []
    xLabel = mode
    cpar_ind = find_cpar_ind(par_nm, mode)
    for i in range(len(pars1)):
        par1, par2 = pars1[i], pars2[i]
        par_dealing1, par_dealing2 = par1[:], par2[:]
        f2rtail1, f2rtail2 = "MC", "MC"
        label1, label2 = "N%d_"%pars1[i][0], "N%d_"%pars2[i][0]
        for j in range(len(par1)):
            if j == cpar_ind:
                f2rtail1 += "_" + par_nm[j] + "s"
                f2rtail2 += "_" + par_nm[j] + "s"
            else:
                f2rtail1 += "_" + par_nm[j] + "%.*f" % (par_dg[j], par_dealing1[j])
                f2rtail2 += "_" + par_nm[j] + "%.*f" % (par_dg[j], par_dealing2[j])
                if par_dealing1[j] != par_dealing2[j] or par_dealing1[j] != pars1[0][j] or par_dealing1[j] != pars1[-1][j]:
                    label1 += par_nm[j] + "%.*f," % (par_dg[j], par_dealing1[j])
                    label2 += par_nm[j] + "%.*f," % (par_dg[j], par_dealing2[j])
        f2rtail1 += "_ana.csv"
        f2rtail2 += "_ana.csv"
        head = "/O_"
        filename1 = foldername + head + f2rtail1
        filename2 = foldername + head + f2rtail2
        data1.append(np.loadtxt(filename1, skiprows=1, delimiter=",", unpack=True))
        data2.append(np.loadtxt(filename2, skiprows=1, delimiter=",", unpack=True))
        O_label1.append(label1)
        O_label2.append(label2)
    data1 = np.transpose(np.array(data1), axes=(1, 0, 2))
    data2 = np.transpose(np.array(data2), axes=(1, 0, 2))
    print("len(data1)", len(data1))
    # total energy
    cpar1, E_ave1 = data1[0], data1[1]
    cpar2, E_ave2 = data2[0], data2[1]
    # perimeter
    Le_ave_bf = []
    for e in range(pars1[0][2]):
        Le_ave_bf.append(data1[4 + 3 * e])
    Le_ave1 = np.sum(Le_ave_bf,axis=0)
    Le_ave_bf = []
    for e in range(pars2[0][2]):
        Le_ave_bf.append(data2[4 + 3 * e])
    Le_ave2 = np.sum(Le_ave_bf,axis=0)
    # mean curvature
    I2H2_ave1, p2uu_ave1, uuc_ave1, un2_ave1 = data1[10 + 3 * (pars1[0][2])], data1[19 + 3 * (pars1[0][2])], data1[22 + 3 * (pars1[0][2])], data1[25 + 3 * (pars1[0][2])]
    I2H2_ave2, p2uu_ave2, uuc_ave2, un2_ave2 = data2[10 + 3 * (pars2[0][2])], data2[19 + 3 * (pars2[0][2])], data2[22 + 3 * (pars2[0][2])], data2[25 + 3 * (pars2[0][2])]

    ppi = 72
    plt.figure()
    #plt.rc("text", usetex=True)
    fig, axs = plt.subplots(len(E_ave1), 7, figsize=(246 / ppi * 7, 246 / ppi * len(E_ave1) * 0.8), sharex=True)  # , sharex=True
    for i in range(len(E_ave1)):
        axs[i, 0].plot(cpar1[i], E_ave1[i], linestyle="--", label=O_label1[i])
        axs[i, 0].plot(cpar2[i], E_ave2[i], linestyle=":", label=O_label2[i])
        axs[i, 0].legend()
        axs[i, 1].plot(cpar1[i], E_ave1[i] - E_ave2[i], label="E_1-E_2")
        axs[i, 1].plot(cpar1[i], np.zeros(len(cpar1[i])))
        axs[i, 2].plot(cpar1[i], Le_ave1[i] - Le_ave2[i], label="Le_1-Le_2")
        axs[i, 2].plot(cpar1[i], np.zeros(len(cpar1[i])))
        axs[i, 3].plot(cpar1[i], I2H2_ave1[i] - I2H2_ave2[i], label="I2H2_1-I2H2_2")
        axs[i, 3].plot(cpar1[i], np.zeros(len(cpar1[i])))
        axs[i, 4].plot(cpar1[i], p2uu_ave1[i] - p2uu_ave2[i], label="p2uu_1-p2uu_2")
        axs[i, 4].plot(cpar1[i], np.zeros(len(cpar1[i])))
        axs[i, 5].plot(cpar1[i], uuc_ave1[i] - uuc_ave2[i], label="uuc_1-uuc_2")
        axs[i, 5].plot(cpar1[i], np.zeros(len(cpar1[i])))
        axs[i, 6].plot(cpar1[i], un2_ave1[i] - un2_ave2[i], label="un2_1-un2_2")
        axs[i, 6].plot(cpar1[i], np.zeros(len(cpar1[i])))
    axs[0, 0].set_ylabel("E/N")
    axs[0, 2].set_ylabel(r"$\int ds$")
    axs[0, 3].set_ylabel(r"$\int (2H)^2 dA$")
    axs[0, 4].set_ylabel(r"$\left<1.5 (u_i\cdot u_j)^2-0.5\right>_{(i,j)}$")
    axs[0, 5].set_ylabel(r"$u_c=\left<(u_i\times u_j)\cdot\hat{r}_{ij} (u_i\cdot u_j)\right>_{(i,j)}$")
    axs[0, 6].set_ylabel(r"$u_n=\left<(u_i\cdot n_i)^2\right>_{i}$")
    axs[-1, 0].set_xlabel(mode)
    axs[0, 1].legend()
    axs[-1, 1].set_xlabel(mode)

    plt.tight_layout()
    plt.savefig(foldername + "/E_compare_" + mode + ".pdf", format="pdf", transparent=True)
    plt.close()



def twist_calc(R1,R2,u1,u2):
    R1,R2,u1,u2 = np.array((R1,R2,u1,u2))
    r12 = R2-R1
    l12 = np.sqrt(np.inner(r12,r12))
    uuc = np.inner(np.cross(u1,u2),r12)
    uuc = uuc*np.inner(u1,u2)
    return uuc

# analyze tilt u.n and twist <uuc> correlation
def tilt_twist_calc(config_filename):
    # reutn un and uc
    data = np.loadtxt(config_filename, skiprows=6, delimiter=",", unpack=True)
    x, y, z, ux, uy, uz, nx, ny, nz, dA, d2H, ds, dAK, un2, enum, en0, en1 = data[:17]
    ns = np.transpose(data[17:])
    print("ns[0]:",ns[0])

    # recenter
    #x, y, z = x - np.average(x), y - np.average(y), z - np.average(z)
    uc = np.zeros(len(x))
    #find average uc per bead with it's neib

    # for each bead
    for i in range(len(x)):
        # for each neighbour
        ncount = 0
        for j in range(len(ns[i])):
            ni = int(ns[i][j])
            print("ni %d \n"%ni)
            if(ni!=-1):
                ncount+=1
                R1 = [x[i],y[i],z[i]]
                R2 = [x[ni],y[ni],z[ni]]
                u1 = [ux[i],uy[i],uz[i]]
                u2 = [ux[ni],uy[ni],uz[ni]]
                uc[i]+=0.5*twist_calc(R1,R2,u1,u2)
        uc[i]/=ncount
    return (un2,uc)

def simple_leastsq_linearreg(x,y):
    x_mean = np.mean(x)
    y_mean = np.mean(y)

    # calculate slope (m) and y-intercept (c) using formula: y = mx + c
    numerator = np.sum((x - x_mean) * (y - y_mean))
    denominator = np.sum((x - x_mean)**2)
    m = numerator / denominator
    c = y_mean - m * x_mean

    # predict the y-values using the linear model
    y_pred = m * x + c

    return(m,c,y_pred)

def tilt_twist_plot():

    lf = 25.0
    Kd = 4.0
    Cn = 4.0
    qs = [0, 1, 2.5, 3.5]

    un2s, ucs = [],[]
    for i in range(len(qs)):
        un2,uc = tilt_twist_calc("../data/Ne2/May12_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q%.1f_Cn4.0.csv"%qs[i])
        un2s.append(un2)
        ucs.append(uc)

    plt.figure()
    for i in range(len(un2s)):
        plt.scatter(un2s[i],ucs[i])
        m,c,uc_pred = simple_leastsq_linearreg(un2s[i],ucs[i])
        plt.plot(un2s[i],uc_pred, label=r"$u_c = %f (u\cdot n)^2 + %f ; q=%.f, corr=%.1f$"%(m,c,qs[i],np.corrcoef(un2s[i], ucs[i])[0, 1]))
    plt.legend()
    plt.show()
    plt.close()





