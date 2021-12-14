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
