import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
from autocorrelation import *
from analyze import *
from scipy.optimize import curve_fit
from scipy import stats
import scipy.signal
import scipy.fft


def Os_pars_plot(foldername, pars, par_nm, par_dg, mode):
    colors, alphas = None, None
    Alldata, O_label = [], []
    xLabel = mode
    cpar_ind = find_cpar_ind(par_nm, mode)
    for i in range(len(pars)):
        par = pars[i]
        par_dealing = par[:]
        # par_dealing[cpar_ind] = par[cpar_ind][i]
        f2rtail = "MC"
        label = ""
        for j in range(len(par)):
            if j == cpar_ind:
                f2rtail += "_" + par_nm[j] + "s"
            else:
                f2rtail += "_" + par_nm[j] + "%.*f" % (par_dg[j], par_dealing[j])
                label += par_nm[j] + "%.*f," % (par_dg[j], par_dealing[j])
        f2rtail += "_ana.csv"
        head = "/O_"
        filename = foldername + head + f2rtail
        if(os.path.exists(filename)):
            print(filename,"exists")
            data = np.loadtxt(filename, skiprows=1, delimiter=",", unpack=True)
            if(len(Alldata)==0):
                Alldata = [[] for i in range(len(data))]
            for i in range(len(data)):
                Alldata[i].append(data[i])
            O_label.append(label)

    #data = np.transpose(np.array(data), axes=(1, 0, 2))
    Alldata = np.array(Alldata)
    cpar, E_ave, E_tau, E_err = Alldata[:4]
    for i in range(len(cpar)):
        print("np.shape(cpar[%d])"%i,np.shape(cpar[i]))
        print("np.shape(E_ave[%d])"%i,np.shape(E_ave[i]))
    print("E_ave", E_ave)
    Les_ave, Les_tau, Les_err = [], [], []
    Ik2s_ave, Ik2s_tau, Ik2s_err = [], [], []
    Ne = pars[0][2]
    for e in range(Ne):
        Les_ave.append(Alldata[4 + 3 * e])
        Les_tau.append(Alldata[5 + 3 * e])
        Les_err.append(Alldata[6 + 3 * e])
    Les_ave = np.array(Les_ave)
    Le_ave = np.sum(Les_ave, axis=0)
    Les_err = np.array(Les_err)
    #Le_err = np.sqrt(np.sum(np.power(Les_err, 2)))
    Le_err = np.zeros(np.shape(Le_ave))

    #IdA_ave, IdA_tau, IdA_err, I2H_ave, I2H_tau, I2H_err, I2H2_ave, I2H2_tau, I2H2_err, I2H2dis_ave, I2H2dis_tau, I2H2dis_err, IK_ave, IK_tau, IK_err, p2uu_ave, p2uu_tau, p2uu_err, uuc_ave, uuc_tau, uuc_err, un2_ave, un2_tau, un2_err, uz2_ave, uz2_tau, uz2_err, lb_ave, lb_tau, lb_err, Eubias_ave, Eubias_tau, Eubias_err = Alldata[7 + 3 * (Ne - 1) : 40 + 3 * (Ne - 1)]

    # got rid of Eu, added uz_abs (Aug 11, 2022)
    IdA_ave, IdA_tau, IdA_err, I2H_ave, I2H_tau, I2H_err, I2H2_ave, I2H2_tau, I2H2_err, I2H2dis_ave, I2H2dis_tau, I2H2dis_err, IK_ave, IK_tau, IK_err, p2uu_ave, p2uu_tau, p2uu_err, uuc_ave, uuc_tau, uuc_err, un2_ave, un2_tau, un2_err, uz2_ave, uz2_tau, uz2_err,  uz_ave, uz_tau, uz_err, lb_ave, lb_tau, lb_err = Alldata[7 + 3 * (Ne - 1) : 40 + 3 * (Ne - 1)]
    Lasym_ave, Lasym_tau, Lasym_err = Alldata[40 + 3 * (Ne - 1) :]
    uuc_grad, uuc_grad_err = [], []
    un2_grad, un2_grad_err = [], []
    cpar_tsqN = []
    # some special cases
    F_ave, F_err = [], []  # just for calculating the pulling force
    F4_ave, F4_err = [], []  # just for calculating the pulling force starting 4th data
    qn = []
    for i in range(len(pars)):
        if mode == "q":
            print(np.shape(uuc_ave),np.shape(cpar),np.shape(pars))
            print(np.shape(uuc_ave[i]),np.shape(cpar[i]))
            uuc_grad.append(np.gradient(uuc_ave[i], cpar[i]))
            uuc_grad_err.append(np.zeros(len(cpar[i])))
            un2_grad.append(np.gradient(un2_ave[i], cpar[i]))
            un2_grad_err.append(np.zeros(len(cpar[i])))

            # TODO: renormalized q related
            # normalized q
            print("pars[-1]=",pars[i][-1])
            qn.append(cpar[i]/pars[i][-1])
            print("qn[-1]\n")
            print(qn[-1])


        if mode == "lf":
            print("lf: cpar[i]",cpar[i])
            fa, fe = Chi2_gradient(cpar[i], E_ave[i], E_err[i], k=2)  # use near 2k+1 points
            F_ave.append(fa)
            F_err.append(fe)
            fa4, fe4 = Chi2_gradient(cpar[i][4:], E_ave[i][4:], E_err[i][4:], k=2)  # use near 2k+1 points
            F4_ave.append(fa4)
            F4_err.append(fe4)
    ppi = 72
    # LineWidth, FontSize, LabelSize = 1, 9, 8
    plt.figure()
    plt.rc("text", usetex=True)
    na = 13
    fig, axs = plt.subplots(na, 2, figsize=(246 / ppi * 2, 246 / ppi * 0.5 * na))#, sharex=True)  # , sharex=True
    # cpar_aj = cpar-np.outer([2.8, 2.0, 1.5, 0.8, 0], np.ones(len(cpar[0])))
    O_cpar_plot(axs[0, 0], E_ave, E_err, O_label, "E", r"$E/N$", cpar, colors, alphas)
    O_cpar_plot(axs[1, 0], Le_ave, Le_err, O_label, "Le", r"$\int ds$", cpar, colors, alphas)

    #O_cpar_plot(axs[0, 1], Lasym_ave, Lasym_err, O_label, "Lasym", r"assym of L", cpar, colors, alphas)
    axs[0,1].set_ylim(-0.1,1.1)
    if Ne == 2:
        if mode == "lf":
            O_cpar_plot(axs[1, 1], F_ave, F_err, O_label, "F", r"$\partial{E/N}/\partial{l_f}$", cpar, colors, alphas)
            axs[1, 1].set_ylim(-0.1, 0.5)
            #O_cpar_plot(axs[2, 1], F4_ave[:,4:], F4_err[:,4:], O_label, "F", r"$\partial{E/N}/\partial{l_f}$", cpar[:,4:], colors, alphas)
            axs[1, 1].set_ylim(-0.1, 0.5)

        Le_ave_diff = np.abs(Les_ave[1] - Les_ave[0])
        #Le_err_diff = np.sqrt(np.power(Les_err[1], 2) + np.power(Les_err[0], 2))
        #O_cpar_plot(axs[1, 1], Le_ave_diff, Le_err_diff, O_label, "Le_diff", r"$|\left<\int_0 ds\right>-\left<\int_1 ds\right>|$", cpar, colors, alphas)

        #O_cpar_plot(axs[3, 1], Le_ave_diff / Le_ave, Le_err_diff * 0, O_label, "Le_diff/Le_ave", r"$|\left<\int_0 ds\right>-\left<\int_1 ds\right>|/sum$", cpar, colors, alphas)
        #O_cpar_plot(axs[4, 1], Le_ave_diff / Le_ave, Le_ave_diff * 0, O_label, "Le_diff/Le_ave'", r"$\left<|\int_0 ds-\int_1 ds|\right>/sum$", cpar, colors, alphas)
        axs[4, 1].set_xlabel(xLabel)
        #Lelong_ave = (Le_ave + Le_ave_diff) / 2
        #Leshort_ave = (Le_ave - Le_ave_diff) / 2
        #O_cpar_plot(axs[5, 1], Lelong_ave, Lelong_ave * 0, O_label, "max(L1,L2)", r"$max(L1,L2)$", cpar, colors, alphas)
        axs[5, 1].set_xlabel(xLabel)
        #O_cpar_plot(axs[6, 1], Leshort_ave, Leshort_ave * 0, O_label, "min(L1,L2)", r"$min(L1,L2)$", cpar, colors, alphas)
        axs[6, 1].set_xlabel(xLabel)

    O_cpar_plot(axs[2, 0], IdA_ave, IdA_err, O_label, "IdA", r"$\int dA$", cpar, colors, alphas)
    O_cpar_plot(axs[3, 0], I2H_ave, I2H_err, O_label, "I2H", r"$\int dA (2H)$", cpar, colors, alphas)

    O_cpar_plot(axs[4, 0], I2H2_ave / (16 * np.pi), I2H2_err / (16 * np.pi), O_label, "I2H2", r"$\int dA (2H)^2/(16\pi)$", cpar, colors, alphas)
    O_cpar_plot(axs[5, 0], I2H2dis_ave, I2H2dis_err, O_label, "I2H2dis", r"$\int dA (2H-C_0)^2$", cpar, colors, alphas)

    O_cpar_plot(axs[6, 0], IK_ave, IK_err, O_label, "IK", r"$\int dA K$", cpar, colors, alphas)
    O_cpar_plot(axs[7, 0], p2uu_ave, p2uu_err, O_label, "p2uu", r"$\left<1.5 (u_i\cdot u_j)^2-0.5\right>_{(i,j)}$", cpar, colors, alphas)
    O_cpar_plot(axs[8, 0], uuc_ave, uuc_err, O_label, "uuc", r"$u_c=\left<(u_i\times u_j)\cdot\hat{r}_{ij} (u_i\cdot u_j)\right>_{(i,j)}$", cpar, colors, alphas)
    # plot the differential if mode is q
    O_cpar_plot(axs[9, 0], un2_ave, un2_err, O_label, "un2", r"$u_n^2=\left<(u_i\cdot n_i)^2\right>_{i}$", cpar, colors, alphas)
    O_cpar_plot(axs[10, 0], uz2_ave, uz2_err, O_label, "uz2", r"$u_z^2=\left<(u_i\cdot \hat{z})^2\right>_{i}$", cpar, colors, alphas)
    O_cpar_plot(axs[11, 0], uz_ave, uz_err, O_label, "uz2", r"$|u_z|=\left<|u_i\cdot \hat{z}|\right>_{i}$", cpar, colors, alphas)
    O_cpar_plot(axs[12, 0], lb_ave, lb_err, O_label, "lb", r"$\left<l_{ij}\right>$", cpar, colors, alphas)
    #O_cpar_plot(axs[12, 0], Eubias_ave, Eubias_err, O_label, "Eubias", r"$\left<Eu\right>_{w}$", cpar, colors, alphas)

    if mode == "q":
        #O_cpar_plot(axs[7, 1], uuc_grad, uuc_grad_err, O_label, "uuc_grad", r"$\partial u_c /\partial q$", cpar, colors, alphas)
        #O_cpar_plot(axs[8, 1], un2_grad, un2_grad_err, O_label, "un2_grad", r"$\partial u_n /\partial q$", cpar, colors, alphas)
        O_cpar_plot(axs[8, 1], uuc_ave, uuc_err, O_label, "uuc", r"$u_c$", qn, colors, alphas)
        O_cpar_plot(axs[9, 1], un2_ave, un2_err, O_label, "un2", r"$u_n^2$", qn, colors, alphas)
    # O_cpar_plot(axs[10,0], IKun2_ave, IKun2_err, O_label, "IKun2", r"$\left<\int dA K (u_i\cdot n_i)\right>$",cp ar, colors, alphas)
    # O_cpar_plot(axs[9,0], Itau2_ave, Itau2_err, O_label, "Itau2", r"$\left<\int ds \tau^2\right>$",cpar, colors, alphas)
    axs[10, 0].set_xlabel(xLabel)
    # axs[0,0].xaxis.set_major_locator(MultipleLocator(2))
    axs[0, 0].xaxis.set_minor_locator(MultipleLocator(1))
    lgd = axs[0, 0].legend(loc="upper center", bbox_to_anchor=(0.5, 0.25 * len(axs)))
    plt.tight_layout(pad=0)
    plt.savefig(foldername + "/Os_" + mode + ".pdf", format="pdf", bbox_extra_artists=(lgd,), bbox_inches="tight", transparent=True)
    plt.close()


def Geig_pars_plot(foldername, pars, par_nm, par_dg, mode):
    colors, alphas = None, None
    data, O_label = [], []
    xLabel = mode
    cpar_ind = find_cpar_ind(par_nm, mode)
    for i in range(len(pars)):
        par = pars[i]
        par_dealing = par[:]
        f2rtail = "/Gij"
        label = ""
        for j in range(len(par)):
            if j == cpar_ind:
                f2rtail += "_" + par_nm[j] + "s"
            else:
                f2rtail += "_" + par_nm[j] + "%.*f" % (par_dg[j], par_dealing[j])
                label += par_nm[j] + "%.*f," % (par_dg[j], par_dealing[j])
        f2rtail += "_ana.csv"
        filename = foldername + f2rtail
        data.append(np.loadtxt(filename, skiprows=1, delimiter=",", unpack=True))
        O_label.append(label)
    data = np.transpose(np.array(data), axes=(1, 0, 2))
    cpar, Dedge_ave, Dedge_tau, Dedge_err, Gxx_ave, Gxx_tau, Gxx_err, Gyy_ave, Gyy_tau, Gyy_err, Gzz_ave, Gzz_tau, Gzz_err, Geig0_ave, Geig0_tau, Geig0_err, Geig1_ave, Geig1_tau, Geig1_err, Geig2_ave, Geig2_tau, Geig2_err, GRg2_ave, GRg2_tau, GRg2_err, Gb_ave, Gb_tau, Gb_err, Gkap2_ave, Gkap2_tau, Gkap2_err = data

    ppi = 72
    plt.figure()
    plt.rc("text", usetex=True)
    fig, axs = plt.subplots(4, 3, figsize=(246 / ppi * 3, 246 / ppi * 0.6 * 4), sharex=True, sharey=False)

    O_cpar_plot(axs[0, 0], Gxx_ave, Gxx_err, O_label, "Gxx", r"$G_{xx}$", cpar, colors, alphas)
    O_cpar_plot(axs[1, 0], Gyy_ave, Gyy_err, O_label, "Gyy", r"$G_{yy}$", cpar, colors, alphas)
    O_cpar_plot(axs[2, 0], Gzz_ave, Gzz_err, O_label, "Gzz", r"$G_{zz}$", cpar, colors, alphas)

    O_cpar_plot(axs[0, 1], Geig0_ave, Geig0_err, O_label, "Geig0", r"$\lambda_0^G$", cpar, colors, alphas)
    O_cpar_plot(axs[1, 1], Geig1_ave, Geig1_err, O_label, "Geig1", r"$\lambda_1^G$", cpar, colors, alphas)
    O_cpar_plot(axs[2, 1], Geig2_ave, Geig2_err, O_label, "Geig2", r"$\lambda_2^G$", cpar, colors, alphas)
    O_cpar_plot(axs[3, 1], Dedge_ave, Dedge_err, O_label, "Dedge", r"$D_{e}$", cpar, colors, alphas)

    O_cpar_plot(axs[0, 2], GRg2_ave, GRg2_err, O_label, "GRg2", r"$R_g^2$", cpar, colors, alphas)
    O_cpar_plot(axs[1, 2], Gb_ave, Gb_err, O_label, "Gb", r"$asphericity:b_{G}$", cpar, colors, alphas)
    O_cpar_plot(axs[2, 2], Gkap2_ave, Gkap2_err, O_label, "Gkap2", r"$anisotropy:\kappa^2_{G}$", cpar, colors, alphas)

    axs[2, 1].set_xlabel(xLabel)
    # axs[0].set_ylim = (ylim_min,ylim_max)
    # axs[0].xaxis.set_minor_locator(MultipleLocator(0.5))
    lgd = axs[0, 0].legend(loc="upper center", bbox_to_anchor=(0.5, 0.7 + 0.25 * 4))
    plt.tight_layout(pad=0.5)
    plt.savefig(foldername + "/Gij_" + mode + ".pdf", format="pdf", bbox_extra_artists=(lgd,), bbox_inches="tight", transparent=True)
    plt.close()


def rhor_ave_plot(foldername, par, par_nm, par_dg, mode, del_r, tag="", leg_num=5):
    # plot density distribution from the center of mass
    rhor_all = []
    cpar_ind = find_cpar_ind(par_nm, mode)
    print("par,cpar_ind")
    print(par, cpar_ind)
    cpar = par[cpar_ind]
    for i in range(len(cpar)):
        par_dealing = par[:]
        par_dealing[cpar_ind] = par[cpar_ind][i]
        f2rtail = "/rhor"
        for j in range(0, len(par_dealing)):
            f2rtail += "_" + par_nm[j] + "%.*f" % (par_dg[j], par_dealing[j])
        f2rtail += ".csv"
        file2read = foldername + f2rtail
        rhodata = np.loadtxt(file2read, skiprows=0, delimiter=",", unpack=True)
        rhor_all.append(np.average(rhodata, axis=1))
    # save result to file
    f2stail = "/rhor"
    for j in range(0, len(par)):
        if j == cpar_ind:
            f2stail += "_" + par_nm[j] + "s"
        else:
            f2stail += "_" + par_nm[j] + "%.*f" % (par_dg[j], par_dealing[j])
    f2stail += "_ana.csv"
    savefile = foldername + f2stail

    if leg_num > len(cpar):
        leg_num = cpar
    leg_ind = np.linspace(0, len(cpar) - 1, leg_num, dtype=np.int)
    rhor_shift = np.linspace(0, 2, len(cpar))  ## this is negotiable

    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.8))
    ax = fig.add_subplot(111)
    bin_num = len(rhor_all[0])
    rplot = np.linspace(0.5 * del_r, bin_num * del_r - 0.5 * del_r, bin_num)

    for i in range(len(cpar)):
        cpar_plot = np.ones(bin_num) * cpar[i]
        LS = ":"
        Lb = None
        if i in leg_ind:
            LS = "-"
            Lb = mode + "=%.*f" % (par_dg[cpar_ind], cpar[i])

        ylabel = r"$\rho(r)/4\pi r^2\delta r$+ shift"
        ax.plot(rplot, rhor_all[i] / (4 * np.pi * rplot * rplot * del_r) + rhor_shift[i], linestyle=LS, label=Lb)

    ax.set_xlabel(r"r")

    ax.set_ylabel(ylabel)
    ax.legend(title=tag)
    plt.tight_layout()
    plt.savefig(savefile[:-4] + ".pdf", format="pdf")
    plt.close()


def O_cpar_plot(ax, O, O_err, O_label, O_name, O_latex, cpar, Colors, Alphas, Xscale="linear", Yscale="linear", ylim=None, Ms=3):
    print("plotting %s vs cpar" % O_name)
    print(len(O), len(cpar), len(O_err))
    for i in range(len(O)):
        if Colors:
            ax.errorbar(cpar[i], O[i], yerr=O_err[i], linestyle="None", color=Colors[i], alpha=Alphas[i], label=O_label[i])
        else:
            ax.errorbar(cpar[i], O[i], yerr=O_err[i], linestyle="None", marker="o", mfc="None", ms=Ms, label=O_label[i])
    # ax.xlabel("$\lambda\sigma_0/k_B T$", fontsize=FontSize)
    ax.set_ylabel(O_latex)
    ax.set_xscale(Xscale)
    ax.set_yscale(Yscale)
    if ylim:
        ax.set_ylim(ylim)
    ax.tick_params(direction="in", top="on", right="on")


def Os_Nes_pars_plot(foldername, Ne1pars, Ne2pars, mode):
    Ne1data, Ne1O_label, Ne1xLabel = Os_pars_data(foldername, Ne1pars, mode)
    Ne2data, Ne2O_label, Ne2xLabel = Os_pars_data(foldername, Ne2pars, mode)
    Ne1data = np.transpose(np.array(Ne1data), axes=(1, 0, 2))
    Ne1cpar, Ne1E_ave, Ne1E_tau, Ne1E_err = Ne1data[:4]
    Ne1Les_ave, Ne1Les_tau, Ne1Les_err = [], [], []
    for e in range(1):
        Ne1Les_ave.append(Ne1data[4 + 3 * e])
        Ne1Les_tau.append(Ne1data[5 + 3 * e])
        Ne1Les_err.append(Ne1data[6 + 3 * e])
    Ne1Le_ave = np.sum(Ne1Les_ave, axis=0)
    Ne1Le_err = np.sqrt(np.sum(np.power(Ne1Les_err, 2), axis=0))
    Ne2data = np.transpose(np.array(Ne2data), axes=(1, 0, 2))
    Ne2cpar, Ne2E_ave, Ne2E_tau, Ne2E_err = Ne2data[:4]
    Ne2Les_ave, Ne2Les_tau, Ne2Les_err = [], [], []
    for e in range(2):
        Ne2Les_ave.append(Ne2data[4 + 3 * e])
        Ne2Les_tau.append(Ne2data[5 + 3 * e])
        Ne2Les_err.append(Ne2data[6 + 3 * e])
    Ne2Le_ave = np.sum(Ne2Les_ave, axis=0)
    Ne2Le_err = np.sqrt(np.sum(np.power(Ne2Les_err, 2), axis=0))

    fig_row = len(Ne1cpar)
    ppi = 72
    # LineWidth, FontSize, LabelSize = 1, 9, 8
    plt.figure()
    plt.rc("text", usetex=True)
    fig, axs = plt.subplots(fig_row, 2, figsize=(246 / ppi, 246 / ppi * fig_row * 0.55))
    for i in range(len(axs)):
        E_ave = [Ne1E_ave[i], Ne2E_ave[i]]
        E_err = [Ne1E_err[i], Ne2E_err[i]]
        Le_ave = [Ne1Le_ave[i], Ne2Le_ave[i]]
        Le_err = [Ne1Le_err[i], Ne2Le_err[i]]
        O_label = [Ne1O_label[i], Ne2O_label[i]]
        cpar = [Ne1cpar[i], Ne2cpar[i]]
        O_cpar_plot(axs[i, 0], E_ave, E_err, O_label, "E", r"$\left<E\right>$", cpar, None, None)
        O_cpar_plot(axs[i, 1], Le_ave, Le_err, O_label, "Le", r"$\left<\int d s\right>$", cpar, None, None)
        axs[i, 0].xaxis.set_minor_locator(MultipleLocator(1))
        lgd = axs[i, 0].legend(loc="lower center")
    axs[-1, 0].set_xlabel(Ne1xLabel)
    axs[-1, 1].set_xlabel(Ne1xLabel)

    plt.tight_layout(pad=0)
    plt.savefig(foldername + "/Os_Nes_" + mode + ".pdf", format="pdf", bbox_extra_artists=(lgd,), bbox_inches="tight", transparent=True)
    plt.close()


def E_parts_cpar_plot(data, para, cpar_label, savename):
    N, kar, lam, karg, B, Bc, tau0 = para
    cpar, E_ave, E_tau, E_err, Le_ave, Le_tau, Le_err, I2H2_ave, I2H2_tau, I2H2_err, Ikg_ave, Ikg_tau, Ikg_err, Ik2_ave, Ik2_tau, Ik2_err, Itau_ave, Itau_tau, Itau_err, Itau2_ave, Itau2_tau, Itau2_err = data
    print("plotting energy by parts")
    ppi = 72
    plt.figure(figsize=(2 * 246 / ppi * 1, 246 / ppi * 0.8))
    plt.rc("text", usetex=True)
    plt.plot(cpar, -karg * Ikg_ave, "--", alpha=0.6, label=r"$-\bar{\kappa}\left<\int k_g ds\right>$")
    E_acum = -karg * Ikg_ave
    plt.fill_between(cpar, y1=E_acum, y2=lam * Le_ave + E_acum, alpha=0.6, label=r"$\lambda \left<\int ds\right>$")
    E_acum += lam * Le_ave
    plt.fill_between(cpar, y1=E_acum, y2=0.5 * kar * I2H2_ave + E_acum, alpha=0.6, label=r"$\frac{\kappa}{2} \left<\int (2H)^2 dA\right>$")
    E_acum += 0.5 * kar * I2H2_ave
    plt.fill_between(cpar, y1=E_acum, y2=0.5 * B * Ik2_ave + E_acum, alpha=0.6, label=r"$\frac{B}{2}\left<\int  k^2 ds\right>$")
    E_acum += 0.5 * B * Ik2_ave
    plt.fill_between(cpar, y1=E_acum, y2=0.5 * Bc * (Itau2_ave - 2 * tau0 * Itau_ave + tau0 * tau0 * Le_ave) + E_acum, alpha=0.6, label=r"$\frac{B'}{2}\left<\int (\tau-\tau_0)^2 ds\right>$")
    E_acum += 0.5 * Bc * (Itau2_ave - 2 * tau0 * Itau_ave + tau0 * tau0 * Le_ave)
    plt.plot(cpar, E_ave, "--", alpha=0.6, label=r"$\left<E \right>$")
    plt.xlabel(cpar_label)
    plt.ylabel("$E_i$")
    plt.ylim(-50, 2 * N)
    plt.tick_params(direction="in", top="on", right="on")
    plt.title(r"$N=%.0f,\kappa=%.0f,\bar{\kappa}=%.1f,B=%.0f,B'=%.0f,\tau_g^*=%.2f$" % (N, kar, karg, B, Bc, tau0))
    plt.legend()
    plt.tight_layout(pad=0)
    plt.savefig(savename, format="pdf", transparant=True)
    plt.close()


def exp_fit(x, a, b):
    return a * np.power(x, b)


def sqrt_fit(x, a):
    return a * np.sqrt(x)


def sqrt_fit_neg(x, a):
    return a / np.sqrt(x)


def a_bx(x, a, b):
    return a + b * x


def lamp_pars_plot(foldername, pars, par_nm, par_dg, mode, head):
    xLabel = mode
    Kd_ind = find_cpar_ind(par_nm, "Kd")
    Cn_ind = find_cpar_ind(par_nm, "Cn")
    Kd, Cn = [], []
    cpar_ind = find_cpar_ind(par_nm, mode)
    data, O_label = [], []
    for i in range(len(pars)):
        par = pars[i]
        Kd.append(par[Kd_ind])
        Cn.append(par[Cn_ind])
        par_dealing = par[:]
        f2rtail = "MC"
        label = ""
        for j in range(len(par)):
            if j == cpar_ind:
                f2rtail += "_" + par_nm[j] + "s"
            else:
                f2rtail += "_" + par_nm[j] + "%.*f" % (par_dg[j], par_dealing[j])
                label += par_nm[j] + "%.*f," % (par_dg[j], par_dealing[j])
        f2rtail += "_ana.txt"
        filename = foldername + "/" + head + "_" + f2rtail
        data.append(np.loadtxt(filename, skiprows=1, delimiter=",", unpack=True))
        O_label.append(label)
    data = np.transpose(np.array(data), axes=(1, 0, 2))
    cpar, lamp, lamperr, r, r_err, r_std, r_std_err = data

    ppi = 72
    plt.figure()
    plt.rc("text", usetex=True)
    fig, axs = plt.subplots(4, 2, figsize=(246 / ppi * 2, 246 / ppi * 4 * 0.8))
    Ylim = (0, 0.4)
    # original data
    O_cpar_plot(axs[0, 0], lamp, lamperr, O_label, "lamp", r"$\lambda_p$", cpar, None, None, ylim=Ylim, Ms=3)
    O_cpar_plot(axs[1, 0], r, r_err, O_label, "r", r"$\overline{r}$", cpar, None, None, Ms=3)
    O_cpar_plot(axs[2, 0], r_std, r_std_err, O_label, "r_std", r"$\delta\overline{r}$", cpar, None, None, Ms=3)

    # fit data (exponent) if cpar is Cn
    if mode == "Cn":
        None_label = ["" for i in range(len(cpar))]
        Kd = np.array(Kd)
        cpar_pKd = cpar / Kd[:, np.newaxis]
        O_cpar_plot(axs[0, 1], lamp, lamperr, None_label, "lamp", r"$\lambda_p$", cpar_pKd, None, None, Ms=0)
        cpar_pKd_all = cpar_pKd.flatten()
        lamp_all = lamp.flatten()
        lamperr_all = lamperr.flatten()
        cpmin = cpar_pKd_all.min()
        cpmax = cpar_pKd_all.max()
        cpar_plt = np.linspace(cpmin, cpmax, 100)

        popt, pcov = curve_fit(sqrt_fit_neg, cpar_pKd_all, lamp_all, sigma=lamperr_all, absolute_sigma=True)
        popterr = np.diag(pcov) ** 0.5
        print("[0,1] popt,popterr", popt, popterr)
        # chi2 = np.sum(np.power((sqrt_fit_neg(cpar_pKd_all,*popt) - lamp_all)/lamperr_all,2))
        chi2 = np.sum(np.power((sqrt_fit_neg(cpar_pKd_all, *popt) - lamp_all), 2) / sqrt_fit_neg(cpar_pKd_all, *popt))
        chi2nu = chi2 / (len(cpar_pKd_all) - 1)
        p = 1 - stats.chi2.cdf(chi2, df=len(cpar_pKd_all) - 1)
        print(p, 1 - stats.chi2.cdf(chi2nu * (len(cpar_pKd_all) - 1), (len(cpar_pKd_all) - 1)))
        print("popt", popt)
        axs[0, 1].plot(cpar_plt, sqrt_fit_neg(cpar_plt, *popt), alpha=0.8, label=r"$\lambda_p=%.2f/\sqrt{\frac{C_n}{K_d}}$" % popt[0])
        # axs[0,1].fill_between(cpar_plt,sqrt_fit_neg(cpar_plt,*popt-popterr),sqrt_fit_neg(cpar_plt,*popt+popterr),alpha=0.4,label=r"$\lambda_p=%.2f/\sqrt{\frac{C_n}{K_d}}$"%popt[0])
        axs[0, 1].legend()

        """
        popt,pcov=curve_fit(a_bx,np.log(cpar_pKd_all),np.log(lamp_all),sigma=lamperr_all/lamp_all,absolute_sigma=True)
        popterr = np.diag(pcov)**0.5
        print("line log fit:",popt,popterr)
        axs[1,1].plot(cpar_plt,np.exp(a_bx(np.log(cpar_plt),*popt)),alpha=0.8,label=r"$\lambda_p=%.2f(C_n/K_d)^{%.2f\pm%.2f},\chi^2_\nu=%.2f,p=%.3f$"%(np.exp(popt[0]),popt[1],popterr[1],chi2nu,p))
        #axs[1,1].fill_between(cpar_plt,np.exp(a_bx(np.log(cpar_plt),*popt-popterr)),np.exp(a_bx(np.log(cpar_plt),*popt+popterr)),alpha=0.4,label=r"$\lambda_p=%.2f(C_n/K_d)^{%.2f\pm%.2f},\chi^2_\nu=%.2f,p=%.3f$"%(np.exp(popt[0]),popt[1],popterr[1],chi2nu,p))
        """
        popt, pcov = curve_fit(exp_fit, cpar_pKd_all, lamp_all, sigma=lamperr_all, p0=[1.0, -0.5], absolute_sigma=False)
        popterr = np.diag(pcov) ** 0.5
        print("[1,1] popt,popterr", popt, popterr)
        chi2nu = np.sum(np.power((exp_fit(cpar_pKd_all, *popt) - lamp_all), 2) / exp_fit(cpar_pKd_all, *popt)) / (len(cpar_pKd_all) - 1)
        p = 1 - stats.chi2.cdf(chi2nu * (len(cpar_pKd_all) - 1), (len(cpar_pKd_all) - 1))
        print("chi2nu,p,", chi2nu, p)
        axs[1, 1].fill_between(cpar_plt, exp_fit(cpar_plt, *popt - popterr), exp_fit(cpar_plt, *popt + popterr), alpha=0.4, label=r"$\lambda_p=%.2f(C_n/K_d)^{%.2f\pm%.2f},\chi^2_\nu=%.2f,p=%.3f$" % (popt[0], popt[1], popterr[1], chi2nu, p))
        print(popt)

        axs[1, 1].legend()
        O_cpar_plot(axs[1, 1], lamp, lamperr, O_label, "lamp", r"$\lambda_p$", cpar_pKd, None, None, Ms=0)

        axs[1, 1].set_xscale("log")
        axs[1, 1].set_yscale("log")
        axs[1, 1].set_xlabel(r"$C_n/K_d$")
        # select data for lamp<0.2
        lp_max = 0.15
        lp_min = 0.0
        boollamp = np.logical_and(lp_min < lamp_all, lamp_all < lp_max)
        cpar_pKd_select = cpar_pKd_all[boollamp]
        lamp_select = lamp_all[boollamp]
        lamperr_select = lamperr_all[boollamp]
        popt, pcov = curve_fit(sqrt_fit_neg, cpar_pKd_select, lamp_select, sigma=lamperr_select, absolute_sigma=True)
        popterr = np.diag(pcov) ** 0.5
        print("[2,1] (select) popt,popterr", popt, popterr)
        chi2nu = np.sum(np.power((sqrt_fit_neg(cpar_pKd_select, *popt) - lamp_select), 2) / sqrt_fit_neg(cpar_pKd_select, *popt)) / (len(cpar_pKd_select) - 1)
        p = 1 - stats.chi2.cdf(chi2nu * (len(cpar_pKd_all) - 1), (len(cpar_pKd_all) - 1))
        print("chi2nu,p", chi2nu, p)
        axs[2, 1].fill_between(cpar_plt, sqrt_fit_neg(cpar_plt, *popt - popterr), sqrt_fit_neg(cpar_plt, *popt + popterr), alpha=0.4, label=r"$\lambda_p=%.2f/\sqrt{\frac{C_n}{K_d}}$" % popt[0])
        O_cpar_plot(axs[2, 1], [lamp_select], [lamperr_select], [""], "lamp", r"$\lambda_p(<%.1f)$ only" % lp_max, [cpar_pKd_select], None, None, Ms=0)
        axs[2, 1].legend()
        """
        popt,pcov=curve_fit(a_bx,np.log(cpar_pKd_select),np.log(lamp_select))
        popterr = np.diag(pcov)**0.5
        print("line log fit:",popt,popterr)
        axs[3,1].fill_between(cpar_plt,np.exp(a_bx(np.log(cpar_plt),*popt-popterr)),np.exp(a_bx(np.log(cpar_plt),*popt+popterr)),alpha=0.4,label=r"$\lambda_p=%.2f(C_n/K_d)^{%.2f\pm%.2f},\chi^2_\nu=%.2f,p=%.3f$"%(np.exp(popt[0]),popt[1],popterr[1],chi2nu,p))
        """
        popt, pcov = curve_fit(exp_fit, cpar_pKd_select, lamp_select, sigma=lamperr_select, p0=[1.0, -0.5], absolute_sigma=True)
        popterr = np.diag(pcov) ** 0.5
        print("[3,1] popt,popterr", popt, popterr)
        chi2nu = np.sum(np.power((exp_fit(cpar_pKd_select, *popt) - lamp_select), 2) / exp_fit(cpar_pKd_select, *popt)) / (len(cpar_pKd_select) - 1)
        print("chi2nu (select),", chi2nu)
        axs[3, 1].fill_between(cpar_plt, exp_fit(cpar_plt, *popt - popterr), exp_fit(cpar_plt, *popt + popterr), alpha=0.4, label=r"$\lambda_p=%.2f(C_n/K_d)^{%.2f\pm%.2f}$" % (popt[0], popt[1], popterr[1]))

        axs[3, 1].legend()
        O_cpar_plot(axs[3, 1], [lamp_select], [lamperr_select], [""], "lamp", r"$\lambda_p(<%.1f)$ only" % lp_max, [cpar_pKd_select], None, None, Ms=0)
        axs[3, 1].set_xscale("log")
        axs[3, 1].set_yscale("log")
        axs[3, 1].set_xlabel(r"$C_n/K_d$")

    elif mode == "Kd":
        Cn = np.array(Cn)
        cpar_pCn = cpar / Cn[:, np.newaxis]
        O_cpar_plot(axs[0, 1], lamp, lamperr, O_label, "lamp", r"$\lambda_p$", cpar_pCn, None, None, ylim=Ylim, Ms=5)
        popt, pcov = curve_fit(exp_fit, cpar_pCn.flatten(), lamp.flatten(), sigma=lamperr.flatten())
        popterr = np.diag(pcov) ** 0.5
        cpmin = cpar_pCn.flatten().min()
        cpmax = cpar_pCn.flatten().max()
        cpar_plt = np.linspace(cpmin, cpmax, 100)
        axs[1, 1].fill_between(cpar_plt, exp_fit(cpar_plt, *popt - popterr), exp_fit(cpar_plt, *popt + popterr), alpha=0.4, label=r"$\sim (K_d/C_n)^{%.2f \pm %.2f}$" % (popt[1], popterr[1]))
        axs[1, 1].legend()
        O_cpar_plot(axs[1, 1], lamp, lamperr, O_label, "lamp", r"$\lambda_p$", cpar_pCn, None, None, Ms=5)
        axs[1, 1].set_xscale("log")
        axs[1, 1].set_yscale("log")
        axs[1, 1].set_xlabel(r"$K_d/C_n$")
    axs[1, 0].set_xlabel(xLabel)
    lgd = axs[0, 0].legend(loc="upper center", bbox_to_anchor=(0.5, 0.75 * len(axs)))
    plt.tight_layout(pad=0.5)
    plt.savefig(foldername + "/O" + head + "_" + mode + ".pdf", format="pdf", bbox_extra_artists=(lgd,), bbox_inches="tight", transparent=True)
    plt.close()


def un2below_pars_plot(foldername, pars, par_nm, par_dg, mode, head):
    print("plotting un2below stats versus par (q?)")
    xLabel = mode  # let's say it's q
    cpar_ind = find_cpar_ind(par_nm, mode)
    data, O_label = [], []
    for i in range(len(pars)):
        par = pars[i]
        par_dealing = par[:]
        f2rtail = "MC"
        label = ""
        for j in range(len(par)):
            if j == cpar_ind:
                f2rtail += "_" + par_nm[j] + "s"
            else:
                f2rtail += "_" + par_nm[j] + "%.*f" % (par_dg[j], par_dealing[j])
                label += par_nm[j] + "%.*f," % (par_dg[j], par_dealing[j])
        f2rtail += "_ana.txt"
        filename = foldername + "/" + head + "_" + f2rtail
        data.append(np.loadtxt(filename, skiprows=1, delimiter=",", unpack=True))
        O_label.append(label)
    data = np.transpose(np.array(data), axes=(1, 0, 2))
    cpar, un2below, un2below_err = data

    # plot it
    ppi = 72
    plt.figure()
    plt.rc("text", usetex=True)
    fig, ax = plt.subplots(1, 1, figsize=(246 / ppi * 1, 246 / ppi * 1.5 * 0.8))
    # original data
    O_cpar_plot(ax, un2below, un2below_err, O_label, "un2below", r"$P(un2<1/2)$", cpar, None, None, Ms=3)
    ax.set_xlabel(xLabel)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.75))
    plt.tight_layout(pad=0.5)
    plt.savefig(foldername + "/O" + head + "_" + mode + ".pdf", format="pdf")
    plt.close()


def Chi2_gradient(x, y, yerr, k):

    # input x,y (n,1) array
    # output b,berr (n-2,1) array
    # fit near 2k+1 point using chi2 method

    b = np.zeros(len(x))
    berr = np.zeros(len(x))

    # first k points

    for i in range(k):
        x_t = np.array(x[: i + k + 1])
        y_t = np.array(y[: i + k + 1])
        sig_t = np.array(yerr[: i + k + 1])
        Delta = np.sum(1 / sig_t / sig_t) * np.sum(np.power(x_t / sig_t, 2)) - np.power(np.sum(x_t / sig_t / sig_t), 2)
        b[i] = 1 / Delta * (np.sum(np.power(sig_t, -2)) * np.sum(x_t * y_t / sig_t / sig_t) - np.sum(x_t / sig_t / sig_t) * np.sum(y_t / sig_t / sig_t))
        berr[i] = np.sqrt(1 / Delta * np.sum(np.power(sig_t, -2)))

    for i in range(len(x) - k)[k:]:
        x_t = np.array(x[i - k : i + k + 1])
        y_t = np.array(y[i - k : i + k + 1])
        sig_t = np.array(yerr[i - k : i + k + 1])
        Delta = np.sum(np.power(sig_t, -2)) * np.sum(np.power(x_t / sig_t, 2)) - np.power(np.sum(x_t / sig_t / sig_t), 2)
        # print("Delta: ",Delta)
        b[i] = 1 / Delta * (np.sum(1 / sig_t / sig_t) * np.sum(x_t * y_t / sig_t / sig_t) - np.sum(x_t / sig_t / sig_t) * np.sum(y_t / sig_t / sig_t))
        # print("b[%.0f]="%i,b[i])
        berr[i] = np.sqrt(1 / Delta * np.sum(np.power(sig_t, -2)))

    # last k points
    for i in range(k):
        i = -i - 1
        x_t = np.array(x[i - k :])
        y_t = np.array(y[i - k :])
        sig_t = np.array(yerr[i - k :])
        Delta = np.sum(1 / sig_t / sig_t) * np.sum(np.power(x_t / sig_t, 2)) - np.power(np.sum(x_t / sig_t / sig_t), 2)
        b[i] = 1 / Delta * (np.sum(np.power(sig_t, -2)) * np.sum(x_t * y_t / sig_t / sig_t) - np.sum(x_t / sig_t / sig_t) * np.sum(y_t / sig_t / sig_t))
        # print("b[%.0f]="%i,b[i])
        berr[i] = np.sqrt(1 / Delta * np.sum(np.power(sig_t, -2)))

    return b, berr


def plot_O_MCsteps(filename):
    data = np.loadtxt(filename, skiprows=13, delimiter=",", unpack=True)
    Ne = 1
    E, Le, IdA, I2H, I2H2, IK, Tp2uu, Tuuc, Bond_num, Tun2, IKun2 = data
    p2uu = Tp2uu / Bond_num
    uuc = Tuuc / Bond_num
    # N in file name is not real N
    print("OMCsteps")
    ppi = 72
    # LineWidth, FontSize, LabelSize = 1, 9, 8
    plt.figure()
    plt.rc("text", usetex=True)
    fig, axs = plt.subplots(5, 1, figsize=(246 / ppi, 246 / ppi * 0.6 * 5), sharex=True)
    # cpar_aj = cpar-np.outer([2.8, 2.0, 1.5, 0.8, 0], np.ones(len(cpar[0])))
    MCstep = np.array(range(len(E))) + 1
    axs[0].plot(MCstep, E, "-")
    axs[0].set_ylabel(r"$E$")
    axs[1].plot(MCstep, I2H2, "-")
    axs[1].set_ylabel(r"$\int (2H)^2 dA$")
    axs[2].plot(MCstep, Le, "-")
    axs[2].set_ylabel(r"$\oint ds$")
    axs[3].plot(MCstep, IK, "-")
    axs[3].set_ylabel(r"$\int K dA$")
    axs[4].plot(MCstep, Tp2uu, "-")
    axs[4].set_ylabel(r"$\sum_{(i,j)}P_2(u\cdot n)$")
    axs[4].set_xlabel(r"MC steps")
    plt.tight_layout(pad=0.5)
    plt.savefig(filename[:-4] + ".pdf")


def plot_vs_rescaled_q(foldername, pars, par_nm, par_dg, mode):
    pass

def un2_distribution_plot(filename,tag):
    data = np.loadtxt(filename, delimiter = ",", unpack = True)

    bin_num = len(data)
    un2pdf = np.average(data,axis=1)*bin_num
    un2x = np.linspace(0.5/bin_num,1-0.5/bin_num,bin_num)
    un2x_ave = np.sum(un2pdf*un2x)/bin_num

    plt.figure()
    plt.plot(un2x,un2pdf)
    plt.plot([un2x_ave,un2x_ave],[0,6],"-")
    plt.xlabel(r"$(u\cdot n)^2$")
    plt.ylabel(r"$p((u\cdot n)^2)$")
    plt.ylim(0,6)
    plt.legend(title=tag)
    plt.savefig(filename[:-4]+".png")
    plt.close()

def dA2H2_distribution_plot(filename,tag):
    data = np.loadtxt(filename, delimiter = ",", unpack = True)

    bin_num = len(data)
    dA2H2pdf = np.average(data,axis=1)*bin_num
    dA2H2x = np.linspace(0.5/bin_num,1-0.5/bin_num,bin_num)
    dA2H2x_ave = np.sum(dA2H2pdf*dA2H2x)/bin_num
    print("hello")
    plt.figure()
    plt.plot(dA2H2x,dA2H2pdf)
    #plt.plot([dA2H2x_ave,dA2H2x_ave],[0,6],"-")
    plt.xlabel(r"$dA (2H)^2$")
    plt.ylabel(r"$p(dA (2H)^2)$")
    #plt.ylim(0,6)
    plt.legend(title=tag)
    plt.savefig(filename[:-4]+".png")
    plt.close()
