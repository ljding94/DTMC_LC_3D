import numpy as np
import matplotlib.pyplot as plt
from config_plot import *


## FIGURE: external force related

def force_normalized_kar_pull_data_get():
    # get data  C0=0 for Kd=Cn=q=0 puilling case
    foldername = "../data/Ne2/Oct18_2021"
    kar0=40
    kars=[20,40,60,80]
    datas, labels, colors, markers = [], [], [], []
    colors = ["red", "green", "blue", "royalblue"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(kars)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lfs_kar%.0f_C00.0_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0_ana.csv" % kars[i]
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))

    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    lfs, E_aves, E_errs = datas[0], datas[1], datas[3]
    f_kar_norm_aves, f_kar_norm_errs = [], []
    for i in range(len(E_aves)):
        fa, fe = Chi2_gradient(lfs[i], E_aves[i], E_errs[i], k=2)
        f_kar_norm_aves.append(fa * kar0 / kars[i])
        f_kar_norm_errs.append(fe * kar0 / kars[i])
    labels = list(map(str, kars))
    legendtitle=r"$\kappa_0=%.0f,\kappa$"%kar0
    return [lfs, f_kar_norm_aves, f_kar_norm_errs, labels, colors, markers, legendtitle]

def force_pull_config_data_get():
    foldername = "../data/Ne2/Oct18_2021"
    kar0=40
    lfs = [10.0,20.0,30.0]
    fnames= []
    for i in range(len(lfs)):
        fnames.append(foldername + "/State_N300_imod3_Ne2_lf%.1f_kar%.0f_C00.0_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0.csv" % (lfs[i],kar0))

    return [lfs,fnames]


def force_pull_C0_pull_data_get():
    foldername = "../data/Ne2/Oct18_2021"
    C0s = [0.0,0.1,0.2,0.3]
    kar0=40
    datas, labels, colors, markers = [], [], [], []
    colors = ["red", "green", "blue", "royalblue"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(C0s)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lfs_kar%.0f_C0%.1f_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0_ana.csv" % (kar0,C0s[i])
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))

    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    lfs, E_aves, E_errs = datas[0], datas[1], datas[3]
    f_kar_norm_aves, f_kar_norm_errs = [], []
    for i in range(len(E_aves)):
        fa, fe = Chi2_gradient(lfs[i], E_aves[i], E_errs[i], k=2)
        f_kar_norm_aves.append(fa)
        f_kar_norm_errs.append(fe)
    labels = list(map(str, C0s))
    legendtitle = r"$\kappa=%.0f, C_0$"%kar0
    return [lfs, f_kar_norm_aves, f_kar_norm_errs, labels, colors, markers, legendtitle]

def force_pull_c0_config_data_get():
    foldername = "../data/Ne2/Oct18_2021"
    kar0=40
    lf=20.0
    C0s = [0.0,0.1,0.2]
    fnames= []
    for i in range(len(C0s)):
        fnames.append( foldername + "/State_N300_imod3_Ne2_lf%.1f_kar%.0f_C0%.1f_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0.csv" % (lf,kar0,C0s[i]))
    return [lf,C0s,fnames]

def force_pull_plot(LineWidth, FontSize, LabelSize):
    # kappa and C0, lambda only, no LC interaction
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.8))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axfkar = plt.subplot2grid((2, 2), (0, 0), colspan=1)
    axfkar_cfg = plt.subplot2grid((2, 2), (0, 1), colspan=1)
    axfC0 = plt.subplot2grid((2, 2), (1, 0), colspan=1, sharex=axfkar)
    axfC0_cfg = plt.subplot2grid((2, 2), (1, 1), colspan=1)
    msize = 3
    n=2 # plot data inteval

    ## f normalized versus lf for different kar at C0=0
    lfs, f_kar_norm_aves, f_kar_norm_errs, labels, colors, markers, legendtitle = force_normalized_kar_pull_data_get()
    for i in range(len(lfs)):
        axfkar.errorbar(lfs[i][::n], f_kar_norm_aves[i][::n], f_kar_norm_errs[i][::n], ls="None", color=colors[i], mfc= "None", marker = markers[i], ms=msize, label=labels[i])

    axfkar.tick_params(which="both",direction="in", top="on", right="on",labelbottom=False, labelleft=True,labelsize=LabelSize)
    axfkar.set_ylabel(r"$F \kappa_0/\kappa$", fontsize=FontSize)
    axfkar.set_ylim(0.0,0.35)
    axfkar.legend(title=legendtitle,ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)

    ## fkar_cfg, configuration for different lf with kar0
    #lfplt = [15.0, 25.0, 35.0]
    lfs,fnames = force_pull_config_data_get()
    for i in range(len(lfs)):
        ax_config_plot_xyz(axfkar_cfg, fnames[i], "gray", LineWidth, pov="zx",xshift=0.5*lfs[i]-5,yshift=-13*i, mesh=1, bead=0)
        axfkar_cfg.text(lfs[i]-2,-13*i,r"$l_f=%.0f$"%lfs[i],fontsize=FontSize)
    axfkar_cfg.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)

    ## axfC0
    lfs, f_C0_aves, f_C0_errs, labels, colors, markers, legendtitle = force_pull_C0_pull_data_get()
    for i in range(len(lfs)):
        axfC0.errorbar(lfs[i][::n], f_C0_aves[i][::n], f_C0_errs[i][::n], ls="None", color=colors[i], mfc= "None", marker = markers[i], ms=msize, label=labels[i])

    axfC0.tick_params(which="both",direction="in", top="off", right="on",labelbottom=True, labelleft=True, labelsize=LabelSize)
    axfC0.set_ylabel(r"$F$", fontsize=FontSize)
    axfC0.set_ylim(0.0,0.35)
    axfC0.legend(title=legendtitle,ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)
    axfC0.set_xlabel(r"$l_f$",fontsize=FontSize)

    ## axfC0_cfg
    lf,C0s,fnames = force_pull_c0_config_data_get()
    for i in range(len(C0s)):
        ax_config_plot_xyz(axfC0_cfg, fnames[i], "gray", LineWidth, pov="zx",xshift=0.5*lf-5,yshift=-13*i, mesh=1, bead=0)
        axfC0_cfg.text(lf-2,-13*i,r"$C_0=%.1f$"%C0s[i],fontsize=FontSize)
    axfC0_cfg.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)


    plt.tight_layout(pad=0.1)
    plt.savefig("force_pull.pdf",format="pdf")



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
