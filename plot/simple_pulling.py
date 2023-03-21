import numpy as np
import matplotlib.pyplot as plt
from config_plot import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator

## FIGURE: external force related


def force_normalized_kar_pull_data_get(gni=2):
    # get data  C0=0 for Kd=Cn=q=0 puilling case
    # foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    #foldername = "../data/Ne2/Apr_2022/Apr2_2022"
    foldername = "../data/Ne2/May13_2022"
    #kars = [20, 30, 40, 50, 60]
    kars = [20, 30, 50, 80]
    N = 300
    datas, labels, colors, markers = [], [], [], []
    colors = ["red", "green", "blue", "royalblue", "purple"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(kars)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lfs_kar%.0f_C00.0_karg0.0_lam6.0_Kd0.0_q0.0_Cn0.0_ana.csv" % kars[i]
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))
    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    lfs, E_aves, E_errs = datas[0], N * datas[1], N * datas[3]
    f_kar_norm_aves, f_kar_norm_errs = [], []
    for i in range(len(E_aves)):
        fa, fe = Chi2_gradient(lfs[i][gni:], E_aves[i][gni:], E_errs[i][gni:], k=3)
        f_kar_norm_aves.append(fa / kars[i])
        f_kar_norm_errs.append(fe / kars[i])
    labels = list(map(str, kars))
    legendtitle = r"$\kappa$"
    return [lfs[:,gni:], f_kar_norm_aves, f_kar_norm_errs, labels, colors, markers, legendtitle]


def force_pull_config_data_get():
    # foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    #foldername = "../data/Ne2/Apr_2022/Apr2_2022"
    #foldername = "../data/Ne2/May12_2022"
    foldername = "../data/Ne2/May13_2022"
    karp = 80
    lfs = [10, 20, 30]
    fnames = []
    for i in range(len(lfs)):
        fnames.append(foldername + "/State_N300_imod3_Ne2_lf%.1f_kar%.0f_C00.0_karg0.0_lam6.0_Kd0.0_q0.0_Cn0.0.csv" % (lfs[i], karp))
    return [lfs, fnames]


def force_q_pull_data_get(gni=2):
    foldername = "../data/Ne2/May13_2022"
    N = 300
    kar=50
    Kd = 4.0
    qs = [0,1,3]
    datas, labels, colors, markers = [], [], [], []
    colors = ["red", "green", "blue", "royalblue", "purple"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(qs)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lfs_kar%.0f_C00.0_karg0.0_lam6.0_Kd4.0_q%.1f_Cn4.0_ana.csv" % (kar,qs[i])
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))
    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    lfs, E_aves, E_errs = datas[0], N * datas[1], N * datas[3]
    f_kar_norm_aves, f_kar_norm_errs = [], []
    for i in range(len(E_aves)):
        fa, fe = Chi2_gradient(lfs[i][gni:], E_aves[i][gni:], E_errs[i][gni:], k=3)
        f_kar_norm_aves.append(fa / kar)
        f_kar_norm_errs.append(fe / kar)
    labels = list(map(str, qs))
    legendtitle = r"$\epsilon_{LL}=C=%.0f,k_c$"%Kd
    return [lfs[:,gni:], f_kar_norm_aves, f_kar_norm_errs, labels, colors, markers, legendtitle]

def force_q_Kd_pull_data_get(gni=2):
    foldername = "../data/Ne2/May13_2022"
    N = 300
    kar=50
    q = 1
    Kds = [2,4,6,8]
    datas, labels, colors, markers = [], [], [], []
    colors = ["red", "green", "blue", "royalblue", "purple"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(Kds)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lfs_kar50_C00.0_karg0.0_lam6.0_Kd%.1f_q1.0_Cn%.1f_ana.csv" % (Kds[i],Kds[i])
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))
    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    lfs, E_aves, E_errs = datas[0], N * datas[1], N * datas[3]
    f_kar_norm_aves, f_kar_norm_errs = [], []
    for i in range(len(E_aves)):
        fa, fe = Chi2_gradient(lfs[i][gni:], E_aves[i][gni:], E_errs[i][gni:], k=3)
        f_kar_norm_aves.append(fa / kar)
        f_kar_norm_errs.append(fe / kar)
    labels = list(map(str, Kds))
    legendtitle = r"$k_c=%.0f,\epsilon_{LL}=C$"%q
    return [lfs[:,gni:], f_kar_norm_aves, f_kar_norm_errs, labels, colors, markers, legendtitle]



def force_pull_plot(LineWidth, FontSize, LabelSize):
    # kappa and C0, lambda only, no LC interaction
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi, 246 / ppi * 1))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axfkar_cfg = plt.subplot2grid((2, 2), (0, 0), colspan=1)
    axfkar = plt.subplot2grid((2, 2), (0, 1), colspan=1)
    axfq = plt.subplot2grid((2, 2), (1, 0), colspan=1)
    axfKd = plt.subplot2grid((2, 2), (1, 1), colspan=1, sharey=axfq)

    msize = 4
    ni = 1  # initial point
    n = 2  # plot data inteval

    ## fkar_cfg, configuration for different lf with kar0
    # lfplt = [15.0, 25.0, 35.0]
    lfs, fnames = force_pull_config_data_get()
    for i in range(len(lfs)):
        ax_config_plot_xyz(axfkar_cfg, fnames[i], "gray", LineWidth, pov="zx", xshift=0.5 * lfs[i], yshift=13 * i, mesh=1, bead=0)
        axfkar_cfg.text(lfs[i] + 5, 13 * i, r"$%.0f$" % lfs[i], fontsize=FontSize)
    axfkar_cfg.text(lfs[-1], 13 * (len(lfs) - 1) + 5, r"$l_f=$", fontsize=FontSize)

    axfkar_cfg.tick_params(which="both", direction="in", bottom="off", top="off", right="off", left="off", labelbottom=False, labelleft=False, labelsize=LabelSize)
    x1, y1 = -0.25, 0.15
    axfkar_cfg.text(x1, y1, r"(a)", fontsize=FontSize, transform=axfkar.transAxes)
    # axfkar_cfg.text(axfkar_cfg.get_xlim()[1]*x1+axfkar_cfg.get_xlim()[0]* (1-x1),  axfkar_cfg.get_ylim()[1]*y1+axfkar_cfg.get_ylim()[0]* (1-y1), r"(a)", fontsize=FontSize,transform=axfkar_cfg.transAxes)


    # fkar
    ## f normalized versus lf for different kar at C0=0
    lfs, f_kar_norm_aves, f_kar_norm_errs, labels, colors, markers, legendtitle = force_normalized_kar_pull_data_get(gni=2)
    print("len(lfs)", len(lfs), len(f_kar_norm_aves))
    for i in range(len(lfs)):
        axfkar.errorbar(lfs[i][ni::n], f_kar_norm_aves[i][ni::n], f_kar_norm_errs[i][ni::n], ls="None", color=colors[i], mfc="None", marker=markers[i], ms=msize, label=labels[i])
    axfkar.yaxis.set_label_position("right")
    axfkar.yaxis.tick_right()
    axfkar.tick_params(which="both", direction="in", left="on", bottom="on",top="on", right="on", labelbottom=True, labelleft=False, labelsize=LabelSize)
    axfkar.set_ylabel(r"$F/\kappa$", fontsize=FontSize)
    axfkar.set_ylim(-0.25, 2.75)
    axfkar.xaxis.set_major_locator(MultipleLocator(5))
    axfkar.xaxis.set_minor_locator(MultipleLocator(2.5))
    axfkar.yaxis.set_major_locator(MultipleLocator(0.5))
    axfkar.yaxis.set_minor_locator(MultipleLocator(0.25))
    axfkar.legend(title=legendtitle, loc="upper left", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)
    axfkar.set_xlabel(r"$l_f$", fontsize=FontSize)
    x1, y1 = 0.85, 0.15
    # axfkar.text(axfkar.get_xlim()[1]*x1+axfkar.get_xlim()[0]* (1-x1),  axfkar.get_ylim()[1]*y1+axfkar.get_ylim()[0]* (1-y1), r"(b)", fontsize=FontSize,transform=axfkar.transAxes)
    axfkar.text(x1, y1, r"(b)", fontsize=FontSize, transform=axfkar.transAxes)


    # fq
    ni = 1  # initial point
    nf = 9
    n = 1  # plot data inteval
    ## f normalized versus lf for different q
    lfs, f_kar_norm_aves, f_kar_norm_errs, labels, colors, markers, legendtitle = force_q_pull_data_get()
    print("len(lfs)", len(lfs), len(f_kar_norm_aves))

    for i in range(len(lfs)):
        axfq.errorbar(lfs[i][ni:nf:n], f_kar_norm_aves[i][ni:nf:n], f_kar_norm_errs[i][ni:nf:n], ls="None", color=colors[i], mfc="None", marker=markers[i], ms=msize, label=labels[i])
    axfq.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=LabelSize)
    axfq.set_ylabel(r"$F/\kappa$", fontsize=FontSize)
    axfq.set_ylim(-0.05, 1.05)
    axfq.xaxis.set_major_locator(MultipleLocator(2))
    axfq.xaxis.set_minor_locator(MultipleLocator(1))
    axfq.yaxis.set_major_locator(MultipleLocator(0.2))
    axfq.yaxis.set_minor_locator(MultipleLocator(0.1))
    axfq.legend(title=legendtitle, loc="upper left", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)
    axfq.set_xlabel(r"$l_f$", fontsize=FontSize)
    x1, y1 = 0.85, 0.15
    # axfq.text(axfq.get_xlim()[1]*x1+axfq.get_xlim()[0]* (1-x1),  axfq.get_ylim()[1]*y1+axfq.get_ylim()[0]* (1-y1), r"(b)", fontsize=FontSize,transform=axfq.transAxes)
    axfq.text(x1, y1, r"(c)", fontsize=FontSize, transform=axfq.transAxes)



    # fq Kd

    lfs, f_kar_norm_aves, f_kar_norm_errs, labels, colors, markers, legendtitle = force_q_Kd_pull_data_get()
    print("len(lfs)", len(lfs), len(f_kar_norm_aves))

    for i in range(len(lfs)):
        axfKd.errorbar(lfs[i][ni:nf:n], f_kar_norm_aves[i][ni:nf:n], f_kar_norm_errs[i][ni:nf:n], ls="None", color=colors[i], mfc="None", marker=markers[i], ms=msize, label=labels[i])
    axfKd.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=False, labelsize=LabelSize)
    #axfKd.set_ylabel(r"$F/\kappa$", fontsize=FontSize)
    #axfKd.set_ylim(-10, 100)
    axfKd.xaxis.set_major_locator(MultipleLocator(2))
    axfKd.xaxis.set_minor_locator(MultipleLocator(1))
    #axfKd.yaxis.set_major_locator(MultipleLocator(20))
    #axfKd.yaxis.set_minor_locator(MultipleLocator(10))
    axfKd.legend(title=legendtitle, loc="upper left", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)
    axfKd.set_xlabel(r"$l_f$", fontsize=FontSize)
    x1, y1 = 0.85, 0.15

    axfKd.text(x1, y1, r"(d)", fontsize=FontSize, transform=axfKd.transAxes)

    plt.tight_layout(pad=0.1)
    plt.savefig("figures/force_pull.pdf", format="pdf")



def force_q_pull_plot(LineWidth, FontSize, LabelSize):
    # kappa and C0, lambda only, no LC interaction
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi, 246 / ppi * 0.5))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axfq = plt.subplot2grid((1, 2), (0, 0), colspan=1)
    axfKd = plt.subplot2grid((1, 2), (0, 1), colspan=1,sharey=axfq)
    # axfC0 = plt.subplot2grid((2, 2), (1, 0), colspan=1, sharex=axfkar)
    # axfC0_cfg = plt.subplot2grid((2, 2), (1, 1), colspan=1)
    msize = 4
    ni = 3  # initial point
    nf = 10
    n = 1  # plot data inteval

    # fq
    ## f normalized versus lf for different q
    lfs, f_kar_norm_aves, f_kar_norm_errs, labels, colors, markers, legendtitle = force_q_pull_data_get()
    print("len(lfs)", len(lfs), len(f_kar_norm_aves))

    for i in range(len(lfs)):
        axfq.errorbar(lfs[i][ni:nf:n], f_kar_norm_aves[i][ni:nf:n], f_kar_norm_errs[i][ni:nf:n], ls="None", color=colors[i], mfc="None", marker=markers[i], ms=msize, label=labels[i])
    axfq.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=LabelSize)
    axfq.set_ylabel(r"$F/\kappa$", fontsize=FontSize)
    axfq.set_ylim(-0.05, 0.75)
    axfq.xaxis.set_major_locator(MultipleLocator(2))
    axfq.xaxis.set_minor_locator(MultipleLocator(1))
    axfq.yaxis.set_major_locator(MultipleLocator(0.2))
    axfq.yaxis.set_minor_locator(MultipleLocator(0.1))
    axfq.legend(title=legendtitle, loc="upper left", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)
    axfq.set_xlabel(r"$l_f$", fontsize=FontSize)
    x1, y1 = 0.85, 0.15
    # axfq.text(axfq.get_xlim()[1]*x1+axfq.get_xlim()[0]* (1-x1),  axfq.get_ylim()[1]*y1+axfq.get_ylim()[0]* (1-y1), r"(b)", fontsize=FontSize,transform=axfq.transAxes)
    axfq.text(x1, y1, r"(a)", fontsize=FontSize, transform=axfq.transAxes)



    # fq Kd
    ni = 3  # initial point
    nf = 10
    n = 1  # plot data inteval
    lfs, f_kar_norm_aves, f_kar_norm_errs, labels, colors, markers, legendtitle = force_q_Kd_pull_data_get()
    print("len(lfs)", len(lfs), len(f_kar_norm_aves))

    for i in range(len(lfs)):
        axfKd.errorbar(lfs[i][ni:nf:n], f_kar_norm_aves[i][ni:nf:n], f_kar_norm_errs[i][ni:nf:n], ls="None", color=colors[i], mfc="None", marker=markers[i], ms=msize, label=labels[i])
    axfKd.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=False, labelsize=LabelSize)
    #axfKd.set_ylabel(r"$F/\kappa$", fontsize=FontSize)
    #axfKd.set_ylim(-10, 100)
    axfKd.xaxis.set_major_locator(MultipleLocator(2))
    axfKd.xaxis.set_minor_locator(MultipleLocator(1))
    #axfKd.yaxis.set_major_locator(MultipleLocator(20))
    #axfKd.yaxis.set_minor_locator(MultipleLocator(10))
    axfKd.legend(title=legendtitle, loc="upper left", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)
    axfKd.set_xlabel(r"$l_f$", fontsize=FontSize)
    x1, y1 = 0.85, 0.15

    axfKd.text(x1, y1, r"(b)", fontsize=FontSize, transform=axfKd.transAxes)
    plt.tight_layout(pad=0.1)
    plt.savefig("figures/force_q_pull.pdf", format="pdf")

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


def force_pull_C0_pull_data_get():
    foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    C0s = [0.0, 0.1, 0.2, 0.3]
    kar0 = 40
    N = 300
    datas, labels, colors, markers = [], [], [], []
    colors = ["red", "green", "blue", "royalblue"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(C0s)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lfs_kar%.0f_C0%.1f_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0_ana.csv" % (kar0, C0s[i])
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))

    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    lfs, E_aves, E_errs = datas[0], N * datas[1], N * datas[3]
    f_kar_norm_aves, f_kar_norm_errs = [], []
    for i in range(len(E_aves)):
        fa, fe = Chi2_gradient(lfs[i], E_aves[i], E_errs[i], k=2)
        f_kar_norm_aves.append(fa)
        f_kar_norm_errs.append(fe)
    labels = list(map(str, C0s))
    legendtitle = r"$\kappa=%.0f, C_0$" % kar0
    return [lfs, f_kar_norm_aves, f_kar_norm_errs, labels, colors, markers, legendtitle]


def force_pull_c0_config_data_get():
    foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    kar0 = 40
    lf = 20.0
    C0s = [0.0, 0.1, 0.2]
    fnames = []
    for i in range(len(C0s)):
        fnames.append(foldername + "/State_N300_imod3_Ne2_lf%.1f_kar%.0f_C0%.1f_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0.csv" % (lf, kar0, C0s[i]))
    return [lf, C0s, fnames]
