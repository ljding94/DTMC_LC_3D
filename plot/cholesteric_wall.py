import numpy as np
import matplotlib.pyplot as plt
from config_plot import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator

# similar to nematic wall, plot smectic to cholesteric phase transition


def tilt_qs_Cn_data_get():
    foldername = "../data/Ne2/Nov30_2022"
    lf = 25
    Kd = 2.0
    Cns = [6, 9, 12]
    datas, labels, colors, markers = [], [], [], []
    colors = ["red", "green", "blue", "royalblue"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(Cns)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lf%.1f_kar50_C00.0_karg0.0_lam6.0_Kd2.0_qs_Cn%.1f_ana.csv" % (lf, Cns[i])
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))

    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    qs, un2_aves, un2_errs = datas[0], datas[31], datas[33]
    # normalize Kd by Cn

    labels = list(map(str, Cns))
    legendtitle = r"$l_f=%.0f, C$" % lf

    return [qs, un2_aves, un2_errs, labels, colors, markers, legendtitle]


def tilt_qs_config_data_get():
    # foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    # foldername = "../data/Ne2/Apr_2022/Apr9_2022"
    foldername = "../data/Ne2/Nov30_2022"
    lf = 25.0
    Cn = 6.0
    qs = [0.0, 1.5]
    qp = []
    fnames, povs, rotxyzs, xyshift, zslice = [], [], [], [], []
    for i in range(len(qs)):
        for j in range(3):
            fnames.append(foldername + "/State_N300_imod3_Ne2_lf%.1f_kar50_C00.0_karg0.0_lam6.0_Kd2.0_q%.1f_Cn%.1f.csv" % (lf, qs[i], Cn))
            qp.append(qs[i])
            povs.append("xy")
            rotxyzs.append([0, 0, 0])
            xyshift.append((10 * j, 10 * i))
        # 2 crossection each
        zslice += [(7, 13), (-2.5, 2.5), None]
        xyshift[3 * i + 2] = (25, 10 * i)
        povs[3 * i + 2] = "zx"
        rotxyzs[3 * i + 2] = [-np.pi * 3 / 8, np.pi * 3 / 8, 0]
    # fnames.append(fnames[-1])
    # povs.append("zx")
    # Kdp.append(Kds[-1])
    # rotxyzs.append([-np.pi*3/8,np.pi*3/8,0])
    # xyshift.append((10*2,0))
    # zslice.append(None)
    return [qp, fnames, povs, rotxyzs, xyshift, zslice]


def tilt_qs_lf_data_get():
    foldername = "../data/Ne2/Nov30_2022"
    Cn = 6.0
    lfs = [20, 25, 30]
    datas, labels, colors, markers = [], [], [], []
    colors = ["red", "green", "blue", "royalblue"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(lfs)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lf%.1f_kar50_C00.0_karg0.0_lam6.0_Kd2.0_qs_Cn%.1f_ana.csv" % (lfs[i], Cn)
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))

    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    qs, un2_aves, un2_errs = datas[0], datas[31], datas[33]
    # normalize Kd by lf

    labels = list(map(str, lfs))
    legendtitle = r"$C=%.0f,l_f$" % Cn

    return [qs, un2_aves, un2_errs, labels, colors, markers, legendtitle]

#TODO: change config to one per each
#TODO: change figure size
#TODO: fix labeling
def cholesteric_q_plot(LineWidth, FontSize, LabelSize):
    print("👌交给我吧")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axcfg = plt.subplot2grid((2, 2), (0, 0), colspan=2)
    axCns = plt.subplot2grid((2, 2), (1, 0), rowspan=1)
    axlfs = plt.subplot2grid((2, 2), (1, 1), rowspan=1, sharey=axCns, sharex=axCns)

    msize = 4

    ## configuration shows the tilt wall formation
    qs, fnames, povs, rotxyzs, xyshifts, zslices = tilt_qs_config_data_get()
    print(xyshifts)
    for i in range(len(qs)):
        ax_config_plot_xyz(axcfg, fnames[i], "gray", LineWidth, pov=povs[i], rotxyz=rotxyzs[i], xshift=xyshifts[i][0], yshift=xyshifts[i][1], zslice=zslices[i], mesh=1, bead=1, rod=1, d=1)
        if i % 3 == 0:
            # axcfg.text(xyshifts[i][0]-10,xyshifts[i][1]-5.5,r"$\epsilon_{LL}=%.1f$"%Kds[i],fontsize=FontSize)
            print("qs[i]", qs[i])
            #axcfg.text(xyshifts[i][0] - 10, xyshifts[i][1] - 5.5, r"$k_c=%.0f$" % qs[i], fontsize=FontSize)
    #axcfg.text(-10, 11, r"$C=6$", fontsize=FontSize)
    axcfg.tick_params(which="both", direction="in", bottom="off", top="off", right="off", left="off", labelbottom=False, labelleft=False, labelsize=LabelSize)
    x1, y1 = 0.85, 1.2
    axcfg.text(x1, y1, r"(a)", fontsize=FontSize, transform=axlfs.transAxes)

    ni = 0  # initial point
    n = 2  # date inteval
    ## tilt drops as Kd increases
    qs, un2_aves, un2_errs, labels, colors, markers, legendtitle = tilt_qs_Cn_data_get()
    for i in range(len(qs)):
        axCns.errorbar(qs[i][ni::n], un2_aves[i][ni::n], un2_errs[i][ni::n], ls=":", color=colors[i], mfc="None", marker=markers[i], ms=msize, label=labels[i])
    axCns.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=LabelSize)
    axCns.set_ylabel(r"$(\vu{u}\cdot\vu{n})^2$", fontsize=FontSize)
    axCns.set_ylim(0.4, 1.0)
    axCns.xaxis.set_major_locator(MultipleLocator(1))
    axCns.xaxis.set_minor_locator(MultipleLocator(0.5))
    axCns.yaxis.set_major_locator(MultipleLocator(0.1))
    axCns.yaxis.set_minor_locator(MultipleLocator(0.05))
    axCns.set_xlabel(r"$q$", fontsize=FontSize)
    axCns.legend(title=legendtitle, ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)
    x1, y1 = 0.85, 0.1
    axCns.text(x1, y1, r"(b)", fontsize=FontSize, transform=axCns.transAxes)

    ## critical q depends on the edge seperation distance
    qs, un2_aves, un2_errs, labels, colors, markers, legendtitle = tilt_qs_lf_data_get()
    for i in range(len(qs)):
        axlfs.errorbar(qs[i][ni::n], un2_aves[i][ni::n], un2_errs[i][ni::n], ls=":", color=colors[i], mfc="None", marker=markers[i], ms=msize, label=labels[i])
    axlfs.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=False, labelsize=LabelSize)
    # axlfs.set_ylabel(r"$(\vu{u}\cdot\vu{n})^2$", fontsize=FontSize)
    # axlfs.set_ylim(0.35,1.0)
    axlfs.set_xlabel(r"$q$", fontsize=FontSize)
    axlfs.legend(title=legendtitle, loc="upper right", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)
    x1, y1 = 0.85, 0.1
    axlfs.text(x1, y1, r"(c)", fontsize=FontSize, transform=axlfs.transAxes)

    plt.tight_layout(pad=0.1)
    plt.savefig("figures/cholesteric_wall.pdf", format="pdf")
