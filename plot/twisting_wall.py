import numpy as np
import matplotlib.pyplot as plt
from config_plot import *
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

def twist_q_data_get():
    foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    lf=20.0
    Kds=[3.0,5.0,7.0]
    Cn=5.0
    datas, labels, colors, markers = [], [], [], []
    colors = ["red", "green", "blue", "royalblue"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(Kds)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lf%.1f_kar40_C00.0_karg0.0_lam5.0_Kd%.1f_qs_Cn%.1f_ana.csv" % (lf,Kds[i],Cn)
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))
    datas = np.transpose(np.array(datas), axes=(1, 0, 2))

    qs, uc_aves, uc_errs, uz2_aves,uz2_errs = datas[0], datas[28], datas[30],datas[34], datas[36]
    labels = list(map(str, Kds))
    legendtitle=r"$l_f=%.0f,\epsilon_{LL}$"%lf
    return [qs,uc_aves, uc_errs, uz2_aves,uz2_errs,labels,colors,markers,legendtitle]


def twist_q_config_data_get():
    foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    lf=20.0
    Kd=5.0
    Cn=5.0
    qs = [0.0,0.4,0.8,1.2]
    fnames,povs,rotxyzs=[],[],[]
    for i in range(len(qs)):
        fnames.append(foldername + "/State_N300_imod3_Ne2_lf%.1f_kar40_C00.0_karg0.0_lam5.0_Kd%.1f_q%.1f_Cn%.1f.csv" % (lf,Kd,qs[i],Cn))
        povs.append("zx")
        rotxyzs.append([np.pi/3,0,np.pi/2])
    return [qs,fnames,povs,rotxyzs]

def twist_q_plot(LineWidth, FontSize, LabelSize):
    print("👌接着交给我吧")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1.5))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axcfg = plt.subplot2grid((3, 1), (0, 0))
    axuc = plt.subplot2grid((3, 1), (2, 0),rowspan=1)
    axuz = plt.subplot2grid((3, 1), (1, 0),rowspan=1,sharex=axuc)
    msize=4

    # configuration plot
    qs,fnames,povs,rotxyzs = twist_q_config_data_get()
    for i in range(len(qs)):
        ax_config_plot_xyz(axcfg, fnames[i], "gray", LineWidth, pov=povs[i], rotxyz=rotxyzs[i],xshift=14*i, mesh=1, bead=0,rod=1, d=1)
        #axcfg.text(14*i-3,-14,r"$k_c=%.1f$"%qs[i],fontsize=FontSize)
        axcfg.text(14*i-3,-13,r"$k_c=%.1f$"%qs[i],fontsize=FontSize)
    axcfg.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)
    x1, y1 = 0.9, 1.06
    axcfg.text(x1,y1, r"(a)", fontsize=FontSize,transform=axuz.transAxes)

    # o vs q plot

    qs,uc_aves, uc_errs, uz2_aves,uz2_errs,labels,colors,markers,legendtitle =twist_q_data_get()

    for i in range(len(qs)):
        axuc.errorbar(qs[i],uc_aves[i],uc_errs[i], ls="None", color=colors[i],mfc="None",marker=markers[i],ms=msize,label=labels[i])

    axuc.tick_params(which="both",direction="in", top="on", right="on",labelbottom=False, labelleft=True,labelsize=LabelSize)
    axuc.set_ylabel(r"$\left<(\vu{u}_i\cross\vu{u}_j)\cdot \vu{r}_{ij} (\vu{u}_i\cdot\vu{u}_j)\right>_{(i,j)}$", fontsize=FontSize)
    axuc.set_ylim(0.0,0.22)
    axuc.xaxis.set_major_locator(MultipleLocator(0.2))
    axuc.xaxis.set_minor_locator(MultipleLocator(0.1))
    axuc.yaxis.set_major_locator(MultipleLocator(0.05))
    axuc.yaxis.set_minor_locator(MultipleLocator(0.025))
    #axuc.set_xlabel(r"$k_c$",fontsize=FontSize)
    axuc.legend(title=legendtitle,ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)
    x1, y1 = 0.9, 0.1
    axuc.text(x1,y1, r"(b)", fontsize=FontSize,transform=axuc.transAxes)


    for i in range(len(qs)):
        axuz.errorbar(qs[i],uz2_aves[i],uz2_errs[i], ls="None", color=colors[i],mfc="None",marker=markers[i],ms=msize,label=labels[i])

    axuz.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelleft=True,labelsize=LabelSize)
    axuz.set_ylabel(r"$(\vu{u}\cdot \vu{z})^2$", fontsize=FontSize)
    axuz.set_ylim(0.0,0.2)
    axuz.xaxis.set_major_locator(MultipleLocator(0.2))
    axuz.xaxis.set_minor_locator(MultipleLocator(0.1))
    axuz.yaxis.set_major_locator(MultipleLocator(0.05))
    axuz.yaxis.set_minor_locator(MultipleLocator(0.025))
    axuz.set_xlabel(r"$k_c$",fontsize=FontSize)
    axuz.legend(title=legendtitle,ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)
    x1, y1 = 0.9, 0.1
    axuz.text(x1,y1, r"(c)", fontsize=FontSize,transform=axuz.transAxes)

    plt.tight_layout(pad=0.1)
    plt.savefig("figures/twisting_wall.pdf",format="pdf")