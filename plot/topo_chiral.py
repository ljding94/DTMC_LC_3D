import numpy as np
import matplotlib.pyplot as plt
from config_plot import *

def topo_data_get(Ne=2):
    foldername = "../data/Ne2/Feb3_2022"
    if(Ne==2):
        fname = foldername+"/O_MC_N200_imod3_Ne2_lf0.0_kar20_C00.0_karg0.0_lam5.0_Kd3.0_qs_Cn3.0_ana.csv"
        data = np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True)
        q, Lasym_ave, Lasym_err = data[0],data[-3],data[-1]
    elif(Ne==3):
        fname = foldername+"/O_MC_N200_imod3_Ne3_lf0.0_kar20_C00.0_karg0.0_lam5.0_Kd3.0_qs_Cn3.0_ana.csv"
        data = np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True)
        q, Lasym_ave, Lasym_err = data[0],data[-3],data[-1]
    return (q,Lasym_ave, Lasym_err)


def topo_config_data_get(Ne=2):
    foldername = "../data/Ne2/Feb3_2022"
    if(Ne==2):
        qs = [0.6,1.8]
        rots = [(0.0,0.5,-0.8),(-1.0,0.0,0)]
        fnames = []
        for i in range(len(qs)):
            fnames.append(foldername+"/State_N200_imod3_Ne2_lf0.0_kar20_C00.0_karg0.0_lam5.0_Kd3.0_q%.1f_Cn3.0.csv"%qs[i])
        qps = [0.9,1.85]

    elif(Ne==3):
        qs = [0.6,1.8,2.7]
        rots = [(-0.4,-0.2,0.3),(0.0,0.0,0),(0.0,0.0,0.0)]
        fnames = []
        for i in range(len(qs)):
            fnames.append(foldername+"/State_N200_imod3_Ne3_lf0.0_kar20_C00.0_karg0.0_lam5.0_Kd3.0_q%.1f_Cn3.0.csv"%qs[i])
        qps = [0.6,1.85,2.6]
    return (qps,rots,fnames)


def topo_change_plot(LineWidth, FontSize, LabelSize):
    print("👌交给我吧")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1.4))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axNe2 = plt.subplot2grid((2, 1), (0, 0))
    axNe3 = plt.subplot2grid((2, 1), (1, 0), sharex=axNe2)
    msize = 3

    # Ne=2
    q,Lasym_ave, Lasym_err = topo_data_get(Ne=2)
    axNe2.errorbar(q,Lasym_ave,yerr=Lasym_err,linestyle="-",marker="o",mfc="None",ms=msize,label=r"$N_e=2$")

    axNe2.tick_params(which="both",direction="in", top="on", right="on",labelbottom=False, labelleft=True,labelsize=LabelSize)
    axNe2.set_ylabel(r"$|L_1-L_2|/(L_1+L_2)$", fontsize=FontSize)
    axNe2.set_ylim(0,1.0)
    axNe2.legend(ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)

    # config
    axNe2c = []
    qcs,rots,fnamecs = topo_config_data_get(Ne=2)
    for i in range(len(qcs)):
        axNe2c.append(fig.add_axes([qcs[i]/q[-1], 0.6, 0.15, 0.3]))
        ax_config_plot_xyz(axNe2c[i], fnamecs[i], "gray", LineWidth, rotxyz=rots[i],mesh=1, bead=0,rod=0,d=1)
        ax_config_plot_xyz(axNe2c[i], fnamecs[i], "gray", LineWidth, rotxyz=rots[i],yshift=-15,mesh=1, bead=0,rod=1,d=1)

    # Ne=3
    q,Lasym_ave, Lasym_err = topo_data_get(Ne=3)
    axNe3.errorbar(q,Lasym_ave,yerr=Lasym_err,linestyle="-",marker="o",mfc="None",ms=msize,label=r"$N_e=3$")

    axNe3.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelleft=True,labelsize=LabelSize)
    axNe3.set_ylabel(r"$|\sum L_k e^{i\pi k 2/3}/\sum L_k|$",fontsize=FontSize)
    axNe3.set_xlabel(r"$k_c$", fontsize=FontSize)
    axNe3.set_ylim(0,1.0)
    axNe3.legend(ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)

    # Ne3 config
    axNe3c = []
    qcs,rots,fnamecs = topo_config_data_get(Ne=3)
    qcs[0]=0.9
    for i in range(len(qcs)):
        axNe3c.append(fig.add_axes([qcs[i]/q[-1]-0.08, 0.1, 0.15, 0.3]))
        ax_config_plot_xyz(axNe3c[i], fnamecs[i], "gray", LineWidth, rotxyz=rots[i],mesh=1, bead=0,rod=0,d=1)
        ax_config_plot_xyz(axNe3c[i], fnamecs[i], "gray", LineWidth, rotxyz=rots[i],yshift=-15, mesh=1, bead=0,rod=1,d=1)

    '''

    # plot config
    Con_lams = [2.0, 5.0, 8.0]
    # configuration files
    for i in range(len(Con_lams)):
        #fConnames.append("../data/Ne1/Feb20_2020/State_N200_Ne1_L0_kar10_lam%.1f_karg0.0_B0.0_Bc0_tau00.00.txt" % Con_lams[i])
        fConnames.append(
            "../data/Ne1/Jun16_2020/State_N200_Ne1_L0_kar10_lam%.1f_karg0.0_B0.0_Bc0_tau00.00.txt" % Con_lams[i])
    axc = []
    axc.append(fig.add_axes([Con_lams[0]/10-0.022, 0.73, 0.13, 0.13]))
    config_xy_plot(axc[0], fConnames[0], "blue",LineWidth/4, multi=True)
    axc.append(fig.add_axes([Con_lams[1]/10, 0.73, 0.13, 0.13]))
    config_xy_plot(axc[1], fConnames[1], "blue",LineWidth/4,rotation=(0,1.6), multi=True)
    axc.append(fig.add_axes([Con_lams[2]/10, 0.73, 0.13, 0.13]))
    config_xy_plot(axc[2], fConnames[2], "blue",LineWidth/4, multi=True)
    axc[0].set_title("BP", y=-0.4, fontsize=LabelSize)
    axc[1].set_title("disk", y=-0.4, fontsize=LabelSize)
    axc[2].set_title("vesicle", y=-0.4, fontsize=LabelSize)
    '''
    plt.tight_layout(pad=0.5)
    plt.savefig("figures/topo_chirality.pdf",format="pdf")
