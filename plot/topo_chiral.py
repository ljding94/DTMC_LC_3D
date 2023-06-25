import numpy as np
import matplotlib.pyplot as plt
from config_plot import *
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

def topo_id_select(Ne=2):
    id_select = []
    if Ne==2:
        id_select = [4]*31
        id_select[5] = 0
        id_select[10] = 2
        id_select[11] = 2
        id_select[26] = 2
    elif Ne==3:
        id_select = [5]*31
        id_select[6] = 2
        id_select[15],id_select[16] = 1,3
    return id_select

def Ts_data_get(Ne=2):
    #foldername = "../data/Ne2/Feb_2022/Feb28_2022"
    foldername = "../data/Ne2/May16_2023"
    id_select = topo_id_select(Ne)
    max_id_num = np.max(id_select)
    q_per_id, uc_ave_per_id,uc_err_per_id = [0]*(max_id_num+1),[0]*(max_id_num+1),[0]*(max_id_num+1)
    for id in range(max_id_num+1):
        fname = foldername+"/O_MC_N300_imod3_Ne%d_lf0.0_kar30_C00.0_karg0.0_lam6.0_Kd4.0_qs_Cn4.0_id%d_ana.csv"%(Ne,id)
        print("fname:",fname)
        data = np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True)
        q, uc_ave, uc_err = data[0], data[3*Ne+22], data[3*Ne+24]
        q_per_id[id] = q
        uc_ave_per_id[id] = uc_ave
        uc_err_per_id[id] = uc_err
    select_q,select_uc_ave,select_uc_err = [],[],[]
    for i in range(len(q)):
        select_q.append(q_per_id[id_select[i]][i])
        select_uc_ave.append(uc_ave_per_id[id_select[i]][i])
        select_uc_err.append(uc_err_per_id[id_select[i]][i])
    return (select_q,select_uc_ave,select_uc_err)

def topo_data_get(Ne=2):
    foldername = "../data/Ne2/May16_2023"
    id_select = topo_id_select(Ne)
    max_id_num = np.max(id_select)
    q_per_id, Lasym_ave_per_id,Lasym_err_per_id = [0]*(max_id_num+1),[0]*(max_id_num+1),[0]*(max_id_num+1)
    for id in range(max_id_num+1):
        fname = foldername+"/O_MC_N300_imod3_Ne%d_lf0.0_kar30_C00.0_karg0.0_lam6.0_Kd4.0_qs_Cn4.0_id%d_ana.csv"%(Ne,id)
        data = np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True)
        q, Lasym_ave, Lasym_err = data[0],data[-3],data[-1]
        q_per_id[id], Lasym_ave_per_id[id],Lasym_err_per_id[id] = q, Lasym_ave, Lasym_err
    select_q, select_Lasym_ave, select_Lasym_err = [],[],[]
    for i in range(len(q)):
        select_q.append(q_per_id[id_select[i]][i])
        select_Lasym_ave.append(Lasym_ave_per_id[id_select[i]][i])
        select_Lasym_err.append(Lasym_err_per_id[id_select[i]][i])

    return (select_q, select_Lasym_ave, select_Lasym_err)

def topo_config_data_get(Ne=2):
    #TODO: find the config map from q_select
    id_select = topo_id_select(Ne)
    q_all = np.arange(0.0,3.01,0.1)
    print("q_all",q_all)
    if(Ne==2):
        #foldername = "../data/Ne2/Feb_2022/Feb28_2022"
        foldername = "../data/Ne2/May16_2023"
        qs = [0.2,2.0]
        qps = [0.4,1.6]
        #rots = [(0.0,0.5,-0.8),(-1.0,0.0,0)]
        rots = [(1.4,0.0,0.0),(1.0,0.0,0)]
        fnames = []
        for i in range(len(qs)):
            print("qs[i]",qs[i])
            #print("np.where(q_all==qs[i]),",np.where(q_all==qs[i]))
            #print("np.where(q_all==qs[i])[0][0],",np.where(q_all==qs[i])[0][0])
            id = id_select[int(qs[i]/0.1)]
            fnames.append(foldername+"/State_N300_imod3_Ne2_lf0.0_kar30_C00.0_karg0.0_lam6.0_Kd4.0_q%.1f_Cn4.0_id%d.csv"%(qs[i],id))


    elif(Ne==3):
        #foldername = "../data/Ne2/Feb_2022/Feb28_2022"
        foldername = "../data/Ne2/May16_2023"
        qs = [0.2,1.4,2.2]
        qps = [0.4,1.2,1.8]
        rots = [(0.5,0.0,0.0),(0.0,0.0,0.5),(0.0,-0.5,0.0)]
        fnames = []
        for i in range(len(qs)):
            id = id_select[int(qs[i]/0.1)]
            #id = id_select[np.where(q_all==qs[i])[0][0]]
            fnames.append(foldername+"/State_N300_imod3_Ne3_lf0.0_kar30_C00.0_karg0.0_lam6.0_Kd4.0_q%.1f_Cn4.0_id%d.csv"%(qs[i],id))

    return (qps,rots,fnames)

def topo_change_plot(LineWidth, FontSize, LabelSize):
    print("👌交给我吧")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1.2))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    axNe2 = plt.subplot2grid((3, 1), (0, 0))
    axNe3 = plt.subplot2grid((3, 1), (1, 0),sharex=axNe2,sharey=axNe2)
    axTs = plt.subplot2grid((3, 1), (2, 0),sharex=axNe2)
    msize = 4
    start = 0
    cut = -3
    skip = 2


    # Ne=2
    q,Lasym_ave, Lasym_err = topo_data_get(Ne=2)
    q,Lasym_ave, Lasym_err = q[start:cut:skip],Lasym_ave[start:cut:skip], Lasym_err[start:cut:skip]
    axNe2.errorbar(q,Lasym_ave,yerr=Lasym_err,linestyle="-",marker="o",mfc="None",ms=msize,label=r"$N_e=2$")


    axNe2.tick_params(which="both",direction="in", top="on", right="on",labelbottom=False, labelleft=True,labelsize=LabelSize)
    axNe2.set_ylabel(r"$|L_1-L_2|/(L_1+L_2)$", fontsize=FontSize)
    axNe2.set_ylim(0,0.95)
    #axNe2.set_xlabel(r"$k_c$", fontsize=FontSize)
    axNe2.xaxis.set_major_locator(MultipleLocator(0.5))
    axNe2.xaxis.set_minor_locator(MultipleLocator(0.25))
    axNe2.yaxis.set_major_locator(MultipleLocator(0.2))
    axNe2.yaxis.set_minor_locator(MultipleLocator(0.1))
    axNe2.legend(ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)
    x1, y1 = 0.9, 0.15
    axNe2.text(x1,y1, r"(a)", fontsize=FontSize,transform=axNe2.transAxes)

    # config
    axNe2c = []
    qcs,rots,fnamecs = topo_config_data_get(Ne=2)
    for i in range(len(qcs)):
        pass
        axNe2c.append(fig.add_axes([qcs[i]/q[cut], 0.3*1/3+2/3, 0.15, 0.3*2/3]))
        ax_config_plot_xyz(axNe2c[i], fnamecs[i], "gray", LineWidth, rotxyz=rots[i],mesh=1, bead=0,rod=0,d=1,pwlim=0)
        ax_config_plot_xyz(axNe2c[i], fnamecs[i], "gray", LineWidth, rotxyz=rots[i],yshift=-20,mesh=1, bead=0,rod=0,d=0.8,pwlim=np.pi/3)

    # Ne=3
    q,Lasym_ave, Lasym_err = topo_data_get(Ne=3)
    q,Lasym_ave, Lasym_err = q[start:cut:skip],Lasym_ave[start:cut:skip], Lasym_err[start:cut:skip]
    axNe3.errorbar(q,Lasym_ave,yerr=Lasym_err,linestyle="-",marker="o",mfc="None",ms=msize,label=r"$N_e=3$")

    axNe3.tick_params(which="both",direction="in", top="on", right="on",labelbottom=False, labelleft=True,labelsize=LabelSize)
    axNe3.set_ylabel(r"$|\sum L_k e^{i\pi k 2/3}/\sum L_k|$",fontsize=FontSize)
    #axNe3.set_xlabel(r"$k_c$", fontsize=FontSize)
    #axNe3.set_ylim(0,1.05)
    axNe3.xaxis.set_major_locator(MultipleLocator(0.5))
    axNe3.xaxis.set_minor_locator(MultipleLocator(0.25))
    axNe3.yaxis.set_major_locator(MultipleLocator(0.2))
    axNe3.yaxis.set_minor_locator(MultipleLocator(0.1))
    axNe3.legend(ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)
    x1, y1 = 0.9, 0.15
    axNe3.text(x1,y1, r"(b)", fontsize=FontSize,transform=axNe3.transAxes)

    # Ne3 config
    axNe3c = []
    qcs,rots,fnamecs = topo_config_data_get(Ne=3)
    for i in range(len(qcs)):
        pass
        axNe3c.append(fig.add_axes([qcs[i]/q[cut], 0.3*1/3+1/3, 0.15, 0.3*2/3]))
        ax_config_plot_xyz(axNe3c[i], fnamecs[i], "gray", LineWidth, rotxyz=rots[i],mesh=1, bead=0,rod=0,d=0.8,pwlim=0)
        ax_config_plot_xyz(axNe3c[i], fnamecs[i], "gray", LineWidth, rotxyz=rots[i],yshift=-20, mesh=1, bead=1,rod=0,d=0.8,pwlim=np.pi/3)


    # Tc
    for Ne in [2,3]:
        q, uc_ave, uc_err = Ts_data_get(Ne)
        if Ne == 3:
            pass
            #q = np.delete(q,[3])
            #uc_ave = np.delete(uc_ave,[3])
            #uc_err = np.delete(uc_err,[3])
        axTs.errorbar(q[start:cut:skip],uc_ave[start:cut:skip],yerr=uc_err[start:cut:skip],linestyle="-",marker="o",mfc="None",ms=msize,label=r"$N_e=%d$"%Ne)
    axTs.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelleft=True,labelsize=LabelSize)
    axTs.set_ylabel(r"$T_s$", fontsize=FontSize)
    axTs.set_xlabel(r"$k_c$", fontsize=FontSize)
    axTs.xaxis.set_major_locator(MultipleLocator(0.5))
    axTs.xaxis.set_minor_locator(MultipleLocator(0.25))
    axTs.yaxis.set_major_locator(MultipleLocator(0.1))
    axTs.yaxis.set_minor_locator(MultipleLocator(0.05))
    axTs.legend(ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)
    x1, y1 = 0.9, 0.15
    axTs.text(x1,y1, r"(c)", fontsize=FontSize,transform=axTs.transAxes)


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
    plt.tight_layout(pad=0.1)
    plt.savefig("figures/topo_chirality.pdf",format="pdf")



def dis2H_data_get(Ne=2,shape="open edge"):
    sprow = 200
    id_select = topo_id_select(Ne)
    #foldername = "../data/Ne2/May21_2023"
    foldername = "../data/Ne2/May26_2023"
    if(shape=="disk"):
        fname = foldername+"/2Hdis_N300_imod3_Ne3_lf0.0_kar30_C00.0_karg0.0_lam6.0_Kd4.0_q0.0_Cn4.0_id0.csv"
        data = np.loadtxt(fname, skiprows=1+sprow, delimiter=",", unpack=True)
    elif(shape=="vesicle"):
        #foldername = "../data/Ne2/May14_2023"
        fname = foldername+"/2Hdis_N300_imod3_Ne3_lf0.0_kar30_C00.0_karg0.0_lam50.0_Kd4.0_q0.0_Cn4.0_id0.csv"
        data = np.loadtxt(fname, skiprows=1+sprow, delimiter=",", unpack=True)
    elif(Ne==2):
        q = 1.3
        #id = id_select[int(q/0.1)]
        fname = foldername+"/2Hdis_N300_imod3_Ne3_lf0.0_kar30_C00.0_karg0.0_lam6.0_Kd4.0_q%.1f_Cn4.0_id0.csv"%(q)
        data = np.loadtxt(fname, skiprows=1+sprow, delimiter=",", unpack=True)
    elif(Ne==3):
        #foldername = "../data/Ne2/May14_2023"
        q = 2.3
        #sprow = 150
        #id = id_select[int(q/0.1)]
        fname = foldername+"/2Hdis_N300_imod3_Ne3_lf0.0_kar30_C00.0_karg0.0_lam6.0_Kd4.0_q%.1f_Cn4.0_id0.csv"%(q)
        data = np.loadtxt(fname, skiprows=1+sprow, delimiter=",", unpack=True)
    bin_num = len(data)
    twoHpdf = np.average(data,axis=1)*bin_num
    twoHx = np.linspace(-1+1/bin_num,1-1/bin_num,bin_num)
    return(twoHx,twoHpdf)


def topo_Hdis_plot(LineWidth, FontSize, LabelSize):
    print("👌交给我吧")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.5))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    ax2H = plt.subplot2grid((1,1),(0,0))
    #axcfg = plt.subplot2grid((1,2),(0,1))

    msize = 4
    colors = ["blue", "orange", "purple", "red"]
    # 2H distribution
    twoHx, twoHpdf = dis2H_data_get(2,shape="disk")
    Hx, Hpdf = twoHx/2,  twoHpdf*2
    ax2H.plot(Hx,Hpdf,label="disk", color = colors[0])

    twoHx, twoHpdf = dis2H_data_get(2,shape="vesicle")
    Hx, Hpdf = twoHx/2,  twoHpdf*2
    ax2H.plot(Hx,Hpdf,label="vesicle", color = colors[1])


    twoHx, twoHpdf = dis2H_data_get(2,shape="open edge")
    Hx, Hpdf = twoHx/2,  twoHpdf*2
    ax2H.plot(Hx,Hpdf,label="catenoid", color = colors[2])

    twoHx, twoHpdf = dis2H_data_get(3,shape="open edge")
    Hx, Hpdf = twoHx/2,  twoHpdf*2
    ax2H.plot(Hx,Hpdf,label="trinoid", color = colors[3])

    ax2H.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelleft=True,labelsize=LabelSize)
    ax2H.set_xlim(-0.3,0.5)
    ax2H.set_xlabel(r"$H_i$")
    ax2H.set_ylabel(r"$P(H_i)$")
    ax2H.legend(ncol=1,columnspacing=0.5,handlelength=0.5,handletextpad=0.5,frameon=False,fontsize=FontSize)
    x1, y1 = 0.9, 0.1
    #ax2H.text(x1,y1, r"(d)", fontsize=FontSize,transform=ax2H.transAxes)


    plt.tight_layout(pad=0.1)
    plt.savefig("figures/topo_Hdis.pdf",format="pdf")