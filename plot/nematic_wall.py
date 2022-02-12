import numpy as np
import matplotlib.pyplot as plt
from config_plot import *



def tilt_Kds_Cn_data_get():
    foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    lf = 20.0
    Cns = [1.0,3.0,5.0,7.0]
    datas, labels, colors, markers = [], [], [], []
    colors = ["red", "green", "blue", "royalblue"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(Cns)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lf%.1f_kar40_C00.0_karg0.0_lam5.0_Kds_q0.0_Cn%.1f_ana.csv" % (lf,Cns[i])
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))

    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    Kds, un2_aves, un2_errs = datas[0], datas[31], datas[33]
    # normalize Kd by Cn
    KdCns=[]
    for i in range(len(Cns)):
        KdCns.append(Kds[i]/Cns[i])

    labels = list(map(str, Cns))
    legendtitle=r"$l_f=%.0f, C_n$"%lf

    return [Kds, un2_aves,un2_errs,labels, colors, markers,legendtitle]


def tilt_Kds_config_data_get():
    foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    lf=20.0
    Cn = 5.0
    Kds = [4.2,1.8]
    Kdp = []
    fnames,povs,rotxyzs, xyshift, zslice = [],[],[],[],[]
    for i in range(len(Kds)):
        for j in range(3):
            fnames.append( foldername + "/State_N300_imod3_Ne2_lf%.1f_kar40_C00.0_karg0.0_lam5.0_Kd%.1f_q0.0_Cn%.1f.csv" % (lf,Kds[i],Cn))
            Kdp.append(Kds[i])
            povs.append("xy")
            rotxyzs.append(None)
            xyshift.append((10*j,10*i))
        # 2 crossection each
        zslice+=[(7,12),(-2.5,2.5),None]
        xyshift[3*i+2]=(25,10*i)
        povs[3*i+2]="zx"
        rotxyzs[3*i+2]=[-np.pi*3/8,np.pi*3/8,0]
    #fnames.append(fnames[-1])
    #povs.append("zx")
    #Kdp.append(Kds[-1])
    #rotxyzs.append([-np.pi*3/8,np.pi*3/8,0])
    #xyshift.append((10*2,0))
    #zslice.append(None)
    return [Kdp,fnames,povs,rotxyzs, xyshift, zslice]

def tilt_Kds_lf_data_get():
    foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    Cn=5.0
    lfs = [10.0,20.0,30.0]
    lf0 = 20.0
    datas, labels, colors, markers = [], [], [], []
    colors = ["red", "green", "blue", "royalblue"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(lfs)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lf%.1f_kar40_C00.0_karg0.0_lam5.0_Kds_q0.0_Cn%.1f_ana.csv" % (lfs[i],Cn)
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))

    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    Kds, un2_aves, un2_errs = datas[0], datas[31], datas[33]
    # normalize Kd by lf
    Kdlfs=[]
    for i in range(len(lfs)):
        Kdlfs.append(Kds[i]/lfs[i])

    labels = list(map(str, lfs))
    legendtitle=r"$C=%.0f,l_f$"%Cn

    return [Kds, un2_aves,un2_errs,labels, colors, markers,legendtitle]


def tilt_Kd_plot(LineWidth, FontSize, LabelSize):
    print("👌交给我吧")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1.4))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axcfg = plt.subplot2grid((3, 1), (0, 0))
    axCns = plt.subplot2grid((3, 1), (1, 0),rowspan=1)
    axlfs = plt.subplot2grid((3, 1), (2, 0),rowspan=1,sharex=axCns)

    msize = 3

    ## configuration shows the tilt wall formation
    Kds, fnames, povs, rotxyzs, xyshifts, zslices = tilt_Kds_config_data_get()
    print(xyshifts)
    for i in range(len(Kds)):
        ax_config_plot_xyz(axcfg, fnames[i], "gray", LineWidth, pov=povs[i], rotxyz=rotxyzs[i],xshift=xyshifts[i][0],yshift=xyshifts[i][1], zslice=zslices[i], mesh=1, bead=0,rod=1,d=1)
        if(i%3==0):
            #axcfg.text(xyshifts[i][0]-10,xyshifts[i][1]-5.5,r"$\epsilon_{LL}=%.1f$"%Kds[i],fontsize=FontSize)
            axcfg.text(xyshifts[i][0]-10,xyshifts[i][1]-5.5,r"$K=%.1f$"%Kds[i],fontsize=FontSize)
    axcfg.text(-10,11,r"$C=5$",fontsize=FontSize)
    axcfg.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)

    ## tilt drops as Kd increases
    Kds, un2_aves,un2_errs,labels, colors, markers,legendtitle = tilt_Kds_Cn_data_get()
    for i in range(len(Kds)):
        axCns.errorbar(Kds[i],un2_aves[i],un2_errs[i], ls="None", color=colors[i],mfc="None",marker=markers[i],ms=msize,label=labels[i])
    axCns.tick_params(which="both",direction="in", top="on", right="on",labelbottom=False, labelleft=True,labelsize=LabelSize)
    axCns.set_ylabel(r"$(\vu{u}\cdot\vu{n})^2$", fontsize=FontSize)
    axCns.set_ylim(0.4,1.0)
    axCns.legend(title=legendtitle,ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)

    ## critical Kd depends on the edge seperation distance
    Kds, un2_aves,un2_errs,labels, colors, markers,legendtitle = tilt_Kds_lf_data_get()
    for i in range(len(Kds)):
        axlfs.errorbar(Kds[i],un2_aves[i],un2_errs[i], ls="None", color=colors[i],mfc="None",marker=markers[i],ms=msize,label=labels[i])
    axlfs.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelleft=True,labelsize=LabelSize)
    axlfs.set_ylabel(r"$(\vu{u}\cdot\vu{n})^2$", fontsize=FontSize)
    axlfs.set_ylim(0.4,1.0)
    axlfs.set_xlabel(r"$\epsilon_{LL}$",fontsize=FontSize)
    axlfs.legend(title=legendtitle,loc="upper right",ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)

    plt.tight_layout(pad=0.1)
    plt.savefig("figures/nematic_wall.pdf",format="pdf")





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
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1.4))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axcfg = plt.subplot2grid((3, 1), (0, 0))
    axuc = plt.subplot2grid((3, 1), (1, 0),rowspan=1)
    axuz = plt.subplot2grid((3, 1), (2, 0),rowspan=1,sharex=axuc)
    msize=3



    # configuration plot
    qs,fnames,povs,rotxyzs = twist_q_config_data_get()
    for i in range(len(qs)):
        ax_config_plot_xyz(axcfg, fnames[i], "gray", LineWidth, pov=povs[i], rotxyz=rotxyzs[i],xshift=14*i, mesh=1, bead=0,rod=1, d=1)
        #axcfg.text(14*i-3,-14,r"$k_c=%.1f$"%qs[i],fontsize=FontSize)
        axcfg.text(14*i-3,-14,r"$q=%.1f$"%qs[i],fontsize=FontSize)
    axcfg.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)

    # o vs q plot

    qs,uc_aves, uc_errs, uz2_aves,uz2_errs,labels,colors,markers,legendtitle =twist_q_data_get()

    for i in range(len(qs)):
        axuc.errorbar(qs[i],uc_aves[i],uc_errs[i], ls="None", color=colors[i],mfc="None",marker=markers[i],ms=msize,label=labels[i])

    axuc.tick_params(which="both",direction="in", top="on", right="on",labelbottom=False, labelleft=True,labelsize=LabelSize)
    axuc.set_ylabel(r"$\left<(\vu{u}_i\cross\vu{u}_j)\cdot \vu{r}_{ij} (\vu{u}_i\cdot\vu{u}_j)\right>_{(i,j)}$", fontsize=FontSize)
    axuc.set_ylim(0.0,0.22)
    axuc.legend(title=legendtitle,ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)


    for i in range(len(qs)):
        axuz.errorbar(qs[i],uz2_aves[i],uz2_errs[i], ls="None", color=colors[i],mfc="None",marker=markers[i],ms=msize,label=labels[i])

    axuz.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelleft=True,labelsize=LabelSize)
    axuz.set_ylabel(r"$(\vu{u}\cdot \vu{z})^2$", fontsize=FontSize)
    axuz.set_ylim(0.0,0.2)
    axuz.set_xlabel(r"$k_c$",fontsize=FontSize)
    axuz.legend(title=legendtitle,ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)

    plt.tight_layout(pad=0.1)
    plt.savefig("twisting_wall.pdf",format="pdf")