import numpy as np
import matplotlib.pyplot as plt
from config_plot import *
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

def twist_q_data_get():
    #foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    #foldername = "../data/Ne2/Mar_2022/Mar23_2022"
    foldername = "../data/Ne2/May12_2022"
    lfs=[15,25,35]
    Kd=4.0
    Cn=4.0
    datas, labels, colors, markers = [], [], [], []
    colors = ["red", "green", "blue", "royalblue","purple"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(lfs)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lf%.1f_kar50_C00.0_karg0.0_lam6.0_Kd%.1f_qs_Cn%.1f_ana.csv" % (lfs[i],Kd,Cn)
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))
    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    qs, uc_aves, uc_errs, uz2_aves,uz2_errs = datas[0], datas[28], datas[30],datas[34], datas[36]
    labels = list(map(str, lfs))
    legendtitle=r"$l_f$"
    return [qs,uc_aves, uc_errs, uz2_aves,uz2_errs,labels,colors,markers,legendtitle]


def twist_q_config_data_get():
    #foldername = "../data/Ne2/Mar_2022/Mar23_2022"
    foldername = "../data/Ne2/May12_2022"
    lf=25.0
    Kd=4.0
    Cn=4.0
    lfs=[25,25,25,35,35,35]
    qs = [0,1,2.5,0,1,2.5]
    fnames,povs,rotxyzs,xysfts=[],[],[],[]
    for i in range(len(qs)):
        fnames.append(foldername + "/State_N300_imod3_Ne2_lf%.1f_kar50_C00.0_karg0.0_lam6.0_Kd%.1f_q%.1f_Cn%.1f.csv" % (lfs[i],Kd,qs[i],Cn))
        povs.append("zx")
        #rotxyzs.append([np.pi/3,0,np.pi/2])
        rotxyzs.append([np.pi/3,0,0])
    rotxyzs[3]=[np.pi/2-0.1,0,0]
    xysfts = [[0,0],[0,13],[0,26],[35,0],[35,13],[35,26]]
    return [lfs,qs,fnames,povs,rotxyzs,xysfts]

def twist_q_plot(LineWidth, FontSize, LabelSize):
    print("👌")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 2, 246 / ppi * 1))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axcfg = plt.subplot2grid((2, 10), (0, 3),colspan=7,rowspan=2)
    #axcfgt = plt.subplot2grid((2, 5), (1, 3),colspan=2)
    axuc = plt.subplot2grid((2, 10), (0, 0),colspan=3)
    axuz = plt.subplot2grid((2, 10), (1, 0),colspan=3,sharex=axuc)
    msize=4


    '''
    # twist bond color
    for i in range(len(qs)):
        pass
        ax_config_plot_xyz(axcfgt, fnames[i], "gray", LineWidth, pov=povs[i], rotxyz=rotxyzs[i],xshift=0,yshift=15*i, mesh=1, bead=0,rod=0,twistcolor=1)
        axcfgt.text(-15, 15*i-7,r"$k_c=%.0f$"%qs[i],fontsize=FontSize)
    cbar=plt.colorbar(cm.ScalarMappable(norm=Normalize(vmin=0,vmax=0.25*np.pi), cmap=cm.get_cmap("ocean")),ax=axcfgt,ticks=[0,0.25*np.pi,0.5*np.pi],orientation="horizontal")
    cbar.ax.set_xticklabels([r"$0$",r"$\pi/4$",r"$\pi/2$"],fontsize=FontSize)
    cbar.ax.tick_params(direction="in",labelsize=LabelSize)
    cbar.ax.set_title(r"$T_s$",fontsize=FontSize)


    axcfgt.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)
    x1, y1 = 0.8, -0.0
    axcfgt.text(x1,y1, r"(d)", fontsize=FontSize,transform=axcfgt.transAxes)
    '''

    # o vs q plot

    qs,uc_aves, uc_errs, uz2_aves,uz2_errs,labels,colors,markers,legendtitle =twist_q_data_get()
    nf=30
    n = 2
    for i in range(len(qs)):
        axuc.errorbar(qs[i][:nf:n],uc_aves[i][:nf:n],uc_errs[i][:nf:n], ls="None", color=colors[i],mfc="None",marker=markers[i],ms=msize,label=labels[i])

    axuc.tick_params(which="both",direction="in", top="on", right="on",labelbottom=False, labelleft=True,labelsize=LabelSize)
    #axuc.set_ylabel(r"$\left<(\vu{u}_i\cross\vu{u}_j)\cdot \vu{r}_{ij} (\vu{u}_i\cdot\vu{u}_j)\right>_{(i,j)}$", fontsize=FontSize)
    axuc.set_ylabel(r"$T_s$", fontsize=FontSize)
    #axuc.set_ylim(0.0,0.22)
    axuc.xaxis.set_major_locator(MultipleLocator(1.0))
    axuc.xaxis.set_minor_locator(MultipleLocator(0.5))
    axuc.yaxis.set_major_locator(MultipleLocator(0.05))
    axuc.yaxis.set_minor_locator(MultipleLocator(0.025))
    #axuc.set_xlabel(r"$k_c$",fontsize=FontSize)
    axuc.legend(loc="upper left",title=legendtitle,ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)
    x1, y1 = 0.8, 0.1
    axuc.text(x1,y1, r"(a)", fontsize=FontSize,transform=axuc.transAxes)



    for i in range(len(qs)):
        axuz.errorbar(qs[i][:nf:n],uz2_aves[i][:nf:n],uz2_errs[i][:nf:n], ls="None", color=colors[i],mfc="None",marker=markers[i],ms=msize,label=labels[i])

    axuz.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelleft=True,labelsize=LabelSize)
    axuz.set_ylabel(r"$(\vu{u}\cdot \vu{z})^2$", fontsize=FontSize)
    #axuz.set_ylim(0.0,0.2)
    axuz.xaxis.set_major_locator(MultipleLocator(0.5))
    axuz.xaxis.set_minor_locator(MultipleLocator(0.25))
    axuz.set_xlim(-0.2,3.2)
    axuz.yaxis.set_major_locator(MultipleLocator(0.05))
    axuz.yaxis.set_minor_locator(MultipleLocator(0.025))
    axuz.set_xlabel(r"$k_c$",fontsize=FontSize)
    axuz.legend(loc="upper left",title=legendtitle,ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)
    x1, y1 = 0.8, 0.1
    axuz.text(x1,y1, r"(b)", fontsize=FontSize,transform=axuz.transAxes)


    # configuration plot
    lfs,qs,fnames,povs,rotxyzs,xysfts = twist_q_config_data_get()
    for i in range(len(qs)):
        pass
        print(xysfts[i])
        ax_config_plot_xyz(axcfg, fnames[i], "gray", LineWidth, pov=povs[i], rotxyz=rotxyzs[i],xshift=xysfts[i][0],yshift=xysfts[i][1], mesh=1, bead=0,rod=0, d=1,pwlim=np.pi/3)
        axcfg.text(xysfts[i][0]-lfs[i]/2, xysfts[i][1]-6,r"$k_c=%.1f,l_f=%.0f$"%(qs[i],lfs[i]),fontsize=FontSize)

        #ax_config_plot_xyz(axcfg, fnames[i], "gray", LineWidth, xshift=10,yshift=25*i-5,zslice=(8,12), mesh=1, bead=0,rod=0,pwlim=0.8, d=1)
    axcfg.text(xysfts[2][0]-5, xysfts[2][1]+lfs[i]/2-1,r"$k_c=%.1f$"%qs[-1],fontsize=FontSize)
    axcfg.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)
    x1, y1 = 3.1, 0.1
    axcfg.text(x1,y1, r"(c)", fontsize=FontSize,transform=axuz.transAxes)
    #axcfg.text(x1,y1, r"(c)", fontsize=FontSize,transform=axcfg.transAxes)

    #fig.tight_layout(pad=0.1)
    plt.tight_layout(pad=0.1)
    plt.savefig("figures/twisting_wall.pdf",format="pdf")

    # extra plot just to see the trend agains 1/2arctan(kc/3)
    '''
    plt.figure()
    for i in range(len(qs)):
        axuz.errorbar(np.arctan(qs[i][:nf:n]/3)/2,np.sqrt(uz2_aves[i][:nf:n]),uz2_errs[i][:nf:n], ls="None", color=colors[i],mfc="None",marker=markers[i],ms=msize,label=labels[i])
    plt.show()
    plt.close()
    '''