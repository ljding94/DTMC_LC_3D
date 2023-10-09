import numpy as np
import matplotlib.pyplot as plt
from config_plot import *
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)


def tilt_Kds_Cn_data_get():
    #foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    #foldername = "../data/Ne2/Apr_2022/Apr9_2022"
    foldername = "../data/Ne2/data_2022/May12_2022"
    lf = 25
    Cns = [2,4,6,8]
    datas, labels, colors, markers = [], [], [], []
    #colors = ["red", "purple", "blue", "royalblue"]
    colors = ["blue", "orange", "purple", "red"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(Cns)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lf%.1f_kar50_C00.0_karg0.0_lam6.0_Kds_q0.0_Cn%.1f_ana.csv" % (lf,Cns[i])
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))

    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    #Kds, un2_aves, un2_errs = datas[0], datas[31], datas[33] # test using p2uu
    Kds, p2uu_aves, p2uu_errs = datas[0], datas[25], datas[27] # test using p2uu
    # normalize Kd by Cn
    KdCns=[]
    for i in range(len(Cns)):
        KdCns.append(Kds[i]/Cns[i])

    labels = list(map(str, Cns))
    legendtitle=r"$l_f=%.0f, C$"%lf

    #return [Kds, un2_aves,un2_errs,labels, colors, markers,legendtitle]
    return [Kds, p2uu_aves,p2uu_errs,labels, colors, markers,legendtitle] # trying p2uu


def tilt_Kds_config_data_get():
    #foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    #foldername = "../data/Ne2/Apr_2022/Apr9_2022"
    foldername = "../data/Ne2/data_2022/May12_2022"
    lf=25.0
    Cn = 6.0
    Kds = [4.0,2.0]
    Kdp = []
    fnames,povs,rotxyzs, xyshift, zslice = [],[],[],[],[]
    for i in range(len(Kds)):
        for j in range(3):
            fnames.append( foldername + "/State_N300_imod3_Ne2_lf%.1f_kar50_C00.0_karg0.0_lam6.0_Kd%.1f_q0.0_Cn%.1f.csv" % (lf,Kds[i],Cn))
            Kdp.append(Kds[i])
            povs.append("xy")
            rotxyzs.append([0,0,0])
            xyshift.append((10*j,10*i))
        # 2 crossection each
        zslice+=[(7,13),(-2.5,2.5),None]
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
    #foldername = "../data/Ne2/Oct_2021/Oct18_2021"
    #foldername = "../data/Ne2/Apr_2022/Apr9_2022"
    foldername = "../data/Ne2/data_2022/May12_2022"
    Cn=4.0
    lfs = [15,25,35]
    datas, labels, colors, markers = [], [], [], []
    #colors = ["red", "green", "blue", "royalblue"]
    colors = ["blue", "orange", "purple", "red"]
    markers = ["v", "s", "p", "h", "o"]
    for i in range(len(lfs)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lf%.1f_kar50_C00.0_karg0.0_lam6.0_Kds_q0.0_Cn%.1f_ana.csv" % (lfs[i],Cn)
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))

    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    Kds, un2_aves, un2_errs = datas[0], datas[31], datas[33]
    Kds, p2uu_aves, p2uu_errs = datas[0], datas[25], datas[27] # test using p2uu
    # normalize Kd by lf
    Kdlfs=[]
    for i in range(len(lfs)):
        Kdlfs.append(Kds[i]/lfs[i])

    labels = list(map(str, lfs))
    legendtitle=r"$C=%.0f,l_f$"%Cn

    #return [Kds, un2_aves,un2_errs,labels, colors, markers,legendtitle]
    return [Kds, p2uu_aves,p2uu_errs,labels, colors, markers,legendtitle] # trying p2uu


def nematic_Kd_plot(LineWidth, FontSize, LabelSize):
    print("👌交给我吧")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axcfg = plt.subplot2grid((2, 2), (0, 0),colspan=2)
    axCns = plt.subplot2grid((2, 2), (1, 0),rowspan=1)
    axlfs = plt.subplot2grid((2, 2), (1, 1),rowspan=1,sharey=axCns,sharex=axCns)

    msize = 4

    ## configuration shows the tilt wall formation
    Kds, fnames, povs, rotxyzs, xyshifts, zslices = tilt_Kds_config_data_get()
    print(xyshifts)
    for i in range(len(Kds)):
        ax_config_plot_xyz(axcfg, fnames[i], "gray", LineWidth, pov=povs[i], rotxyz=rotxyzs[i],xshift=xyshifts[i][0],yshift=xyshifts[i][1], zslice=zslices[i], mesh=1, bead=1,rod=1,d=1)
        if(i%3==0):
            #axcfg.text(xyshifts[i][0]-10,xyshifts[i][1]-5.5,r"$\epsilon_{LL}=%.1f$"%Kds[i],fontsize=FontSize)
            print("Kds[i]",Kds[i])
            axcfg.text(xyshifts[i][0]-10,xyshifts[i][1]-5.5,r"$\epsilon_{LL}=%.0f$"%Kds[i],fontsize=FontSize)
    axcfg.text(-10,11,r"$C=6$",fontsize=FontSize)
    axcfg.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)
    x1, y1 = 0.85, 1.2
    axcfg.text(x1,y1, r"(a)", fontsize=FontSize,transform=axlfs.transAxes)


    ni = 0 # initial point
    n = 2 # date inteval
    ## tilt drops as Kd increases
    Kds, un2_aves,un2_errs,labels, colors, markers,legendtitle = tilt_Kds_Cn_data_get()
    for i in range(len(Kds)):
        axCns.errorbar(Kds[i][ni::n],un2_aves[i][ni::n],un2_errs[i][ni::n], ls=":", color=colors[i],mfc="None",marker=markers[i],ms=msize,label=labels[i])
    axCns.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelleft=True,labelsize=LabelSize)
    axCns.set_ylabel(r"$(\vu{u}\cdot\vu{n})^2$", fontsize=FontSize)
    axCns.set_ylim(0.4,1.0)
    axCns.xaxis.set_major_locator(MultipleLocator(1))
    axCns.xaxis.set_minor_locator(MultipleLocator(0.5))
    axCns.yaxis.set_major_locator(MultipleLocator(0.1))
    axCns.yaxis.set_minor_locator(MultipleLocator(0.05))
    axCns.set_xlabel(r"$\epsilon_{LL}$",fontsize=FontSize)
    axCns.legend(title=legendtitle,ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)
    x1, y1 = 0.85, 0.1
    axCns.text(x1,y1, r"(b)", fontsize=FontSize,transform=axCns.transAxes)

    ## critical Kd depends on the edge seperation distance
    Kds, un2_aves,un2_errs,labels, colors, markers,legendtitle = tilt_Kds_lf_data_get()
    for i in range(len(Kds)):
        axlfs.errorbar(Kds[i][ni::n],un2_aves[i][ni::n],un2_errs[i][ni::n], ls=":", color=colors[i],mfc="None",marker=markers[i],ms=msize,label=labels[i])
    axlfs.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelleft=False,labelsize=LabelSize)
    #axlfs.set_ylabel(r"$(\vu{u}\cdot\vu{n})^2$", fontsize=FontSize)
    #axlfs.set_ylim(0.35,1.0)
    axlfs.set_xlabel(r"$\epsilon_{LL}$",fontsize=FontSize)
    axlfs.legend(title=legendtitle,loc="upper right",ncol=1,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)
    x1, y1 = 0.85, 0.1
    axlfs.text(x1,y1, r"(c)", fontsize=FontSize,transform=axlfs.transAxes)

    plt.tight_layout(pad=0.1)
    plt.savefig("figures/nematic_wall.pdf",format="pdf")


def smectic_to_walls_config_data_get():
    lf=25.0
    Cn = 6.0
    fns = ["../data/Ne2/data_2022/May12_2022","../data/Ne2/data_2022/May12_2022","../data/Ne2/data_2022/Nov30_2022"]
    Kds = [2.0,4.0,2.0]
    qs = [0.0,0.0,1.5]
    Kps = []
    qps = []
    fnames,povs,rotxyzs, xyshift, zslice, hides = [],[],[],[],[],[]
    rotate_per_j = [[0,0,0],[-np.pi/2,0,0],[-np.pi*3/8,np.pi*3/8,0]]
    for i in range(len(Kds)):
        for j in range(3):
            fnames.append( fns[i] + "/State_N300_imod3_Ne2_lf%.1f_kar50_C00.0_karg0.0_lam6.0_Kd%.1f_q%.1f_Cn%.1f.csv" % (lf,Kds[i],qs[i],Cn))
            povs.append("xy")
            rotxyzs.append(rotate_per_j[j])
            xyshift.append((10*j,-10*i))
            Kps.append(Kds[i])
            qps.append(qs[i])
        # 2 crossection each
        #zslice+=[(7,13),(-2.5,2.5),None]
        zslice+=[(10,13),(7,13),None] # usd same slice for the two crosssection shown
        hides+=[None,0.5,None]
        xyshift[3*i+2]=(25,-10*i)
        povs[3*i+2]="zx"
        #rotxyzs[3*i+2]=[-np.pi*3/8,np.pi*3/8,0]

    return [Kps,qps,fnames,povs,rotxyzs, xyshift, zslice, hides]

from cholesteric_wall import *
# merge nematic wall and cholesteric wall plot
def walls_Cn_lf_vs_Kd_q(LineWidth, FontSize, LabelSize):
    print("👌交给我吧")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1.6))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axcfg = plt.subplot2grid((7, 2), (0, 0),colspan=2,rowspan=3)
    #smectic to nematic
    axCn_Kd = plt.subplot2grid((7, 2), (3, 0),rowspan=2)
    axlf_Kd = plt.subplot2grid((7, 2), (5, 0),rowspan=2,sharey=axCn_Kd)
    #smectic to cholesteric
    axCn_q = plt.subplot2grid((7, 2), (3, 1),rowspan=2,sharey=axCn_Kd)
    axlf_q = plt.subplot2grid((7, 2), (5, 1),rowspan=2,sharey=axlf_Kd,sharex=axCn_q)

    msize = 4

    ## configuration shows the tilt wall formation
    Kds, qs, fnames, povs, rotxyzs, xyshifts, zslices, hides = smectic_to_walls_config_data_get()
    print(xyshifts)

    for i in range(len(Kds)):
        ax_config_plot_xyz(axcfg, fnames[i], "gray", LineWidth, pov=povs[i], rotxyz=rotxyzs[i],xshift=xyshifts[i][0],yshift=xyshifts[i][1], zslice=zslices[i], hide_behind=hides[i], mesh=1, bead=1,rod=1,d=1)
        if(i%3==0):
            #axcfg.text(xyshifts[i][0]-10,xyshifts[i][1]-5.5,r"$\epsilon_{LL}=%.1f$"%Kds[i],fontsize=FontSize)
            print("Kds[i]",Kds[i])
            axcfg.text(xyshifts[i][0]-10,xyshifts[i][1]-5.5,r"$(\epsilon_{LL},k_c)=(%.0f,%.1f)$"%(Kds[i],qs[i]),fontsize=FontSize)

    axcfg.tick_params(which="both",direction="in", bottom="off",top="off", right="off",left="off",labelbottom=False,labelleft=False, labelsize=LabelSize)
    x1, y1 = 0.85, 1.2
    axcfg.text(x1,y1, r"(a)", fontsize=FontSize,transform=axCn_q.transAxes)


    ni = 0 # initial point
    n = 2 # date inteval
    # smectic to nematic
    ## tilt drops as Kd increases
    Kds, un2_aves,un2_errs,labels, colors, markers,legendtitle = tilt_Kds_Cn_data_get()
    for i in range(len(Kds)):
        axCn_Kd.errorbar(Kds[i][ni::n],un2_aves[i][ni::n],un2_errs[i][ni::n], ls=":", color=colors[i],mfc="None",marker=markers[i],ms=msize,label=labels[i])
    axCn_Kd.tick_params(which="both",direction="in", top="on", right="on",labelbottom=False, labelleft=True,labelsize=LabelSize)
    #axCn_Kd.set_ylabel(r"$\left<(\vu{u}\cdot\vu{n})^2\right>$", fontsize=FontSize)
    axCn_Kd.set_ylabel(r"$\left<\frac{3}{2}(\vu{u}\cdot\vu{u})-\frac{1}{2}\right>$", fontsize=FontSize)
    axCn_Kd.set_ylim(0.42,0.98)
    axCn_Kd.xaxis.set_major_locator(MultipleLocator(1))
    axCn_Kd.xaxis.set_minor_locator(MultipleLocator(0.5))
    axCn_Kd.yaxis.set_major_locator(MultipleLocator(0.1))
    axCn_Kd.yaxis.set_minor_locator(MultipleLocator(0.05))
    #axCn_Kd.set_xlabel(r"$\epsilon_{LL}$",fontsize=FontSize)
    axCn_Kd.legend(title=legendtitle,ncol=2,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)
    x1, y1 = 0.85, 0.1
    axCn_Kd.text(x1,y1, r"(b)", fontsize=FontSize,transform=axCn_Kd.transAxes)

    ## critical Kd depends on the edge seperation distance
    Kds, un2_aves,un2_errs,labels, colors, markers,legendtitle = tilt_Kds_lf_data_get()
    for i in range(len(Kds)):
        axlf_Kd.errorbar(Kds[i][ni::n],un2_aves[i][ni::n],un2_errs[i][ni::n], ls=":", color=colors[i],mfc="None",marker=markers[i],ms=msize,label=labels[i])
    axlf_Kd.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelleft=True,labelsize=LabelSize)
    axlf_Kd.set_ylabel(r"$\left<(\vu{u}\cdot\vu{n})^2\right>$", fontsize=FontSize)
    axlf_Kd.set_ylabel(r"$\left<\frac{3}{2}(\vu{u}\cdot\vu{u})-\frac{1}{2}\right>$", fontsize=FontSize)
    #axlf_Kd.set_ylim(0.35,1.0)
    axlf_Kd.xaxis.set_major_locator(MultipleLocator(1))
    axlf_Kd.xaxis.set_minor_locator(MultipleLocator(0.5))
    axlf_Kd.yaxis.set_major_locator(MultipleLocator(0.1))
    axlf_Kd.yaxis.set_minor_locator(MultipleLocator(0.05))
    axlf_Kd.set_xlabel(r"$\epsilon_{LL}$",fontsize=FontSize)
    axlf_Kd.legend(title=legendtitle,loc="upper right",ncol=1,columnspacing=0.5,handlelength=0.5,handletextpad=0.1,frameon=False,fontsize=FontSize)
    x1, y1 = 0.85, 0.1
    axlf_Kd.text(x1,y1, r"(c)", fontsize=FontSize,transform=axlf_Kd.transAxes)

    # smectic to cholesteric
    qs, un2_aves, un2_errs, labels, colors, markers, legendtitle = tilt_qs_Cn_data_get()
    for i in range(len(qs)):
        axCn_q.errorbar(qs[i][ni::n], un2_aves[i][ni::n], un2_errs[i][ni::n], ls=":", color=colors[i], mfc="None", marker=markers[i], ms=msize, label=labels[i])
    axCn_q.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=False, labelsize=LabelSize)
    #axCn_q.set_ylabel(r"$(\vu{u}\cdot\vu{n})^2$", fontsize=FontSize)
    axCn_q.set_ylim(0.42, 0.98)
    axCn_q.xaxis.set_major_locator(MultipleLocator(0.5))
    axCn_q.xaxis.set_minor_locator(MultipleLocator(0.25))
    #axCn_q.yaxis.set_major_locator(MultipleLocator(0.1))
    #axCn_q.yaxis.set_minor_locator(MultipleLocator(0.05))
    #axCn_q.set_xlabel(r"$q$", fontsize=FontSize)
    axCn_q.legend(title=legendtitle, ncol=3, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)
    x1, y1 = 0.85, 0.1
    axCn_q.text(x1, y1, r"(d)", fontsize=FontSize, transform=axCn_q.transAxes)

    ## critical q depends on the edge seperation distance
    qs, un2_aves, un2_errs, labels, colors, markers, legendtitle = tilt_qs_lf_data_get()
    for i in range(len(qs)):
        axlf_q.errorbar(qs[i][ni::n], un2_aves[i][ni::n], un2_errs[i][ni::n], ls=":", color=colors[i], mfc="None", marker=markers[i], ms=msize, label=labels[i])
    axlf_q.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=False, labelsize=LabelSize)
    # axlf_q.set_ylabel(r"$(\vu{u}\cdot\vu{n})^2$", fontsize=FontSize)
    # axlf_q.set_ylim(0.35,1.0)
    axlf_q.xaxis.set_major_locator(MultipleLocator(0.5))
    axlf_q.xaxis.set_minor_locator(MultipleLocator(0.25))
    axlf_q.set_xlabel(r"$k_c$", fontsize=FontSize)
    axlf_q.legend(title=legendtitle, loc="upper right", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=FontSize)
    x1, y1 = 0.85, 0.1
    axlf_q.text(x1, y1, r"(e)", fontsize=FontSize, transform=axlf_q.transAxes)

    plt.tight_layout(pad=0.1)
    plt.savefig("figures/smectic_to_walls.pdf",format="pdf")


def walls_membrane_shape_slice(LineWidth, FontSize, LabelSize):
    print("see if (cyindrical) cross section shape twist with pi wall")
    ppi = 72
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.4))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axcfg = plt.subplot2grid((1, 1), (0, 0),colspan=1,rowspan=1)
    fname = "../data/Ne2/data_2022/Nov30_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd2.0_q1.5_Cn6.0.csv"
    zs = np.arange(-14,14,4)
    for i in range(len(zs)-1):
        ax_config_plot_xyz(axcfg, fname, "gray", LineWidth, pov="xy",xshift = 8*i,zslice = [zs[i],zs[i+1]], mesh=1,bead=1,rod=1,d=1)
    plt.show()




def walls_3phase_diagram_data_get():
    foldername = "../data/Ne2/Aug8_2023"
    lf=15.0
    Kd=2.0
    Cns = np.arange(1.0,13.1,1.0)
    #Cns = np.arange(15.0,34.1,1.0)

    datas = []
    for i in range(len(Cns)):
        fname = foldername + "/O_MC_N300_imod3_Ne2_lf%.1f_kar50_C00.0_karg0.0_lam6.0_Kd2.0_qs_Cn%.1f_id0_ana.csv" % (lf,Cns[i])
        #datas.append(np.loadtxt(fname, skiprows=1,max_rows=15, delimiter=",", unpack=True))
        datas.append(np.loadtxt(fname, skiprows=1, delimiter=",", unpack=True))
    datas = np.transpose(np.array(datas), axes=(1, 0, 2))
    qs, uc_aves, uc_errs, un_aves, un_errs,uz2_aves, uz2_errs = datas[0], datas[28], datas[30], datas[31], datas[33], datas[34], datas[36]

    return [qs[0],Cns,uc_aves]


def walls_3phase_diagram_config_data_get():
    foldername="../data/Ne2/Aug8_2023"
    qCns = [(0.5,7.0),(1.0,2.0),(2.2,4.0)]
    ypos = 14.1
    arrowpos=[(0.6,ypos),(1.8,ypos),(2.8,ypos)]
    fnames = []
    povs = []
    rotxyzs = [(0,0.7,0),(-0.3,0.7,0),(0.2,0.7,0)]
    for i in range(len(qCns)):
        fnames.append(foldername + "/State_N300_imod3_Ne2_lf15.0_kar50_C00.0_karg0.0_lam6.0_Kd2.0_q%.1f_Cn%.1f_id0.csv" % (qCns[i]))
        povs.append("zx")
        #rotxyzs.append([0,0.7,0])

    return [qCns,arrowpos,fnames,povs,rotxyzs]


def walls_3phase_diagram(LineWidth, FontSize, LabelSize):
    print("PRINT 3-PHASE phase diagram plot for m=0,2,3")

    qs,Cns,uc_aves = walls_3phase_diagram_data_get()
    qCns,arrowpos,fnames,povs,rotxyzs = walls_3phase_diagram_config_data_get()

    ppi = 72
    #fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1.05))
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.8))
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{physics}')

    ax = plt.subplot2grid((1,1), (0, 0))
    #axcfg = ax.inset_axes([-0.1, 0.98, 1.3, 0.4])
    axcfg = ax.inset_axes([-0.1, 0.98, 1.3, 0.6])

    qp,Cnp = np.meshgrid(qs,Cns)
    im = ax.pcolormesh(qp,Cnp,uc_aves,shading="auto",cmap=cm.get_cmap("rainbow"))
    cbar=fig.colorbar(im,ax=ax)
    cbar.ax.tick_params(direction="in",labelsize=LabelSize)
    #cbar.set_label(r"$T_s$") # this set on the right side
    cbar.ax.set_title(r"$T_s$", fontsize=FontSize) # this set on the top
    ax.tick_params(which="both",direction="in", top="on", right="on",labelbottom=True, labelleft=True,labelsize=LabelSize)
    ax.set_xlabel(r"$k_c$", fontsize=FontSize)
    ax.set_ylabel(r"$C$", fontsize=FontSize)

    anotext = ["smectic-A", "2-wall", "3-wall"]

    for i in range(len(fnames)):
        ax_config_plot_xyz(axcfg, fnames[i], "gray", LineWidth, pov=povs[i], rotxyz=rotxyzs[i],xshift=i*20,yshift=0, mesh=1, bead=1,rod=1,d=1)
        ax.text(arrowpos[i][0],arrowpos[i][1]+0.2,anotext[i],ha="center",fontsize=LabelSize)
        ax.annotate("",xy=qCns[i], xytext=arrowpos[i],arrowprops=dict(arrowstyle="->",lw=LineWidth,color="gray"),annotation_clip=False)
        #axconf.text(i/rn+0.1,0, r"$%.0f$"%lams[i], fontsize=FontSize)

    x1, y1 = -0.035, 0.8
    #axcfg.text(axcfg.get_xlim()[1]*x1+axcfg.get_xlim()[0]* (1-x1),axcfg.get_ylim()[1]*y1+axcfg.get_ylim()[0]* (1-y1), r"(a)", fontsize=FontSize)

    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_minor_locator(MultipleLocator(1))

    plt.tight_layout(0.05)
    #plt.subplots_adjust(wspace=0.05,hspace=0.3,left=0.12,right=0.97,bottom=0.08,top=1.37)
    plt.savefig("figures/phase_diagram_walls.pdf",format="pdf")