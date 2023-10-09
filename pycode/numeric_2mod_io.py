import numpy as np
from nematic_cholesteric_2mod_cylinder_cal import *
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
# single set of parameter
def Ftot_par_run(par, bn_phi, bn_z, method):
    K, C, q, m = par
    opt = opt_alpha_gamma_fun(obj_Ftot, K, C, q, m, bn_phi, bn_z, method)
    savename = "optFtot_K%.1f_C%.1f_m%0f_q%.1f.csv" % (K, C, m, q)
    with open(savename, "w") as f:
        f.write("optFtot,optalpha,optgamma\n")
        f.write("%f,%f,%f" % (opt.fun, opt.x[0], opt.x[1]))


# various q
def Ftot_qs_cal(folder, K, C, qs, m, bn_phi, bn_z, R, method):
    optFtots, optalphas, optgammas = [], [], []
    for i in range(len(qs)):
        q = qs[i]
        print(q, "/", qs[-1])
        opt = opt_alpha_gamma_fun(obj_Ftot, K, C, q, m, bn_phi, bn_z, R,method)
        optFtots.append(opt.fun)
        optalphas.append(opt.x[0])
        optgammas.append(opt.x[1])
    savename = folder + "/optFtot_K%.2f_C%.1f_m%.0f_R%.1f_qs.csv" % (K, C, m,R)
    with open(savename, "w") as f:
        f.write("q,optFtot,optalpha,optgamma\n")
        for i in range(len(qs)):
            f.write("%f,%f,%f,%f\n" % (qs[i],optFtots[i],optalphas[i], optgammas[i]))

def del_Ftot_Ks_qs_plot(folder,K,Cs,R,mi=0,mj=2):
    LineWidth, FontSize, LabelSize = 1,9,8
    optFtot_m0, optalpha_m0, optgamma_m0 = [], [], []
    optFtot_m2, optalpha_m2, optgamma_m2 = [], [], []
    for C in Cs:
        filename_m0 = folder + "/optFtot_K%.2f_C%.1f_m%d_R%.1f_qs.csv" % (K, C,mi, R)
        filename_m2 = folder + "/optFtot_K%.2f_C%.1f_m%d_R%.1f_qs.csv" % (K, C,mj, R)

        qs,optFtot,optalpha,optgamma = np.loadtxt(filename_m0, skiprows=1, delimiter=",", unpack=True)
        optFtot_m0.append(optFtot)
        optalpha_m0.append(optalpha)
        optgamma_m0.append(optgamma)

        qs,optFtot,optalpha,optgamma = np.loadtxt(filename_m2, skiprows=1, delimiter=",", unpack=True)
        optFtot_m2.append(optFtot)
        optalpha_m2.append(optalpha)
        optgamma_m2.append(optgamma)

    #qmesh,Kmesh = np.meshgrid(qs,Ks)
    qmesh,Kmesh = np.meshgrid(qs,Cs) # testing CR^2/K
    optFtot_diff = np.array(optFtot_m2) - np.array(optFtot_m0)
    absFdiffmax = np.max(np.abs(optFtot_diff))
    print("absFdiffmax=",absFdiffmax)
    ppi = 72
    plt.figure()
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.7))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axFdiff = plt.subplot2grid((2, 3), (0, 0),rowspan=2, colspan=2)
    axalpha = plt.subplot2grid((2, 3), (0, 2),rowspan=1)
    axgamma = plt.subplot2grid((2, 3), (1, 2),rowspan=1,sharex=axalpha)
    msize = 4

    cmap = cm.get_cmap("bwr")
    #cf = axs[0].pcolormesh(qmesh,Kmesh,optFtot_diff, shading="nearest", cmap = cmap,
                            #norm = colors.DivergingNorm(vmin=np.min(optFtot_diff), vcenter = 0,vmax = np.max(optFtot_diff) ))
    Fdiffmin,Fdiffmax = np.min(optFtot_diff),np.max(optFtot_diff)
    print(Fdiffmin,Fdiffmax)
    #Fdiffmin,Fdiffmax = round(Fdiffmin,0),round(Fdiffmax,2)
    levels = [Fdiffmin,Fdiffmin*3/4,Fdiffmin*2/4,Fdiffmin/4,0,Fdiffmax/4,Fdiffmax*2/4,Fdiffmax*3/4,Fdiffmax]
    print("qmesh.shape(),Kmesh.shape,optFtot.shape()",np.shape(qmesh),np.shape(Kmesh),np.shape(optFtot_diff))
    #cf = axFdiff.pcolormesh(qs,Cs,optFtot_diff,shading="gouraud",cmap=cmap)
    #cf = axFdiff.pcolormesh(qmesh,Kmesh,optFtot_diff,shading="gouraud",cmap=cmap, norm=colors.TwoSlopeNorm(vmin=Fdiffmin,vcenter=0,vmax=Fdiffmax))
    cf = axFdiff.pcolormesh(qmesh,Kmesh,optFtot_diff, shading="nearest",cmap=cmap, norm=colors.TwoSlopeNorm(vmin=Fdiffmin,vcenter=0,vmax=Fdiffmax))

    cfc = axFdiff.contourf(qmesh,Kmesh,optFtot_diff,levels,cmap=cmap,alpha=0, norm=colors.TwoSlopeNorm(vmin=Fdiffmin,vcenter=0,vmax=Fdiffmax))
    cs = axFdiff.contour(cfc,levels=[0],colors=("k"),linestyles="--",linewidths=(LineWidth))

    divider = make_axes_locatable(axFdiff)
    cax = divider.append_axes("top", size="3%", pad=0.02)
    #cbar=plt.colorbar(cf, cax = cax,ax=axFdiff,orientation="vertical",ticks=levels)
    cbar=plt.colorbar(cf, cax = cax,ax=axFdiff,orientation="horizontal",ticks=levels)
    #cbar.add_lines(cs)
    #cbar.ax.set_yticklabels([str(round(Fdiffmin,1)),"",str(round(Fdiffmin*2/4,1)),"","0","",str(round(Fdiffmax*2/4,2)),"",str(round(Fdiffmax,2))],fontsize=LabelSize)
    cbar.ax.set_xticklabels([str(round(Fdiffmin,1)),"",str(round(Fdiffmin*2/4,1)),"","0","",str(round(Fdiffmax*2/4,1)),"",str(round(Fdiffmax,1))],fontsize=LabelSize)
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.tick_params(direction="in")

    #cbar.ax.set_yticklabels([str(round(Fdiffmin,1)),"",str(round(Fdiffmin*2/4,1)),"","0","",str(round(Fdiffmax*2/4,2)),"",str(round(Fdiffmax,2))],fontsize=LabelSize)
    axFdiff.set_xlabel(r"$qR$",fontsize=LabelSize)
    axFdiff.set_ylabel(r"$CR^2/K$",fontsize=LabelSize)
    axFdiff.tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=True, labelleft=True, labelsize=LabelSize)
    cbar.ax.set_title(r"$\Delta E'R^2/K$",fontsize=LabelSize)
    axFdiff.xaxis.set_major_locator(MultipleLocator(1))
    axFdiff.xaxis.set_minor_locator(MultipleLocator(0.2))
    axFdiff.yaxis.set_major_locator(MultipleLocator(1))
    axFdiff.yaxis.set_minor_locator(MultipleLocator(0.5))

    if(mi==0 and mj==2):
        axFdiff.text(0.5,6,r"$m=%d$"%mi)
        axFdiff.text(1.5,2,r"$m=%d$"%mj)

    x1, y1 = 0.8, 0.1
    axFdiff.text(x1, y1, r"(a)", fontsize=FontSize, transform=axFdiff.transAxes)


    # plot tan alpha
    nlegend = 4
    gap=len(optalpha_m2)//nlegend
    if not gap:
        gap=1
    print("gap",gap)
    for i in range(len(optalpha_m2))[::gap]:
        select = optFtot_diff[i]<0
        #select = optFtot_diff[i]<100000
        axalpha.plot(qs[select][:],np.tan(optalpha_m2[i][select][:]),"-",ms=msize,mfc="None",lw=LineWidth,label="%.1f"%Cs[i])
        #print("K[i]=",Ks[i])
        print("optalpha_m2[i][select]=",optalpha_m2[i][select])
    #axalpha.set_xlabel(r"$qR$",fontsize=LabelSize)
    axalpha.legend(loc="upper left",ncol=1,columnspacing=0.1,handlelength=0.5,handletextpad=0.5,labelspacing=0.1,markerscale=1,frameon=False,fontsize=LabelSize)
    axalpha.set_ylabel(r"$\tan\alpha$",fontsize=LabelSize)
    axalpha.yaxis.set_label_position("right")
    axalpha.yaxis.tick_right()
    axalpha.tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=False, labelleft=False, labelsize=LabelSize)
    axalpha.xaxis.set_major_locator(MultipleLocator(1))
    axalpha.xaxis.set_minor_locator(MultipleLocator(0.5))
    axalpha.yaxis.set_major_locator(MultipleLocator(0.5))
    axalpha.yaxis.set_minor_locator(MultipleLocator(0.25))
    x1, y1 = 0.7, 0.1
    axalpha.text(x1, y1, r"(b)", fontsize=FontSize, transform=axalpha.transAxes)


    # plot gamma
    for i in range(len(optalpha_m2))[::gap]:
        select = optFtot_diff[i]<0
        axgamma.plot(qs[select],optgamma_m2[i][select],"-",ms=msize,mfc="None",lw=LineWidth)
    axgamma.set_xlabel(r"$qR$",fontsize=LabelSize)
    axgamma.set_ylabel(r"$\gamma$",fontsize=LabelSize)
    axgamma.yaxis.set_label_position("right")
    axgamma.yaxis.tick_right()
    axgamma.tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=True, labelleft=False, labelsize=LabelSize)

    axgamma.xaxis.set_major_locator(MultipleLocator(1))
    axgamma.xaxis.set_minor_locator(MultipleLocator(0.5))

    axgamma.yaxis.set_major_locator(MultipleLocator(0.2))
    axgamma.yaxis.set_minor_locator(MultipleLocator(0.1))

    x1, y1 = 0.7, 0.1
    axgamma.text(x1, y1, r"(c)", fontsize=FontSize, transform=axgamma.transAxes)


    plt.tight_layout(pad=0.01)
    plt.savefig(folder+"/two_mod_del_E_m%d_m%d.pdf"%(mi,mj),format="pdf")
    #plt.show()
    plt.close()


# see how these various with director field
def F_compot_param(folder, m, alphas, gammas,  bn_phi, bn_z, R):
    intSS,intTT,intT,intBB,intC = [],[],[],[],[]
    alpha_write,gamma_write = [],[]
    for alpha in alphas:
        for gamma in gammas:
            intSS.append(intSS_unit_length(m,alpha,gamma,bn_phi,bn_z,R))
            intTT.append(intTT_unit_length(m,alpha,gamma,bn_phi,bn_z,R))
            intT.append(intT_unit_length(m,alpha,gamma,bn_phi,bn_z,R))
            intBB.append(intBB_unit_length(m,alpha,gamma,bn_phi,bn_z,R))
            intC.append(intC_unit_length(m,alpha,gamma,bn_phi,bn_z,R))
            alpha_write.append(alpha)
            gamma_write.append(gamma)

    savename = folder + "/F_compot_m%.0f_R%.1f_alphas_gammas.csv" % (m,R)
    with open(savename, "w") as f:
        f.write("alpha,gamma,intSS,intTT,intT,intBB,intC\n")
        for i in range(len(intSS)):
            f.write("%f,%f,%f,%f,%f,%f,%f\n" % (alpha_write[i],gamma_write[i],intSS[i],intTT[i],intT[i],intBB[i],intC[i]))

from mpl_toolkits.axes_grid1 import make_axes_locatable
def F_compot_param_plot(filename,m,ms=12):
    data = np.loadtxt(filename,delimiter=",",skiprows=1,unpack=True)
    alpha,gamma,intSS,intTT,intT,intBB,intC = data
    cmap = cm.get_cmap("jet_r")
    plt.figure()
    fig, axs = plt.subplots(5,1,figsize=(4,3*4))
    axs[0].set_title("m=%d"%m)
    im0 = axs[0].scatter(alpha,gamma,c=intSS,cmap=cmap,s=ms)
    select = intSS==0
    axs[0].scatter(alpha[select],gamma[select],marker="+",color="white")
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(im0, cax=cax, orientation='vertical')
    cbar.ax.set_ylabel(r"$\int S^2$")

    im1 = axs[1].scatter(alpha,gamma,c=intTT,cmap=cmap,s=ms)
    select = intTT==0
    axs[1].scatter(alpha[select],gamma[select],marker="+",color="white")
    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(im1, cax=cax, orientation='vertical')
    cbar.ax.set_ylabel(r"$\int T^2$")

    im1 = axs[2].scatter(alpha,gamma,c=intT,cmap=cmap,s=ms)
    select = intT==0
    axs[2].scatter(alpha[select],gamma[select],marker="+",color="white")
    divider = make_axes_locatable(axs[2])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(im1, cax=cax, orientation='vertical')
    cbar.ax.set_ylabel(r"$\int T$")

    im2 = axs[3].scatter(alpha,gamma,c=intBB,cmap=cmap,s=ms)
    select = intBB==0
    axs[3].scatter(alpha[select],gamma[select],marker="+",color="white")
    divider = make_axes_locatable(axs[3])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(im2, cax=cax, orientation='vertical')
    cbar.ax.set_ylabel(r"$\int B^2$")

    im3 = axs[4].scatter(alpha,gamma,c=intC,cmap=cmap,s=ms)
    select = intC==0
    axs[4].scatter(alpha[select],gamma[select],marker="+",color="white")
    divider = make_axes_locatable(axs[4])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(im3, cax=cax, orientation='vertical')
    cbar.ax.set_ylabel(r"$\int 1- (u\cdot n)^2$")

    for ax in axs:
        ax.set_xlabel(r"$\alpha$")
        ax.set_ylabel(r"$\gamma$")
    #cbar=plt.colorbar(ax=axs[0],orientation="horizontal")
    #plt.show()
    plt.tight_layout()
    plt.savefig(filename[:-4]+".pdf",format="pdf")
    plt.close()


def Ftot_per_q_gamma_plot_draft(C, m, bn_phi,bn_z,R):
    qs = np.arange(0.0,5.0,0.2)
    gammas = np.arange(0.0,1.0,0.05)
    Ftots = np.zeros((len(qs),len(gammas)))
    print("np.shape(Ftots),len(Ftots),len(Ftots[0])",np.shape(Ftots),len(Ftots),len(Ftots[0]))
    print("np.shape(qs)",np.shape(qs))
    print("np.shape(gammas)",np.shape(gammas))
    #return 0
    for i in range(len(qs)):
        print("calculating q=%.1f"%qs[i])
        alpha = np.arctan(qs[i]) # set tan alpha to q
        print("alpha",alpha)
        for j in range(len(gammas)):
            Ftot = Ftot_unit_length(1, C, qs[i], m, alpha, gammas[j], bn_phi, bn_z, R)
            print("i,j",i,j)
            print("q,gamma",qs[i],gammas[j])
            print("Ftot,", Ftot)
            Ftots[i][j] = Ftot #- np.pi*qs[i]**2

    ppi = 72
    plt.figure()
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 1.7))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axFtot = plt.subplot2grid((2, 1), (0, 0))
    axFtotq = plt.subplot2grid((2, 1), (1, 0),sharex=axFtot)
    msize = 4
    gamma_mesh, q_mesh  = np.meshgrid(gammas,qs)
    #axFtot.contourf(q_mesh,gamma_mesh,Ftot)
    cfc = axFtot.contourf(gamma_mesh,q_mesh,Ftots,alpha=1)
    #cs = axFdiff.contour(cfc,levels=[0],colors=("k"),linestyles="--",linewidths=(LineWidth))
    for i in range(len(qs))[::5]:
        axFtotq.plot(gammas,Ftots[i],label="q=%.1f"%qs[i])
    #axFtot.set_xlabel(r"$\gamma$")
    axFtot.set_ylabel(r"$q$")

    axFtotq.set_xlabel(r"$\gamma$")
    axFtotq.legend()
    #plt.show()
    plt.savefig("Ftot_per_q_gamma_C%.1f_m%.0f.pdf"%(C,m))


def Ftot_per_q_gamma_plot( m, bn_phi,bn_z,R):
    Cs = [1,4,8]
    qs = np.arange(0.0,3.1,0.5)
    gammas = np.arange(0.0,1.0,0.05)
    #return 0

    ppi = 72
    plt.figure()
    fig = plt.figure(figsize=(246 / ppi * 0.4*len(qs), 246 / ppi * 0.8*len(Cs)))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axss = []
    for i in range(len(Cs)):
        axFtotq = plt.subplot2grid((3, 1), (len(Cs)-i-1, 0))
        divider = make_axes_locatable(axFtotq)
        axs = [axFtotq]
        for i in range(len(qs)-1):
            axs.append(divider.append_axes("right",size="100%",pad=0.03))
        axss.append(axs)

    msize = 4
    for i in range(len(Cs)):
        for j in range(len(qs)):
            print("plotting C%f,q%f"%(Cs[i],qs[j]))
            alpha = np.arctan(qs[j]) # set tan alpha to q
            Ftots = []
            for k in range(len(gammas)):
                Ftot = Ftot_unit_length(1, Cs[i], qs[j], m, alpha, gammas[k], bn_phi, bn_z, R)
                Ftots.append(Ftot)

            axss[i][j].plot(Ftots,gammas,label="q=%.1f"%qs[j])
            axss[i][j].set_title(r"$(C,q)=(%.0f,%.1f)$"%(Cs[i],qs[j]),fontsize=8)
            argmin_gamma = np.argmin(Ftots)
            axss[i][j].plot(Ftots[argmin_gamma],gammas[argmin_gamma],"o")
        axss[i][0].set_ylabel(r"$\gamma$")

    axss[0][len(qs)//2].set_xlabel(r"$E(\tan\alpha=q),m=%d$"%m)

    plt.tight_layout()
    #plt.show()
    plt.savefig("Ftot_per_q_gamma_m%.0f.pdf"%m)

def min_m_plot(folder,K,Cs,R,ms):
    optFms = []
    optalphas,optgammas=[],[]
    for i in range(len(ms)):
        optFms.append([])
        optalphas.append([])
        optgammas.append([])
        for C in Cs:
            filename_m = folder + "/optFtot_K%.2f_C%.1f_m%d_R%.1f_qs.csv" % (K, C, ms[i], R)
            qs,optFtot,optalpha,optgamma = np.loadtxt(filename_m, skiprows=1, delimiter=",", unpack=True)
            optFms[i].append(optFtot)
            optalphas[i].append(optalpha)
            optgammas[i].append(optgamma)

    qmesh,Cmesh = np.meshgrid(qs,Cs)
    optFms = np.array(optFms)
    #Ftot_min = np.min(allFtot,axes)
    optFms = np.moveaxis(optFms,[0,1,2],[2,0,1])
    m_of_minFtot = np.argmin(optFms,axis=2)
    print(m_of_minFtot)


    plt.figure()
    axF = plt.subplot2grid((2, 3), (0, 0),rowspan=2, colspan=2)
    axalpha = plt.subplot2grid((2, 3), (0, 2),rowspan=1)
    axgamma = plt.subplot2grid((2, 3), (1, 2),rowspan=1,sharex=axalpha)
    msize = 4
    # tell the colorbar to tick at integers
    axF.pcolormesh(qmesh,Cmesh,m_of_minFtot,shading="auto",cmap=cm.get_cmap("rainbow"),vmin=-0.5,vmax=3.5)
    #axF.colorbar(ticks=np.arange(0,4))
    axF.set_xlabel(r"$qR$")
    axF.set_ylabel(r"$CR^2/K$")

    ccolors = ["gold","red", "blue", "purple", "black"]
    lss = [":","--","-"]
    c = 0
    #C_to_plot = [1.0,1.5,2.5,4]
    C_to_plot = [1.0,2.5,4,9]

    for C in C_to_plot:
        i = int(C*10)-1
        for m in [1,2,3]:
            select = m_of_minFtot[i]==m
            #select = optFtot_diff[i]<100000
            if m==2:
                optalpha =optalphas[m][i]
                axalpha.plot(qs[select][:],np.tan(optalpha[select][:]),ls=lss[m-1],ms=msize,mfc="None",label="%.1f"%Cs[i],color=ccolors[c])

                optgamma =optgammas[m][i]
                axgamma.plot(qs[select][:],np.tan(optgamma[select][:]),ls=lss[m-1],ms=msize,mfc="None",label="%.1f"%Cs[i],color=ccolors[c])

            else:
                optalpha =optalphas[m][i]
                axalpha.plot(qs[select][:],np.tan(optalpha[select][:]),ls=lss[m-1],ms=msize,mfc="None",color=ccolors[c])

                optgamma =optgammas[m][i]
                axgamma.plot(qs[select][:],np.tan(optgamma[select][:]),ls=lss[m-1],ms=msize,mfc="None",color=ccolors[c])

        c=c+1

    axalpha.legend(ncol=1,columnspacing=0.1,handlelength=0.5,handletextpad=0.25,labelspacing=0.1,markerscale=1,frameon=False)
    axalpha.set_ylabel(r"$\tan\alpha$")
    axalpha.yaxis.set_label_position("right")
    axalpha.yaxis.tick_right()
    axgamma.set_ylabel(r"$\gamma$")
    axgamma.yaxis.set_label_position("right")
    axgamma.yaxis.tick_right()

    #plt.show()

    plt.savefig(folder+"/min_Ftot_m.pdf",format="pdf")

    return m_of_minFtot
