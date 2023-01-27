import numpy as np
from nematic_cholesteric_2mod_cylinder_cal import *
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

def del_Ftot_Ks_qs_plot(folder,K,Cs,R):
    LineWidth, FontSize, LabelSize = 1,9,8
    optFtot_m0, optalpha_m0, optgamma_m0 = [], [], []
    optFtot_m2, optalpha_m2, optgamma_m2 = [], [], []
    for C in Cs:
        filename_m0 = folder + "/optFtot_K%.2f_C%.1f_m0_R%.1f_qs.csv" % (K, C, R)
        filename_m2 = folder + "/optFtot_K%.2f_C%.1f_m2_R%.1f_qs.csv" % (K, C, R)

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
    fig = plt.figure(figsize=(246 / ppi * 1, 246 / ppi * 0.8))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    axFdiff = plt.subplot2grid((2, 3), (0, 0),rowspan=2, colspan=2)
    axalpha = plt.subplot2grid((2, 3), (0, 2),rowspan=1)
    axgamma = plt.subplot2grid((2, 3), (1, 2),rowspan=1)
    msize = 4

    cmap = cm.get_cmap("bwr")
    #cf = axs[0].pcolormesh(qmesh,Kmesh,optFtot_diff, shading="nearest", cmap = cmap,
                            #norm = colors.DivergingNorm(vmin=np.min(optFtot_diff), vcenter = 0,vmax = np.max(optFtot_diff) ))
    Fdiffmin,Fdiffmax = np.min(optFtot_diff),np.max(optFtot_diff)
    print(Fdiffmin,Fdiffmax)
    #Fdiffmin,Fdiffmax = round(Fdiffmin,0),round(Fdiffmax,2)
    levels = [Fdiffmin,Fdiffmin*3/4,Fdiffmin*2/4,Fdiffmin/4,0,Fdiffmax/4,Fdiffmax*2/4,Fdiffmax*3/4,Fdiffmax]
    cf = axFdiff.contourf(qmesh,Kmesh,optFtot_diff,levels,cmap=cmap, norm=colors.TwoSlopeNorm(vmin=Fdiffmin,vcenter=0,vmax=Fdiffmax))
    cs = axFdiff.contour(cf,levels=[0],colors=("k"),linestyles="--",linewidths=(LineWidth))

    divider = make_axes_locatable(axFdiff)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    #cbar=plt.colorbar(cf, cax = cax,ax=axFdiff,orientation="vertical",ticks=levels)
    cbar=plt.colorbar(cf, cax = cax,ax=axFdiff,orientation="horizontal",ticks=levels)
    cbar.add_lines(cs)
    #cbar.ax.set_yticklabels([str(round(Fdiffmin,1)),"",str(round(Fdiffmin*2/4,1)),"","0","",str(round(Fdiffmax*2/4,2)),"",str(round(Fdiffmax,2))],fontsize=LabelSize)
    cbar.ax.set_xticklabels([str(round(Fdiffmin,1)),"",str(round(Fdiffmin*2/4,1)),"","0","",str(round(Fdiffmax*2/4,2)),"",str(round(Fdiffmax,2))],fontsize=LabelSize)
    cbar.ax.xaxis.set_ticks_position('top')

    #cbar.ax.set_yticklabels([str(round(Fdiffmin,1)),"",str(round(Fdiffmin*2/4,1)),"","0","",str(round(Fdiffmax*2/4,2)),"",str(round(Fdiffmax,2))],fontsize=LabelSize)
    axFdiff.set_xlabel(r"$qR$",fontsize=LabelSize)
    axFdiff.set_ylabel(r"$CR^2/K$",fontsize=LabelSize)
    axFdiff.tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=True, labelleft=True, labelsize=LabelSize)
    cbar.ax.set_title(r"$\Delta ER^2/K$",fontsize=LabelSize)

    # plot tan alpha
    for i in range(len(optalpha_m2))[::1]:
        select = optFtot_diff[i]<0
        axalpha.plot(qs[select][:],np.tan(optalpha_m2[i][select][:]),"-",lw=LineWidth)
        #print("K[i]=",Ks[i])
        print("optalpha_m2[i][select]=",optalpha_m2[i][select])
    axalpha.set_xlabel(r"$qR$",fontsize=LabelSize)
    axalpha.set_ylabel(r"$\tan(\alpha)$",fontsize=LabelSize)
    axalpha.yaxis.set_label_position("right")
    axalpha.yaxis.tick_right()
    axalpha.tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=True, labelleft=False, labelsize=LabelSize)

    # plot gamma
    for i in range(len(optalpha_m2)):
        select = optFtot_diff[i]<0
        axgamma.plot(qs[select],optgamma_m2[i][select],"-",lw=LineWidth)
    axgamma.set_xlabel(r"$qR$",fontsize=LabelSize)
    axgamma.set_ylabel(r"$\gamma$",fontsize=LabelSize)
    axgamma.yaxis.set_label_position("right")
    axgamma.yaxis.tick_right()
    axgamma.tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=True, labelleft=False, labelsize=LabelSize)


    plt.tight_layout()
    plt.savefig(folder+"/two_mod_del_E.pdf",format="pdf")
    #plt.show()
    plt.close()