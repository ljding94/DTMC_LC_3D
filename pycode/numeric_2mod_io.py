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

def del_Ftot_Ks_qs_plot(folder,Ks,C,qs,R):
    LineWidth, FontSize, LabelSize = 1,9,8
    optFtot_m0, optalpha_m0, optgamma_m0 = [], [], []
    optFtot_m2, optalpha_m2, optgamma_m2 = [], [], []
    for K in Ks:
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

    qmesh,Kmesh = np.meshgrid(qs,Ks)
    optFtot_diff = np.array(optFtot_m2) - np.array(optFtot_m0)
    absFdiffmax = np.max(np.abs(optFtot_diff))
    print("absFdiffmax=",absFdiffmax)
    ppi = 72
    plt.figure()
    fig, axs = plt.subplots(3, 2, figsize=(246 / ppi * 1, 246 / ppi * 1.6))
    cmap = cm.get_cmap("bwr")
    #cf = axs[0].pcolormesh(qmesh,Kmesh,optFtot_diff, shading="nearest", cmap = cmap,
                            #norm = colors.DivergingNorm(vmin=np.min(optFtot_diff), vcenter = 0,vmax = np.max(optFtot_diff) ))
    Fdiffmin,Fdiffmax = np.min(optFtot_diff),np.max(optFtot_diff)
    print(Fdiffmin,Fdiffmax)
    #Fdiffmin,Fdiffmax = round(Fdiffmin,0),round(Fdiffmax,2)
    levels = [Fdiffmin,Fdiffmin*3/4,Fdiffmin*2/4,Fdiffmin/4,0,Fdiffmax/4,Fdiffmax*2/4,Fdiffmax*3/4,Fdiffmax]
    cf = axs[0].contourf(qmesh,Kmesh,optFtot_diff,levels,cmap=cmap, norm=colors.TwoSlopeNorm(vmin=Fdiffmin,vcenter=0,vmax=Fdiffmax))
    cs = axs[0].contour(cf,levels=[0],colors=("k"),linestyles="--",linewidths=(LineWidth))

    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="3%", pad=0.05)
    cbar=plt.colorbar(cf, cax = cax,ax=axs[0],orientation="vertical",ticks=levels)
    cbar.add_lines(cs)
    cbar.ax.set_yticklabels([str(round(Fdiffmin,1)),"",str(round(Fdiffmin*2/4,1)),"","0","",str(round(Fdiffmax*2/4,2)),"",str(round(Fdiffmax,2))],fontsize=LabelSize)
    #cbar.ax.xaxis.set_ticks_position('top')
    #cbar.ax.set_yticklabels([str(round(Fdiffmin,1)),"",str(round(Fdiffmin*2/4,1)),"","0","",str(round(Fdiffmax*2/4,2)),"",str(round(Fdiffmax,2))],fontsize=LabelSize)
    axs[0].set_xlabel(r"$qR$",fontsize=LabelSize)
    axs[0].set_ylabel(r"$K/(CR^2)$",fontsize=LabelSize)
    axs[0].tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=True, labelleft=True, labelsize=LabelSize)
    cbar.ax.set_title(r"$\Delta E/CR$",fontsize=LabelSize)


    # plot tan alpha
    for i in range(len(optalpha_m2)):
        select = optFtot_diff[i]<0
        axs[1].plot(qs[select],optalpha_m2[i][select],"-",lw=LineWidth)
    axs[1].set_xlabel(r"$qR$",fontsize=LabelSize)
    axs[1].set_ylabel(r"$\tan(\alpha$",fontsize=LabelSize)
    axs[1].yaxis.set_label_position("right")
    axs[1].yaxis.tick_right()
    axs[1].tick_params(which="both", direction="in", bottom="on", top="off", right="on", left="off", labelbottom=True, labelleft=False, labelsize=LabelSize)
    plt.tight_layout()
    plt.savefig(folder+"/two_mod_del_E.pdf",format="pdf")
    #plt.show()
    plt.close()