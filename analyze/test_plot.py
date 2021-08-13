import numpy as np
import matplotlib.pyplot as plt
from Oplot import *

def E_compare(foldername, pars1, pars2,par_nm,par_dg,mode):
    colors, alphas = None, None
    data1,data2,O_label1, O_label2 = [], [],[],[]
    xLabel = mode
    cpar_ind = find_cpar_ind(par_nm,mode)
    for i in range(len(pars1)):
        par1,par2 = pars1[i],pars2[i]
        par_dealing1, par_dealing2 = par1[:],par2[:]
        f2rtail1,f2rtail2 = "MC","MC"
        label1, label2 = "",""
        for j in range(len(par1)):
            if(j==cpar_ind):
                f2rtail1+="_"+par_nm[j]+"s"
                f2rtail2+="_"+par_nm[j]+"s"
            else:
                f2rtail1+="_"+par_nm[j]+"%.*f"%(par_dg[j],par_dealing1[j])
                f2rtail2+="_"+par_nm[j]+"%.*f"%(par_dg[j],par_dealing2[j])
                if (par_dealing1[j]!=par_dealing2[j] or par_dealing1[j]!=pars1[0][j]or par_dealing1[j]!=pars1[-1][j]):
                    label1+=par_nm[j]+"%.*f,"%(par_dg[j],par_dealing1[j])
                    label2+=par_nm[j]+"%.*f,"%(par_dg[j],par_dealing2[j])
        f2rtail1+="_ana.csv"
        f2rtail2+="_ana.csv"
        head = "/O_"
        filename1 = foldername + head+f2rtail1
        filename2 = foldername + head+f2rtail2
        data1.append(np.loadtxt(filename1, skiprows=1,delimiter=",", unpack=True))
        data2.append(np.loadtxt(filename2, skiprows=1,delimiter=",", unpack=True))
        O_label1.append(label1)
        O_label2.append(label2)
    data1 = np.transpose(np.array(data1), axes=(1, 0, 2))
    data2 = np.transpose(np.array(data2), axes=(1, 0, 2))
    print("len(data1)",len(data1))
    cpar1, E_ave1, E_tau1, E_err1 = data1[:4]
    cpar2, E_ave2, E_tau2, E_err2 = data2[:4]
    print("len(cpar1)",len(cpar1))
    print("len(E_ave1)",len(E_ave1))

    ppi = 72
    plt.figure()
    plt.rc('text', usetex=True)
    fig, axs = plt.subplots(len(E_ave1), 2, figsize=(
        246 / ppi*2, 246 / ppi * len(E_ave1)*0.8), sharex=True)  # , sharex=True
    for i in range(len(E_ave1)):
        axs[i,0].errorbar(cpar1[i],E_ave1[i],yerr=E_err1[i],linestyle="--",label=O_label1[i])
        axs[i,0].errorbar(cpar2[i],E_ave2[i],yerr=E_err2[i],linestyle=":",label=O_label2[i])
        axs[i,1].plot(cpar1[i],E_ave1[i]-E_ave2[i],label="E1-E2")
        axs[i,0].legend()
    axs[0,0].set_ylabel("E/N")
    axs[-1,0].set_xlabel(mode)
    axs[0,1].legend()
    axs[-1,1].set_xlabel(mode)

    plt.tight_layout()
    plt.savefig(foldername +"/E_compare_"+mode+".pdf",format="pdf", transparent=True)
    plt.close()
