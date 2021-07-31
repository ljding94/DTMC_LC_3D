#!/opt/homebrew/bin/python3

import time
import numpy as np
import matplotlib.pyplot as plt
from plot import *
from Oplot import *
from analyze import *


def main():
    print("hello! dtmc_lc analysis")

    config_plot3D("../data/scratch_local/State_N200_imod1_Ne1_kar10.0_C00.2_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0.csv",mesh=1,rod=0)
    return 0

    #foldername = "../data/Ne2/Apr29_2021" # kinda switch to Ne2 simulation for next paper
    #foldername = "../data/Ne1/May4_2021" # mobius strip~,
    foldername = "../data/Ne2/Jul27_2021"
    print("analyzing "+foldername)
    N = 300
    imod=1 # 1 for rhombus, 2 disk, 3 cylinder, 4 for mobius strip
    Ne=2
    kars = np.arange(4.0,20.1,4.0)
    kargs = np.arange(2.0,20.1,2.0)
    lam=10.0
    lams=np.arange(3.0,11.1,2.0)
    Kd= 12.0
    q = 0.7
    #qs=np.arange(0.2,1.1,0.2)
    Cn=0.0
    pars = []
    for kar in kars:
        Cn = Kd
        pars.append([N, imod, Ne, kar, kargs, lam, Kd, q, Cn])
    par_nm = ["N","imod", "Ne", "kar","karg", "lam","Kd","q", "Cn"]
    par_dg = [0,0,0,1,1,1,1,1,1] # number of digit for each
    mod="karg"
    for i in range(len(pars)):
        print("analyzing",pars[i])
        N, imod, Ne, kar,karg, lam, Kd, q, Cn = pars[i]
        O_stat_ana(foldername,pars[i],par_nm,par_dg, mode=mod, CnequalsKc=0, tau_c=6)
        if(1):
            pass
            Gij_stat_ana(foldername,pars[i],par_nm,par_dg,mode=mod, tau_c=6)

        for karg in kargs[::1]:
            if(i%1==0):
                pass
                filename = foldername + "/State_N%.0f_imod%.0f_Ne%.0f_kar%.1f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f.csv" % (N, imod,Ne,kar,karg,lam,Kd,q,Cn)
                #config_plot_xyz(filename,mesh=0,rod=1,tag=r"$\bar{\kappa}=%.1f,\lambda=%.1f$" % (karg,lam),Format="png")
                config_plot_xyz(filename,mesh=1,rod=0,tag=r"$\bar{\kappa}=%.1f,\lambda=%.1f$" % (karg,lam),Format="png")

    colors = None
    alphas = None

    Geig_pars_plot(foldername, pars,par_nm,par_dg, mode=mod)
    Os_pars_plot(foldername, pars,par_nm,par_dg,mode=mod)

    # additional test

if __name__ == '__main__':
    main()

##!/usr/local/bin/python3