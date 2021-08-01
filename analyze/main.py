#!/opt/homebrew/bin/python3

import time
import numpy as np
import matplotlib.pyplot as plt
from plot import *
from Oplot import *
from analyze import *


def main():
    print("hello! dtmc_lc analysis")

    #config_plot3D("../data/scratch_local/State_N200_imod1_Ne1_kar10.0_C00.2_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0.csv",mesh=1,rod=0)
    #return 0

    #foldername = "../data/Ne2/Apr29_2021" # kinda switch to Ne2 simulation for next paper
    #foldername = "../data/Ne1/May4_2021" # mobius strip~,
    foldername = "../data/Ne1/Aug1_2021"
    print("analyzing "+foldername)
    N = 200
    imod=1 # 1 for rhombus, 2 disk, 3 cylinder, 4 for mobius strip
    Ne=1
    kar = 10.0
    C0s = np.arange(0.0,0.91,0.1)
    karg = 0.0
    lams=np.arange(2.0,10.1,0.5)
    Kd= 0.0
    q = 0.0
    #qs=np.arange(0.2,1.1,0.2)
    Cn=0.0
    pars = []
    for C0 in C0s:
        Cn = Kd
        pars.append([N, imod, Ne, kar, C0, karg, lams, Kd, q, Cn])
    par_nm = ["N","imod", "Ne", "kar", "C0","karg", "lam","Kd","q", "Cn"]
    par_dg = [0,0,0,1,1,1,1,1,1,1] # number of digit for each
    mod="lam"
    for i in range(len(pars)):
        print("analyzing",pars[i])
        N, imod, Ne, kar, C0, karg, lam, Kd, q, Cn = pars[i]
        O_stat_ana(foldername,pars[i],par_nm,par_dg, mode=mod, CnequalsKc=0, tau_c=6)
        if(1):
            pass
            Gij_stat_ana(foldername,pars[i],par_nm,par_dg,mode=mod, tau_c=6)

        for lam in lams[::1]:
            if(i%1==0):
                pass
                filename = foldername + "/State_N%.0f_imod%.0f_Ne%.0f_kar%.1f_C0%.1f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f.csv" % (N, imod,Ne,kar,C0,karg,lam,Kd,q,Cn)
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