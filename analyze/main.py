#!/opt/homebrew/bin/python3

import time
import numpy as np
import matplotlib.pyplot as plt
from plot import *
from Oplot import *
from analyze import *
from test_plot import *


def main():
    print("hello! dtmc_lc analysis")

    #config_plot3D("../data/scratch_local/State_N100_imod3_Ne2_kar5_J0.0_C00.0_karg0.0_lam5.0_B1.0_Kd0.0_q0.0_Cn0.0.csv",mesh=1,rod=1,piwall=0,fnormal=0)
    #config_plot3D("../data/scratch_local/State_N100_imod3_Ne2_kar5_J0.0_C00.0_karg0.0_lam5.0_B1.0_Kd0.0_q0.0_Cn0.0.csv",mesh=1,rod=0,piwall=1,fnormal=0)
    #config_plot3D("../data/scratch_local/State_N200_imod3_Ne2_kar20_J0.0_C00.3_karg5.0_lam5.0_Kd4.0_q1.0_Cn4.0.csv",mesh=1,rod=0,fnormal=0)
    #config_plot3D("../data/scratch_local/State_N200_imod3_Ne2_kar40_J0.0_C00.1_karg2.0_lam5.0_B5.0_Kd0.0_q0.0_Cn0.0.csv",mesh=1,rod=0,fnormal=0)
    config_plot3D("../data/scratch_local/State_N200_imod3_Ne2_kar40_J0.0_C00.0_karg0.0_lam5.0_B0.0_Kd1.0_q2.0_Cn10.0.csv",mesh=1,rod=1,fnormal=0)
    #config_plot3D("../data/Ne2/Aug21_2021/State_N300_imod1_Ne2_kar40_J0.40_C00.3_karg2.0_lam5.0_Kd4.0_q1.5_Cn4.0.csv",mesh=1,rod=1,fnormal=0)
    #config_plot3D("../data/Ne2/Aug29_2021_1/State_N500_imod3_Ne2_kar40_J0.00_C00.0_karg3.0_lam5.0_B0.0_Kd16.0_q0.8_Cn40.0.csv",mesh=1,rod=1,piwall=0,fnormal=0)
    #config_plot3D("../data/Ne2/Aug28_2021_1/State_N300_imod3_Ne2_kar40_J0.00_C00.0_karg0.0_lam5.0_B0.0_Kd5.0_q1.20_Cn10.0.csv",mesh=1,rod=1,piwall=0,fnormal=0)
    return 0

    #foldername = "../data/Ne2/Apr29_2021" # kinda switch to Ne2 simulation for next paper
    #foldername = "../data/Ne1/May4_2021" # mobius strip~,
    foldername = "../data/Ne2/Aug28_2021_1"
    print("analyzing "+foldername)
    N = 300
    imod=3 # 1 for rhombus, 2 disk, 3 cylinder, 4 for mobius strip
    Ne=2
    kar = 40
    J = 0.00
    #C0s = np.arange(0.1,0.51,0.1)
    C0 = 0.0
    karg = 0.0
    lam = 5.0
    B = 0.0
    Kd = 1.0
    qs = np.arange(0.80,1.60,0.05)
    Cns = np.arange(4.0,20.1,4.0)
    Cns = [10.0]
    Kd = 5.0
    pars = []
    pars1,pars2 = [],[]
    for Cn in Cns:
        #Cn = 2*Kd
        pars.append([N, imod, Ne, kar, J, C0, karg, lam, B, Kd, qs, Cn])
        #pars1.append([N, 1, 1, kar, C0, karg, lam, Kd, qs, Cn])
        #pars2.append([N, 3, 2, kar, C0, karg, lam, Kd, qs, Cn])
    par_nm = ["N","imod", "Ne", "kar", "J" ,"C0","karg", "lam", "B", "Kd", "q", "Cn"]
    par_dg = [0,0,0,0,2,1,1,1,1,1,2,1] # number of digit for each
    mod="q"

    # run test plot
    #E_compare(foldername, pars1, pars2, par_nm, par_dg, mod)
    #return 0

    for i in range(len(pars)):
        print("analyzing",pars[i])
        N, imod, Ne, kar, J, C0, karg, lam, B, Kd, q, Cn = pars[i]
        #O_stat_ana(foldername,pars[i],par_nm,par_dg, mode=mod, CnequalsKc=0, tau_c=6)
        if(1):
            pass
            #Gij_stat_ana(foldername,pars[i],par_nm,par_dg,mode=mod, tau_c=6)

        for q in qs[::1]:
            if(i%1==0):
                pass
                filename = foldername + "/State_N%.0f_imod%.0f_Ne%.0f_kar%.0f_J%.2f_C0%.1f_karg%.1f_lam%.1f_B%.1f_Kd%.1f_q%.2f_Cn%.1f.csv" % (N, imod,Ne,kar,J,C0,karg,lam,B,Kd,q,Cn)
                ctag = r"$\bar{\kappa}=%.1f,C_0=%.1f,B=%.1f,\lambda=%.1f$" %(karg,C0,B,lam)
                config_plot_xyz(filename,mesh=1,rod=0,piwall=1,tag=ctag,Format="png")
                #config_plot_xyz(filename,mesh=0,rod=1,tag=ctag,Format="png")

    colors = None
    alphas = None

    #Geig_pars_plot(foldername, pars,par_nm,par_dg, mode=mod)
    #Os_pars_plot(foldername, pars,par_nm,par_dg,mode=mod)
    # additional test

if __name__ == '__main__':
    main()

##!/usr/local/bin/python3