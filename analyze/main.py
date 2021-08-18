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

    #config_plot3D("../data/scratch_local/State_N200_imod1_Ne1_kar20_J2.0_C00.5_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0_init.csv",mesh=1,rod=0,fnormal=1)
    #config_plot3D("../data/scratch_local/State_N200_imod1_Ne1_kar20_J2.0_C00.5_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0.csv",mesh=1,rod=0,fnormal=0)
    #config_plot3D("../data/Ne2/Aug4_2021_1/State_N300_imod3_Ne2_kar50_C00.0_karg0.0_lam8.0_Kd10.0_q2.0_Cn20.0.csv",mesh=1,rod=0)

    #return 0

    #foldername = "../data/Ne2/Apr29_2021" # kinda switch to Ne2 simulation for next paper
    #foldername = "../data/Ne1/May4_2021" # mobius strip~,
    foldername = "../data/Ne2/Aug17_2021"
    print("analyzing "+foldername)
    N = 300
    imod=1 # 1 for rhombus, 2 disk, 3 cylinder, 4 for mobius strip
    Ne=1
    kar = 100
    Js = np.arange(0.03,0.61,0.03)
    C0s = np.arange(0.0,0.91,0.1)
    karg = 0.0
    lam = 5.0
    Kd = 0.0
    q = 0.0
    Cn = 0.0
    pars = []
    pars1,pars2 = [],[]
    for C0 in C0s:
        Cn = Kd
        pars.append([N, imod, Ne, kar, Js, C0, karg, lam, Kd, q, Cn])
        #pars1.append([N, 1, 1, kar, C0, karg, lam, Kd, qs, Cn])
        #pars2.append([N, 3, 2, kar, C0, karg, lam, Kd, qs, Cn])
    par_nm = ["N","imod", "Ne", "kar", "J" ,"C0","karg", "lam","Kd","q", "Cn"]
    par_dg = [0,0,0,0,2,1,1,1,1,1,1] # number of digit for each
    mod="J"

    # run test plot
    #E_compare(foldername, pars1, pars2, par_nm, par_dg, mod)
    #return 0

    for i in range(len(pars)):
        print("analyzing",pars[i])
        N, imod, Ne, kar, J, C0, karg, lam, Kd, q, Cn = pars[i]
        O_stat_ana(foldername,pars[i],par_nm,par_dg, mode=mod, CnequalsKc=0, tau_c=6)
        if(1):
            pass
            Gij_stat_ana(foldername,pars[i],par_nm,par_dg,mode=mod, tau_c=6)

        for J in Js[::1]:
            if(i%1==0):
                pass
                filename = foldername + "/State_N%.0f_imod%.0f_Ne%.0f_kar%.0f_J%.2f_C0%.1f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f.csv" % (N, imod,Ne,kar,J,C0,karg,lam,Kd,q,Cn)
                ctag = r"$J=%.2f,C_0=%.1f$" %(J,C0)
                #config_plot_xyz(filename,mesh=0,rod=1,tag=ctag,Format="png")
                #config_plot_xyz(filename,mesh=1,rod=0,tag=ctag,Format="png")

    colors = None
    alphas = None

    Geig_pars_plot(foldername, pars,par_nm,par_dg, mode=mod)
    Os_pars_plot(foldername, pars,par_nm,par_dg,mode=mod)
    # additional test

if __name__ == '__main__':
    main()

##!/usr/local/bin/python3