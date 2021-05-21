#!/opt/homebrew/bin/python3

import time
import numpy as np
import matplotlib.pyplot as plt
from plot import *
from Oplot import *
from analyze import *


def main():
    print("hello! dtmc_lc analysis")
    #config_plot3D("../data/Ne2/Apr18_2021/State_N300_imod3_Ne3_kar100_lam4.0_Kd5.0_q1.4_Cn5.0_kard0.0.txt",mesh=1,rod=0)
    #config_plot3D("../data/scratch_local/State_N200_imod1_Ne1_kar20_lam0.0_Kd3.0_q0.0_Cn3.0_kard0.0_lamd5.0_init.csv",mesh=1,rod=1)
    #config_plot3D("../data/scratch_local/State_N200_imod1_Ne1_kar20_lam0.0_Kd3.0_q0.0_Cn3.0_kard0.0_lamd5.0.csv",mesh=1,rod=1)
    #config_plot3D("../data/Ne1/May18_2021/State_N200_imod1_Ne1_kar100_lam5.0_Kd6.0_q0.0_Cn6.0_kard0.0_lamd7.5.csv",mesh=1,rod=1)
    #return 0

    #foldername = "../data/Ne2/Apr29_2021" # kinda switch to Ne2 simulation for next paper
    #foldername = "../data/Ne1/May4_2021" # mobius strip~,
    foldername = "../data/Ne1/May20_2021"
    print("analyzing "+foldername)
    N = 200
    imod=1 # 1 for rhombus, 2 disk, 3 cylinder, 4 for mobius strip
    Ne=1
    kar = 100
    karg=0.0
    lam=5.0
    lams=np.arange(4.0,10.1,1.0)
    Kd=6.0
    qs=np.arange(0.0,2.01,0.5)
    kard=0.0
    lamds=np.arange(0.0,10.01,0.5)
    pars = []
    for q in qs:
        Cn=Kd
        pars.append([N, imod,Ne, kar, lam, Kd, q, Cn, kard,lamds])
    par_nm = ["N","imod", "Ne", "kar","lam","Kd","q", "Cn", "kard","lamd"]
    par_dg = [0,0,0,0,1,1,1,1,1,1] # nsumber of digit for each
    mod="lamd"
    for i in range(len(pars)):
        print("analyzing",pars[i])
        N, imod, Ne, kar, lam, Kd, q,Cn,kard,lamd = pars[i]
        O_stat_ana(foldername,pars[i],par_nm,par_dg, mode=mod, tau_c=6)
        if(1):
            pass
            Gij_stat_ana(foldername,pars[i],par_nm,par_dg,mode=mod, tau_c=6)

        for lamd in lamds[::1]:
            if(i%1==0):
                pass
                filename = foldername + "/State_N%.0f_imod%.0f_Ne%.0f_kar%.0f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_kard%.1f_lamd%.1f.csv" % (N, imod,Ne, kar,lam,Kd,q,Cn,kard,lamd)
                #config_plot_xyz(filename,mesh=1,rod=0,tag=r"$K_d=%.1f,q=%.1f,C_n=%.1f$" % (Kd,q,Cn),Format="png")
                config_plot_xyz(filename,mesh=0,rod=1,tag=r"$K_d=%.1f,q=%.1f,C_n=%.1f$" % (Kd,q,Cn),Format="png")

    colors = None
    alphas = None

    Geig_pars_plot(foldername, pars,par_nm,par_dg, mode=mod)
    Os_pars_plot(foldername, pars,par_nm,par_dg,mode=mod)

    # additional test

if __name__ == '__main__':
    main()

##!/usr/local/bin/python3