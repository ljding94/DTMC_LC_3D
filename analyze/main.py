#!/usr/local/bin/python3
import time
import numpy as np
import matplotlib.pyplot as plt
from plot import *
from Oplot import *
from analyze import *


def main():
    print("hello! dtmc_lc analysis")
    #config_plot3D("../data/Ne1/Mar23_2021/State_N500_Ne1_L-1_kar100_karg0.0_lam6.0_Kd7.0_q2.0_Cn7.0_kard0.0.txt",mesh=1,rod=0)
    #return 0

    foldername = "../data/Ne2/Apr15_2021" # kinda switch to Ne2 simulation for next paper
    print("analyzing "+foldername)
    N = 400
    imod=3 # 1 for rhombus, 2 disk, 3 cylinder
    Ne=2
    kar = 100
    karg=0.0
    lam=5.0
    Kd=7.0
    qs=np.arange(0.0,2.01,0.1)
    Cns = np.arange(3.0,7.1,2.0)
    Cn=5.0
    kard=0.0
    pars = []
    for Cn in Cns:
        pars.append([N, imod,Ne, kar, lam, Kd, qs, Cn, kard])
    par_nm = ["N","imod" ,"Ne", "kar","lam","Kd","q", "Cn", "kard"]
    par_dg = [0,0,0,0,1,1,1,1,1] # nsumber of digit for each
    mod="q"
    for i in range(len(pars)):
        print("analyzing",pars[i])
        N, imod, Ne, kar, lam, Kd, q,Cn,kard = pars[i]
        #O_stat_ana(foldername,pars[i],par_nm,par_dg, mode=mod, tau_c=6)
        if(1):
            pass
            #Gij_stat_ana(foldername,pars[i],par_nm,par_dg,mode=mod, tau_c=6)

        for q in qs[::2]:
            if(i%1==0):
                pass
                filename = foldername + "/State_N%.0f_imod%.0f_Ne%.0f_kar%.0f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_kard%.1f.txt" % (N, imod,Ne, kar,lam,Kd,q,Cn, kard)
                #config_plot_xyz(filename,mesh=0,rod=1,tag=r"$K_d=%.1f,q=%.1f,C_n=%.1f$" % (Kd,q,Cn),Format="png")

    colors = None
    alphas = None

    Geig_pars_plot(foldername, pars,par_nm,par_dg, mode=mod)
    #Os_pars_plot(foldername, pars,par_nm,par_dg,mode=mod)

    # additional test

if __name__ == '__main__':
    main()