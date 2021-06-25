#!/opt/homebrew/bin/python3

import time
import numpy as np
import matplotlib.pyplot as plt
from plot import *
from Oplot import *
from analyze import *


def main():
    print("hello! dtmc_lc analysis")

    #config_plot3D("../data/Ne1/May23_2021/State_N300_imod1_Ne1_kar100_lam5.0_Kd3.0_q0.0_Cn3.0_kard0.0_lamd4.0.csv",mesh=1,rod=1)
    #return 0

    #foldername = "../data/Ne2/Apr29_2021" # kinda switch to Ne2 simulation for next paper
    #foldername = "../data/Ne1/May4_2021" # mobius strip~,
    foldername = "../data/Ne2/Jun23_2021"
    print("analyzing "+foldername)
    N = 300
    imod=3 # 1 for rhombus, 2 disk, 3 cylinder, 4 for mobius strip
    Ne=2
    kar = 20
    kargs=np.arange(0.5,5.1,0.5)
    lam=5.0
    lams=np.arange(3.0,9.1,2.0)
    Kd=3.0
    q=1.0
    Cn=3.0
    pars = []
    for lam in lams:
        pars.append([N, imod, Ne, kar, kargs,lam, Kd, q, Cn])
    par_nm = ["N","imod", "Ne", "kar","karg", "lam","Kd","q", "Cn"]
    par_dg = [0,0,0,0,1,1,1,1,1] # nsumber of digit for each
    mod="karg"
    for i in range(len(pars)):
        print("analyzing",pars[i])
        N, imod, Ne, kar,karg, lam, Kd, q,Cn = pars[i]
        O_stat_ana(foldername,pars[i],par_nm,par_dg, mode=mod, tau_c=6)
        if(1):
            pass
            Gij_stat_ana(foldername,pars[i],par_nm,par_dg,mode=mod, tau_c=6)

        for karg in kargs[::1]:
            if(i%1==0):
                pass
                filename = foldername + "/State_N%.0f_imod%.0f_Ne%.0f_kar%.0f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f.csv" % (N, imod,Ne,kar,karg,lam,Kd,q,Cn)

                config_plot_xyz(filename,mesh=0,rod=1,tag=r"$\bar{\kappa}=%.1f,\lambda=%.1f$" % (karg,lam),Format="png")

    colors = None
    alphas = None

    Geig_pars_plot(foldername, pars,par_nm,par_dg, mode=mod)
    Os_pars_plot(foldername, pars,par_nm,par_dg,mode=mod)

    # additional test

if __name__ == '__main__':
    main()

##!/usr/local/bin/python3