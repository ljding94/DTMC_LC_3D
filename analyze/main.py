#!/usr/local/bin/python3
#TODO: can be replaced with /opt/homebrew/bin/python3 once scipy can be installed for apple silicon
import time
import numpy as np
import matplotlib.pyplot as plt
from plot import *
from Oplot import *
from analyze import *


def main():
    print("hello! dtmc_lc analysis")
    #config_plot3D("../data/Ne2/Apr18_2021/State_N300_imod3_Ne3_kar100_lam4.0_Kd5.0_q1.4_Cn5.0_kard0.0.txt",mesh=1,rod=0)
    config_plot3D("../data/scratch_local/State_N198_imod4_Ne1_kar5_lam3.0_Kd5.0_q1.0_Cn5.0_Cnp0.0_rCnp0.0_init.csv",mesh=1)
    config_plot3D("../data/scratch_local/State_N198_imod4_Ne1_kar5_lam3.0_Kd5.0_q1.0_Cn5.0_Cnp0.0_rCnp0.0.csv",mesh=1)
    #config_plot_xyz("../data/scratch_local/State_N200_imod1_Ne2_kar100_lam6.0_Kd5.0_q1.0_Cn5.0_kard0.0_init.txt",mesh=1)
    return 0

    foldername = "../data/Ne2/Apr29_2021" # kinda switch to Ne2 simulation for next paper
    foldername = "../data/Ne1/May2_2021" # mobius strip~,
    print("analyzing "+foldername)
    N = 198
    imod=4 # 1 for rhombus, 2 disk, 3 cylinder, 4 for mobius strip
    Ne=1
    kars = [2,5,10,15]
    karg=0.0
    lam=5.0
    lams=np.arange(3.0,10.1,0.5)
    Kd=3.0
    Cn=Kd
    q=2.0
    Cnp=0
    rCnp=0
    Cnps=np.arange(1.0,10.1,1.0)
    rCnps=np.arange(0.1,0.91,0.2)
    pars = []
    for kar in kars:
        pars.append([N, imod,Ne, kar, lams, Kd, q, Cn, Cnp, rCnp])
    par_nm = ["N","imod", "Ne", "kar","lam","Kd","q", "Cn", "Cnp", "rCnp"]
    par_dg = [0,0,0,0,1,1,1,1,1,1] # nsumber of digit for each
    mod="lam"
    for i in range(len(pars)):
        print("analyzing",pars[i])
        N, imod, Ne, kar, lam, Kd, q,Cn,Cnp,rCnp = pars[i]
        O_stat_ana(foldername,pars[i],par_nm,par_dg, mode=mod, tau_c=6)
        if(1):
            pass
            Gij_stat_ana(foldername,pars[i],par_nm,par_dg,mode=mod, tau_c=6)

        for lam in lams[::2]:
            if(i%1==0):
                pass
                filename = foldername + "/State_N%.0f_imod%.0f_Ne%.0f_kar%.0f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_Cnp%.1f_rCnp%.1f.csv" % (N, imod,Ne, kar,lam,Kd,q,Cn,Cnp,rCnp)
                config_plot_xyz(filename,mesh=1,rod=0,tag=r"$K_d=%.1f,q=%.1f,C_n=%.1f$" % (Kd,q,Cn),Format="png")
                config_plot_xyz(filename,mesh=1,rod=1,tag=r"$K_d=%.1f,q=%.1f,C_n=%.1f$" % (Kd,q,Cn),Format="png")

    colors = None
    alphas = None

    Geig_pars_plot(foldername, pars,par_nm,par_dg, mode=mod)
    Os_pars_plot(foldername, pars,par_nm,par_dg,mode=mod)

    # additional test

if __name__ == '__main__':
    main()