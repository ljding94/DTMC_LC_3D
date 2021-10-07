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
    config_plot3D("../data/scratch_local/State_N200_imod1_Ne1_lf0.0_kar50_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0_g3.0.csv",mesh=1,rod=0,fnormal=0)
    #config_plot3D("../data/Ne1/Sep8_2021/State_N500_imod1_Ne1_kar100_J0.00_C00.0_karg0.0_lam6.0_B0.0_Kd5.0_q3.0_Cn40.0.csv",mesh=1,rod=1,fnormal=0,piwall=1)
    return 0

    foldername = "../data/Ne1/Oct6_2021"
    print("analyzing "+foldername)
    N = 400
    imod=1 # 1 for rhombus, 2 disk, 3 cylinder, 4 for mobius strip
    Ne=1
    #lfs = np.arange(20.0,50.1,2.0)
    lf = 0.0
    kars = np.arange(10,61,10)
    karg = 0.0
    #lams = np.arange(5.0,10.1,1.0)
    lam = 11.0
    Kd = 0.0
    #qs = np.arange(0.0, 0.4, 0.1)
    q = 0.0
    Cn = 0.0
    gs = np.arange(0.0,4.01,0.5)
    pars = []
    for kar in kars:
        Cn = Kd
        pars.append([N, imod, Ne, lf, kar, karg, lam, Kd, q, Cn, gs])
    par_nm = ["N","imod", "Ne", "lf", "kar","karg", "lam", "Kd", "q", "Cn","g"]
    par_dg = [0,0,0,1,0,1,1,1,1,1,1,1,1] # number of digit for each
    mod="g"

    # run test plot
    #E_compare(foldername, pars1, pars2, par_nm, par_dg, mod)
    #return 0

    for i in range(len(pars)):
        print("analyzing",pars[i])
        N, imod, Ne, lf, kar, karg, lam, Kd, q, Cn, g = pars[i]
        O_stat_ana(foldername,pars[i],par_nm,par_dg, mode=mod, CnequalsKc=0, tau_c=6)
        if(1):
            pass
            Gij_stat_ana(foldername,pars[i],par_nm,par_dg,mode=mod, tau_c=6)
            rhor_ave_plot(foldername, pars[i], par_nm, par_dg, mode=mod, del_r=0.2, tag="",leg_num=5)

        for g in gs[::1]:
            if(i%1==0):
                pass
                filename = foldername + "/State_N%.0f_imod%.0f_Ne%.0f_lf%.1f_kar%.0f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_g%.1f.csv" % (N, imod,Ne,lf,kar,karg,lam,Kd,q,Cn,g)
                ctag = r"$l_f=%.1f,\bar{\kappa}=%.1f,\lambda=%.1f,g=%.1f$" %(lf,karg,lam,g)
                config_plot_xyz(filename,mesh=1,rod=0,piwall=0,tag=ctag,Format="png")
                #config_plot_xyz(filename,mesh=0,rod=1,tag=ctag,Format="png")

    colors = None
    alphas = None

    Geig_pars_plot(foldername, pars,par_nm,par_dg, mode=mod)
    Os_pars_plot(foldername, pars,par_nm,par_dg,mode=mod)
    # additional test

if __name__ == '__main__':
    main()
