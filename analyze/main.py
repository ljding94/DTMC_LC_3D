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
    config_plot3D("../data/scratch_local/State_N200_imod3_Ne2_lf15.0_kar50_C02.0_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0_init.csv",mesh=1,rod=0,fnormal=0)
    config_plot3D("../data/scratch_local/State_N200_imod3_Ne2_lf15.0_kar50_C02.0_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0.csv",mesh=1,rod=0,fnormal=0)
    # config_plot3D("../data/Ne2/Oct12_2021/State_N400_imod3_Ne2_lf5.0_kar100_karg0.0_lam5.0_Kd7.0_q0.9_Cn7.0_g0.0.csv",mesh=1,rod=1,fnormal=0,piwall=0)
    return 0

    foldername = "../data/Ne2/Oct13_2021"
    print("analyzing " + foldername)
    N = 400
    imod = 3  # 1 for rhombus, 2 disk, 3 cylinder, 4 for mobius strip
    Ne = 2
    lfs = np.arange(5.0, 14.1, 1.0)
    # lf = 0.0
    kar = 80
    karg = 0.0
    # lams = np.arange(5.0,10.1,1.0)
    lam = 5.0
    Kd = 6.0
    qs = np.arange(0.0, 1.4, 0.2)
    # q = 0.0
    # Cns = [5.0,10.0,15.0,20.0]
    Cn = 6.0
    # gs = np.arange(0.2,2.01,0.2)
    g = 0.0
    pars = []
    for q in qs:
        Cn = Kd
        pars.append([N, imod, Ne, lfs, kar, karg, lam, Kd, q, Cn, g])
    par_nm = ["N", "imod", "Ne", "lf", "kar", "karg", "lam", "Kd", "q", "Cn", "g"]
    par_dg = [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1]  # number of digit for each
    mod = "lf"

    # run test plot
    # E_compare(foldername, pars1, pars2, par_nm, par_dg, mod)
    # return 0

    for i in range(len(pars)):
        print("analyzing", pars[i])
        N, imod, Ne, lf, kar, karg, lam, Kd, q, Cn, g = pars[i]
        O_stat_ana(foldername, pars[i], par_nm, par_dg, mode=mod, CnequalsKc=0, tau_c=6)
        if 1:
            pass
            Gij_stat_ana(foldername, pars[i], par_nm, par_dg, mode=mod, tau_c=6)
            rhor_ave_plot(foldername, pars[i], par_nm, par_dg, mode=mod, del_r=0.1, tag="", leg_num=5)

        for lf in lfs[::1]:
            if i % 1 == 0:
                filename = foldername + "/State_N%.0f_imod%.0f_Ne%.0f_lf%.1f_kar%.0f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_g%.1f.csv" % (N, imod, Ne, lf, kar, karg, lam, Kd, q, Cn, g)
                ctag = r"$l_f=%.1f,\bar{\kappa}=%.1f,\lambda=%.1f,g=%.1f$" % (lf, karg, lam, g)
                config_plot_xyz(filename, mesh=1, rod=0, piwall=1, tag=ctag, Format="png")
                # config_plot_xyz(filename,mesh=0,rod=1,tag=ctag,Format="png")

    colors = None
    alphas = None
    Geig_pars_plot(foldername, pars, par_nm, par_dg, mode=mod)
    Os_pars_plot(foldername, pars, par_nm, par_dg, mode=mod)
    # additional test


if __name__ == "__main__":
    main()
