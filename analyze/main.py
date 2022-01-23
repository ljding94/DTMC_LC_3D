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
    # config_plot3D("../data/Ne2/Nov14_2021_1/State_N100_imod1_Ne2_lf0.0_kar15_C00.0_karg0.0_lam2.0_Kd6.0_q2.0_Cn9.0.csv", mesh=1, rod=0, piwall=1, fnormal=0)
    config_plot3D("../data/scratch_local/State_N300_imod3_Ne2_lf0.0_kar50_C00.0_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0_therm.csv",mesh=1,rod=0,piwall=0,fnormal=0)
    #config_plot3D("../data/scratch_local/State_N200_imod3_Ne2_lf0.0_kar10_C00.0_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0_therm.csv",mesh=1,rod=0,piwall=0,fnormal=0)

    return 0

    foldername = "../data/Ne2/Jan22_2022"
    print("analyzing " + foldername)
    N = 300
    Ns = [100, 200, 300, 400, 500]
    imod = 3  # 1 for rhombus, 2 disk, 3 cylinder, 4 for mobius strip
    Ne = 4
    lfs = np.arange(6.0, 35.1, 1.0)
    lf = 0.0
    kars = [10, 40]
    kar = 40
    C0s = [0.1, 0.2, 0.3]
    C0 = 0.0
    karg = 0.0
    kargs = np.arange(0.2, 2.01, 0.2)
    lams = np.arange(2.0, 4.01, 0.4)
    lam = 5.0
    lams = [5.0,10.0]
    Kd = 4.0
    Kds = np.arange(4.0, 8.1, 2.0)
    qs = np.arange(0.0, 3.01, 0.2)
    q = 2.0
    Cns = [3.0, 4.0, 5.0, 7.0, 9.0, 11.0]
    Cn = 4.0
    ku = 0
    pars = []
    pars1, pars2 = [], []
    for Kd in Kds[:]:
        Cn = Kd
        pars.append([N, imod, Ne, lf, kar, C0, karg, lam, Kd, qs, Cn])
        #pars1.append([N, 1, 2, lf, kar, C0, karg, lam, Kd, qs, Cn])
        #pars2.append([N, 3, 2, lf, kar, C0, karg, lam, Kd, qs, Cn])
    par_nm = ["N", "imod", "Ne", "lf", "kar", "C0", "karg", "lam", "Kd", "q", "Cn"]
    par_dg = [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1]  # number of digit for each
    mod = "q"

    # run test plot
    #E_compare(foldername, pars1, pars2, par_nm, par_dg, mod)
    #return 0

    for i in range(len(pars)):
        print("analyzing", pars[i])
        N, imod, Ne, lf, kar, C0, karg, lam, Kd, q, Cn = pars[i]
        O_stat_ana(foldername, pars[i], par_nm, par_dg, mode=mod, thermN=2000,CnequalsKc=0, tau_c=6)
        if 1:
            pass
            Gij_stat_ana(foldername, pars[i], par_nm, par_dg, mode=mod, tau_c=6)
            #rhor_ave_plot(foldername, pars[i], par_nm, par_dg, mode=mod, del_r=0.1, tag="", leg_num=5)

        for q in qs[::1]:
            if i % 1 == 0:
                filename = foldername + "/State_N%.0f_imod%.0f_Ne%.0f_lf%.1f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f.csv" % (N, imod, Ne, lf, kar, C0, karg, lam, Kd, q, Cn)
                ctag = r"$l_f=%.1f,C_0=%.1f,\lambda=%.1f$" % (lf, C0, lam)
                config_plot_xyz(filename[:-4]+"_therm.csv", mesh=1, rod=0, piwall=1, tag=ctag, Format="png")
                config_plot_xyz(filename, mesh=1, rod=0, piwall=1, tag=ctag, Format="png")

                Ofilename = foldername + "/O_MC_N%.0f_imod%.0f_Ne%.0f_lf%.1f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f.csv" % (N, imod, Ne, lf, kar, C0, karg, lam, Kd, q, Cn)
                O_MCstep_plot(Ofilename, Ne)
                #config_plot_xyz(filename,mesh=0,rod=1,tag=ctag,Format="png")

    colors = None
    alphas = None
    Geig_pars_plot(foldername, pars, par_nm, par_dg, mode=mod)
    Os_pars_plot(foldername, pars, par_nm, par_dg, mode=mod)
    # additional test


if __name__ == "__main__":
    main()
