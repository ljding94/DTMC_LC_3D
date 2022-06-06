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
    #config_plot3D("../data/Ne2/May12_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q0.0_Cn4.0.csv", mesh=1, rod=0, piwall=1, fnormal=0)

    #return 0

    # sequence config plotting
    for i in range(1):
        filename = "../data/scratch_local/State_N200_imod1_Ne2_lf0.0_kar30_C00.0_karg0.0_lam5.0_Kd4.0_q3.0_Cn4.0_seq%d.csv" % i
        ctag = "seq=%d" % i
        if os.path.exists(filename):
            pass
            # config_plot_xyz(filename, mesh=1, rod=0, piwall=0, tag=ctag, Format="png")

    # return 0

    foldername = "../data/Ne2/May24_2022"
    print("analyzing " + foldername)
    N = 300
    Ns = [100, 200, 300, 400, 500]
    imod = 3  # 1 for rhombus, 2 disk, 3 cylinder, 4 for mobius strip
    Ne = 2
    lfs = np.arange(8.0, 37.1, 1.0)
    #lfs = np.arange(8.0, 27.1, 1.0)
    # lfs = [15.0,25.0,35.0]
    lf = 0.0

    # kars = [20,30,50,80,120]
    kars = [20, 30, 40, 60, 80, 100]
    kar = 40
    C0 = 0.0
    karg = 0.0
    lams = np.arange(2.0, 4.01, 0.4)
    lam = 6.0
    lams = [2.0, 4.0, 6.0]

    # Kds = np.arange(1.0, 7.01, 0.2)
    # Kds = [1.0,2.0,3.0,4.0,5.0]
    Kds = [2.0,4.0,6.0]
    #Kds = [2.0,4.0,6.0,8.0]
    # Kds = [0.0,2.0,4.0,6.0]
    Kd = 4.0
    qs = np.arange(0.0, 4.01, 0.1)
    qs = [0.0, 1.0, 3.0, 4.0]
    q = 0.0
    Cns = [2.0, 4.0, 6.0, 8.0]
    Cn = 4.0
    pars = []
    pars1, pars2 = [], []
    for Kd in Kds[:]:
        #Cn = Kd
        pars.append([N, imod, Ne, lfs, kar, C0, karg, lam, Kd, q, Cn])
        # pars1.append([N, 1, 2, lf, kar, C0, karg, lam, Kd, qs, Cn])
        # pars2.append([N, 3, 2, lf, kar, C0, karg, lam, Kd, qs, Cn])
    par_nm = ["N", "imod", "Ne", "lf", "kar", "C0", "karg", "lam", "Kd", "q", "Cn"]
    par_dg = [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1]  # number of digit for each
    mod = "lf"

    # run test plot
    # E_compare(foldername, pars1, pars2, par_nm, par_dg, mod)
    # return 0

    for i in range(len(pars)):
        print("analyzing", pars[i])
        N, imod, Ne, lf, kar, C0, karg, lam, Kd, q, Cn = pars[i]
        O_stat_ana(foldername, pars[i], par_nm, par_dg, mode=mod, thermN=5000, CnequalsKc=0, tau_c=6)
        if 1:
            pass
            # Gij_stat_ana(foldername, pars[i], par_nm, par_dg, mode=mod, tau_c=6)
            # rhor_ave_plot(foldername, pars[i], par_nm, par_dg, mode=mod, del_r=0.1, tag="", leg_num=5)

        if 0:
            for lf in lfs[::1]:
                if kar == 20:
                    filename = foldername + "/State_N%.0f_imod%.0f_Ne%.0f_lf%.1f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f.csv" % (N, imod, Ne, lf, kar, C0, karg, lam, Kd, q, Cn)
                    ctag = r"$l_f=%.1f,C_0=%.1f,\lambda=%.1f$" % (lf, C0, lam)
                    # config_plot_xyz(filename[:-4]+"_therm.csv", mesh=1, rod=0, piwall=1, tag=ctag, Format="png")
                    if os.path.exists(filename):
                        pass
                        config_plot_xyz(filename, mesh=1, rod=1, piwall=0, tag=ctag, Format="png")

                    Ofilename = foldername + "/O_MC_N%.0f_imod%.0f_Ne%.0f_lf%.1f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f.csv" % (N, imod, Ne, lf, kar, C0, karg, lam, Kd, q, Cn)
                    if os.path.exists(Ofilename) and lf < 10:
                        pass
                        O_MCstep_plot(Ofilename, 1000, Ne)

    colors = None
    alphas = None
    # Geig_pars_plot(foldername, pars, par_nm, par_dg, mode=mod)
    Os_pars_plot(foldername, pars, par_nm, par_dg, mode=mod)
    # additional test


if __name__ == "__main__":
    main()
