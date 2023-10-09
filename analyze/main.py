#!/usr/local/bin/python3

import time
import numpy as np
import matplotlib.pyplot as plt
from plot import *
from Oplot import *
from analyze import *
from test_plot import *


def main():
    print("hello! dtmc_lc analysis")



    # test tilt twist cal:
    #tilt_twist_calc("../data/Ne2/May12_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q1.0_Cn4.0.csv")
    #tilt_twist_plot()
    #return 0
    #config_plot3D("../data/Ne2/May3_2023/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd2.0_q4.0_Cn10.0_id0.csv", mesh=1, rod=1, piwall=1, fnormal=0)
    #config_plot3D("../data/Ne2/May12_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q2.5_Cn4.0.csv", mesh=1, rod=1, piwall=0, fnormal=0)
    #uz_distribution_plot("../data/Ne2/May12_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q2.5_Cn4.0.csv")
    #tan_alpha_distribution_plot("../data/Ne2/May12_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q1.0_Cn4.0.csv")
    #return 0

    tanalpha = []
    alpha_u = []
    qs = np.arange(0.0,3.01,0.1)

    #qs = [2.3]
    if (0):
        for q in qs:
        #for q in [1.9]:
            print("q=",q)
            #filename = "../data/Ne2/Dec7_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q%.1f_Cn2.0_id0.csv"%q
            #filename = "../data/Ne2/Apr30_2023/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q%.1f_Cn2.0_id0.csv"%q
            #filename = "../data/Ne2/data_2022/May12_2022/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q%.1f_Cn4.0.csv"%q
            #filename = "../data/Ne2/Jun29_2023/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q%.1f_Cn2.0.csv"%q
            filename = "../data/Ne2/Jul1_2023/data_use/State_N300_imod3_Ne2_lf25.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q%.1f_Cn10.0_id1.csv"%q
            m = 2
            if q>2.3:
                m=3
            ans = tilt_slice_distri_plot(filename,m)
            tanalpha.append(ans)
            ans_u = wall_director_alpha_calc(filename,pwlim = np.pi/3)
            alpha_u.append(ans_u)
        print("alpha_u",alpha_u, "\n")
        #alpha_u_plot()
        print("tanalpha",tanalpha,"\n")
        #tanalpha_q_plot()
        return 0


    # sequence config plotting
    for i in range(1):
        filename = "../data/scratch_local/State_N200_imod1_Ne2_lf0.0_kar30_C00.0_karg0.0_lam5.0_Kd4.0_q3.0_Cn4.0_seq%d.csv" % i
        ctag = "seq=%d" % i
        if os.path.exists(filename):
            pass
            # config_plot_xyz(filename, mesh=1, rod=0, piwall=0, tag=ctag, Format="png")

    # return 0

    foldername = "../data/Ne2/Oct8_2023"
    print("analyzing " + foldername)
    N = 300
    #Ns = [100, 200, 300, 400, 500]
    imod = 3  # 1 for rhombus, 2 disk, 3 cylinder, 4 for mobius strip
    Ne = 2
    lf = 25.0
    lfs = [15.0,20.0,25.0,30.0,35.0]

    kar = 50
    C0 = 0.0
    karg = 0.0

    lam = 6.0
    Kd = 2.0
    Kds = [2.0]
    qs = np.arange(0.3, 6.01, 0.3)
    #Cns = np.arange(15.0,34.1,1.0)

    qs = np.arange(0.1, 3.01, 0.1)
    Cns = np.arange(1.0,13.1,1.0)
    Cns = [2,4,6,8]
    Cn = 4.0
    lfs = [15,25,35]
    Kds = np.arange(1.0,7.01,0.2)
    q = 0
    #ids = [0,1,2,3,4,5]
    id = 0
    #Cns = np.arange(6.0,12.1,2.0)
    pars = []
    pars1, pars2 = [], []
    #for lf in lfs[:]:
    for Cn in Cns[:]:
        #Cn = Kd
        pars.append([N, imod, Ne, lf, kar, C0, karg, lam, Kds, q, Cn, id])
        # pars1.append([N, 1, 2, lf, kar, C0, karg, lam, Kd, qs, Cn])
        # pars2.append([N, 3, 2, lf, kar, C0, karg, lam, Kd, qs, Cn])
    par_nm = ["N", "imod", "Ne", "lf", "kar", "C0", "karg", "lam", "Kd", "q", "Cn", "id"]
    par_dg = [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0]  # number of digit for each
    mod = "Kd"

    #test Q ana
    Qij_stat_ana("../data/scratch_local", [101,1,1,0.0,50,0.0,0.0,6.0,10.0,[0.0],4.0,0], par_nm, par_dg, "q", 1,tau_c=6)
    #config_plot_xyz("../data/scratch_local/"+"State_N91_imod2_Ne2_lf0.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q0.0_Cn4.0_id0_init.csv", mesh=1, rod=1, piwall=0, tag=ctag, Format="png")
    #config_plot_xyz("../data/scratch_local/"+"State_N91_imod2_Ne2_lf0.0_kar50_C00.0_karg0.0_lam6.0_Kd4.0_q0.0_Cn4.0_id0.csv", mesh=1, rod=1, piwall=0, tag=ctag, Format="png")
    return 0

    # run test plot
    # E_compare(foldername, pars1, pars2, par_nm, par_dg, mod)
    #return 0

    for i in range(len(pars)):
        print("analyzing", pars[i])
        N, imod, Ne, lf, kar, C0, karg, lam, Kd, q, Cn, id = pars[i]
        #O_stat_ana(foldername, pars[i], par_nm, par_dg, mode=mod, thermN=5000, CnequalsKc=0, tau_c=6)
        #O_stat_ana(foldername, pars[i], par_nm, par_dg,mode=mod, thermN=5000, CnequalsKc=0, tau_c=6)
        #Os_pars_plot(foldername, [pars[i]], par_nm, par_dg, mode=mod, subfix="Cn%d"%Cn)
        #O_stat_ana_multi(foldername, pars[i], par_nm, par_dg, mode=mod, nmulti=13, thermN=5000, CnequalsKc=0, tau_c=6)

        if 1:
            Qij_stat_ana(foldername, pars[i], par_nm, par_dg, mode=mod, thermN=100, tau_c=6)
            #Gij_stat_ana(foldername, pars[i], par_nm, par_dg, mode=mod, tau_c=6)
            # rhor_ave_plot(foldername, pars[i], par_nm, par_dg, mode=mod, del_r=0.1, tag="", leg_num=5)

        if 0:
            for q in qs[::2]:
                if 0:
                    print("plotting O vs MC steps")
                    filename = foldername + "/O_MC_N%.0f_imod%.0f_Ne%.0f_lf%.1f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_id%d.csv" %(N, imod, Ne, lf, kar, C0, karg, lam, Kd, q, Cn, id)
                    if os.path.exists(filename):
                        O_MCstep_plot(filename, 0, Ne, 1000)

                if 0:
                    filename = foldername + "/State_N%.0f_imod%.0f_Ne%.0f_lf%.1f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_id%d.csv" % (N, imod, Ne, lf, kar, C0, karg, lam, Kd, q, Cn, id)
                    ctag = r"$l_f=%.1f,C_0=%.1f,\lambda=%.1f,q=%.1f, Kd=%.1f$" % (lf, C0, lam, q, Kd)
                    #config_plot_xyz(filename[:-4]+"_therm.csv", mesh=1, rod=0, piwall=1, tag=ctag, Format="png")
                    if os.path.exists(filename):
                        pass
                        config_plot_xyz(filename, mesh=1, rod=0, piwall=1, tag=ctag, Format="png")

                if 0:
                    print("plotting un2 distribution")
                    un2filename = foldername + "/un2dis_N%.0f_imod%.0f_Ne%.0f_lf%.1f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_Kd%.0f_q%.1f_Cn%.1f_id0.csv" % (N, imod, Ne, lf, kar, C0, karg, lam, Kd, q, Cn)
                    ctag = r"$l_f=%.1f,q=%.1f,C=%.1f$" % (lf, q, Cn)
                    if os.path.exists(un2filename):
                        pass
                        un2_distribution_plot(un2filename,ctag)

                if 0:
                    print("plotting un2theta distribution")
                    un2thetafilename = foldername + "/un2thetadis_N%.0f_imod%.0f_Ne%.0f_lf%.1f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_id0.csv" % (N, imod, Ne, lf, kar, C0, karg, lam, Kd, q, Cn)
                    ctag = r"$l_f=%.1f,q=%.1f,C=%.1f$" % (lf, q, Cn)
                    if os.path.exists(un2thetafilename):
                        pass
                        un2theta_distribution_plot(un2thetafilename,ctag)

                if 0:
                    print("plotting dA2H2 distribution")
                    dA2H2filename = foldername + "/dA2H2dis_N%.0f_imod%.0f_Ne%.0f_lf%.1f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_Kd%.0f_q%.1f_Cn%.1f_id0.csv" % (N, imod, Ne, lf, kar, C0, karg, lam, Kd, q, Cn)
                    ctag = r"$l_f=%.1f,q=%.1f,C=%.1f$" % (lf, q, Cn)
                    if os.path.exists(dA2H2filename):
                        pass
                        dA2H2_distribution_plot(dA2H2filename,ctag)
                    else:
                        print(dA2H2filename + "doen'st exist")

                if 0:
                    print("plotting 2H distribution")
                    twoHfilename = foldername + "/2Hdis_N%.0f_imod%.0f_Ne%.0f_lf%.1f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_id%d.csv" % (N, imod, Ne, lf, kar, C0, karg, lam, Kd, q, Cn, id)
                    ctag = r"$l_f=%.1f,Kd=%.1f,q=%.1f,C=%.1f$" % (lf, Kd, q, Cn)
                    if os.path.exists(twoHfilename):
                        pass
                        twoH_distribution_plot(twoHfilename,ctag)
                    else:
                        print(twoHfilename + "doen'st exist")
                if 0:
                    print("plotting dA distribution")
                    dAfilename = foldername + "/dAdis_N%.0f_imod%.0f_Ne%.0f_lf%.1f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f_id%d.csv" % (N, imod, Ne, lf, kar, C0, karg, lam, Kd, q, Cn, id)
                    ctag = r"$l_f=%.1f,Kd=%.1f,q=%.1f,C=%.1f$" % (lf, Kd, q, Cn)
                    if os.path.exists(dAfilename):
                        pass
                        dA_distribution_plot(dAfilename,ctag)
                    else:
                        print(dAfilename + "doen'st exist")

                if 0:
                    Ofilename = foldername + "/O_MC_N%.0f_imod%.0f_Ne%.0f_lf%.1f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f.csv" % (N, imod, Ne, lf, kar, C0, karg, lam, Kd, q, Cn)
                    if os.path.exists(Ofilename) and lsf < 10:
                        pass
                        #O_MCstep_plot(Ofilename, 1000, Ne)

    colors = None
    alphas = None
    #Geig_pars_plot(foldername, pars, par_nm, par_dg, mode=mod)
    Qeig_pars_plot(foldername, pars, par_nm, par_dg, mode=mod)
    Os_pars_plot(foldername, pars, par_nm, par_dg, mode=mod)
    # additional test


if __name__ == "__main__":
    main()
