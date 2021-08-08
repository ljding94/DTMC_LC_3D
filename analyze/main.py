#!/opt/homebrew/bin/python3

import time
import numpy as np
import matplotlib.pyplot as plt
from plot import *
from Oplot import *
from analyze import *


def main():
    print("hello! dtmc_lc analysis")

    config_plot3D("../data/scratch_local/State_N200_imod1_Ne2_kar100_C00.0_karg0.0_lam5.0_Kd0.0_q0.0_Cn0.0.csv",mesh=1,rod=0)
    #config_plot3D("../data/Ne2/Aug4_2021_1/State_N300_imod3_Ne2_kar50_C00.0_karg0.0_lam8.0_Kd10.0_q2.0_Cn20.0.csv",mesh=1,rod=0)

    return 0

    #foldername = "../data/Ne2/Apr29_2021" # kinda switch to Ne2 simulation for next paper
    #foldername = "../data/Ne1/May4_2021" # mobius strip~,
    foldername = "../data/Ne2/Aug4_2021_1"
    print("analyzing "+foldername)
    N = 300
    imod=3 # 1 for rhombus, 2 disk, 3 cylinder, 4 for mobius strip
    Ne=2
    kar = 50
    C0 = 0.0
    karg = 0.0
    lam = 8.0
    Kds = np.arange(6.0,10.1,1.0)
    qs = np.arange(0.2,2.1,0.2)
    #qs=np.arange(0.2,1.1,0.2)
    Cn =5.0
    pars = []
    for Kd in Kds:
        Cn = 2*Kd
        pars.append([N, imod, Ne, kar, C0, karg, lam, Kd, qs, Cn])
    par_nm = ["N","imod", "Ne", "kar", "C0","karg", "lam","Kd","q", "Cn"]
    par_dg = [0,0,0,0,1,1,1,1,1,1] # number of digit for each
    mod="q"
    for i in range(len(pars)):
        print("analyzing",pars[i])
        N, imod, Ne, kar, C0, karg, lam, Kd, q, Cn = pars[i]
        O_stat_ana(foldername,pars[i],par_nm,par_dg, mode=mod, CnequalsKc=0, tau_c=6)
        if(1):
            pass
            Gij_stat_ana(foldername,pars[i],par_nm,par_dg,mode=mod, tau_c=6)

        for q in qs[::1]:
            if(i%1==0):
                pass
                filename = foldername + "/State_N%.0f_imod%.0f_Ne%.0f_kar%.0f_C0%.1f_karg%.1f_lam%.1f_Kd%.1f_q%.1f_Cn%.1f.csv" % (N, imod,Ne,kar,C0,karg,lam,Kd,q,Cn)
                ctag = r"$K_d=%.1f,q=%.1f$" %(Kd,q)
                config_plot_xyz(filename,mesh=0,rod=1,tag=ctag,Format="png")
                config_plot_xyz(filename,mesh=1,rod=0,tag=ctag,Format="png")

    colors = None
    alphas = None

    Geig_pars_plot(foldername, pars,par_nm,par_dg, mode=mod)
    Os_pars_plot(foldername, pars,par_nm,par_dg,mode=mod)
    '''
    E_imod3=[[-4606.085 ,-4647.325205 ,-4644.769535 ,-4700.72541  ,-4698.958655,-4725.382595 ,-4765.38003  ,-4888.81189  ,-4972.98778  ,-5121.578275],[-5604.926245 ,-5655.661895 ,-5650.947515 ,-5678.887585 ,-5690.480095,-5685.379475 ,-5703.39637  ,-5836.721035 ,-5837.82122  ,-5929.458295],[-6561.96697  ,-6555.43771  ,-6643.53209  ,-6648.061655 ,-6649.78069,-6672.08846  ,-6710.733325 ,-6751.80299  ,-6450.921065 ,-6836.391095],[-7532.169695 ,-7534.372125 ,-7515.95098  ,-7455.09314  ,-7291.337255,-7603.09252  ,-7723.733425 ,-7718.725995 ,-7758.155295 ,-7823.584525],[-8530.323535 ,-8552.12995  ,-8578.01737  ,-8355.815315 ,-8411.933615,-8625.944025 ,-8504.524975 ,-8683.217    ,-8529.88371  ,-8259.345755]]
    print(E_imod3)
    E_imod1=[[-4588.418815, -4582.10823 , -4573.67018 , -4575.93939 , -4582.125405, -4704.31282 , -4717.63663 , -4841.98612 , -4943.28134 , -5125.876185],[-5498.52288 , -5511.715175, -5497.737305, -5494.51312 , -5484.942885, -5597.180385, -5541.3707  , -5579.209585, -5774.83976 , -5895.537035],[-6423.88423 , -6424.71934 , -6404.916825, -6414.56019 , -6428.70782, -6425.968365, -6535.61825 , -6481.4472  , -6551.668025, -6746.7734  ],[-7343.1547  , -7336.488365, -7333.746245, -7336.00785 , -7337.12866 ,-7367.95394 , -7414.20379 , -7381.28256 , -7470.60073 , -7527.70312 ] ,[-8258.326575, -8248.327985, -8253.014915, -8255.37931 , -8261.47933, -8268.9012  , -8282.048735, -8292.111465, -8515.325225, -8433.62854 ]]
    print(E_imod1)
    print("E_imod3-E_imod1",np.array(E_imod3)-np.array(E_imod1))
    '''
    # additional test

if __name__ == '__main__':
    main()

##!/usr/local/bin/python3