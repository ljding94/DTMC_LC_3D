#!#!/usr/local/bin/python3
#from plot import *
#from catelinder_cal import *
#from catelinder_plot import *
#from unduloid_cal import *
#from unduloid_plot import *
#from catenoid_cal import *

from numeric_2mod_io import *
def main():
    print("python code for numerical calculation of pulling experiment o/n unduloid")
    #cylinder_2mod_plot("./pydata/optFtot_K1.0_C10.0_m2_qs.csv")
    #cylinder_2mod_plot("./pydata/optFtot_K10.0_C1.0_m2_qs.csv")
    F_compot_param_plot("../data/pydata/local/Mar18_2023_params_test/F_compot_m2_R1.0_alphas_gammas.csv",2)
    F_compot_param_plot("../data/pydata/local/Mar18_2023_params_test/F_compot_m1_R1.0_alphas_gammas.csv",1)
    F_compot_param_plot("../data/pydata/local/Mar18_2023_params_test/F_compot_m0_R1.0_alphas_gammas.csv",0)
    return 0
    K = 2 #[0.1,0.5,2.0,4.0]
    qs = np.arange(1.0,3.5,0.1)
    m = 2
    R = 5
    Rs =np.arange(2,10.1,1)
    C = 6
    Cs =np.arange(2,10.1,1)


    #m,alpha,gamma,bn_phi,bn_z=0,0,1,100,100
    #director_field_plot(m, alpha, gamma, bn_phi, bn_z)
    #return 0

    pars = []
    for R in Rs:
        pars.append([K,C,m,R])
    foldername = "./pydata/Nov30_2022"
    #cylinder_2mod_plot(foldername, pars)
    #cylinder_2mod_plot_normlq(foldername, pars, normqlim=0.5)
    #cylinder_c_q_hist(foldername, pars)
    #cylinder_critical_q_C(foldername, pars)

    #director_field_plot(2, 0.5, 0.8, 100, 30)
    #director_field_plot(2, 0, 0.5, 100, 30)
    if(0):
        Ks = np.arange(0.02,0.421,0.1)
        Ks = np.arange(0.4,1.01,0.1)
        Cs = np.arange(1,6.1,1)
        K = 1
        del_Ftot_Ks_qs_plot("../data/pydata/local/Jan26_2023_R1",K,Cs,1)
        return 0

    #Ks = np.arange(0.01,0.101,0.01)
    #qs = np.arange(0.0,8.1,0.5)
    #Ks = np.arange(0.02,0.401,0.02)
    Ks = np.arange(0.01,0.801,0.01)
    #qs = np.arange(0.0,8.1,0.4)
    K = 1
    C, R = 1,1
    Cs = np.arange(0.5,10.1,0.5)

    del_Ftot_Ks_qs_plot("../data/pydata/Jan26_2023",K,Cs,1)

if __name__ == "__main__":
    main()

