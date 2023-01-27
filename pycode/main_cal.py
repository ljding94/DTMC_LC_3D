#!/usr/local/bin/python3
from cholesteric_cylinder_cal import *
from numeric_2mod_io import *
import time
import sys

def main():
    #print("python code for numerical calculation cholesteric on cylinder")
    #find_alpha()
    #find_alpha_opt()
    if(len(sys.argv)>2):
        # on ccv
        folder = "/users/lding3/scratch/dtmc_lc_3d"
        K= float(sys.argv[1])
        C = int(sys.argv[2])
        m = float(sys.argv[3])
        bn_phi = int(sys.argv[4])
        bn_z = int(sys.argv[5])
        R = 1
        qs = np.arange(0.00,3.01,0.05)
        Ftot_qs_cal(folder, K, C, qs, m, bn_phi, bn_z, R, "Powell")
    else:
        print("python code for numerical calculation cholesteric (2 mod) on cylinder")
        #Ftot_par_run(par, bn_phi, bn_z, method)
        #K, C, q, m = par
        #Ftot_par_run((1,1,0.5,2), 50, 50, "Powell")


        #print(intSS_unit_length(m, alpha, gamma, bn_phi, bn_z))
        m, alpha, gamma, bn_phi, bn_z = 2, 0.0, 0, 50, 50
        R = 1
        print("SS: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intSS_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("TT: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intTT_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("T: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intT_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("BB: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intBB_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("C: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intC_unit_length(m, alpha, gamma, bn_phi, bn_z,R))

        gamma=1
        print("for gamma = 1")
        print("SS: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intSS_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("TT: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intTT_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("T: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intT_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("BB: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intBB_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("C: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intC_unit_length(m, alpha, gamma, bn_phi, bn_z,R))

        gamma = 0
        R=2
        print("for R = 2")
        print("SS: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intSS_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("TT: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intTT_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("T: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intT_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("BB: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intBB_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("C: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intC_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        '''
        print("C=",intC_unit_length(2, 0, 0, 50, 50))
        print("C=",intC_unit_length(2, 0, 1, 50, 50))
        print("C=",intC_unit_length(3, 0, 0, 50, 50))
        print("C=",intC_unit_length(3, 0, 1, 50, 50))
        '''
        m = 0
        print("for m = 0")
        print("SS: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intSS_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("TT: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intTT_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("T: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intT_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("BB: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intBB_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("C: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intC_unit_length(m, alpha, gamma, bn_phi, bn_z,R))

        R = 2
        print("for m = 0, R=2")
        print("SS: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intSS_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("TT: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intTT_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("T: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intT_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("BB: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intBB_unit_length(m, alpha, gamma, bn_phi, bn_z,R))
        print("C: (m=%.0f,alpha=%.1f,gamma=%.1f,)"%(m,alpha,gamma),intC_unit_length(m, alpha, gamma, bn_phi, bn_z,R))

        bn_phi,bn_z = 200, 200
        #E_qs_cal(folder, K, C, qs, m, bn_phi, bn_z, method):

        #Ftot_qs_cal("./pydata", 10, 1, qs, 2, 50, 50, "Powell")
        #
        #qs = np.arange(0.01,0.1,0.002)
        Cs = np.arange(1,6.1,1)
        qs = np.arange(0.0,3.1,0.2)
        m = 0
        R = 1
        K = 1
        alpha = 0
        gamma = 0

        for m in [0,2]:
            print("m=",m)
            for C in Cs:
                print("C=",C)
                for q in qs:
                    pass
                    #print("R=%f,q=%f \n"%(R,q))
                    #print("Ftot_unit_length/R= ",Ftot_unit_length(K, C, q, m, alpha, gamma, bn_phi, bn_z, R))
                Ftot_qs_cal("../data/pydata/local/Jan26_2023_R1", K, C, qs, m, bn_phi, bn_z, R, "Powell")


if __name__ == "__main__":
    main()

