import numpy as np
import matplotlib.pyplot as plt
from plot import *
from autocorrelation import *
from scipy.optimize import curve_fit
from scipy import odr
import scipy.fft
import os


def find_cpar_ind(par_nm, mode):
    cpar_ind = -1
    for i in range(len(par_nm)):
        if par_nm[i] == mode:
            cpar_ind = i
            break
    return cpar_ind


def O_stat_cal_weight(O, tau_c, wnu=1, wnu_ave=1, wnu_err=0):
    # for each O, return O_ave,O_tau and O_err using O and wnu of bias
    if wnu_ave != 1:
        # with bias wnu
        Obias_ave = np.average(O * wnu)
        O_ave = Obias_ave / wnu_ave
        rho, cov0 = autocorrelation_function_fft(O * wnu)
        O_tau, tau_err = tau_int_cal_rho(rho, tau_c)
        Obias_err = np.sqrt(2 * O_tau / len(O) * cov0)
        O_err = O_ave * np.sqrt(np.power(Obias_err / Obias_ave, 2) + np.power(wnu_err / wnu_ave, 2))
    else:
        O_ave = np.average(O)
        rho, cov0 = autocorrelation_function_fft(O)
        O_tau, tau_err = tau_int_cal_rho(rho, tau_c)
        O_err = np.sqrt(2 * O_tau / len(O) * cov0)
    return (O_ave, O_tau, O_err)


def O_stat_cal(O, tau_c):
    O_ave = np.average(O)
    rho, cov0 = autocorrelation_function_fft(O)
    O_tau, tau_err = tau_int_cal_rho(rho, tau_c)
    O_err = np.sqrt(2 * O_tau / len(O) * cov0)
    return (O_ave, O_tau, O_err)

def O_stat_cal_sq(O, tau_c):
    O_ave = np.average(O)
    rho, cov0 = autocorrelation_function_fft(O)
    O_tau, tau_err = tau_int_cal_rho(rho, tau_c)
    O_err_sq = 2 * O_tau / len(O) * cov0
    return (O_ave, O_tau, O_err_sq)



def eff_Le_sym_cal(Les):
    # calculate the assymetry of edge length
    # for Ne=2: Asym = |(L0-L1)/sum(Ls),0|
    # Ne=3: Asym = |L0+L1*e^{i/3}+L2*e^{i*2/3},0|
    # Ne=4 Asym = symmetric direction in 3d, weighted by Les
    # Les = np.array(Les)
    Ne = len(Les)
    if Ne == 1:
        asym = 1
    elif Ne == 2:
        asym = np.abs(Les[0] - Les[1]) / (Les[0] + Les[1])
    elif Ne == 3:
        asymx = Les[0] + np.cos(np.pi * 2 / 3) * (Les[1] + Les[2])
        asymx /= Les[0] + Les[1] + Les[2]
        asymy = np.sin(np.pi * 2 / 3) * (Les[1] - Les[2])
        asymy /= Les[0] + Les[1] + Les[2]
        asym = np.sqrt(np.power(asymx, 2) + np.power(asymy, 2))
    elif Ne == 4:
        vs = []
        vs.append([np.sqrt(8.0 / 9), 0, -1.0 / 3])
        vs.append([-np.sqrt(2.0 / 9), np.sqrt(2.0 / 3), -1.0 / 3])
        vs.append([-np.sqrt(2.0 / 9), -np.sqrt(2.0 / 3), -1.0 / 3])
        vs.append([0, 0, 1.0])
        vs = np.array(vs)
        vasym = [0, 0, 0]
        Lesum = Les[0] + Les[1] + Les[2] + Les[3]

        for k in range(3):
            for e in range(Ne):
                vasym[k] += Les[e] * vs[e][k]
            vasym[k] /= Lesum
        asym = np.sqrt(vasym[0] * vasym[0] + vasym[1] * vasym[1] + vasym[2] * vasym[2])

    return asym


def O_stat_ana(foldername, par, par_nm, par_dg, mode, thermN=0, CnequalsKc=0, tau_c=6):
    cpar_valid = []
    E_ave, E_tau, E_err = [], [], []
    Ne = par[find_cpar_ind(par_nm, "Ne")]

    Les_ave, Les_tau, Les_err = [[] for i in range(Ne)], [[] for i in range(Ne)], [[] for i in range(Ne)]
    Lasym_ave, Lasym_tau, Lasym_err = [], [], []
    IdA_ave, IdA_tau, IdA_err = [], [], []
    I2H_ave, I2H_tau, I2H_err = [], [], []
    I2H2_ave, I2H2_tau, I2H2_err = [], [], []
    I2H2dis_ave, I2H2dis_tau, I2H2dis_err = [], [], []
    IK_ave, IK_tau, IK_err = [], [], []
    p2uu_ave, p2uu_tau, p2uu_err = [], [], []
    uuc_ave, uuc_tau, uuc_err = [], [], []
    # may need to add uuc2 to study the spontaneous symmetry breaking
    un2_ave, un2_tau, un2_err = [], [], []
    un2p_ave, un2p_tau, un2p_err = [], [], []
    uz2_ave, uz2_tau, uz2_err = [], [], []
    uz_ave, uz_tau, uz_err = [], [], [] # abs(uz)
    lb_ave, lb_tau, lb_err = [], [], []  # average length of bond
    #Eu_ave, Eu_tau, Eu_err = [], [], []  # average length of bond
    #wnu_ave, wnu_tau, wnu_err = [], [], []  # <wu^-1> = <1/wu>=<exp(Eu)>

    cpar_ind = find_cpar_ind(par_nm, mode)
    cpar = par[cpar_ind]
    Cn_ind = -1
    if CnequalsKc and mode == "Kd":
        Cn_ind = find_cpar_ind(par_nm, "Cn")

    Ncpar_valid = len(cpar)
    for i in range(len(cpar)):
        par_dealing = par[:]
        par_dealing[cpar_ind] = par[cpar_ind][i]
        if CnequalsKc and mode == "Kd":
            par_dealing[Cn_ind] = par[cpar_ind][i]  # adjust the ratio accordingly
        f2rtail = "MC"
        for j in range(len(par_dealing)):
            f2rtail += "_" + par_nm[j] + "%.*f" % (par_dg[j], par_dealing[j])
        f2rtail += "_id0.csv"
        # print("f2rtail",f2rtail)
        file2read = foldername + "/O_" + f2rtail
        print("file2read", file2read)
        # check for file existence
        if not os.path.exists(file2read):
            print(file2read, "not exist")
            continue
        data = np.loadtxt(file2read, skiprows=14 + thermN, delimiter=",", unpack=True)
        cpar_valid.append(cpar[i])
        N = par_dealing[0]
        E = data[0] / N
        Les = data[1 : 1 + Ne]
        IdA, I2H, I2H2, I2H2dis, IK, Tp2uu, Tuuc, Bond_num, Tun2, Tuz2, Tuz_abs, Tlb = data[1 + Ne : 13 + Ne]
        Eu = data[12 + Ne]
        Lasym = eff_Le_sym_cal(Les)
        p2uu = Tp2uu / Bond_num
        uuc = Tuuc / Bond_num
        lb = Tlb / Bond_num
        N = par[find_cpar_ind(par_nm, "N")]
        print("N=", N)
        # rCnp = par[find_cpar_ind(par_nm,"rCnp")]
        # Np = int(N*rCnp)
        # un2=Tun2/(N-Np)
        un2 = Tun2 / N
        print(file2read, "Tun2[-1]=%.1f,Tun2=%.1f,un2=%.1f" % (Tun2[-1], np.average(Tun2), np.average(un2)))
        uz2 = Tuz2 / N
        uz = Tuz_abs/N

        #Euave = np.average(Eu)
        #wnu = np.exp(Eu - Euave)  # biased weight function, assumed beta=1, normalize to Eu_ave
        # un2p=Tun2p/Np if Np>0 else Tun2p
        # Ne2 case, need Ledif for additional info
        # was used for checking energy
        # again, it's correct
        """
        print("energy slicing E",E)
        print("cpar",cpar[i])
        print("Les[0]+Les[1]",Les[0]+Les[1])
        Et = E-0.5*par[find_cpar_ind(par_nm,"kar")]*I2H2
        print("E-0.5kar*I2H2",Et)
        Et = Et-par[find_cpar_ind(par_nm,"lam")]*(Les[0]+Les[1])
        print("E-0.5kar*I2H2-lam*L",Et)
        Et=Et+par[find_cpar_ind(par_nm,"Kd")]*Tp2uu
        print("E-0.5kar*I2H2-lam*L+Kd*Tp2uu",Et)
        Et=Et+par[find_cpar_ind(par_nm,"Kd")]*par[find_cpar_ind(par_nm,"q")][i]*Tuuc
        print("E-0.5kar*I2H2-lam*L+Kd*Tp2uu+Kd*q*Tuuc",Et)
        Et=Et+0.5*par[find_cpar_ind(par_nm,"Cn")]*(Tun2-N)
        print("E-0.5kar*I2H2-lam*L+Kd*Tp2uu+Kd*q*Tuuc+0.5Cn*(Tun2-N)",Et)
        """

        # wnu = 1 when Enu=0

        # wnu_avei, wnu_taui, wnu_erri = O_stat_cal(wnu, tau_c)
        # wnu_ave.append(wnu_avei)
        # wnu_tau.append(wnu_taui)
        # wnu_err.append(wnu_erri)
        # print("wnu_avei, wnu_taui, wnu_erri", wnu_avei, wnu_taui, wnu_erri)

        # E
        # E_avei, E_taui, E_erri = O_stat_cal_weight(E, tau_c, wnu, wnu_avei, wnu_erri)
        E_avei, E_taui, E_erri = O_stat_cal(E, tau_c)
        E_ave.append(E_avei)
        E_tau.append(E_taui)
        E_err.append(E_erri)

        # Le and Ik2s
        for e in range(Ne):
            # Le_avei, Le_taui, Le_erri = O_stat_cal(Les[e], tau_c, wnu, wnu_avei, wnu_erri)
            Le_avei, Le_taui, Le_erri = O_stat_cal(Les[e], tau_c)
            Les_ave[e].append(Le_avei)
            Les_tau[e].append(Le_taui)
            Les_err[e].append(Le_erri)

        # Lasym_avei, Lasym_taui, Lasym_erri = O_stat_cal(Lasym, tau_c, wnu, wnu_avei, wnu_erri)
        Lasym_avei, Lasym_taui, Lasym_erri = O_stat_cal(Lasym, tau_c)
        Lasym_ave.append(Lasym_avei)
        Lasym_tau.append(Lasym_taui)
        Lasym_err.append(Lasym_erri)

        # IdA
        IdA_avei, IdA_taui, IdA_erri = O_stat_cal(IdA, tau_c)
        IdA_ave.append(IdA_avei)
        IdA_tau.append(IdA_taui)
        IdA_err.append(IdA_erri)

        # I2H
        I2H_avei, I2H_taui, I2H_erri = O_stat_cal(I2H, tau_c)
        I2H_ave.append(I2H_avei)
        I2H_tau.append(I2H_taui)
        I2H_err.append(I2H_erri)

        # I2H2
        I2H2_avei, I2H2_taui, I2H2_erri = O_stat_cal(I2H2, tau_c)
        I2H2_ave.append(I2H2_avei)
        I2H2_tau.append(I2H2_taui)
        I2H2_err.append(I2H2_erri)

        # I2H2dis
        I2H2dis_avei, I2H2dis_taui, I2H2dis_erri = O_stat_cal(I2H2dis, tau_c)
        I2H2dis_ave.append(I2H2dis_avei)
        I2H2dis_tau.append(I2H2dis_taui)
        I2H2dis_err.append(I2H2dis_erri)
        # IK
        IK_avei, IK_taui, IK_erri = O_stat_cal(IK, tau_c)
        IK_ave.append(IK_avei)
        IK_tau.append(IK_taui)
        IK_err.append(IK_erri)

        # p2uu
        p2uu_avei, p2uu_taui, p2uu_erri = O_stat_cal(p2uu, tau_c)
        p2uu_ave.append(p2uu_avei)
        p2uu_tau.append(p2uu_taui)
        p2uu_err.append(p2uu_erri)

        # uuc
        uuc_avei, uuc_taui, uuc_erri = O_stat_cal(uuc, tau_c)
        uuc_ave.append(uuc_avei)
        uuc_tau.append(uuc_taui)
        uuc_err.append(uuc_erri)

        # un2
        un2_avei, un2_taui, un2_erri = O_stat_cal(un2, tau_c)
        un2_ave.append(un2_avei)
        un2_tau.append(un2_taui)
        un2_err.append(un2_erri)

        # uz2
        uz2_avei, uz2_taui, uz2_erri = O_stat_cal(uz2, tau_c)
        uz2_ave.append(uz2_avei)
        uz2_tau.append(uz2_taui)
        uz2_err.append(uz2_erri)

        # uz_abs
        uz_avei, uz_taui, uz_erri = O_stat_cal(uz2, tau_c)
        uz_ave.append(uz2_avei)
        uz_tau.append(uz2_taui)
        uz_err.append(uz2_erri)

        # lb
        lb_avei, lb_taui, lb_erri = O_stat_cal(lb, tau_c)
        lb_ave.append(lb_avei)
        lb_tau.append(lb_taui)
        lb_err.append(lb_erri)

        # Eu_bias
        '''
        Eu_avei, Eu_taui, Eu_erri = O_stat_cal(Eu, tau_c)
        Eu_ave.append(Eu_avei)
        Eu_tau.append(Eu_taui)
        Eu_err.append(Eu_erri)
        '''
    # generalize using par_nm list
    f2stail = "MC"
    for j in range(len(par)):
        if j == cpar_ind:
            f2stail += "_" + par_nm[j] + "s"
        else:
            f2stail += "_" + par_nm[j] + "%.*f" % (par_dg[j], par_dealing[j])
    f2stail += "_ana.csv"
    # only changed "lam" and "B" mode here, others waiting for further decision
    savefile = foldername + "/O_" + f2stail

    with open(savefile, "w") as f:
        f.write(mode + ",Ep_ave,Ep_tau,Ep_err")
        # p stands for per bead
        for e in range(Ne):
            f.write(",Les_ave[%d],Les_tau[%d],Les_err[%d]" % (e, e, e))

        f.write(",IdA_ave,IdA_tau,IdA_err,I2H_ave,I2H_tau,I2H_err,I2H2_ave,I2H2_tau,I2H2_err,I2H2dis_ave,I2H2dis_tau,I2H2dis_err,IK_ave,IK_tau,IK_err,p2uu_ave,p2uu_tau,p2uu_err,uuc_ave,uuc_tau,uuc_err,un2_ave,un2_tau,un2_err,uz2_ave,uz2_tau,uz2_err,uz_ave,uz_tau,uz_err,lb_ave,lb_tau,lb_err")
        f.write(",Lasym_ave,Lasym_tau,Lasym_err")
        f.write("\n")
        for i in range(len(cpar_valid)):
            f.write("%f,%f,%f,%f" % (cpar_valid[i], E_ave[i], E_tau[i], E_err[i]))
            for e in range(Ne):
                f.write(",%f,%f,%f" % (Les_ave[e][i], Les_tau[e][i], Les_err[e][i]))
            f.write(
                ",%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f"
                % (IdA_ave[i], IdA_tau[i], IdA_err[i], I2H_ave[i], I2H_tau[i], I2H_err[i], I2H2_ave[i], I2H2_tau[i], I2H2_err[i], I2H2dis_ave[i], I2H2dis_tau[i], I2H2dis_err[i], IK_ave[i], IK_tau[i], IK_err[i], p2uu_ave[i], p2uu_tau[i], p2uu_err[i], uuc_ave[i], uuc_tau[i], uuc_err[i], un2_ave[i], un2_tau[i], un2_err[i], uz2_ave[i], uz2_tau[i], uz2_err[i], uz_ave[i], uz_tau[i], uz_err[i], lb_ave[i], lb_tau[i], lb_err[i])
            )
            f.write(",%f,%f,%f" % (Lasym_ave[i], Lasym_tau[i], Lasym_err[i]))
            f.write("\n")


def O_stat_ana_multi(foldername, par, par_nm, par_dg, mode, nmulti=1, thermN=0, CnequalsKc=0, tau_c=6):
    #TODO: optimize these code for different observables using dictionary, dataframe
    cpar_valid = []
    E_ave, E_tau, E_err = [], [], []
    Ne = par[find_cpar_ind(par_nm, "Ne")]
    Les_ave, Les_tau, Les_err = [[] for i in range(Ne)], [[] for i in range(Ne)], [[] for i in range(Ne)]
    Lasym_ave, Lasym_tau, Lasym_err = [], [], []
    IdA_ave, IdA_tau, IdA_err = [], [], []
    I2H_ave, I2H_tau, I2H_err = [], [], []
    I2H2_ave, I2H2_tau, I2H2_err = [], [], []
    I2H2dis_ave, I2H2dis_tau, I2H2dis_err = [], [], []
    IK_ave, IK_tau, IK_err = [], [], []
    p2uu_ave, p2uu_tau, p2uu_err = [], [], []
    uuc_ave, uuc_tau, uuc_err = [], [], []
    # may need to add uuc2 to study the spontaneous symmetry breaking
    un2_ave, un2_tau, un2_err = [], [], []
    un2p_ave, un2p_tau, un2p_err = [], [], []
    uz2_ave, uz2_tau, uz2_err = [], [], []
    lb_ave, lb_tau, lb_err = [], [], []  # average length of bond
    Eu_ave, Eu_tau, Eu_err = [], [], []  # average length of bond
    wnu_ave, wnu_tau, wnu_err = [], [], []  # <wu^-1> = <1/wu>=<exp(Eu)>

    cpar_ind = find_cpar_ind(par_nm, mode)
    cpar = par[cpar_ind]
    Cn_ind = -1
    if CnequalsKc and mode == "Kd":
        Cn_ind = find_cpar_ind(par_nm, "Cn")

    Ncpar_valid = len(cpar)
    for i in range(len(cpar)):
        par_dealing = par[:]
        par_dealing[cpar_ind] = par[cpar_ind][i]
        if CnequalsKc and mode == "Kd":
            par_dealing[Cn_ind] = par[cpar_ind][i]  # adjust the ratio accordingly
        f2rtail = "MC"
        for j in range(len(par_dealing)):
            f2rtail += "_" + par_nm[j] + "%.*f" % (par_dg[j], par_dealing[j])
        ## add multi run data
        E_ave_multi, E_tau_multi, E_err_multi = [], [], []
        Les_ave_multi, Les_tau_multi, Les_err_multi = [[] for i in range(Ne)], [[] for i in range(Ne)], [[] for i in range(Ne)]
        Lasym_ave_multi, Lasym_tau_multi, Lasym_err_multi = [], [], []
        IdA_ave_multi, IdA_tau_multi, IdA_err_multi = [], [], []
        I2H_ave_multi, I2H_tau_multi, I2H_err_multi = [], [], []
        I2H2_ave_multi, I2H2_tau_multi, I2H2_err_multi = [], [], []
        I2H2dis_ave_multi, I2H2dis_tau_multi, I2H2dis_err_multi = [], [], []
        IK_ave_multi, IK_tau_multi, IK_err_multi = [], [], []
        p2uu_ave_multi, p2uu_tau_multi, p2uu_err_multi = [], [], []
        uuc_ave_multi, uuc_tau_multi, uuc_err_multi = [], [], []
        # may need to add uuc2 to study the spontaneous symmetry breaking
        un2_ave_multi, un2_tau_multi, un2_err_multi = [], [], []
        un2p_ave_multi, un2p_tau_multi, un2p_err_multi = [], [], []
        uz2_ave_multi, uz2_tau_multi, uz2_err_multi = [], [], []
        lb_ave_multi, lb_tau_multi, lb_err_multi = [], [], []  # average length of bond
        for k in range(nmulti):
            f2rtail_multi = f2rtail + "_id" + "%d" % k + ".csv"
            file2read = foldername + "/O_" + f2rtail_multi
            print("file2read", file2read)
            # check for file existence
            if not os.path.exists(file2read):
                print(file2read, "not exist")
                continue
            data = np.loadtxt(file2read, skiprows=14 + thermN, delimiter=",", unpack=True)
            if(k==0):
                cpar_valid.append(cpar[i])
            N = par_dealing[0]
            E = data[0] / N
            Les = data[1 : 1 + Ne]
            IdA, I2H, I2H2, I2H2dis, IK, Tp2uu, Tuuc, Bond_num, Tun2, Tuz2, Tlb = data[1 + Ne : 12 + Ne]
            Eu = data[12 + Ne]
            Lasym = eff_Le_sym_cal(Les)
            p2uu = Tp2uu / Bond_num
            uuc = Tuuc / Bond_num
            lb = Tlb / Bond_num
            N = par[find_cpar_ind(par_nm, "N")]
            un2 = Tun2 / N
            uz2 = Tuz2 / N

            E_avei, E_taui, E_erri_sq = O_stat_cal_sq(E, tau_c)
            E_ave_multi.append(E_avei)
            E_tau_multi.append(E_taui)
            E_err_multi.append(E_erri_sq)

            # Le and Ik2s
            for e in range(Ne):
                Le_avei, Le_taui, Le_erri_sq = O_stat_cal_sq(Les[e], tau_c)
                Les_ave_multi[e].append(Le_avei)
                Les_tau_multi[e].append(Le_taui)
                Les_err_multi[e].append(Le_erri_sq)

            Lasym_avei, Lasym_taui, Lasym_erri_sq = O_stat_cal_sq(Lasym, tau_c)
            Lasym_ave_multi.append(Lasym_avei)
            Lasym_tau_multi.append(Lasym_taui)
            Lasym_err_multi.append(Lasym_erri_sq)

            # IdA
            IdA_avei, IdA_taui, IdA_erri_sq = O_stat_cal_sq(IdA, tau_c)
            IdA_ave_multi.append(IdA_avei)
            IdA_tau_multi.append(IdA_taui)
            IdA_err_multi.append(IdA_erri_sq)

            # I2H
            I2H_avei, I2H_taui, I2H_erri_sq = O_stat_cal_sq(I2H, tau_c)
            I2H_ave_multi.append(I2H_avei)
            I2H_tau_multi.append(I2H_taui)
            I2H_err_multi.append(I2H_erri_sq)

            # I2H2
            I2H2_avei, I2H2_taui, I2H2_erri_sq = O_stat_cal_sq(I2H2, tau_c)
            I2H2_ave_multi.append(I2H2_avei)
            I2H2_tau_multi.append(I2H2_taui)
            I2H2_err_multi.append(I2H2_erri_sq)

            # I2H2dis
            I2H2dis_avei, I2H2dis_taui, I2H2dis_erri_sq = O_stat_cal_sq(I2H2dis, tau_c)
            I2H2dis_ave_multi.append(I2H2dis_avei)
            I2H2dis_tau_multi.append(I2H2dis_taui)
            I2H2dis_err_multi.append(I2H2dis_erri_sq)
            # IK
            IK_avei, IK_taui, IK_erri_sq = O_stat_cal_sq(IK, tau_c)
            IK_ave_multi.append(IK_avei)
            IK_tau_multi.append(IK_taui)
            IK_err_multi.append(IK_erri_sq)

            # p2uu
            p2uu_avei, p2uu_taui, p2uu_erri_sq = O_stat_cal_sq(p2uu, tau_c)
            p2uu_ave_multi.append(p2uu_avei)
            p2uu_tau_multi.append(p2uu_taui)
            p2uu_err_multi.append(p2uu_erri_sq)

            # uuc
            uuc_avei, uuc_taui, uuc_erri_sq = O_stat_cal_sq(uuc, tau_c)
            uuc_ave_multi.append(uuc_avei)
            uuc_tau_multi.append(uuc_taui)
            uuc_err_multi.append(uuc_erri_sq)

            # un2
            un2_avei, un2_taui, un2_erri_sq = O_stat_cal_sq(un2, tau_c)
            un2_ave_multi.append(un2_avei)
            un2_tau_multi.append(un2_taui)
            un2_err_multi.append(un2_erri_sq)

            # uz2
            uz2_avei, uz2_taui, uz2_erri_sq = O_stat_cal_sq(uz2, tau_c)
            uz2_ave_multi.append(uz2_avei)
            uz2_tau_multi.append(uz2_taui)
            uz2_err_multi.append(uz2_erri_sq)

            # lb
            lb_avei, lb_taui, lb_erri_sq = O_stat_cal_sq(lb, tau_c)
            lb_ave_multi.append(lb_avei)
            lb_tau_multi.append(lb_taui)
            lb_err_multi.append(lb_erri_sq)

        # calculate average based on independent run
        # \mu = \sum(x/\sigma^2)/\sum(1/\sigma^2))
        # \sigma^2 = 1/\sum(1/\sigma^2))
        def ob_mean_err_cal(x_ave,x_tau,x_err_sq):
            x_ave,x_tau,x_err_sq=np.array(x_ave),np.array(x_tau),np.array(x_err_sq)
            ob_err_sq = 1/np.sum(1/x_err_sq)
            ob_mean = np.sum(x_ave/x_err_sq)*ob_err_sq
            return (ob_mean,x_tau.mean(),np.sqrt(ob_err_sq+np.var(x_ave)/len(x_ave)))

        # E
        mn,ta,er = ob_mean_err_cal(E_ave_multi,E_tau_multi,E_err_multi)
        E_ave.append(mn)
        E_tau.append(ta)
        E_err.append(er)

        # Le and Ik2s
        for e in range(Ne):
            mn,ta,er = ob_mean_err_cal(Les_ave_multi[e],Les_tau_multi[e],Les_err_multi[e])
            Les_ave[e].append(mn)
            Les_tau[e].append(ta)
            Les_err[e].append(er)

        # Lasym
        mn,ta,er = ob_mean_err_cal(Lasym_ave_multi,Lasym_tau_multi,Lasym_err_multi)
        Lasym_ave.append(mn)
        Lasym_tau.append(ta)
        Lasym_err.append(er)

        # IdA
        mn,ta,er = ob_mean_err_cal(IdA_ave_multi,IdA_tau_multi,IdA_err_multi)
        IdA_ave.append(mn)
        IdA_tau.append(ta)
        IdA_err.append(er)

        # I2H
        mn,ta,er = ob_mean_err_cal(I2H_ave_multi,I2H_tau_multi,I2H_err_multi)
        I2H_ave.append(mn)
        I2H_tau.append(ta)
        I2H_err.append(er)

        # I2H2
        mn,ta,er = ob_mean_err_cal(I2H2_ave_multi,I2H2_tau_multi,I2H2_err_multi)
        I2H2_ave.append(mn)
        I2H2_tau.append(ta)
        I2H2_err.append(er)

        # I2H2dis
        mn,ta,er = ob_mean_err_cal(I2H2dis_ave_multi,I2H2dis_tau_multi,I2H2dis_err_multi)
        I2H2dis_ave.append(mn)
        I2H2dis_tau.append(ta)
        I2H2dis_err.append(er)
        # IK
        mn,ta,er = ob_mean_err_cal(IK_ave_multi,IK_tau_multi,IK_err_multi)
        IK_ave.append(mn)
        IK_tau.append(ta)
        IK_err.append(er)

        # p2uu
        mn,ta,er = ob_mean_err_cal(p2uu_ave_multi,p2uu_tau_multi,p2uu_err_multi)
        p2uu_ave.append(mn)
        p2uu_tau.append(ta)
        p2uu_err.append(er)

        # uuc
        mn,ta,er = ob_mean_err_cal(uuc_ave_multi,uuc_tau_multi,uuc_err_multi)
        uuc_ave.append(mn)
        uuc_tau.append(ta)
        uuc_err.append(er)

        # un2
        mn,ta,er = ob_mean_err_cal(un2_ave_multi,un2_tau_multi,un2_err_multi)
        un2_ave.append(mn)
        un2_tau.append(ta)
        un2_err.append(er)

        # uz2
        mn,ta,er = ob_mean_err_cal(uz2_ave_multi,uz2_tau_multi,uz2_err_multi)
        uz2_ave.append(mn)
        uz2_tau.append(ta)
        uz2_err.append(er)

        # lb
        mn,ta,er = ob_mean_err_cal(lb_ave_multi,lb_tau_multi,lb_err_multi)
        lb_ave.append(mn)
        lb_tau.append(ta)
        lb_err.append(er)

    # generalize using par_nm list
    f2stail = "MC"
    for j in range(len(par)):
        if j == cpar_ind:
            f2stail += "_" + par_nm[j] + "s"
        else:
            f2stail += "_" + par_nm[j] + "%.*f" % (par_dg[j], par_dealing[j])
    f2stail += "_ana.csv"
    # only changed "lam" and "B" mode here, others waiting for further decision
    savefile = foldername + "/O_" + f2stail
    print("len(E_ave)",len(E_ave))
    print("len(cpar)",len(cpar))
    print("len(cpar_valid)",len(cpar_valid))

    with open(savefile, "w") as f:
        f.write(mode + ",Ep_ave,Ep_tau,Ep_err")
        # p stands for per bead
        for e in range(Ne):
            f.write(",Les_ave[%d],Les_tau[%d],Les_err[%d]" % (e, e, e))

        f.write(",IdA_ave,IdA_tau,IdA_err,I2H_ave,I2H_tau,I2H_err,I2H2_ave,I2H2_tau,I2H2_err,I2H2dis_ave,I2H2dis_tau,I2H2dis_err,IK_ave,IK_tau,IK_err,p2uu_ave,p2uu_tau,p2uu_err,uuc_ave,uuc_tau,uuc_err,un2_ave,un2_tau,un2_err,uz2_ave,uz2_tau,uz2_err,lb_ave,lb_tau,lb_err,Eubias_ave,Eubias_tau,Eubias_err")
        f.write(",Lasym_ave,Lasym_tau,Lasym_err")
        f.write("\n")
        for i in range(len(cpar_valid)):
            f.write("%f,%f,%f,%f" % (cpar_valid[i], E_ave[i], E_tau[i], E_err[i]))
            for e in range(Ne):
                f.write(",%f,%f,%f" % (Les_ave[e][i], Les_tau[e][i], Les_err[e][i]))
            f.write(
                ",%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f"
                % (IdA_ave[i], IdA_tau[i], IdA_err[i], I2H_ave[i], I2H_tau[i], I2H_err[i], I2H2_ave[i], I2H2_tau[i], I2H2_err[i], I2H2dis_ave[i], I2H2dis_tau[i], I2H2dis_err[i], IK_ave[i], IK_tau[i], IK_err[i], p2uu_ave[i], p2uu_tau[i], p2uu_err[i], uuc_ave[i], uuc_tau[i], uuc_err[i], un2_ave[i], un2_tau[i], un2_err[i], uz2_ave[i], uz2_tau[i], uz2_err[i], lb_ave[i], lb_tau[i], lb_err[i],0, 0,0)
            )
            f.write(",%f,%f,%f" % (Lasym_ave[i], Lasym_tau[i], Lasym_err[i]))
            f.write("\n")


def Gij_stat_ana(foldername, par, par_nm, par_dg, mode, tau_c=6):
    # TODO: add shape descriptors ref: https://en.wikipedia.org/wiki/Gyration_tensor
    # GRg, Gb, Gc, Gkap

    Gxx_ave, Gxx_tau, Gxx_err = [], [], []  # just simplying analyze the xx component
    Gyy_ave, Gyy_tau, Gyy_err = [], [], []
    Gzz_ave, Gzz_tau, Gzz_err = [], [], []
    Geig0_ave, Geig0_tau, Geig0_err = [], [], []
    Geig1_ave, Geig1_tau, Geig1_err = [], [], []
    Geig2_ave, Geig2_tau, Geig2_err = [], [], []
    # Geig02_ave, Geig02_tau, Geig02_err = [], [], []  # eig0/eig2 ratio
    GRg2_ave, GRg2_tau, GRg2_err = [], [], []  # Rg^2 = sum eig_i
    Gb_ave, Gb_tau, Gb_err = [], [], []  # b =  eig2 - (eig0+eig1)/2
    Gkap2_ave, Gkap2_tau, Gkap2_err = [], [], []  # kap^2 = 3/2(sum eigi^2/ (sum eig_i)^2) -1/2
    Dedge_ave, Dedge_tau, Dedge_err = [], [], []
    cpar_ind = find_cpar_ind(par_nm, mode)
    cpar = par[cpar_ind]
    for i in range(len(cpar)):
        par_dealing = par[:]
        par_dealing[cpar_ind] = par[cpar_ind][i]
        f2rtail = "/Gij"
        for j in range(len(par_dealing)):
            # print("par_dealing[j]",j,par_dealing[j])
            f2rtail += "_" + par_nm[j] + "%.*f" % (par_dg[j], par_dealing[j])
        f2rtail += ".csv"
        file2read = foldername + f2rtail
        Gdata = np.loadtxt(file2read, skiprows=1, usecols=range(10), delimiter=",", unpack=True)
        # get edge-edge distance Dedge
        Dedge = Gdata[0]
        Dedge_ave.append(np.average(Dedge))
        rho, cov0 = autocorrelation_function_fft(Dedge)
        tau, tau_err = tau_int_cal_rho(rho, tau_c)
        Dedge_tau.append(tau)
        Dedge_err.append(np.sqrt(2 * tau / len(Dedge) * cov0))

        # get eigenvalues of Gij
        Gijs = np.transpose(Gdata[1:])
        Gxxs, Gyys, Gzzs = Gdata[1], Gdata[5], Gdata[9]
        Geigs = []
        for Gij in Gijs:
            w, v = np.linalg.eig(np.reshape(Gij, (3, 3)))
            Geigs.append(np.sort(w))
        Geigs = np.transpose(Geigs)  # Geigs become np.array() after transpose using numpy
        GRg2s = Geigs[0] + Geigs[1] + Geigs[2]
        Gbs = Geigs[2] - (Geigs[0] + Geigs[1]) / 2
        Gkap2s = 3 / 2 * (Geigs[0] * Geigs[0] + Geigs[1] * Geigs[1] + Geigs[2] * Geigs[2]) / (GRg2s * GRg2s) - 1 / 2

        # put into Observables contaners
        Gxx_ave.append(np.average(Gxxs))
        rho, cov0 = autocorrelation_function_fft(Gxxs)
        tau, tau_err = tau_int_cal_rho(rho, tau_c)
        Gxx_tau.append(tau)
        Gxx_err.append(np.sqrt(2 * tau / len(Gxxs) * cov0))

        Gyy_ave.append(np.average(Gyys))
        rho, cov0 = autocorrelation_function_fft(Gyys)
        tau, tau_err = tau_int_cal_rho(rho, tau_c)
        Gyy_tau.append(tau)
        Gyy_err.append(np.sqrt(2 * tau / len(Gyys) * cov0))

        Gzz_ave.append(np.average(Gzzs))
        rho, cov0 = autocorrelation_function_fft(Gzzs)
        tau, tau_err = tau_int_cal_rho(rho, tau_c)
        Gzz_tau.append(tau)
        Gzz_err.append(np.sqrt(2 * tau / len(Gzzs) * cov0))

        Geig0_ave.append(np.average(Geigs[0]))
        rho, cov0 = autocorrelation_function_fft(Geigs[0])
        tau, tau_err = tau_int_cal_rho(rho, tau_c)
        Geig0_tau.append(tau)
        Geig0_err.append(np.sqrt(2 * tau / len(Geigs[0]) * cov0))

        Geig1_ave.append(np.average(Geigs[1]))
        rho, cov0 = autocorrelation_function_fft(Geigs[1])
        tau, tau_err = tau_int_cal_rho(rho, tau_c)
        Geig1_tau.append(tau)
        Geig1_err.append(np.sqrt(2 * tau / len(Geigs[1]) * cov0))

        Geig2_ave.append(np.average(Geigs[2]))
        rho, cov0 = autocorrelation_function_fft(Geigs[2])
        tau, tau_err = tau_int_cal_rho(rho, tau_c)
        Geig2_tau.append(tau)
        Geig2_err.append(np.sqrt(2 * tau / len(Geigs[2]) * cov0))

        GRg2_ave.append(np.average(GRg2s))
        rho, cov0 = autocorrelation_function_fft(GRg2s)
        tau, tau_err = tau_int_cal_rho(rho, tau_c)
        GRg2_tau.append(tau)
        GRg2_err.append(np.sqrt(2 * tau / len(GRg2s) * cov0))

        Gb_ave.append(np.average(Gbs))
        rho, cov0 = autocorrelation_function_fft(Gbs)
        tau, tau_err = tau_int_cal_rho(rho, tau_c)
        Gb_tau.append(tau)
        Gb_err.append(np.sqrt(2 * tau / len(Gbs) * cov0))

        Gkap2_ave.append(np.average(Gkap2s))
        rho, cov0 = autocorrelation_function_fft(Gkap2s)
        tau, tau_err = tau_int_cal_rho(rho, tau_c)
        Gkap2_tau.append(tau)
        Gkap2_err.append(np.sqrt(2 * tau / len(Gkap2s) * cov0))

    # save result to file
    f2stail = "/Gij"
    for j in range(len(par)):
        if j == cpar_ind:
            f2stail += "_" + par_nm[j] + "s"
        else:
            f2stail += "_" + par_nm[j] + "%.*f" % (par_dg[j], par_dealing[j])
    f2stail += "_ana.csv"
    # only changed "lam" and "B" mode here, others waiting for further decision
    savefile = foldername + f2stail

    # TODO: add additional observables to file
    with open(savefile, "w") as f:
        f.write(mode + ",Dedge_ave,Dedge_tau,Dedge_err,Gxx_ave,Gxx_tau,Gxx_err,Gyy_ave,Gyy_tau,Gyy_err,Gzz_ave,Gzz_tau,Gzz_err,Geig0_ave,Geig0_tau,Geig0_err,Geig1_ave,Geig1_tau,Geig1_err,Geig2_ave,Geig2_tau,Geig2_err,GRg2_ave,GRg2_tau,GRg2_err,Gb_ave,Gb_tau,Gb_err,Gkap2_ave,Gkap2_tau,Gkap2_err\n")
        for i in range(len(cpar)):
            f.write(
                "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n"
                % (cpar[i], Dedge_ave[i], Dedge_tau[i], Dedge_err[i], Gxx_ave[i], Gxx_tau[i], Gxx_err[i], Gyy_ave[i], Gyy_tau[i], Gyy_err[i], Gzz_ave[i], Gzz_tau[i], Gzz_err[i], Geig0_ave[i], Geig0_tau[i], Geig0_err[i], Geig1_ave[i], Geig1_tau[i], Geig1_err[i], Geig2_ave[i], Geig2_tau[i], Geig2_err[i], GRg2_ave[i], GRg2_tau[i], GRg2_err[i], Gb_ave[i], Gb_tau[i], Gb_err[i], Gkap2_ave[i], Gkap2_tau[i], Gkap2_err[i])
            )
