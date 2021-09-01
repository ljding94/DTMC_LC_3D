import numpy as np
import matplotlib.pyplot as plt
from plot import *
from autocorrelation import *
from scipy.optimize import curve_fit
from scipy import odr
import scipy.fft

def find_cpar_ind(par_nm,mode):
    cpar_ind = -1
    for i in range(len(par_nm)):
        if par_nm[i]==mode:
            cpar_ind=i
            break
    return cpar_ind

def O_stat_ana(foldername,par,par_nm,par_dg, mode, CnequalsKc=0, tau_c=6):
    E_ave, E_tau, E_err = [], [], []
    Ne = par[find_cpar_ind(par_nm,"Ne")]

    Les_ave, Les_tau, Les_err = [[] for i in range(Ne)], [[] for i in range(Ne)], [[] for i in range(Ne)]
    #Leuns_ave, Leuns_tau, Leuns_err = [[] for i in range(Ne)], [[] for i in range(Ne)], [[] for i in range(Ne)]
    Ik2s_ave, Ik2s_tau, Ik2s_err = [[] for i in range(Ne)], [[] for i in range(Ne)], [[] for i in range(Ne)]
    IdA_ave, IdA_tau, IdA_err = [], [], []
    I2H_ave, I2H_tau, I2H_err = [], [], []
    I2H2_ave, I2H2_tau, I2H2_err = [], [], []
    phip_ave,phip_err,phip_tau = [], [], [] # phi per bead
    phi2p_ave,phi2p_err, phi2p_tau = [],[],[] # phi interaction per bond
    I2H2dis_ave, I2H2dis_tau, I2H2dis_err = [], [], []
    IK_ave, IK_tau, IK_err = [], [], []
    p2uu_ave, p2uu_tau, p2uu_err = [], [], []
    uuc_ave, uuc_tau, uuc_err = [], [], []
    # may need to add uuc2 to study the spontaneous symmetry breaking
    un2_ave,un2_tau,un2_err = [],[],[]
    un2p_ave,un2p_tau,un2p_err = [],[],[]
    IKun2_ave, IKun2_tau, IKun2_err = [], [], []
    if(Ne==2):
        Ledif_ave,Ledif_tau,Ledif_err=[],[],[]
    cpar_ind = find_cpar_ind(par_nm,mode)
    cpar = par[cpar_ind]
    Cn_ind=-1
    if(CnequalsKc and mode=="Kd"):
        Cn_ind=find_cpar_ind(par_nm,"Cn")

    for i in range(len(cpar)):
        par_dealing = par[:]
        par_dealing[cpar_ind] = par[cpar_ind][i]
        if(CnequalsKc and mode=="Kd"):
            par_dealing[Cn_ind]=par[cpar_ind][i]
        f2rtail = "MC"
        for j in range(len(par_dealing)):
            f2rtail+="_"+par_nm[j]+"%.*f"%(par_dg[j],par_dealing[j])
        f2rtail+=".csv"
        #print("f2rtail",f2rtail)
        file2read = foldername + "/O_"+f2rtail
        data = np.loadtxt(file2read, skiprows=14, delimiter=",", unpack=True)
        N = par_dealing[0]
        E = data[0]/N
        Les = data[1:1+Ne]
        #Leuns = data[1+Ne:1+2*Ne]
        Ik2s = data[1+Ne:1+2*Ne]
        #IdA,I2H,I2H2,phi,Tphi2,I2H2dis,IK,IKphi2,Tp2uu,Tuuc,Bond_num,Tun2 = data[1+Ne:]
        IdA,I2H,I2H2,phi,Tphi2,I2H2dis,IK,IKphi2,Tp2uu,Tuuc,Bond_num,Tun2 = data[1+2*Ne:]
        phip = phi/N
        phi2p = Tphi2/Bond_num
        p2uu = Tp2uu/Bond_num
        uuc = Tuuc/Bond_num
        N =par[find_cpar_ind(par_nm,"N")]
        #rCnp = par[find_cpar_ind(par_nm,"rCnp")]
        #Np = int(N*rCnp)
        #un2=Tun2/(N-Np)
        un2=Tun2/N
        #un2p=Tun2p/Np if Np>0 else Tun2p
        # Ne2 case, need Ledif for additional info
        if(Ne==2):
            Ledif = np.abs(Les[0]-Les[1])
            Ledif_ave.append(np.average(Ledif))
            rho, cov0 = autocorrelation_function_fft(Ledif)
            tau, tau_err = tau_int_cal_rho(rho,tau_c)
            Ledif_tau.append(tau)
            Ledif_err.append(np.sqrt(2 * tau / len(Ledif) * cov0))
        #was used for checking energy
        # again, it's correct
        '''
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
        '''

        # E
        E_ave.append(np.average(E))
        rho, cov0 = autocorrelation_function_fft(E)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        E_tau.append(tau)
        E_err.append(np.sqrt(2 * tau / len(E) * cov0))

        # Le and Ik2s
        for e in range(Ne):
            Les_ave[e].append(np.average(Les[e]))
            rho, cov0 = autocorrelation_function_fft(Les[e])
            tau, tau_err = tau_int_cal_rho(rho,tau_c)
            Les_tau[e].append(tau)
            Les_err[e].append(np.sqrt(2 * tau / len(Les[e]) * cov0))

            Ik2s_ave[e].append(np.average(Ik2s[e]))
            rho, cov0 = autocorrelation_function_fft(Ik2s[e])
            tau, tau_err = tau_int_cal_rho(rho,tau_c)
            Ik2s_tau[e].append(tau)
            Ik2s_err[e].append(np.sqrt(2 * tau / len(Ik2s[e]) * cov0))



        # IdA
        IdA_ave.append(np.average(IdA))
        rho, cov0 = autocorrelation_function_fft(IdA)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        IdA_tau.append(tau)
        IdA_err.append(np.sqrt(2 * tau / len(IdA) * cov0))

        # I2H
        I2H_ave.append(np.average(I2H))
        rho, cov0 = autocorrelation_function_fft(I2H)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        I2H_tau.append(tau)
        I2H_err.append(np.sqrt(2 * tau / len(I2H) * cov0))

        # I2H2
        I2H2_ave.append(np.average(I2H2))
        rho, cov0 = autocorrelation_function_fft(I2H2)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        # autocorrelation_plot(rho, tau, file2read[:-4] + "_autoI2H2.pdf")
        I2H2_tau.append(tau)
        I2H2_err.append(np.sqrt(2 * tau / len(I2H2) * cov0))

        # phi
        phip_ave.append(np.average(phip))
        rho, cov0 = autocorrelation_function_fft(phip)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        phip_tau.append(tau)
        phip_err.append(np.sqrt(2 * tau / len(phip) * cov0))

        # Tphi2
        phi2p_ave.append(np.average(phi2p))
        rho, cov0 = autocorrelation_function_fft(phi2p)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        phi2p_tau.append(tau)
        phi2p_err.append(np.sqrt(2 * tau / len(phi2p) * cov0))

        # I2H2dis
        I2H2dis_ave.append(np.average(I2H2dis))
        rho, cov0 = autocorrelation_function_fft(I2H2dis)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        I2H2dis_tau.append(tau)
        I2H2dis_err.append(np.sqrt(2 * tau / len(I2H2dis) * cov0))

        # IK
        IK_ave.append(np.average(IK))
        rho, cov0 = autocorrelation_function_fft(IK)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        IK_tau.append(tau)
        IK_err.append(np.sqrt(2 * tau / len(IK) * cov0))

        # p2uu
        p2uu_ave.append(np.average(p2uu))
        rho, cov0 = autocorrelation_function_fft(p2uu)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        p2uu_tau.append(tau)
        p2uu_err.append(np.sqrt(2 * tau / len(p2uu) * cov0))

        # uuc
        uuc_ave.append(np.average(uuc))
        rho, cov0 = autocorrelation_function_fft(uuc)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        uuc_tau.append(tau)
        uuc_err.append(np.sqrt(2 * tau / len(uuc) * cov0))

        # un2
        un2_ave.append(np.average(un2))
        rho, cov0 = autocorrelation_function_fft(un2)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        un2_tau.append(tau)
        un2_err.append(np.sqrt(2 * tau / len(un2) * cov0))


    # generalize using par_nm list
    f2stail = "MC"
    for j in range(len(par)):
        if(j==cpar_ind):
            f2stail+="_"+par_nm[j]+"s"
        else:
            f2stail+="_"+par_nm[j]+"%.*f"%(par_dg[j],par_dealing[j])
    f2stail+="_ana.csv"
    # only changed "lam" and "B" mode here, others waiting for further decision
    savefile = foldername + "/O_" + f2stail

    with open(savefile, "w") as f:
        f.write(mode+",Ep_ave,Ep_tau,Ep_err")
        # p stands for per bead
        for e in range(Ne):
            f.write(",Les_ave[%d],Les_tau[%d],Les_err[%d]"%(e,e,e))
        for e in range(Ne):
            f.write(",Ik2s_ave[%d],Ik2s_tau[%d],Ik2s_err[%d]"%(e,e,e))
        #for e in range(Ne):
        #    f.write(",Leuns_ave[%d],Leuns_tau[%d],Leuns_err[%d]"%(e,e,e))

        f.write(",IdA_ave,IdA_tau,IdA_err,I2H_ave,I2H_tau,I2H_err,I2H2_ave,I2H2_tau,I2H2_err,phip_ave,phip_tau,phip_err,phi2p_ave,phi2p_tau,phi2p_err,I2H2dis_ave,I2H2dis_tau,I2H2dis_err,IK_ave,IK_tau,IK_err,p2uu_ave,p2uu_tau,p2uu_err,uuc_ave,uuc_tau,uuc_err,un2_ave,un2_tau,un2_err")
        if(Ne==2):
            f.write(",Ledif_ave,Ledif_tau,Ledif_err")
        f.write("\n")
        for i in range(len(cpar)):
            f.write("%f,%f,%f,%f" % (cpar[i], E_ave[i], E_tau[i], E_err[i]))
            for e in range(Ne):
                f.write(",%f,%f,%f"%(Les_ave[e][i],Les_tau[e][i], Les_err[e][i]))
            for e in range(Ne):
                f.write(",%f,%f,%f"%(Ik2s_ave[e][i],Ik2s_tau[e][i], Ik2s_err[e][i]))
            #for e in range(Ne):
            #    f.write(",%f,%f,%f"%(Leuns_ave[e][i],Leuns_tau[e][i], Leuns_err[e][i]))
            f.write(",%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f"%(IdA_ave[i], IdA_tau[i], IdA_err[i],I2H_ave[i], I2H_tau[i], I2H_err[i],I2H2_ave[i], I2H2_tau[i], I2H2_err[i], phip_ave[i],phip_tau[i],phip_err[i], phi2p_ave[i],phi2p_tau[i],phi2p_err[i],I2H2dis_ave[i],I2H2dis_tau[i],I2H2dis_err[i], IK_ave[i], IK_tau[i], IK_err[i], p2uu_ave[i], p2uu_tau[i], p2uu_err[i], uuc_ave[i], uuc_tau[i], uuc_err[i], un2_ave[i], un2_tau[i], un2_err[i]))
            if(Ne==2):
                f.write(",%f,%f,%f"%(Ledif_ave[i], Ledif_tau[i], Ledif_err[i]))
            f.write("\n")

def Gij_stat_ana(foldername,par,par_nm,par_dg,mode,tau_c=6):
    Geig0_ave, Geig0_tau, Geig0_err = [], [], []
    Geig2_ave, Geig1_tau, Geig1_err = [], [], []
    Geig1_ave, Geig2_tau, Geig2_err = [], [], []
    Dedge_ave, Dedge_tau, Dedge_err= [], [], []
    cpar_ind = find_cpar_ind(par_nm,mode)
    cpar = par[cpar_ind]
    for i in range(len(cpar)):
        par_dealing = par[:]
        par_dealing[cpar_ind] = par[cpar_ind][i]
        f2rtail = "/Gij"
        for j in range(len(par_dealing)):
            #print("par_dealing[j]",j,par_dealing[j])
            f2rtail+="_"+par_nm[j]+"%.*f"%(par_dg[j],par_dealing[j])
        f2rtail+=".csv"
        file2read = foldername +f2rtail
        Gdata = np.loadtxt(file2read, skiprows=1,usecols=range(10), delimiter=",", unpack=True)
        # get edge-edge distance Dedge
        Dedge= Gdata[0]
        Dedge_ave.append(np.average(Dedge))
        rho, cov0 = autocorrelation_function_fft(Dedge)
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        Dedge_tau.append(tau)
        Dedge_err.append(np.sqrt(2 * tau / len(Dedge) * cov0))

        # get eigenvalues of Gij
        Gijs = np.transpose(Gdata[1:])
        Geigs=[]
        for Gij in Gijs:
            w,v = np.linalg.eig(np.reshape(Gij,(3,3)))
            Geigs.append(np.sort(w))
        Geigs = np.transpose(Geigs)
        # put into
        Geig0_ave.append(np.average(Geigs[0]))
        rho, cov0 = autocorrelation_function_fft(Geigs[0])
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        Geig0_tau.append(tau)
        Geig0_err.append(np.sqrt(2 * tau / len(Geigs[0]) * cov0))

        Geig1_ave.append(np.average(Geigs[1]))
        rho, cov0 = autocorrelation_function_fft(Geigs[1])
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        Geig1_tau.append(tau)
        Geig1_err.append(np.sqrt(2 * tau / len(Geigs[1]) * cov0))

        Geig2_ave.append(np.average(Geigs[2]))
        rho, cov0 = autocorrelation_function_fft(Geigs[2])
        tau, tau_err = tau_int_cal_rho(rho,tau_c)
        Geig2_tau.append(tau)
        Geig2_err.append(np.sqrt(2 * tau / len(Geigs[2]) * cov0))

    # save result to file
    f2stail = "/Gij"
    for j in range(len(par)):
        if(j==cpar_ind):
            f2stail+="_"+par_nm[j]+"s"
        else:
            f2stail+="_"+par_nm[j]+"%.*f"%(par_dg[j],par_dealing[j])
    f2stail+="_ana.csv"
    # only changed "lam" and "B" mode here, others waiting for further decision
    savefile = foldername + f2stail

    with open(savefile, "w") as f:
        f.write(mode+",Dedge_ave,Dedge_tau,Dedge_err,Geig0_ave,Geig0_tau,Geig0_err,Geig1_ave,Geig1_tau,Geig1_err,Geig2_ave,Geig2_tau,Geig2_err\n")
        for i in range(len(cpar)):
            f.write("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n" % (cpar[i], Dedge_ave[i],Dedge_tau[i],Dedge_err[i],Geig0_ave[i], Geig0_tau[i], Geig0_err[i], Geig1_ave[i], Geig1_tau[i], Geig1_err[i], Geig2_ave[i], Geig2_tau[i], Geig2_err[i]))


