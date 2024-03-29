#include "dtmc_lc_3d.h"
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

void dtmc_lc::State_write(std::string filename)
{
    std::ofstream f(filename);
    if (f.is_open())
    {
        // find maxim nei.size()
        int max_nei_size = 0;
        for (int i = 0; i < mesh.size(); i++)
        {
            if (mesh[i].nei.size() > max_nei_size)
            {
                max_nei_size = mesh[i].nei.size();
            }
        }
        f << "beta=" << beta << "\n";
        f << "N=" << N << "\n";
        f << "imod=" << imod << "\n";
        f << "l0=" << l0 << "\n";
        f << "max_nei_size=" << max_nei_size << "\n";
        //f << "x,y,z,ux,uy,uz,nx,ny,nz,phi,dA,2H,ds,dAK,un2,edge_num,edge_neibs,neibs";
        f << "x,y,z,ux,uy,uz,nx,ny,nz,dA,2H,ds,dAK,un2,edge_num,edge_neibs,neibs";
        //f << "x,y,z,ux,uy,uz,phi,dA,2H,ds,dAK,un2,edge_num,edge_neibs,neibs";
        for (int i = 0; i < mesh.size(); i++)
        {
            f << "\n"
              << mesh[i].R[0] << "," << mesh[i].R[1] << "," << mesh[i].R[2];
            f << "," << mesh[i].u[0] << "," << mesh[i].u[1] << ","
              << mesh[i].u[2];
            f << "," << mesh[i].n[0] << "," << mesh[i].n[1] << "," << mesh[i].n[2];
            //f << "," << mesh[i].phi;
            f << "," << mesh[i].dAn2H[0] << "," << mesh[i].dAn2H[1] << ","
              << mesh[i].ds << "," << mesh[i].dAK << "," << mesh[i].un2;
            f << "," << mesh[i].edge_num;
            for (int j = 0; j < mesh[i].edge_nei.size(); j++)
            {
                f << "," << mesh[i].edge_nei[j];
            }
            for (int j = 0; j < 2 - mesh[i].edge_nei.size(); j++)
            {
                f << ",-1";
            }
            for (int j = 0; j < mesh[i].nei.size(); j++)
            {
                f << "," << mesh[i].nei[j];
            }
            for (int j = 0; j < max_nei_size - mesh[i].nei.size(); j++)
            {
                f << ",-1";
            }
        }
    }
    f.close();
} // end State_write

void dtmc_lc::State_load(std::string state_file)
{
    // FIXME: this function is indeed broken now
    // need to be re con struct based on new state write code above
    std::ifstream f(state_file);
    std::string buff;
    int max_nei_size;
    while (f.good())
    {
        // skip 8 lines
        for (int i = 0; i < 13; i++)
        {
            std::getline(f, buff);
            // std::cout << buff << "\n";
        }
        // load geometric obsevables
        /*
        std::getline(f, buff, '=');
        std::getline(f, buff);
        Les[0] = std::stof(buff);
        // kinda don't wanna use it any more
        std::getline(f, buff, '=');
        std::getline(f, buff);
        I2H2 = std::stof(buff);
        std::getline(f, buff, '=');
        std::getline(f, buff);
        IK = std::stof(buff);
        std::getline(f, buff, '=');
        std::getline(f, buff);
        */
        max_nei_size = std::stof(buff);
        // skipp another line
        std::getline(f, buff);
        // std::cout << buff << "\n";

        int mcount = 0;

        for (mcount = 0; mcount < mesh.size(); mcount++)
        {
            // load observables
            std::getline(f, buff, ',');
            mesh[mcount].ds = std::stof(buff);
            std::getline(f, buff, ',');
            mesh[mcount].dAn2H[0] = std::stof(buff);
            std::getline(f, buff, ',');
            mesh[mcount].dAn2H[1] = std::stof(buff);
            std::getline(f, buff, ',');
            mesh[mcount].dAK = std::stof(buff);
            // load pos
            for (int k = 0; k < 3; k++)
            {
                std::getline(f, buff, ',');
                mesh[mcount].R[k] = std::stof(buff);
            }
            // load edge_nei
            mesh[mcount].edge_nei.clear();
            for (int k = 0; k < 2; k++)
            {
                std::getline(f, buff, ',');
                if (std::stof(buff) != -1)
                {
                    mesh[mcount].edge_nei.push_back(std::stof(buff));
                }
            }
            // load neibs
            // should read to the end of line, but how?
            mesh[mcount].nei.clear();
            for (int k = 0; k < max_nei_size - 1; k++)
            {
                std::getline(f, buff, ',');
                if (std::stof(buff) != -1)
                {
                    mesh[mcount].nei.push_back(std::stof(buff));
                }
            }
            std::getline(f, buff);
            if (std::stof(buff) != -1)
            {
                mesh[mcount].nei.push_back(std::stof(buff));
            }
        }
    }
    //E = 0;
}

void dtmc_lc::State_write_seq(std::string filename, int MC_sweeps,
                              int step_p_sweep, double delta_s,
                              double delta_theta)
{
    std::string savename;
    for (int sweep_n = 0; sweep_n < MC_sweeps; sweep_n++)
    {
        Thermal(1, step_p_sweep, 1, delta_s, delta_theta);
        savename = filename;
        savename.insert(savename.size() - 4, "_" + std::to_string(sweep_n));
        State_write(savename);
    }
}

void dtmc_lc::Thermal(int MC_sweeps, int step_p_sweep, int beta_steps,
                      double delta_s, double delta_theta)
{
    beta = 0.0;
    for (int n_beta = 0; n_beta < beta_steps; n_beta++)
    {
        beta += 1.0 / beta_steps;
        for (int sweep_n = 0; sweep_n < MC_sweeps / beta_steps; sweep_n++)
        {
            std::cout << sweep_n + n_beta * MC_sweeps / beta_steps << "/" << MC_sweeps << "\n";
            //std::cout << "Eu/sqrt{N}=" << Ob_sys.Eu / std::sqrt(N) << "\n";
            //std::cout << " Les[1]/Les[0]" << Ob_sys.Les[1] / Ob_sys.Les[0] << "\n";
            for (int i = 0; i < step_p_sweep; i++)
            {
                bead_metropolis(delta_s);
                spin_metropolis(delta_theta);
                bond_metropolis();
                bond_metropolis();
                //hop_metropolis();
                if (i % int(std::sqrt(N)) == 0)
                {
                    edge_metropolis();
                    //lifted_edge_metropolis(); // it should works faster
                }
            }
            // std::cout << "thermo, beta=" << beta << "," << sweep_n << "/"<<
            // MC_sweeps << "\n";
        }
    }
}
void dtmc_lc::Thermal_kar1lam1(int MC_sweeps, int step_p_sweep, double kar1, double lam1, double delta_s, double delta_theta)
{
    double kar0, lam0, Kd0;
    kar0 = Epar.kar;
    lam0 = Epar.lam;
    Kd0 = Epar.Kd;
    Epar.kar = kar1;
    Epar.lam = lam1;
    Epar.Kd = 0;
    Ob_sys.E += Ob_sys.I2H2dis * (kar1 - kar0);
    for (int e = 0; e < Ne; e++)
    {
        Ob_sys.E += Ob_sys.Les[e] * (lam1 - lam0);
    }
    Ob_sys.E += (0 - Kd0) * (Ob_sys.Tp2uu + Epar.q * Ob_sys.Tuuc);

    for (int sweep_n = 0; sweep_n < MC_sweeps; sweep_n++)
    {
        std::cout << "thermal_kar1lam1:" << sweep_n << "/" << MC_sweeps << "\n";
        for (int i = 0; i < step_p_sweep; i++)
        {
            bead_metropolis(delta_s);
            spin_metropolis(delta_theta);
            bond_metropolis();
            bond_metropolis();
            //hop_metropolis();
            if (i % int(std::sqrt(N)) == 0)
            {
                edge_metropolis();
                //lifted_edge_metropolis(); // it should works faster
            }
        }
        // std::cout << "thermo, beta=" << beta << "," << sweep_n << "/"<<
        // MC_sweeps << "\n";
    }

    Epar.kar = kar0;
    Epar.lam = lam0;
    Epar.Kd = Kd0;
    Ob_sys.E += Ob_sys.I2H2dis * (kar0 - kar1);
    for (int e = 0; e < Ne; e++)
    {
        Ob_sys.E += Ob_sys.Les[e] * (lam0 - lam1);
    }
    Ob_sys.E += (Kd0 - 0) * (Ob_sys.Tp2uu + Epar.q * Ob_sys.Tuuc);
    // correct energy for new kar0
}

void dtmc_lc::Thermal_pinch(int MC_sweeps, int step_p_sweep, double k_pinch_, double delta_s, double delta_theta)
{
    double lf0;
    lf0 = lf;
    if (lf == 0)
    {
        lf = 5.0;
        edge_zlim.resize(2);
        edge_zlim = {-0.5 * lf, 0.5 * lf};
    }
    k_pinch = k_pinch_;
    for (int sweep_n = 0; sweep_n < MC_sweeps; sweep_n++)
    {
        std::cout << "thermal_pinch:" << sweep_n << "/" << MC_sweeps << "\n";
        for (int i = 0; i < step_p_sweep; i++)
        {
            bead_metropolis(delta_s);
            spin_metropolis(delta_theta);
            bond_metropolis();
            bond_metropolis();
            //hop_metropolis();
            if (i % int(std::sqrt(N)) == 0)
            {
                edge_metropolis();
                //lifted_edge_metropolis(); // it should works faster
            }
        }
        // std::cout << "thermo, beta=" << beta << "," << sweep_n << "/"<<
        // MC_sweeps << "\n";
    }

    k_pinch = 0;
    lf = lf0;
    edge_zlim.clear();
    // correct energy for new kar0
}

void dtmc_lc::O_MC_measure(int MC_sweeps, int sweep_p_G, int step_p_sweep,
                           double delta_s, double delta_theta, double delta_r, double bin_num, std::string folder,
                           std::string finfo, int seq)
{
    std::vector<double> E_all;
    std::vector<double> I2H2_all;
    std::vector<double> I2H2dis_all;
    std::vector<double> IK_all;
    std::vector<std::vector<double>> Les_all;
    Les_all.resize(Ne);
    std::vector<double> Tp2uu_all;
    std::vector<double> Tuuc_all;
    std::vector<double> Tun2_all;
    std::vector<double> IdA_all;
    std::vector<double> I2H_all;
    std::vector<int> Bond_num_all;
    std::vector<double> Tuz2_all;
    std::vector<double> Tuz_abs_all;
    std::vector<double> Tlb_all;
    std::vector<double> Eu_all;

    std::vector<double> D_edge_all;
    //std::vector<std::vector<double>> Gij_all;
    std::vector<std::vector<double>> Qij_all;
    //std::vector<std::vector<double>> rhor_all;
    //std::vector<std::vector<double>> uucdis_all;
    std::vector<std::vector<double>> un2dis_all;
    std::vector<std::vector<double>> un2thetadis_all;
    std::vector<std::vector<double>> dA2H2dis_all;
    std::vector<std::vector<double>> twoHdis_all;
    std::vector<std::vector<double>> dAdis_all;


    double bead_accept = 0;
    double spin_accept = 0;
    double bond_accept = 0;
    double edge_accept = 0;
    double hop_accept = 0;

    int seq_n = 0;
    std::cout << "l0=" << l0 << "\n";
    //std::cout << "Epar.lam=" << Epar.B << "\n";
    //std::cout << "Epar.B=" << Epar.B << "\n";

    std::clock_t c_start = std::clock();
    for (int sweep_n = 0; sweep_n < MC_sweeps; sweep_n++)
    {

        if(seq && (sweep_n*seq % MC_sweeps == 0)){
            seq_n = int(sweep_n*seq/MC_sweeps);
            std::cout<<"seq "<< seq_n<<"\n";
            State_write(folder + "/State_" + finfo + "_seq"+std::to_string(seq_n)+".csv");
        }
        std::cout << sweep_n << "/" << MC_sweeps << "\n";
        //std::cout << "Eu/sqrt{N}=" << Ob_sys.Eu / std::sqrt(N) << "\n";
        E_all.push_back(Ob_sys.E);
        //std::cout << "E=" << Ob_sys.E << "\n";
        I2H2_all.push_back(Ob_sys.I2H2);
        I2H2dis_all.push_back(Ob_sys.I2H2dis);
        for (int e = 0; e < Ne; e++)
        {
            Les_all[e].push_back(Ob_sys.Les[e]);
        }
        Tp2uu_all.push_back(Ob_sys.Tp2uu);
        Tuuc_all.push_back(Ob_sys.Tuuc);
        Tun2_all.push_back(Ob_sys.Tun2);
        IdA_all.push_back(Ob_sys.IdA);
        I2H_all.push_back(Ob_sys.I2H);
        IK_all.push_back(Ob_sys.IK);
        //IKphi2_all.push_back(Ob_sys.IKphi2);
        Bond_num_all.push_back(Ob_sys.Bond_num);
        Tuz2_all.push_back(Ob_sys.Tuz2);
        Tuz_abs_all.push_back(Ob_sys.Tuz_abs);
        Tlb_all.push_back(Ob_sys.Tlb);

        if (sweep_n % sweep_p_G == 0)
        {
            Qij_all.push_back(Qij_m());
            D_edge_all.push_back(D_edge_com_m());
            //rhor_all.push_back(rho_rcom_m(delta_r, bin_num));
            //uucdis_all.push_back(uucdis_m(bin_num));
            //un2dis_all.push_back(un2dis_m(bin_num));
            //un2thetadis_all.push_back(un2thetadis_m(bin_num));
            //dA2H2dis_all.push_back(dA2H2dis_m(bin_num));
            //twoHdis_all.push_back(twoHdis_m(bin_num));
            //dAdis_all.push_back(dAdis_m(bin_num));
        }
        for (int i = 0; i < step_p_sweep; i++)
        {

            bead_accept += bead_metropolis(delta_s);
            spin_accept += spin_metropolis(delta_theta);
            bond_accept += bond_metropolis();
            bond_accept += bond_metropolis();
            //hop_accept += hop_metropolis();
            if (i % int(std::sqrt(N)) == 0)
            {
                edge_accept += edge_metropolis();
                //edge_accept += lifted_edge_metropolis();
            }
        }
    }

    bead_accept /= MC_sweeps * step_p_sweep;
    spin_accept /= MC_sweeps * step_p_sweep;
    bond_accept /= MC_sweeps * step_p_sweep * 2;
    hop_accept /= MC_sweeps * step_p_sweep;
    edge_accept /= MC_sweeps * (step_p_sweep / int(std::sqrt(N)));

    std::clock_t c_end = std::clock();
    std::ofstream f(folder + "/O_MC_" + finfo + ".csv");
    if (f.is_open())
    {
        f << "beta=" << beta << "\n";
        f << "N=" << N << "\n";
        f << "imod=" << imod << "\n";
        f << "Ne=" << Ne << "\n";
        f << "l0=" << l0 << "\n";
        f << "step_p_sweep=" << step_p_sweep << "\n";
        f << "MC_sweeps=" << MC_sweeps << "\n";
        f << "CPUtime=" << double(c_end - c_start) * 1.0 / CLOCKS_PER_SEC
          << "\n";
        f << "bead_accept," << bead_accept << "\n";
        f << "spin_accept," << spin_accept << "\n";
        f << "bond_accept," << bond_accept << "\n";
        f << "hop_accept," << hop_accept << "\n";
        f << "edge_accept," << edge_accept << "\n";
        f << "E,";
        for (int e = 0; e < Ne; e++)
        {
            f << "Les[" << e << "],";
        }
        /*
        for (int e = 0; e < Ne; e++)
        {
            f << "Ik2s[" << e << "],";
        }*/
        //f << "IdA,I2H,I2H2,phi_sum,Tphi2,I2H2dis,IK,IKphi2,Tp2uu,Tuuc,Bond_num,Tun2\n";
        f << "IdA,I2H,I2H2,I2H2dis,IK,Tp2uu,Tuuc,Bond_num,Tun2,Tuz2,Tuz_abs,Tlb\n";
        for (int i = 0; i < E_all.size(); i++)
        {
            f << E_all[i] << ",";
            for (int e = 0; e < Ne; e++)
            {
                f << Les_all[e][i] << ",";
            }
            f << IdA_all[i] << "," << I2H_all[i] << "," << I2H2_all[i] << "," << I2H2dis_all[i] << "," << IK_all[i] << "," << Tp2uu_all[i] << "," << Tuuc_all[i] << "," << Bond_num_all[i] << "," << Tun2_all[i] << "," << Tuz2_all[i] << "," << Tuz_abs_all[i] << "," << Tlb_all[i] << "\n";
        }
    }
    f.close();

    std::ofstream fQ(folder + "/Qij_" + finfo + ".csv");
    if (fQ.is_open())
    {
        fQ << "Qxx,Qxy,Qxz,Qyx,Qyy,Qyz,Qzx,Qzy,Qzz\n";
        for (int i = 0; i < Qij_all.size(); i++)
        {
            fQ << Qij_all[i][0];
            for (int j = 1; j < Qij_all[0].size(); j++)
            {
                fQ << "," << Qij_all[i][j];
            }
            fQ << "\n";
        }
    }
    fQ.close();

    /*
    std::ofstream fG(folder + "/Gij_" + finfo + ".csv");
    if (fG.is_open())
    {
        fG << "D_edge_com,Gxx,Gxy,Gxz,Gyx,Gyy,Gyz,Gzx,Gzy,Gzz\n";
        for (int i = 0; i < Gij_all.size(); i++)
        {
            fG << D_edge_all[i];
            for (int j = 0; j < Gij_all[0].size(); j++)
            {
                fG << "," << Gij_all[i][j];
            }
            fG << "\n";
        }
    }
    fG.close();
*/
/*
    std::ofstream fuuc(folder + "/uucdis_" + finfo + ".csv");
    if (fuuc.is_open())
    {
        for (int i = 0; i < uucdis_all.size(); i++)
        {
            fuuc << uucdis_all[i][0];
            for (int j = 1; j < uucdis_all[i].size(); j++)
            {
                fuuc << "," << uucdis_all[i][j];
            }
            fuuc << "\n";
        }
    }
    fuuc.close();
*/
    /*
    std::ofstream frhor(folder + "/rhor_" + finfo + ".csv");
    if (frhor.is_open())
    {
        for (int i = 0; i < rhor_all.size(); i++)
        {
            frhor << rhor_all[i][0];
            for (int j = 1; j < rhor_all[i].size(); j++)
            {
                frhor << "," << rhor_all[i][j];
            }
            frhor << "\n";
        }
    }
    frhor.close();
    */
   /*
   std::ofstream fun2(folder + "/un2dis_" + finfo + ".csv");
    if (fun2.is_open())
    {
        for (int i = 0; i < un2dis_all.size(); i++)
        {
            fun2 << un2dis_all[i][0];
            for (int j = 1; j < un2dis_all[i].size(); j++)
            {
                fun2 << "," << un2dis_all[i][j];
            }
            fun2 << "\n";
        }
    }
    fun2.close();
    */

   /*
   std::ofstream fun2theta(folder + "/un2thetadis_" + finfo + ".csv");
    if (fun2theta.is_open())
    {
        for (int i = 0; i < un2thetadis_all.size(); i++)
        {
            fun2theta << un2thetadis_all[i][0];
            for (int j = 1; j < un2thetadis_all[i].size(); j++)
            {
                fun2theta << "," << un2thetadis_all[i][j];
            }
            fun2theta << "\n";
        }
    }
    fun2theta.close();
    */
    /*
    std::ofstream fdA2H2(folder + "/dA2H2dis_" + finfo + ".csv");
    if (fdA2H2.is_open())
    {
        for (int i = 0; i < dA2H2dis_all.size(); i++)
        {
            fdA2H2 << dA2H2dis_all[i][0];
            for (int j = 1; j < dA2H2dis_all[i].size(); j++)
            {
                fdA2H2 << "," << dA2H2dis_all[i][j];
            }
            fdA2H2 << "\n";
        }
    }
    fdA2H2.close();
*/
    /*
    std::ofstream f2H(folder + "/2Hdis_" + finfo + ".csv");
    if (f2H.is_open())
    {
        for (int i = 0; i < twoHdis_all.size(); i++)
        {
            f2H << twoHdis_all[i][0];
            for (int j = 1; j < twoHdis_all[i].size(); j++)
            {
                f2H << "," << twoHdis_all[i][j];
            }
            f2H << "\n";
        }
    }
    f2H.close();
    */
   /*
    std::ofstream fdA(folder + "/dAdis_" + finfo + ".csv");
    if (fdA.is_open())
    {
        for (int i = 0; i < dAdis_all.size(); i++)
        {
            fdA << dAdis_all[i][0];
            for (int j = 1; j < dAdis_all[i].size(); j++)
            {
                fdA << "," << dAdis_all[i][j];
            }
            fdA << "\n";
        }
    }
    fdA.close();
    */
}
