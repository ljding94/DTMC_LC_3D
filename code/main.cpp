// Copyright[2021] [Lijie Ding]
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "dtmc_lc_3d.h"

int main(int argc, char const *argv[])
{
    std::clock_t c_start = std::clock();
    double beta = 1;
    int N;        // number of beads
    int imod = 1; //1 for rhombus, 2 disk, 3 cylinder, 4 for mobius
    int Ne;       // number of edges
    double lf;    // height of cylinder init shape, used for fix beads distance at different edges (pulling experiment, when !=0)

    double d0 = 1.4; // initialization distance between neighboring bead (for rhombus shape)
    double l0 = 1.73;
    double delta_s = 0.1;
    double delta_theta = 0.5;
    double delta_r = 0.1;
    double bin_num = 80;
    //int nEu = 3; // number of harmonic oscilattors for umbrella sampling bias energy Eu

    std::string folder;
    N = std::atoi(argv[1]);
    imod = std::atoi(argv[2]); // mode of initialization
    Ne = std::atoi(argv[3]);
    lf = std::atof(argv[4]); //
    E_parameter Epar;
    Epar.kar = std::atof(argv[5]);
    Epar.C0 = std::atof(argv[6]);
    Epar.karg = std::atof(argv[7]);
    Epar.lam = std::atof(argv[8]);
    Epar.Kd = std::atof(argv[9]);
    //Epar.Ksb = std::atof(argv[6]);
    //Epar.Kt = std::atof(argv[7]);
    Epar.q = std::atof(argv[10]);
    Epar.Cn = std::atof(argv[11]);
    //Epar.ku = std::atof(argv[12]);
    //Epar.n_Eu = nEu;
    // Epar.g = std::atof(argv[12]);
    //Epar.kard = std::atof(argv[10]);
    //Epar.lamd = std::atof(argv[11]);

    std::cout << "argc = " << argc << "\n";
    //Epar.ms = std::atof(argv[9]);
    //Epar.mr = std::atof(argv[10]);
    //kard = std::atof(argv[9]);
    dtmc_lc membrane(beta, N, imod, Ne, lf, d0, l0, Epar);
    N = membrane.mesh.size();

    std::string finfo = "N" + std::to_string(N) + "_imod" + std::string(argv[2]) + "_Ne" + std::string(argv[3]) + "_lf" + std::string(argv[4]) + "_kar" + std::string(argv[5]) + "_C0" + std::string(argv[6]) + "_karg" + std::string(argv[7]) + "_lam" + std::string(argv[8]) + "_Kd" + std::string(argv[9]) + "_q" + std::string(argv[10]) + "_Cn" + std::string(argv[11]); // + "_ku" + std::string(argv[12]);

    if (argc == 13)
    {
        // use "prog name par* local" for local running
        // used for local running!
        std::cout << "running on local machine\n";
        folder = "../data/scratch_local";

        membrane.State_write(folder + "/State_" + finfo + "_init.csv");
        //membrane.Thermal_pinch(200, int(N / (delta_s * delta_s)) + 1, 2.0, delta_s, delta_theta);
        if(lf==0){
            membrane.Thermal_kar1lam1(100, int(N / (delta_s * delta_s)) + 1, 10, 20, delta_s, delta_theta);
            membrane.State_write(folder + "/State_" + finfo + "_therm.csv");
        }
        membrane.Thermal(100, int(N / (delta_s * delta_s)) + 1, 10, delta_s, delta_theta);
        membrane.O_MC_measure(100, 50, int(N / (delta_s * delta_s)) + 1, delta_s, delta_theta, delta_r, bin_num, folder, finfo);
        membrane.State_write(folder + "/State_" + finfo + ".csv");

        return 0;
    }
    else if (argc == 12)
    {
        // ccv running
        folder = "/users/lding3/scratch/dtmc_lc_3d";
        // used 2000, 4000 for manuscript
        //membrane.Thermal_pinch(1000, int(N / (delta_s * delta_s)) + 1, 2.0, delta_s, delta_theta);
        //membrane.State_write(folder + "/State_" + finfo + "_therm.csv");
        if(lf==0){
            membrane.Thermal_kar1lam1(2000, int(N / (delta_s * delta_s)) + 1, 10, 20, delta_s, delta_theta);
            membrane.State_write(folder + "/State_" + finfo + "_therm.csv");
        }
        //membrane.Thermal_kar1lam1(2000, int(N / (delta_s * delta_s)) + 1, 10, 20, delta_s, delta_theta);
        membrane.Thermal(5000, int(N / (delta_s * delta_s)), 50, delta_s,delta_theta);
        // see how it evolves for now
        membrane.O_MC_measure(15000, 1000, int(N / (delta_s * delta_s)) + 1, delta_s, delta_theta, delta_r, bin_num, folder, finfo);
        membrane.State_write(folder + "/State_" + finfo + ".csv");

        return 0;
    }
}