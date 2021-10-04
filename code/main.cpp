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
    int N;           // number of beads
    int imod = 1;    //1 for rhombus, 2 disk, 3 cylinder, 4 for mobius
    int Ne;          // number of edges
    double lf;          // height of cylinder init shape, used for fix beads distance at different edges (pulling experiment, when !=0)

    double d0 = 1.4; // initialization distance between neighboring bead (for rhombus shape)
    double l0 = 1.73;
    double delta_s = 0.2;
    double delta_theta = 0.5;
    double delta_r = 0.2;
    double bin_num = 80;

    std::string folder;
    N = std::atoi(argv[1]);
    imod = std::atoi(argv[2]); // mode of initialization
    Ne = std::atoi(argv[3]);
    lf = std::atof(argv[4]); // for N=400, lf=17-36 works
    E_parameter Epar;
    Epar.kar = std::atof(argv[5]);
    //Epar.C0 = std::atof(argv[6]);
    Epar.karg = std::atof(argv[6]);
    Epar.lam = std::atof(argv[7]);
    Epar.Kd = std::atof(argv[8]);
    //Epar.Ksb = std::atof(argv[6]);
    //Epar.Kt = std::atof(argv[7]);
    Epar.q = std::atof(argv[9]);
    Epar.Cn = std::atof(argv[10]);
    //Epar.kard = std::atof(argv[10]);
    //Epar.lamd = std::atof(argv[11]);

    //Epar.ms = std::atof(argv[9]);
    //Epar.mr = std::atof(argv[10]);
    //kard = std::atof(argv[9]);
    dtmc_lc membrane(beta, N, imod, Ne, lf, d0, l0, Epar);
    N = membrane.mesh.size();

    std::string finfo = "N" + std::to_string(N) + "_imod" + std::string(argv[2]) + "_Ne" + std::string(argv[3]) + "_lf" + std::string(argv[4]) + "_kar" + std::string(argv[5]) + "_karg" + std::string(argv[6]) + "_lam" + std::string(argv[7]) + "_Kd" + std::string(argv[8]) + "_q" + std::string(argv[9]) + "_Cn" + std::string(argv[10]);

    if (argc == 12)
    {
        // use "prog name par* local" for local running
        // used for local running!
        std::cout << "running on local machine\n";
        folder = "../data/scratch_local";

        membrane.State_write(folder + "/State_" + finfo + "_init.csv");
        //membrane.Thermal(20, int(N / (delta_s * delta_s)) + 1, 1, delta_s, delta_theta);
        membrane.O_MC_measure(20, 5, int(N / (delta_s * delta_s)) + 1, delta_s, delta_theta, delta_r, bin_num, folder, finfo);
        membrane.State_write(folder + "/State_" + finfo + ".csv");

        return 0;
    }
    else if (argc == 11)
    {
        // ccv running
        folder = "/users/lding3/scratch";
        // used 2000, 4000 for manuscript
        membrane.Thermal(1000, int(N / (delta_s * delta_s)), 100, delta_s,
                         delta_theta);
        membrane.O_MC_measure(2000, 10, int(N / (delta_s * delta_s)) + 1, delta_s,
                              delta_theta, delta_r,bin_num, folder, finfo);
        membrane.State_write(folder + "/State_" + finfo + ".csv");

        return 0;
    }
}