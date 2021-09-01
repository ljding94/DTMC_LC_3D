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
    double d0 = 1.4; // initialization distance between neighboring bead (for rhombus shape)
    double l0 = 1.73;
    double l1 = 1.73;
    //double l0 = 1.65;
    //double l1 = 1.85; // l1<l0*(4-l0^2)
    double delta_s = 0.1;
    double delta_theta = 0.5;
    std::string folder;
    N = std::atoi(argv[1]);
    imod = std::atoi(argv[2]); // mode of initialization
    Ne = std::atoi(argv[3]);
    E_parameter Epar;
    Epar.kar = std::atof(argv[4]);
    Epar.J = std::atof(argv[5]);
    Epar.C0 = std::atof(argv[6]);
    Epar.karg = std::atof(argv[7]);
    Epar.lam = std::atof(argv[8]);
    Epar.B = std::atof(argv[9]);
    Epar.Kd = std::atof(argv[10]);
    //Epar.Ksb = std::atof(argv[6]);
    //Epar.Kt = std::atof(argv[7]);
    Epar.q = std::atof(argv[11]);
    Epar.Cn = std::atof(argv[12]);
    //Epar.kard = std::atof(argv[10]);
    //Epar.lamd = std::atof(argv[11]);

    //Epar.ms = std::atof(argv[9]);
    //Epar.mr = std::atof(argv[10]);
    //kard = std::atof(argv[9]);

    dtmc_lc membrane(beta, N, imod, Ne, d0, l0, l1, Epar);
    N = membrane.mesh.size();
    std::string finfo = "N" + std::to_string(N) + "_imod" + std::string(argv[2]) + "_Ne" + std::string(argv[3]) + "_kar" + std::string(argv[4]) + "_J" + std::string(argv[5]) + "_C0" + std::string(argv[6]) + "_karg" + std::string(argv[7]) + "_lam" + std::string(argv[8]) + "_B" + std::string(argv[9]) + "_Kd" + std::string(argv[10]) + "_q" + std::string(argv[11]) + "_Cn" + std::string(argv[12]);
    if (argc == 14)
    {
        // use "prog name par* local" for local running
        // used for local running!
        std::cout << "running on local machine\n";
        folder = "../data/scratch_local";

        membrane.State_write(folder + "/State_" + finfo + "_init.csv");
        //membrane.Thermal(20, int(N / (delta_s * delta_s)) + 1, 1, delta_s, delta_theta);
        membrane.O_MC_measure(20, 5, int(N / (delta_s * delta_s)) + 1, delta_s, delta_theta, folder, finfo);
        membrane.State_write(folder + "/State_" + finfo + ".csv");

        return 0;
    }
    else if (argc == 13)
    {
        // ccv running
        folder = "/users/lding3/scratch";
        // used 2000, 4000 for manuscript
        membrane.Thermal(1000, int(N / (delta_s * delta_s)), 100, delta_s,
                         delta_theta);
        membrane.O_MC_measure(2000, 10, int(N / (delta_s * delta_s)) + 1, delta_s,
                              delta_theta, folder, finfo);
        membrane.State_write(folder + "/State_" + finfo + ".csv");

        return 0;
    }
}