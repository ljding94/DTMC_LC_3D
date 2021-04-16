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
    int N; // number of beads
    int imod = 0;
    //1 for rhombus, 2 disk, 3 cylinder
    int Ne;          // number of edges
    double d0 = 1.5; // initialization distance between neighboring bead (for rhombus shape)
    double l0 = 1.73;
    double kar;
    double karg;
    double lam;
    double Kd;
    double q; // Kt = Kd*q
    double Kt;
    double Cn;
    double kard;
    // double delta_s = 0.2;
    double delta_s = 0.1;
    double delta_theta = 0.5;

    std::string folder;

    N = std::atoi(argv[1]);
    imod = std::atoi(argv[2]); // mode of initialization
    Ne = std::atoi(argv[3]);
    kar = std::atof(argv[4]);
    lam = std::atof(argv[5]);
    Kd = std::atof(argv[6]);
    q = std::atof(argv[7]);
    Cn = std::atof(argv[8]);
    kard = std::atof(argv[9]);
    Kt = q * Kd;

    std::string finfo;
    finfo = "N" + std::string(argv[1]) + "_imod" + std::string(argv[2]) + "_Ne" + std::string(argv[3]) + "_kar" + std::string(argv[4]) + "_lam" + std::string(argv[5]) + "_Kd" + std::string(argv[6]) + "_q" + std::string(argv[7]) + "_Cn" + std::string(argv[8]) + "_kard" + std::string(argv[9]);
    dtmc_lc membrane(beta, N, imod, Ne, d0, l0, kar, lam, Kd, Kt, Cn, kard);
    N = membrane.mesh.size();
    if (argc == 11)
    {
        // use "triangulation kar lam local" for local running
        // used for local running!
        std::cout << "running on local machine\n";
        folder = "../data/scratch_local";

        membrane.State_write(folder + "/State_" + finfo + "_init.txt");
        membrane.Thermal(0, int(N / (delta_s * delta_s)) + 1, 2, delta_s, delta_theta);
        // membrane.O_MC_measure(5, 1, int(N / (ds * ds)) + 1, ds,
        membrane.O_MC_measure(20, 5, int(N / (delta_s * delta_s)) + 1, delta_s, delta_theta, folder, finfo);
        // membrane.O_MC_measure(2, 1, 0, delta_s, delta_theta, folder,
        // finfo,bin_num_r);
        membrane.State_write(folder + "/State_" + finfo + ".txt");

        return 0;
    }
    else if (argc == 10)
    {
        // ccv running
        folder = "/users/lding3/scratch";
        // membrane.Thermal(500, N / (delta_s * delta_s), 10,
        // delta_s,delta_theta);
        // used 2000, 4000 for manuscript
        membrane.Thermal(1000, N / (delta_s * delta_s), 1, delta_s,
                         delta_theta);
        membrane.O_MC_measure(2000, 50, N / (delta_s * delta_s) + 1, delta_s,
                              delta_theta, folder, finfo);
        membrane.State_write(folder + "/State_" + finfo + ".txt");

        return 0;
    }
}