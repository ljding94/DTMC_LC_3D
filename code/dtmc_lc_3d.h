#ifndef _DTMC_LC_3D_H
#define _DTMC_LC_3D_H
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

// observable
struct observable
{
    // energy
    double E;
    // geometric
    double I2H2;
    std::vector<double> Les;
    // crystalline
    double Tp2uu;
    double Tuuc;
    // coupling
    double Tun2;
    double IKun2;
    // miscellany
    double IdA; // integral of dA
    double I2H; // integral of dA(2H)
    double IK;  // integral of dA(K)

    int Bond_num; // total number of bonds
};
struct vertex
{
    // configuration related
    std::vector<double> R{0, 0, 0}; // position (x,y,z)
    std::vector<double> u{0, 0, 1}; // rod orientation (ux,uy,uz)
    std::vector<double> n{0, 0, 0}; // membrane normal

    std::vector<int> nei; // index of neighbors
    // nei[k+1],nei[k] are connected!!!
    int edge_num; // which edge
    std::vector<int> edge_nei;
    // neighbors form edge with this one (if exist)

    std::vector<int> nei2flip; // index in nei that could be flipped with

    // measurement related
    std::vector<double> dAn2H; // in bulk: (dA, 2H), on edge (0,0)
    // energy related (directly)
    double ds; // beads in bulk: 0 , beads on edge 0.5*(l_{i+}+l_{i-})
    // double dskg; // in bulk: 0, on edge: kg*ds
    double dAK; // in bulk: K*dA, on edge: 0
    double un2; // local un2
};
class dtmc_lc
{
public:
    // eternal parameters
    double beta; // system temperature
    int N;       // number of beads
    int imod;    // mode for initialization shape
    int Ne;      // number of edges
    // also used to set fixed distance betwen two beads
    double l0;   // tether maximum length
    double kar;  // mean curvature bending stiffness
    double lam;  // line tension coefficient
    double Kd;   // liquid crystal interaction moduli
    double Kt;   // liquid crystall twist interaction moduli
    double Cn;   // liquid crystal to membrane normal moduli
    double kard; // depletion-like gaussian curvature bending stiffness

    // system configuration

    std::vector<vertex> mesh; // the mesh
    std::vector<std::vector<int>>
        edge_lists; // list of the edge beads sorted, linked
    std::vector<std::pair<int, int>> bulk_bond_list;
    // bulk_bond_list[i] is a pair of two bead connected in [bulk!}
    // notice only bond in the bulk are included, and there will be on
    // repetition due to symmetry,
    std::vector<int> fixed_beads; // beads can't be moved

    void mesh_bead_info_update(std::vector<int> ind_relate);

    observable Ob_sys;
    observable Ob_m(std::vector<int> ind_relate,
                    std::vector<std::pair<int, int>> bond_relate);
    void Ob_sys_update(observable Ob_new, observable Ob_old);

    // measure observables related to these ind and bond
    double E_m(observable Ob);

    // randomnumber generators
    std::mt19937 gen;
    std::uniform_int_distribution<> rand_pos;  // random position
    std::uniform_real_distribution<> rand_uni; // uniform distribution

    // initialization
    dtmc_lc(double beta_, int N_, int imod_, int Ne_, double d0_, double l0_,
            double kar_, double lam_, double Kd_, double Kt_, double Cn_,
            double kard_);

    // put the beads and bonds in to position accordingly
    void init_rhombus_shape(double d0_);
    void init_disk_shape(double d0_);
    void init_cylinder_shape(double d0_);
    int add_hole_as_edge(int b0, int edgenum); // return 0 if fail to do so with b0

    void reset_config();
    void push_neis_back(int i, std::vector<int> nei_dist);
    void push_eneis_back(int i, std::vector<int> enei_dist);
    void push_bneis_list(int i, std::vector<int> bnei_dist);

    // bond in bulk to bulk_bond_list
    // d0_ initial vertex-vertex distance
    // other parameters are the same as in the class

    // local energy-related measurement
    // _m stand for measurement
    double ds_m(int index);                 // length of the local edge index
    std::vector<double> dAn2H_m(int index); // measure and set dA and |2H|
    double dAK_m(int index);                // K*dA measure the gauss curvature
    double dskg_m(int index);               // kg*ds, for gaussian
    std::vector<double> n_m(int index);     // surface normal
    double p2uu_m(int ind_i, int ind_j);    // spin spin interaction of i and j
    double uuc_m(int ind_i, int ind_j);     // spin twist interaction of i and j
    double un2_m(int index);                // spin normal interaction of i
    double dEgeo_m(int index);              // energy of the local vertex

    // Gyration tensor measurement
    std::vector<double> Gij_m();
    // seperation between 2 edges measurement
    double D_edge_com_m();
    // normalized twist from center of mass measurement
    // (u(r)*nu)^2, how mush director twist about membrane nematic director
    std::vector<double> un2dis_m(int bin_num);
    // distribution of un2 among beads, good indication for pi wall formation

    // useful tools
    double distance2(int ind_1, int ind_2);
    double distance(int ind_1, int ind_2);
    double distancefp(int ind_1, std::vector<double> p);
    double innerproduct(std::vector<double> a, std::vector<double> b);
    void delete_bulk_bond_list(int ind_i, int ind_j);
    // delete bond i-j, including both <i,j> and <j,i> in the bulk bulk_bond_list

    observable get_related_local_observables(std::vector<int> ind_list);

    // Monte Carlo updates
    int bead_metropolis(double delta_s);
    // execute bead move update in [-ds,ds]^3
    // return 1: accepted, 0: rejected

    int spin_metropolis(double delta_theta);
    // execute spin move update along random axis in [-dtheta,dtheta]
    // return 1: accepted, 0: rejected 0

    int bond_metropolis();
    // execute bond switch update
    // return 1: accepted, 0: rejected

    int edge_metropolis();
    // execute bond remove/add update
    // return 1: accepted, 0: rejected

    // experiment
    void State_write(std::string filename);
    // write system state to file
    void State_write_seq(std::string filename, int MC_sweeps, int step_p_sweep,
                         double delta_s, double delta_theta);

    void State_load(std::string state_file);
    // load system state from file
    void Thermal(int MC_sweeps, int step_p_sweep, int beta_steps,
                 double delta_s, double delta_theta);
    // thermalisation of the system, starting from beta=0 (hot start)
    void O_MC_measure(int MC_sweeps, int sweep_p_G, int step_p_sweep,
                      double delta_s, double delta_theta, std::string folder,
                      std::string finfo);
    // measure the obserables
    // energy versus lambda curve testing

    // little tool
    int energy_check(); // check if the energy set are correct, for debug use
    int sort_nei(int index);
    int list_a_nei_b(std::vector<int> a, int b);
    int if_near_edge(int b); // check if b is connected to a edge bead
    int check_nei_connect();
    int check_duplication(int ind_i);
};
#endif
