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
    // total energy
    double E;
    // geometric
    double I2H2; // integral of dA (2H)^2
    //double Iphi;              // sum of phi
    //double Tphi2;             // sum of phi*phi for near sites
    //double I2H2dis;           // integral of dA (2H - 2H0)^2, dis = displacement, where 2H0 = phi*C0
    double IK; // integral of dA(K)
    //double IKphi2;            // integral of dA(K*phi^2)
    std::vector<double> Les; // list of edge length
    // std::vector<double> Ik2s; // list of edge bending
    //std::vector<double> Leuns; // edge length couples with u cdot t (edge tangent)
    // crystalline
    double Tp2uu;
    double Tuuc;
    // coupling
    double Tun2; // tilt coupling
    // miscellany, not directly related to system energy
    double IdA;   // integral of dA
    double I2H;   // integral of dA(2H)
    int Bond_num; // total number of bonds0
};
//hamiltonion parameters
struct E_parameter
{
    double kar; // mean curvature bending stiffness
    //double J;    // side-side phi field Ising-like intereaction
    //double C0;   // spontaneous absolute mean curvature (introduced by the hooping of short rods)
    double karg; // Gaussian curvature bending stiffness
    double lam;  // line tension coefficient
    double Kd;   // liquid crystal interaction moduli
    //double Ksb;  // liquid crystalline interaction for splay and bend
    //double Kt;   //liquid crystalline interaction for the twist only
    double q;  // liquid crystall twist constant
    double Cn; // liquid crystal to membrane normal moduli
};
struct vertex
{
    // configuration related
    std::vector<double> R{0, 0, 0}; // position (x,y,z)
    std::vector<double> u{0, 0, 1}; // rod orientation (ux,uy,uz)
    std::vector<double> n{0, 0, 0}; // membrane normal
    std::vector<int> nei;           // index of neighbors
    // nei[k+1],nei[k] are connected!!!
    int edge_num; // which edge
    std::vector<int> edge_nei;
    // neighbors form edge with this one (if exist)

    // measurement related
    //double phi; // local 2H0 = phi*C0, phi\in (-1,1) // also karg phi
    // interaction among phi field can be added as need
    std::vector<double> dAn2H; // in bulk: (dA, 2H), on edge (0,0)
    // energy related (directly)
    double ds; // beads in bulk: 0 , beads on edge 0.5*(l_{i+}+l_{i-})
    //double dsk2; // edge bending
    // double dskg; // in bulk: 0, on edge: kg*ds
    double dAK; // in bulk: K*dA, on edge: 0
    double un2; // local un2
    //double es;  // energy related strength, 1 for normal rods ms for mixture
    //double dE;  // local energy; sum dE over all beads is the actually total energy.
};
class dtmc_lc
{
public:
    // eternal parameters
    double beta; // system temperature
    int N;       // number of beads
    int imod;    // mode for initialization shape
    int Ne;      // number of edges
    int lf;      //fixed distance along pulling direcion (z for cylinder)
    // also used to set fixed distance betwen two beads
    double l0; // in-bulk tether maximum length

    // system configuration

    std::vector<vertex> mesh; // the mesh
    std::vector<std::vector<int>>
        edge_lists; // list of the edge beads sorted, linked
    std::vector<std::pair<int, int>> bulk_bond_list;
    // bulk_bond_list[i] is a pair of two bead connected in [bulk!}
    // notice only bond in the bulk are included, and there will be on
    // repetition due to symmetry,
    std::vector<int> fixed_beads;   // beads can't be moved
    std::vector<int> fixed_beads_z; // beads can't be moved in the z direction for imod3 cylinder initial shape
    std::vector<double> edge_zlim;
    void mesh_bead_info_update(std::vector<int> ind_relate);

    E_parameter Epar;
    observable Ob_sys, Ob_sys_w; //system observables and mixture weighted observables, later one is directly related to the energy measurement

    observable Ob_m(std::vector<int> ind_relate,
                    std::vector<std::pair<int, int>> bond_relate);
    void Ob_init(observable &Ob);
    void Ob_sys_update(observable Ob_new, observable Ob_old);

    // measure observables related to these ind and bond
    double E_m(observable Ob);

    // randomnumber generators
    std::mt19937 gen;
    std::uniform_int_distribution<> rand_pos;  // random position
    std::uniform_real_distribution<> rand_uni; // uniform distribution

    // initialization
    dtmc_lc(double beta_, int N_, int imod_, int Ne_, double lf_, double d0_, double l0_, E_parameter Epar_);
    //double kar_, double lam_, double Kd_, double Kt_, double Cn_,double kard_);

    // put the beads and bonds in to position accordingly
    void init_rhombus_shape(double d0_);
    void init_disk_shape(double d0_);
    void init_cylinder_shape(double d0_, double lf_); // lf will over write d0_ when lf!=0
    void init_mobius_shape(double d0_);
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
    double ds_m(int index); // length of the local edge bead
    //double dsk2_m(int index);               // edge bending of local edge bead ds(index) needed
    double ut_m(int index);                 // director edge tangent sine angle
    std::vector<double> dAn2H_m(int index); // measure and set dA and |2H|

    double dAK_m(int index);             // K*dA measure the gauss curvature
    double dskg_m(int index);            // kg*ds, for gaussian
    std::vector<double> n_m(int index);  // surface normal
    double p2uu_m(int ind_i, int ind_j); // spin spin interaction of i and j
    double uuc_m(int ind_i, int ind_j);  // spin twist interaction of i and j
    double uusb_m(int ind_i, int ind_j); // splay bend interaction of i and j
    double uut_m(int ind_i, int ind_j);  // twist interaction of i and j
    double un2_m(int index);             // spin normal interaction of i
    double dEgeo_m(int index);           // energy of the local vertex

    // Structure related measurements
    // Gyration tensor measurement
    std::vector<double> Gij_m();
    // seperation between 2 edges measurement
    double D_edge_com_m();
    // normalized twist from center of mass measurement
    // (u(r)*nu)^2, how mush director twist about membrane nematic director
    std::vector<double> un2dis_m(int bin_num);
    // distribution of un2 among beads, good indication for pi wall formation

    std::vector<double> comR_m();                              // center of mass measurement
    std::vector<double> rho_rcom_m(double del_r, int bin_num); // density distribution from center of mass

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

    int hop_metropolis();
    // nearing short - long rods exchange update
    // return 1: accepted, 0: rejected

    //int swap_metropolis(); // reserved for mixture which is not of my concern for now
    // swap two beads with different interaction parameters, like cn
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
                      double delta_s, double delta_theta, double delta_r, double bin_num, std::string folder,
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
