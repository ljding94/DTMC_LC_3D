#include "dtmc_lc_3d.h"
#include <cmath>
#include <iostream>
#include <random>
#define PI 3.14159265358979323846

// initialization
dtmc_lc::dtmc_lc(double beta_, int N_, int imod_, int Ne_, double d0_, double l0_, E_parameter Epar_)
{
    // system related
    beta = beta_;
    N = N_;
    imod = imod_;
    Ne = Ne_;
    l0 = l0_; // sigma0 is always 1;

    // set energy related
    // geometric
    Epar.kar = Epar_.kar;
    Epar.C0 = Epar_.C0;
    Epar.karg = Epar_.karg;
    Epar.lam = Epar_.lam;
    // orientational
    Epar.Kd = Epar_.Kd;
    //Epar.Ksb = Epar_.Ksb;
    //Epar.Kt = Epar_.Kt;
    Epar.q = Epar_.q;
    // coupling
    Epar.Cn = Epar_.Cn;
    //Epar.kard = Epar_.kard;
    //Epar.lamd = Epar_.lamd;
    //Epar.ms = Epar_.ms;
    //Epar.mr = Epar_.mr;

    edge_lists.resize(Ne);
    for (int n = 0; n < Ne; n++)
    {
        edge_lists[n].clear();
    }
    bulk_bond_list.clear();

    // new way, assume (N+2)%L != 0
    // if L=sqrt(N), N can't be 700
    // determine initialization shape based on imod
    int num_edge_exist; // number of edges already exist
    int hole_pos;       // position of the hole to add
    if (imod == 1)
    {
        init_rhombus_shape(d0_);
        num_edge_exist = 1;
    }
    else if (imod == 2)
    {
        init_disk_shape(d0_);
        num_edge_exist = 1;
    }
    else if (imod == 3)
    {
        init_cylinder_shape(d0_);
        num_edge_exist = 2;
    }
    else if (imod == 4)
    {
        init_mobius_shape(d0_);
        num_edge_exist = 1;
    }
    // N is updated based on initialization procedure, especially will decrease if using disk-shape
    hole_pos = N / 4; // add hole at a random position
    while (num_edge_exist < Ne)
    {
        //std::cout << "num_edge_exist," << num_edge_exist << std::endl;

        if (add_hole_as_edge(hole_pos, num_edge_exist) == 1)
        {
            // successfully added a hole >.<
            num_edge_exist += 1;
            std::cout << "successfully added a hole NO." << num_edge_exist << "\n";
            std::cout << "at hole_pos." << hole_pos << std::endl;
        }
        hole_pos += 4;
        //std::cout << "hole_pos=" << hole_pos << std::endl;
    }
    // above shape setting take care of beads position
    // bulk_bond_list (excluding edge bond) and edge_list

    // set inital observable value
    Ob_init(Ob_sys);

    for (int i = 0; i < mesh.size(); i++)
    {
        // set phiH field
        mesh[i].phiH = 1;
        // vertex info measurement
        mesh[i].n = n_m(i);
        mesh[i].dAn2H = dAn2H_m(i);
        mesh[i].ds = ds_m(i);
        mesh[i].un2 = un2_m(i);
        mesh[i].dAK = dAK_m(i);
    }
    /*
    // not study mixture right now
    for (int i = 0; i < Np; i++)
    {
        // cn assignment for first Np rods
        mesh[i].es = Epar.ms;
    }
    */
    // geometric energy related
    for (int i = 0; i < mesh.size(); i++)
    {
        Ob_sys.I2H2 += mesh[i].dAn2H[0] * mesh[i].dAn2H[1] * mesh[i].dAn2H[1];
        Ob_sys.phiH_sum += mesh[i].phiH;
        Ob_sys.I2H2dis += mesh[i].dAn2H[0] * std::pow(mesh[i].dAn2H[1] - mesh[i].phiH * Epar.C0, 2);
        Ob_sys.IK += mesh[i].dAK;
        if (mesh[i].edge_num != -1)
        {
            Ob_sys.Les[mesh[i].edge_num] += mesh[i].ds;
            //Ob_sys.Leuns[mesh[i].edge_num] += mesh[i].ds * std::sqrt(mesh[i].un2);
        }
        // crystalline energy related
        for (int j = 0; j < mesh[i].nei.size(); j++)
        {
            // factor of 1/2 is for double counting
            Ob_sys.Tp2uu += 0.5 * p2uu_m(i, mesh[i].nei[j]);
            Ob_sys.Tuuc += 0.5 * uuc_m(i, mesh[i].nei[j]);
            //Ob_sys_w.Tuuc += 0.5 * (mesh[i].es + mesh[mesh[i].nei[j]].es) * 0.5 * uuc_m(i, mesh[i].nei[j]); // was used for mixture study
            // Ising-like phiH field interaction
            Ob_sys.TphiH2 += 0.5 * mesh[i].phiH*mesh[mesh[i].nei[j]].phiH;
        }
        // coupling related
        Ob_sys.Tun2 += mesh[i].un2;
        //Ob_sys.IKun2 += mesh[i].dAK * mesh[i].un2;
        // miscellany
        Ob_sys.IdA += mesh[i].dAn2H[0];
        Ob_sys.I2H += mesh[i].dAn2H[1] * mesh[i].dAn2H[0];
    }
    //Ob_sys.E = 0.5 * Epar.Cn * (N - Np) + 0.5 * Epar.Cn * Np;
    Ob_sys.E = 0.5 * Epar.Cn * N; // offset tilt energy
    Ob_sys.E += E_m(Ob_sys);

    // set random number generators
    std::random_device rd;
    std::mt19937 gen_set(rd());
    std::uniform_int_distribution<> rand_pos_set(0, mesh.size() - 1); // random pos
    std::uniform_real_distribution<> rand_uni_set(0, 1);
    gen = gen_set;
    rand_pos = rand_pos_set;
    rand_uni = rand_uni_set;

} // end of triangulation

#pragma region : shape initialization
void dtmc_lc::init_rhombus_shape(double d0_)
{
    // works for Ne>=1, will add triangle holes if Ne>1
    int x_n, y_n; // position of the vertex in the two vector coordinate
    mesh.resize(N);
    int L = int(std::sqrt(N));
    if ((N + 2) % L == 1 || (N + 2) % L >= (L - 1))
    {
        std::cout << "invalid N!\n";
    }

    for (int i = 0; i < N; i++)
    {
        // assign position
        if (i < N - N % L - 2)
        {
            x_n = (i + 1) % L;
            y_n = (i + 1) / L;
        }
        else
        {
            // ignore upper-right corner due to #-of-nei limit
            x_n = (i + 2) % L;
            y_n = (i + 2) / L;
        }

        mesh[i].R[0] = d0_ * (x_n + 0.5 * y_n);
        mesh[i].R[1] = d0_ * 0.5 * std::sqrt(3) * y_n;
        mesh[i].R[2] = 0;

        // put bonds
        if (x_n == L - 1)
        {
            // right line
            edge_lists[0].push_back(i);
            mesh[i].edge_num = 0;
            if (y_n == 0)
            {
                // lower-right corner
                push_neis_back(i, {L, L - 1, -1});
                push_eneis_back(i, {L, -1});
                push_bneis_list(i, {L - 1});
            }
            else if (y_n == (N + 1) / L - 2)
            {
                // upper-right corner
                push_neis_back(i, {L - 1, -1, -L});
                push_eneis_back(i, {L - 1, -L});
                push_bneis_list(i, {-1});
            }
            else
            {
                push_neis_back(i, {L, L - 1, -1, -L});
                push_eneis_back(i, {L, -L});
                push_bneis_list(i, {L - 1, -1});
            }
        }
        else if (y_n == (N + 1) / L - 1 && x_n != 0)
        {
            if (x_n > (N + 1) % L)
            {
                edge_lists[0].push_back(i);
                mesh[i].edge_num = 0;
                if (x_n == L - 2)
                {
                    if ((N + 2) % L < (L - 2))
                    {
                        push_neis_back(i, {-1, -L, -L + 1});
                        push_eneis_back(i, {-1, -L + 1});
                        push_bneis_list(i, {-L});
                    }
                    else
                    {
                        push_neis_back(i, {L - 2, -1, -L, -L + 1});
                        push_eneis_back(i, {L - 2, -L + 1});
                        push_bneis_list(i, {-1, -L});
                    }
                }
                else if (x_n == (N + 1) % L + 1)
                {
                    push_neis_back(i, {L - 2, -1, -L, -L + 1, 1});
                    push_eneis_back(i, {L - 2, 1});
                    push_bneis_list(i, {-1, -L, -L + 1});
                }
                else
                {
                    push_neis_back(i, {-1, -L, -L + 1, 1});
                    push_eneis_back(i, {-1, 1});
                    push_bneis_list(i, {-L, -L + 1});
                }
            }
            else
            {
                // in the bulk
                mesh[i].edge_num = -1;
                push_neis_back(i, {1, L - 1, L - 2, -1, -L, -L + 1});
                push_bneis_list(i, {1, L - 1, L - 2, -1, -L, -L + 1});
            }
        }
        else if (y_n == (N + 1) / L)
        {
            edge_lists[0].push_back(i);
            mesh[i].edge_num = 0;
            if (x_n == (N + 1) % L)
            {
                push_neis_back(i, {-1, -L + 1, -L + 2});
                push_eneis_back(i, {-1, -L + 2});
                push_bneis_list(i, {-L + 1});
            }
            else if (x_n == 0)
            {
                push_neis_back(i, {-L + 1, -L + 2, 1});
                push_eneis_back(i, {-L + 1, 1});
                push_bneis_list(i, {-L + 2});
            }
            else
            {
                push_neis_back(i, {-1, -L + 1, -L + 2, 1});
                push_eneis_back(i, {-1, 1});
                push_bneis_list(i, {-L + 1, -L + 2});
            }
        }
        else if (x_n == 0)
        {
            edge_lists[0].push_back(i);
            mesh[i].edge_num = 0;
            if (y_n == 1)
            {
                push_neis_back(i, {-L + 1, 1, L});
                push_eneis_back(i, {-L + 1, L});
                push_bneis_list(i, {1});
            }
            else if (y_n == (N + 1) / L - 1)
            {
                push_neis_back(i, {-L, -L + 1, 1, L - 1});
                push_eneis_back(i, {-L, L - 1});
                push_bneis_list(i, {-L + 1, 1});
            }
            else
            {
                push_neis_back(i, {-L, -L + 1, 1, L});
                push_eneis_back(i, {-L, L});
                push_bneis_list(i, {-L + 1, 1});
            }
        }
        else if (y_n == 0)
        {
            edge_lists[0].push_back(i);
            mesh[i].edge_num = 0;
            if (x_n == 1)
            {
                push_neis_back(i, {1, L, L - 1});
                push_eneis_back(i, {1, L - 1});
                push_bneis_list(i, {L});
            }
            else
            {
                push_neis_back(i, {1, L, L - 1, -1});
                push_eneis_back(i, {1, -1});
                push_bneis_list(i, {L, L - 1});
            }
        }
        else
        {
            mesh[i].edge_num = -1; // in the bulk
            push_neis_back(i, {1, L, L - 1, -1, -L, -L + 1});
            push_bneis_list(i, {1, L, L - 1, -1, -L, -L + 1});
        }
    }
}

void dtmc_lc::init_disk_shape(double d0_)
{
    // works for Ne>=1, will add triangle holes if Ne>1
    int x_n, y_n;
    double disi2c, dise2c;
    int nei_pos = 0;

    int rb = int(std::sqrt(N / 4));
    std::vector<vertex> init_mesh; // initial rhombus-shape mesh
    init_mesh.resize((4 * rb + 1) * (4 * rb + 1));
    std::vector<int> init2after_index(init_mesh.size(), -1);

    // f[# in init_mesh] = # in mesh, -1 as default
    mesh.clear();
    std::vector<int> nei_dist = {1, 4 * rb + 1, 4 * rb, -1, -4 * rb - 1, -4 * rb};
    dise2c = std::pow(d0_ * rb, 2);
    for (int i = 0; i < init_mesh.size(); i++)
    {
        x_n = i % (4 * rb + 1) - 2 * rb;
        y_n = i / (4 * rb + 1) - 2 * rb; //-2r is for centralization
        init_mesh[i].R[0] = d0_ * (x_n + 0.5 * y_n);
        init_mesh[i].R[1] = d0_ * 0.5 * std::sqrt(3) * y_n;
        init_mesh[i].R[2] = 0;
        // put neighbors as if every on is in-bulk
        init_mesh[i].edge_num = -1;
        //ignore bead on init_mesh edge
        for (int j = 0; j < 6; j++)
        {
            nei_pos = i + nei_dist[j];
            if (0 <= nei_pos && nei_pos < init_mesh.size())
            {
                //std::cout << "nei_pos=" << nei_pos << "\n";
                init_mesh[i].nei.push_back(nei_pos);
            }
        }
        // find position in real mesh;

        disi2c = std::pow(init_mesh[i].R[0], 2) + std::pow(init_mesh[i].R[1], 2);
        if (disi2c <= dise2c)
        {
            mesh.push_back(init_mesh[i]);
            init2after_index[i] = mesh.size() - 1;
        }
    }
    //corresting mesh nei index
    N = mesh.size(); // assign the system size parameter
    std::vector<int> new_nei_cache;
    std::vector<int> ind_nei_cache;

    for (int i = 0; i < mesh.size(); i++)
    {
        new_nei_cache.clear();
        ind_nei_cache.clear();
        for (int j = 0; j < mesh[i].nei.size(); j++)
        {
            if (init2after_index[mesh[i].nei[j]] != -1)
            {
                new_nei_cache.push_back(init2after_index[mesh[i].nei[j]]);
                ind_nei_cache.push_back(init2after_index[mesh[i].nei[j]] - i);
            }
        } // end for mesh[i].nei

        mesh[i].nei = new_nei_cache;

        // edge bead?
        if (mesh[i].nei.size() < 6)
        {
            mesh[i].edge_num = 0;
            edge_lists[0].push_back(i);
        }
        else
        { // if in-bulk, add all bond to bond list
            push_bneis_list(i, ind_nei_cache);
        }
    } // end of mesh connection
    // time take care of beads on the edge
    // edge_nei[0] must be nei[0] edge_nei[1] must be nei[-1]
    int nei_p, ind;
    std::vector<int> enei_cache; // neighbors on the edge
    std::vector<int> bnei_cache; // neighbors in the bulk
    for (int i = 0; i < edge_lists[0].size(); i++)
    {
        ind = edge_lists[0][i];
        enei_cache.clear();
        bnei_cache.clear();
        for (int j = 0; j < mesh[ind].nei.size(); j++)
        {
            nei_p = mesh[ind].nei[j];
            if (mesh[nei_p].nei.size() < 6)
            { // nei_p is on the edge
                enei_cache.push_back(nei_p - ind);
            }
            else
            { // nei_p is in the bulk
                bnei_cache.push_back(nei_p - ind);
            }
        } // end for ind nei

        push_eneis_back(ind, enei_cache);
        // to be sorted
        push_bneis_list(ind, bnei_cache);
    }
    int eind, eind_next, eind_cache;
    eind = edge_lists[0][0];

    do
    {
        eind_next = mesh[eind].edge_nei[1];
        eind_cache = mesh[eind_next].edge_nei[0];
        if (eind_cache != eind)
        { //fix edge_nei direction of eind_next
            mesh[eind_next].edge_nei[0] = eind;
            mesh[eind_next].edge_nei[1] = eind_cache;
        }
        sort_nei(eind_next); // sort nei list based on edge-nei
        eind = eind_next;
    } while (eind != edge_lists[0][0]);
}

void dtmc_lc::init_cylinder_shape(double d0_)
{
    // only work for Ne>=2
    int x_n, y_n; // position of the vertex in the two vector coordinate
    // cylinder initial shape
    //use cylinder initialization
    int L = 10; // length of the cylinder
    if (N % L != 0)
    {
        std::cout << "N % L != 0 \n";
    }
    int Lr = N / L; // perimeter of cylinder circular bottom
    double del_theta = 2 * PI / Lr;
    double R = d0_ / (2 * std::sin(del_theta / 2));

    mesh.resize(N);
    for (int i = 0; i < N; i++)
    {
        // assign position
        x_n = i % Lr;
        y_n = i / Lr;
        mesh[i].R[0] = R * std::cos(del_theta * (x_n + 0.5 * y_n));
        mesh[i].R[1] = R * std::sin(del_theta * (x_n + 0.5 * y_n));
        mesh[i].R[2] = d0_ * 0.5 * std::sqrt(3) * y_n;

        // put bonds
        if (y_n == 0)
        {
            // lower edge
            mesh[i].edge_num = 0;
            edge_lists[0].push_back(i);
            if (x_n == 0)
            {
                push_neis_back(i, {1, Lr, 2 * Lr - 1, Lr - 1});
                push_eneis_back(i, {1, Lr - 1});
                push_bneis_list(i, {Lr, 2 * Lr - 1});
            }
            else if (x_n == Lr - 1)
            {
                push_neis_back(i, {-Lr + 1, Lr, Lr - 1, -1});
                push_eneis_back(i, {-Lr + 1, -1});
                push_bneis_list(i, {Lr, Lr - 1});
            }
            else
            {
                push_neis_back(i, {1, Lr, Lr - 1, -1});
                push_eneis_back(i, {1, -1});
                push_bneis_list(i, {Lr, Lr - 1});
            }
        }
        else if (y_n == L - 1)
        {
            // upper edge
            mesh[i].edge_num = 1;
            edge_lists[1].push_back(i);
            if (x_n == 0)
            {
                push_neis_back(i, {Lr - 1, -Lr, -Lr + 1, 1});
                push_eneis_back(i, {Lr - 1, 1});
                push_bneis_list(i, {-Lr, -Lr + 1});
            }
            else if (x_n == Lr - 1)
            {
                push_neis_back(i, {-1, -Lr, -2 * Lr + 1, -Lr + 1});
                push_eneis_back(i, {-1, -Lr + 1});
                push_bneis_list(i, {-Lr, -2 * Lr + 1});
            }
            else
            {
                push_neis_back(i, {-1, -Lr, -Lr + 1, 1});
                push_eneis_back(i, {-1, 1});
                push_bneis_list(i, {-Lr, -Lr + 1});
            }
        }
        else
        {
            // internal
            mesh[i].edge_num = -1;
            if (x_n == 0)
            {
                push_neis_back(i, {1, Lr, 2 * Lr - 1, Lr - 1, -Lr, -Lr + 1});
                push_bneis_list(i, {1, Lr, 2 * Lr - 1, Lr - 1, -Lr, -Lr + 1});
            }
            else if (x_n == Lr - 1)
            {
                push_neis_back(i, {-Lr + 1, Lr, Lr - 1, -1, -Lr, -2 * Lr + 1});
                push_bneis_list(i, {-Lr + 1, Lr, Lr - 1, -1, -Lr, -2 * Lr + 1});
            }
            else
            {
                push_neis_back(i, {1, Lr, Lr - 1, -1, -Lr, -Lr + 1});
                push_bneis_list(i, {1, Lr, Lr - 1, -1, -Lr, -Lr + 1});
            }
        }
    }
}

void dtmc_lc::init_mobius_shape(double d0_)
{
    // initialization of mobius strip
    double w, t; // continues parameter along the width and rotational direction

    //use cylinder initialization
    int W = 5; // width of mobius strip
               // for now, odd only
               //3 is very optimal consider the distance change need to be within (1,l0) when d0=1.5
               // for W=5 N>=300 at least when d0=1.4
    N -= N % W;
    int Lr = N / W;    // perimeter of strip's circular bottom
    int Lw;            // number of bead depending on w_n
    int i, inei_cache; // index in the mesh
    double del_theta = 2 * PI / Lr;
    double R = d0_ / (2 * std::sin(del_theta / 2));
    mesh.resize(N);
    for (int wn = 0; wn < W; wn++)
    {
        Lw = Lr + W / 2 - wn;
        for (int tn = 0; tn < Lw; tn++)
        {
            i = wn * (Lr + W / 2) - (wn - 1) * wn / 2 + tn;
            w = (wn - (W - 1) / 2) * std::sqrt(3) / 2;
            t = tn + 0.5 * wn;
            mesh[i].R[0] = (R + d0_ * w * std::cos(del_theta / 2 * t)) * std::cos(del_theta * t);
            mesh[i].R[1] = (R + d0_ * w * std::cos(del_theta / 2 * t)) * std::sin(del_theta * t);
            mesh[i].R[2] = w * std::sin(del_theta / 2 * t);

            //TODO: implement below putting bond code
            // put bonds
            if (wn == 0)
            {
                // edge~
                mesh[i].edge_num = 0;
                edge_lists[0].push_back(i);
                if (tn == 0)
                {
                    // i==0, connect to i=N-1
                    // left most, connecting to w_n=W
                    inei_cache = N - 1 - (Lr + W / 2 - (W - 1)); // just for index position storage
                    push_neis_back(i, {1, Lw, inei_cache - i, N - 1 - i});
                    push_eneis_back(i, {N - 1 - i, 1});
                    push_bneis_list(i, {Lw, inei_cache - i});
                }
                else if (tn == Lw - 1)
                {
                    // right most, also  connecting to wn=W
                    inei_cache = N - 1 - (Lr + W / 2 - (W - 1)) + 1; // just for index position storage
                    push_neis_back(i, {inei_cache - i, inei_cache - (Lr + W / 2 - (W - 2)) - i, Lw - 1, -1});
                    push_eneis_back(i, {-1, inei_cache - i});
                }
                else
                {
                    // middle
                    push_neis_back(i, {1, Lw, Lw - 1, -1});
                    push_eneis_back(i, {-1, 1});
                    push_bneis_list(i, {Lw, Lw - 1});
                }
            }
            else if (wn == W - 1)
            {
                mesh[i].edge_num = 0;
                edge_lists[0].push_back(i);
                if (tn == 0)
                {
                    // tn=0,w_n=W-1
                    // left most, connecting to w_n=0
                    push_neis_back(i, {1, -Lw, -Lw - 1, Lr + W / 2 - 1 - i});
                    push_eneis_back(i, {Lr + W / 2 - 1 - i, 1});
                    push_bneis_list(i, {-Lw, -Lw - 1});
                }
                else if (tn == Lw - 1)
                {
                    // right most, also  connecting to w_n=0 tn=0
                    push_neis_back(i, {0 - i, -Lw, -Lw - 1, -1});
                    push_eneis_back(i, {-1, 0 - i});
                    push_bneis_list(i, {-Lw, -Lw - 1});
                }
                else
                {
                    // middle
                    push_neis_back(i, {1, -1, -Lw - 1, -Lw});
                    push_eneis_back(i, {-1, 1});
                    push_bneis_list(i, {-Lw - 1, -Lw});
                }
                // edge~
            }
            else
            { // in-bulk~ w_n!=0 and w_n!=W-1
                mesh[i].edge_num = -1;
                if (tn == 0)
                { // left most
                    int wnc = W - 1 - wn;
                    int tnc = Lr + W / 2 - wnc - 1;
                    inei_cache = wnc * (Lr + W / 2) - (wnc - 1) * wnc / 2 + tnc;
                    push_neis_back(i, {1, Lw, inei_cache - tnc - 1 - i, inei_cache - i, -Lw - 1, -Lw});
                    push_bneis_list(i, {1, Lw, inei_cache - tnc - 1 - i, inei_cache - i, -Lw - 1, -Lw});
                }
                else if (tn == Lw - 1)
                { // right most
                    int wnc = W - 1 - wn;
                    inei_cache = wnc * (Lr + W / 2) - (wnc - 1) * wnc / 2;
                    push_neis_back(i, {inei_cache - i,
                                       inei_cache - (Lr + W / 2 - wnc) - 1 - i,
                                       Lw - 1, -1, -Lw - 1, -Lw});
                    push_bneis_list(i, {inei_cache - i,
                                        inei_cache - (Lr + W / 2 - wnc) - 1 - i,
                                        Lw - 1, -1, -Lw - 1, -Lw});
                }
                else
                {
                    //true middle
                    push_neis_back(i, {1, Lw, Lw - 1, -1, -Lw - 1, -Lw});
                    push_bneis_list(i, {1, Lw, Lw - 1, -1, -Lw - 1, -Lw});
                }
            }
        }
    }
}

int dtmc_lc::add_hole_as_edge(int b0, int edgenum)
{
    // add a hole as edge[edge_num], and i0 is one of the edge beads.
    int b1, b2; // another two beads to include for the new triangle edge
    // check if b0 is near the edge, if yes, it can't be on the new edge
    if (if_near_edge(b0))
    {
        //std::cout << b0 << "b0 is near an edge\n";
        return 0;
    }

    int i_nei = 1; // position in mesh[b0].nei
    while (i_nei < mesh[b0].nei.size())
    {
        b1 = mesh[b0].nei[i_nei - 1];
        b2 = mesh[b0].nei[i_nei];
        if (list_a_nei_b(mesh[b1].nei, b2) == -1)
        {
            // b1 b2 somehow not connected?
            std::cout << "b1-b2 not connected, big issue on connection\n";
            continue;
        }
        /* no longer necessary
        if (if_near_edge(b1) || if_near_edge(b2))
        {
            // if b1 or b2 is near another edge, can't use them
            // TODO: this line of code can be optimized
            continue;
        }
        */
        else
        { // found the b1 b2 candidate
            break;
        }
        i_nei += 1;
    }
    // b0-b1-b2-b0 as new edge
    if (i_nei == mesh[b0].nei.size())
    {
        // didn't find appropriate candidate around b0
        return 0;
    }
    // add edge bond, but also keep the b0-b1-b2-b0 order
    mesh[b0].edge_num = edgenum;
    mesh[b0].edge_nei = {b2, b1};
    edge_lists[edgenum].push_back(b0);
    mesh[b1].edge_num = edgenum;
    mesh[b1].edge_nei = {b0, b2};
    edge_lists[edgenum].push_back(b1);
    mesh[b2].edge_num = edgenum;
    mesh[b2].edge_nei = {b1, b0};
    edge_lists[edgenum].push_back(b2);
    // delete these three edge bond from the bulk bulk_bond_list
    delete_bulk_bond_list(b0, b1);
    delete_bulk_bond_list(b1, b2);
    delete_bulk_bond_list(b2, b0);
    return 1;
}

#pragma endregion

#pragma region : useful tools

void dtmc_lc::Ob_init(observable &Ob)
{
    Ob.E = 0;
    Ob.I2H2 = 0;
    Ob.TphiH2 = 0;
    Ob.I2H2dis = 0;
    Ob.IK = 0;
    Ob.Les.resize(Ne);
    for (int e = 0; e < Ne; e++)
    {
        Ob.Les[e] = 0;
    }
    Ob.Tp2uu = 0;
    Ob.Tuuc = 0;
    Ob.Tun2 = 0;
    Ob.IdA = 0;
    Ob.I2H = 0;
    Ob.Bond_num = 0.5 * bulk_bond_list.size();
    for (int n = 0; n < Ne; n++)
    {
        Ob.Bond_num += edge_lists[n].size();
    }
}

int dtmc_lc::if_near_edge(int b)
{
    // check if b is connected to a edge bead
    int bn;
    for (int i = 0; i < mesh[b].nei.size(); i++)
    {
        bn = mesh[b].nei[i];
        if (mesh[bn].edge_num != -1)
        {
            return 1;
        }
    }
    return 0;
}

void dtmc_lc::push_neis_back(int i, std::vector<int> nei_dist)
{
    for (int j = 0; j < nei_dist.size(); j++)
    {
        mesh[i].nei.push_back(i + nei_dist[j]);
    }
}
void dtmc_lc::push_eneis_back(int i, std::vector<int> enei_dist)
{
    for (int j = 0; j < enei_dist.size(); j++)
    {
        mesh[i].edge_nei.push_back(i + enei_dist[j]);
    }
}
void dtmc_lc::push_bneis_list(int i, std::vector<int> bnei_dist)
{
    std::pair<int, int> bond;
    for (int j = 0; j < bnei_dist.size(); j++)
    {
        bond.first = i;
        bond.second = i + bnei_dist[j];
        bulk_bond_list.push_back(bond);
    }
}
void dtmc_lc::delete_bulk_bond_list(int ind_i, int ind_j)
{
    std::pair<int, int> bond0, bond1;
    bond0.first = ind_i;
    bond0.second = ind_j;
    bond1.first = ind_j;
    bond1.second = ind_i;
    for (int i = 0; i < bulk_bond_list.size(); i++)
    {
        if (bulk_bond_list[i] == bond0 || bulk_bond_list[i] == bond1)
        {
            bulk_bond_list.erase(bulk_bond_list.begin() + i);
            i--;
        }
    }
}
#pragma endregion