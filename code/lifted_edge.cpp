#include "dtmc_lc_3d.h"
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>

int dtmc_lc::lifted_edge_metropolis()
{
    // the algorithm
    // 1) find bead of interest, and see if on the longest edge [l]
    // if ([+]&&[l] || [-]&&[s]): res = edge_extend(ind_boi)
    // else res =edge_shrink(ind_boi)
    // if(res==0): rep =- rep
    int ve_ind_boi, ind_boi, num_edge_bead;
    std::vector<int> fedge_list;
    fedge_list.clear();
    for (int i = 0; i < edge_lists.size(); i++)
    {
        fedge_list.insert(fedge_list.end(), edge_lists[i].begin(),
                          edge_lists[i].end());
    } // flatten edge_lists and put in to one int vector
    num_edge_bead = fedge_list.size();
    ve_ind_boi = int(fedge_list.size() * rand_uni(gen));
    ind_boi = fedge_list[ve_ind_boi];

    // found ind_boi on the edge
    bool ledge; // if ind_boi is on the longest edge
    int res;
    int alg = 1;

    if (alg == 0)
    {
        // algorithm 0: replica depends on difference of edge length
        // find if ledge
        ledge = 1;
        for (int i = 0; i < edge_lists.size(); i++)
        {
            if (Ob_sys.Les[mesh[ind_boi].edge_num] > Ob_sys.Les[i])
            {
                ledge = 0;
            }
        }

        if ((ledge && lifted_rep == 1) || (!ledge && lifted_rep == 0))
        {
            res = edge_extend(ind_boi, num_edge_bead);
        }
        else
        {
            res = edge_shrink(ind_boi, num_edge_bead);
        }
        if (res == 0)
        {
            lifted_rep = -lifted_rep;
        }
    }
    else if (alg == 1)
    {
        // algotithm 1: replica depends on edge extension vs shrink
        if (lifted_rep == 1)
        {
            res = edge_extend(ind_boi, num_edge_bead);
        }
        else
        {
            res = edge_shrink(ind_boi, num_edge_bead);
        }
        if (res == 0)
        {
            lifted_rep = -lifted_rep;
        }
    }

    return res;
}

int dtmc_lc::edge_shrink(int ind_boi, int num_edge_bead)
{
    std::vector<int> ind_relate;
    int indi_edge_num_old;
    std::vector<vertex> bead_relate;
    std::vector<std::pair<int, int>> bond_relate_old, bond_relate_new;

    // observables
    observable Ob_relate_old, Ob_relate_new;
    // useful variables
    std::pair<int, int> bond;

#pragma region : [shrink] find update and related beads.
    int e_ind_i, ind_i;
    ind_i = ind_boi;
    e_ind_i = list_a_nei_b(edge_lists[mesh[ind_i].edge_num], ind_i);

    // vive versa for j,k
    int ind_j, ind_k;
    int j_e_nei_i, k_e_nei_i;
    int j_nei_i, k_nei_i;
    int j_nei_k, k_nei_j;
    int j_nei_i_next, j_nei_i_previous;
    int k_nei_i_next, k_nei_i_previous;

    ind_j = mesh[ind_i].edge_nei[0];
    ind_k = mesh[ind_i].edge_nei[1];
    // check number of beads on the edge, need to be greater than 3?
    if (edge_lists[mesh[ind_i].edge_num].size() <= 3)
    {
        // # beads for each edge need to be > 5
        return 0;
    }
    // check j k connection
    if (list_a_nei_b(mesh[ind_j].nei, ind_k) != -1)
    {
        return 0;
    }
    if (distance(ind_j, ind_k) >= l0)
    {
        return 0;
    }
    // check # nei, j and k will be each other's new neighbor, can't be greater than 9
    if (mesh[ind_j].nei.size() >= 9 || mesh[ind_k].nei.size() >= 9)
    {
        return 0;
    }

    // find affected beads
    ind_relate.clear();
    bead_relate.clear();
    bond_relate_old.clear();
    bond_relate_new.clear();
    ind_relate.push_back(ind_i);
    ind_relate.push_back(ind_j);
    ind_relate.push_back(ind_k);
    bead_relate.push_back(mesh[ind_i]);
    bead_relate.push_back(mesh[ind_j]);
    bead_relate.push_back(mesh[ind_k]);
    bond.first = ind_i;
    bond.second = ind_j;
    bond_relate_old.push_back(bond);
    bond_relate_new.push_back(bond);
    bond.first = ind_i;
    bond.second = ind_k;
    bond_relate_old.push_back(bond);
    bond_relate_new.push_back(bond);
    bond.first = ind_j;
    bond.second = ind_k;
    bond_relate_new.push_back(bond);

#pragma endregion

#pragma region : [shrink] store observables of affected beads
    indi_edge_num_old = mesh[ind_i].edge_num;
    Ob_relate_old = Ob_m(ind_relate, bond_relate_old);
#pragma endregion
// neighbors are also sorted during the energy calculation

// [update]
#pragma region : [shrink] edge update
    mesh[ind_i].edge_nei.clear();
    mesh[ind_i].edge_num = -1;
    // no worries, already stored this one

    j_e_nei_i = list_a_nei_b(mesh[ind_j].edge_nei, ind_i);
    if (j_e_nei_i != 1)
    {
        std::cout << "j->i order?\n";
    }
    j_nei_i = list_a_nei_b(mesh[ind_j].nei, ind_i);
    if (list_a_nei_b(mesh[ind_k].nei, ind_j) != -1)
    {
        std::cout << "j-k connected!";
    }
    mesh[ind_j].edge_nei[j_e_nei_i] = ind_k;

    j_nei_i_next = (j_nei_i + 1) % mesh[ind_j].nei.size();
    j_nei_i_previous =
        (j_nei_i - 1 + mesh[ind_j].nei.size()) % mesh[ind_j].nei.size();
    if (mesh[ind_j].nei[j_nei_i_next] ==
        mesh[ind_j].edge_nei[1 - j_e_nei_i])
    {
        // j_nei_i_next is the other edge, can push i forward
        mesh[ind_j].nei.insert(mesh[ind_j].nei.begin() + j_nei_i_next,
                               ind_k);
        j_nei_k = j_nei_i_next;
    }
    else
    {
        mesh[ind_j].nei.insert(mesh[ind_j].nei.begin() + j_nei_i, ind_k);
        j_nei_k = j_nei_i;
    } // okie!

    k_e_nei_i = list_a_nei_b(mesh[ind_k].edge_nei, ind_i);
    if (k_e_nei_i != 0)
    {
        std::cout << "i->k order?\n";
    }
    k_nei_i = list_a_nei_b(mesh[ind_k].nei, ind_i);
    if (list_a_nei_b(mesh[ind_k].nei, ind_j) != -1)
    {
        std::cout << "k-j connected!";
    }
    mesh[ind_k].edge_nei[k_e_nei_i] = ind_j;
    k_nei_i_next = (k_nei_i + 1) % mesh[ind_k].nei.size();
    k_nei_i_previous =
        (k_nei_i - 1 + mesh[ind_k].nei.size()) % mesh[ind_k].nei.size();
    if (mesh[ind_k].nei[k_nei_i_next] ==
        mesh[ind_k].edge_nei[1 - k_e_nei_i])
    {
        // k_nei_i_next in side, can push i forward
        mesh[ind_k].nei.insert(mesh[ind_k].nei.begin() + k_nei_i_next,
                               ind_j);
        k_nei_j = k_nei_i_next;
    }
    else
    {
        // j_nei_i_next is edge, can push i backward
        mesh[ind_k].nei.insert(mesh[ind_k].nei.begin() + k_nei_i, ind_j);
        k_nei_j = k_nei_i;
    }

    if (mesh[ind_i].edge_nei.size() != 0 || mesh[ind_i].edge_num != -1)
    {
        std::cout << "ind_i edge not clear"
                  << "\n";
    }

#pragma endregion

#pragma region : [shrink] get after update observables
    mesh_bead_info_update(ind_relate);
    Ob_relate_new = Ob_m(ind_relate, bond_relate_new);
    // after-update observables
#pragma endregion

#pragma region : [shrink] Metropolis
    // [Metropolis]
    //std::cout<<"[shrink](Ob_relate_new.E - Ob_relate_old.E)"<<(Ob_relate_new.E - Ob_relate_old.E)<<"\n";
    if (rand_uni(gen) <=
        1.0 * num_edge_bead / (num_edge_bead - 1) *
            std::exp(-beta * (Ob_relate_new.E - Ob_relate_old.E)))
    {
        // [accepted]
        // remove ind_i from it's original edge_list
        edge_lists[indi_edge_num_old].erase(
            edge_lists[indi_edge_num_old].begin() + e_ind_i);

        // update system observables
        Ob_sys_update(Ob_relate_new, Ob_relate_old);

        // update bulk_bond_list, add i-j i-k as bulk bond
        std::pair<int, int> bond1, bond2;
        bond1.first = ind_i;
        bond1.second = ind_j;
        bond2.first = ind_j;
        bond2.second = ind_i;
        bulk_bond_list.push_back(bond1);
        bulk_bond_list.push_back(bond2);
        bond1.first = ind_i;
        bond1.second = ind_k;
        bond2.first = ind_k;
        bond2.second = ind_i;
        bulk_bond_list.push_back(bond1);
        bulk_bond_list.push_back(bond2);
        //std::cout<<"edge--\n";
        return 1;
    }
    else
    {
        // [rejected]
        // put everything back (so sad)
        for (int i = 0; i < ind_relate.size(); i++)
        {
            mesh[ind_relate[i]] = bead_relate[i];
        }
        return 0;
    }
#pragma endregion
}

int dtmc_lc::edge_extend(int ind_boi, int num_edge_bead)
{
    std::vector<int> ind_relate;
    int indi_edge_num_old;
    std::vector<vertex> bead_relate;
    std::vector<std::pair<int, int>> bond_relate_old, bond_relate_new;

    // observables
    observable Ob_relate_old, Ob_relate_new;
    // useful variables
    std::pair<int, int> bond;

    // [extend] add i (from the bulk) to edge
    // std::vector<int> i_e_nei, j_e_nei, k_e_nei;
    int ind_i, e_ind_j, ind_j, ind_k;
    int k_nei_j, k_nei_i, j_nei_k;
    int k_e_nei_j, j_e_nei_k;
    int i_nei_j, i_nei_k;

#pragma region : [extend] find update and related beads
    ind_j = ind_boi;
    ind_k = mesh[ind_j].edge_nei[1];
    if (mesh[ind_k].edge_nei[0] != ind_j)
    {
        std::cout << "j-k order wrong here!\n";
    }
    // check # nei, can't be less than 3
    if (mesh[ind_j].nei.size() <= 3 || mesh[ind_k].nei.size() <= 3)
    {
        return 0;
    }

    ind_i = -1;
    int j_next, j_previous;
    int ind_l, l_nei_j;
    for (int l = 0; l < mesh[ind_j].nei.size(); l++)
    {
        ind_l = mesh[ind_j].nei[l];
        if ((list_a_nei_b(mesh[ind_k].nei, ind_l) != -1) &&
            (mesh[ind_l].edge_nei.size() == 0))
        {
            l_nei_j = list_a_nei_b(mesh[ind_l].nei, ind_j);
            j_next = (l_nei_j + 1) % mesh[ind_l].nei.size();
            j_previous = (l_nei_j - 1 + mesh[ind_l].nei.size()) %
                         mesh[ind_l].nei.size();
            if (mesh[ind_l].nei[j_next] == ind_k ||
                mesh[ind_l].nei[j_previous] == ind_k)
            {
                ind_i = ind_l;
                break; // k is the next of j in i_nei
            }
        }
    }

    if (ind_i == -1)
    {
        // std::cout << "couldn't find ind_i!!?\n";
        return 0;
    }
    // check for ind_i, position, need to comply with edge limit in pulling case

    if (lf != 0)
    {
        // pulling experiment for Ne2
        if (mesh[ind_j].edge_num == 0 && mesh[ind_i].R[2] > edge_zlim[0])
        {
            // bottom edge beads has to be z<=0
            return 0;
        }
        if (mesh[ind_j].edge_num == 1 && mesh[ind_i].R[2] < edge_zlim[1])
        {
            // top edge beads has to be z>=lf
            return 0;
        }
    }
    /*
    else if (Epar.g != 0)
    {
        // graviti experiment for edge 0
        if (mesh[ind_j].edge_num == 0 && mesh[ind_i].R[2] > 0)
        {
            // bottom edge beads has to be z<=0
            return 0;
        }
    }
    */
    // check for potential edge bridging bond
    for (int k = 0; k < mesh[ind_i].nei.size(); k++)
    { // check every neighbor of ind_i
        if ((mesh[mesh[ind_i].nei[k]].edge_num != -1) && (mesh[ind_i].nei[k] != ind_j) && (mesh[ind_i].nei[k] != ind_k))
        {
            return 0;
        }
    }

    // order of j k is unclear, unclear now
    // find affected beads
    ind_relate.clear();
    bead_relate.clear();
    bond_relate_old.clear();
    bond_relate_new.clear();
    ind_relate.push_back(ind_i);
    ind_relate.push_back(ind_j);
    ind_relate.push_back(ind_k);
    bead_relate.push_back(mesh[ind_i]);
    bead_relate.push_back(mesh[ind_j]);
    bead_relate.push_back(mesh[ind_k]);
    bond.first = ind_i;
    bond.second = ind_j;
    bond_relate_old.push_back(bond);
    bond_relate_new.push_back(bond);
    bond.first = ind_i;
    bond.second = ind_k;
    bond_relate_old.push_back(bond);
    bond_relate_new.push_back(bond);
    bond.first = ind_j;
    bond.second = ind_k;
    bond_relate_old.push_back(bond);

#pragma endregion

#pragma region : [extend] store observable of affected beads
    Ob_relate_old = Ob_m(ind_relate, bond_relate_old);
    // get pre update info

#pragma endregion

#pragma region : [extend] edge update
    // will remove k_j link and use j-i k_i as new edge
    // need to follow the order of J-K link
    if (mesh[ind_j].edge_nei[0] == ind_k)
    {
        // j<-k
        std::cout << "oh?\n";
        mesh[ind_i].edge_nei.push_back(ind_k);
        mesh[ind_i].edge_nei.push_back(ind_j);
    }
    else if (mesh[ind_j].edge_nei[1] == ind_k)
    {
        // j->k
        mesh[ind_i].edge_nei.push_back(ind_j);
        mesh[ind_i].edge_nei.push_back(ind_k);
    }
    else
    {
        std::cout << "j-k not connected!\n";
    }

    k_nei_j = list_a_nei_b(mesh[ind_k].nei, ind_j);
    k_e_nei_j = list_a_nei_b(mesh[ind_k].edge_nei, ind_j);
    if (k_e_nei_j != 0)
    {
        std::cout << "k_e_nei_j!=0\n";
    }

    if (k_nei_j == -1 || k_e_nei_j == -1)
    {
        std::cout << "k_j not connected here (1)"
                  << "\n";
    }
    // ind_k break up with j
    mesh[ind_k].edge_nei[k_e_nei_j] = ind_i;
    mesh[ind_k].nei.erase(mesh[ind_k].nei.begin() + k_nei_j);

    j_e_nei_k = list_a_nei_b(mesh[ind_j].edge_nei, ind_k);
    j_nei_k = list_a_nei_b(mesh[ind_j].nei, ind_k);
    if (j_nei_k == -1 || j_e_nei_k == -1)
    {
        std::cout << "j_k not connected here (1)"
                  << "\n";
    }

    // ind_j break up with k
    // recall j_e_nei_k = 1s
    mesh[ind_j].edge_nei[j_e_nei_k] = ind_i;
    mesh[ind_j].nei.erase(mesh[ind_j].nei.begin() + j_nei_k);

    mesh[ind_i].edge_num = mesh[ind_j].edge_num;
#pragma endregion

// after-update observables
#pragma region : [extend] get after update observables
    mesh_bead_info_update(ind_relate);
    Ob_relate_new = Ob_m(ind_relate, bond_relate_new);
    // after-update observables

#pragma endregion
// [ Metropolis]
#pragma region : [extend] Metropolis
    //std::cout<<"[extend](Ob_relate_new.E - Ob_relate_old.E)"<<(Ob_relate_new.E - Ob_relate_old.E)<<"\n";
    //std::cout<<"fedge_list.size()"<<fedge_list.size()<<"\n";
    if (rand_uni(gen) <=
        1.0 * num_edge_bead / (num_edge_bead + 1) *
            std::exp(-beta * (Ob_relate_new.E - Ob_relate_old.E)))
    {
        // [accepted]
        if (mesh[ind_i].edge_num == -1)
        {
            std::cout << "somewhat edge_num wrong\n";
        }
        edge_lists[mesh[ind_i].edge_num].push_back(ind_i);
        // I2H_new = 0;
        // i j k are sorted again!
        Ob_sys_update(Ob_relate_new, Ob_relate_old);

        // these are edge bond now
        delete_bulk_bond_list(ind_i, ind_j);
        delete_bulk_bond_list(ind_i, ind_k);
        //std::cout<<"edge++\n";
        return 1;
    }
    else
    {
        // [rejected]
        // rejected put everything back (so sad)
        for (int i = 0; i < ind_relate.size(); i++)
        {
            mesh[ind_relate[i]] = bead_relate[i];
        }
        return 0;
    }
#pragma endregion
}
