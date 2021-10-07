#include "dtmc_lc_3d.h"
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>

int dtmc_lc::bead_metropolis(double delta_s)
{
// change mesh[index].R,
// related energy change: geometric and coupling energy
#pragma region : variable declaration
    // vertex info related, which the update involves
    std::vector<int> ind_relate;
    std::vector<vertex> bead_relate;
    std::vector<std::pair<int, int>> bond_relate;

    // observables
    observable Ob_relate_old, Ob_relate_new;

    // some useful cache variables
    int index, nei_ind;
    std::pair<int, int> bond;
    std::vector<double> delta_pos{0, 0, 0};
    double distance2_cache, uuc_cache; // some cache
    int del_dim;                       // for normal beads it's 3, for fixed z direction it's 2
#pragma endregion

#pragma region : find update, and related beads

    /*do
    {
        index = rand_pos(gen);
    } while (list_a_nei_b(fixed_beads, index) != -1);
    */
    index = rand_pos(gen);
    // TODO: generalize the fixed bead condition, only fix update on the z direction when doing pulling experiment on catenoid, Ne2 membrane

    // checked, it will pass value instead of pointer
    ind_relate.clear();
    bead_relate.clear();
    bond_relate.clear();
    ind_relate.push_back(index);
    bead_relate.push_back(mesh[index]);
    bond.first = index;
    for (int j = 0; j < mesh[index].nei.size(); j++)
    {
        nei_ind = mesh[index].nei[j];

        ind_relate.push_back(nei_ind);
        bead_relate.push_back(mesh[nei_ind]);
        bond.second = nei_ind;
        bond_relate.push_back(bond);
    }

#pragma endregion

#pragma region : store observable of affected beads
    // get pre update info
    Ob_relate_old = Ob_m(ind_relate, bond_relate);

#pragma endregion

#pragma region : bead MC update proposal
    //if (list_a_nei_b(fixed_beads_z, index) != -1)

    for (int k = 0; k < 3; k++)
    {
        delta_pos[k] = 2 * delta_s * rand_uni(gen) - delta_s;
        mesh[index].R[k] += delta_pos[k];
    }

#pragma endregion

#pragma region : hard bead and tether potential between all beads
    // limit on edge beads
    if(lf!=0){
        // pulling for 2 edges
        if(mesh[index].edge_num==0 && mesh[index].R[2]>edge_zlim[0]){
            // edge bead out of range
            mesh[index].R = bead_relate[0].R;
            // return previous position
            return 0;
        }else if(mesh[index].edge_num==1 && mesh[index].R[2]<edge_zlim[1]){
            // edge bead out of range
            mesh[index].R = bead_relate[0].R;
            // return previous position
            return 0;
        }
    }
    else if(Epar.g!=0){
        // pull edge 0 when there is gravity
        if(mesh[index].edge_num==0 && mesh[index].R[2]>0){
            // edge bead out of range
            mesh[index].R = bead_relate[0].R;
            // return previous position
            return 0;
        }
    }

    // hard bead potential
    for (int nei_ind = 0; nei_ind < mesh.size(); nei_ind++)
    {
        if ((std::abs(mesh[index].R[0] - mesh[nei_ind].R[0])) < 1 &&
            (std::abs(mesh[index].R[1] - mesh[nei_ind].R[1])) < 1 &&
            (std::abs(mesh[index].R[2] - mesh[nei_ind].R[2])) < 1)
        {
            distance2_cache = distance2(index, nei_ind);
            if ((distance2_cache - 1 < 0) && (nei_ind != index))
            {
                mesh[index].R = bead_relate[0].R;
                // return to old position
                return 0;
            }
        }
    }
    // tether potential between linked beads
    for (int j = 0; j < mesh[index].nei.size(); j++)
    {
        nei_ind = mesh[index].nei[j];
        distance2_cache = distance2(index, nei_ind);
        /*
        if (distance2_cache - l0 * l0 >= 0)
        {
            mesh[index].R = bead_relate[0].R;
            // return previous position
            return 0;
        }
        */
        // used when there are two tether length limit
        /*
        if (mesh[index].edge_nei.size() != 0 &&
            (nei_ind == mesh[index].edge_nei[0] ||
             nei_ind == mesh[index].edge_nei[1]))
        {
            if (distance2_cache - l1 * l1 >= 0)
            {
                mesh[index].R = bead_relate[0].R;
                // return previous position
                return 0;
            }
        }
        else
        */
        if (distance2_cache - l0 * l0 >= 0)
        {
            mesh[index].R = bead_relate[0].R;
            // return previous position
            return 0;
        }
    }
#pragma endregion

#pragma region : get after - update observables

    mesh_bead_info_update(ind_relate);
    Ob_relate_new = Ob_m(ind_relate, bond_relate);

#pragma endregion

#pragma region : Metropolis
    // std::cout << "Er_new-Er_old=" << Er_new - Er_old << "\n";
    if (rand_uni(gen) <=
        std::exp(-beta * (Ob_relate_new.E - Ob_relate_old.E)))
    {
        // [accepted]
        Ob_sys_update(Ob_relate_new, Ob_relate_old);
        return 1;
    }
    else
    {
        // [rejected]
        for (int i = 0; i < ind_relate.size(); i++)
        {
            mesh[ind_relate[i]] = bead_relate[i];
        }
        return 0;
    }
#pragma endregion
}

int dtmc_lc::spin_metropolis(double delta_theta)
{
    // tracking related beads and bond
    std::vector<int> ind_relate;
    std::vector<vertex> bead_relate;
    std::vector<std::pair<int, int>> bond_relate;

    // observables
    observable Ob_relate_old, Ob_relate_new;
    // some useful variables
    int index, nei_ind; // bead to rotate spin, and cache for it's nei
    std::pair<int, int> bond;
    std::vector<double> ind_u_old;
    int axr0, axr1; // rotate in the axr0-axr1 plane
    double theta2rot, cos_theta,
        sin_theta; // angle of rotation, uniform distribution in
                   // [-delta_theta,delta_theta], store cos and sin for
                   // efficiency
// spin_update
#pragma region : find update, and related beads
    index = rand_pos(gen);
#pragma endregion

#pragma region : store related observables of affected beads
    ind_relate.clear();
    bead_relate.clear();
    bond_relate.clear();
    ind_relate.push_back(index);
    bead_relate.push_back(mesh[index]);
    bond.first = index;
    for (int j = 0; j < mesh[index].nei.size(); j++)
    {
        bond.second = mesh[index].nei[j];
        bond_relate.push_back(bond);
    }
    Ob_relate_old = Ob_m(ind_relate, bond_relate);

#pragma endregion

#pragma region : do spin rotation MC update proposal
    // find rotational plan, randomly
    axr0 = int(3 * rand_uni(gen));
    axr1 = (axr0 + 1) % 3;
    // rotate random angle.
    theta2rot = delta_theta * (2 * rand_uni(gen) - 1);
    cos_theta = std::cos(theta2rot);
    sin_theta = std::sin(theta2rot);
    // carry out rotation
    ind_u_old = mesh[index].u;
    mesh[index].u[axr0] =
        ind_u_old[axr0] * cos_theta - ind_u_old[axr1] * sin_theta;
    mesh[index].u[axr1] =
        ind_u_old[axr0] * sin_theta + ind_u_old[axr1] * cos_theta;
#pragma endregion

#pragma region : get after update observables
    mesh[index].un2 = un2_m(index);
    Ob_relate_new = Ob_m(ind_relate, bond_relate);

#pragma endregion

#pragma region : Metropolis
    if (rand_uni(gen) <=
        std::exp(-beta * (Ob_relate_new.E - Ob_relate_old.E)))
    {
        // [accepted]
        Ob_sys_update(Ob_relate_new, Ob_relate_old);
        return 1;
    }
    else
    {
        //[rejected]
        mesh[index] = bead_relate[0];
        return 0;
    }
#pragma endregion
}

int dtmc_lc::bond_metropolis()
/*
           - (o) -
         -    |       -
       -      |       - (b)
    (a)       |     -
        -     |   -
           - (i)

    to:
              - (o) -
         -            -
       -      - - - - - (b)
    (a) - -         -
        -         -
           - (i)
    will change local n, thus both geometric and spin-normal energy are affected
*/
{
    // vertex info related, which the update involves
    std::vector<int> ind_relate;
    std::vector<vertex> bead_relate;
    std::vector<std::pair<int, int>> bond_relate_old, bond_relate_new;

    // observables
    observable Ob_relate_old, Ob_relate_new;
    // useful variables
    std::pair<int, int> bond;
    int bondlist_ind_i;
    int ind_i, ind_a, ind_b, ind_o;
    int i_nei_o, i_nei_a, i_nei_b; // mesh[ind_i].nei[i_nei_o] = ind_o;
    int o_nei_i, o_nei_a, o_nei_b;
    int a_nei_i, a_nei_o, a_nei_b;
    int b_nei_i, b_nei_o, b_nei_a;

#pragma region : find update and related beads
    if (bulk_bond_list.size() == 0)
    {
        return 0;
    }
    bondlist_ind_i = int(bulk_bond_list.size() * rand_uni(gen));
    ind_i = bulk_bond_list[bondlist_ind_i].first;
    ind_o = bulk_bond_list[bondlist_ind_i].second;

    i_nei_o = list_a_nei_b(mesh[ind_i].nei, ind_o);
    i_nei_a = (i_nei_o + 1) % mesh[ind_i].nei.size();
    i_nei_b = (i_nei_o - 1 + mesh[ind_i].nei.size()) % mesh[ind_i].nei.size();
    ind_a = mesh[ind_i].nei[i_nei_a];
    ind_b = mesh[ind_i].nei[i_nei_b];
    if (list_a_nei_b(mesh[ind_a].nei, ind_b) != -1)
    {
        // a and b can't be connected
        return 0;
    }

    // check the a-b distance, a-b is the new bond
    if (distance(ind_a, ind_b) > l0)
    {
        return 0;
    }
    // check edge_num, can't connect different edges
    // it has became unneseccery since no difference were observed
    // conclusion on mobius strip, this condition is not useful

    // still, helps perserving the topoloty
    // check for in-bulk bond bridging edges
    if ((mesh[ind_a].edge_num != -1) && (mesh[ind_b].edge_num != -1))
    {
        return 0;
    }

    // check # of nei
    // can't be greater than 9, otherwise will have more than 9 after flip
    if (mesh[ind_a].nei.size() >= 9 || mesh[ind_b].nei.size() >= 9)
    {
        return 0;
    }
    // can't be less than 3
    // both has to have at least 4 neighbours, after flip, will have at least
    // 3;
    if (mesh[ind_i].nei.size() <= 3 || mesh[ind_o].nei.size() <= 3)
    {
        return 0;
    }
    // find the rest of bonds
    o_nei_i = list_a_nei_b(mesh[ind_o].nei, ind_i);
    o_nei_a = list_a_nei_b(mesh[ind_o].nei, ind_a);
    o_nei_b = list_a_nei_b(mesh[ind_o].nei, ind_b);
    a_nei_i = list_a_nei_b(mesh[ind_a].nei, ind_i);
    a_nei_o = list_a_nei_b(mesh[ind_a].nei, ind_o);
    b_nei_i = list_a_nei_b(mesh[ind_b].nei, ind_i);
    b_nei_o = list_a_nei_b(mesh[ind_b].nei, ind_o);

    // find related beads
    ind_relate.clear();
    bead_relate.clear();
    bond_relate_old.clear();
    bond_relate_new.clear();
    ind_relate.push_back(ind_i);
    ind_relate.push_back(ind_o);
    ind_relate.push_back(ind_a);
    ind_relate.push_back(ind_b);
    bead_relate.push_back(mesh[ind_i]);
    bead_relate.push_back(mesh[ind_o]);
    bead_relate.push_back(mesh[ind_a]);
    bead_relate.push_back(mesh[ind_b]);
    bond.first = ind_i;
    bond.second = ind_o;
    bond_relate_old.push_back(bond);
    bond.first = ind_a;
    bond.second = ind_b;
    bond_relate_new.push_back(bond);

#pragma endregion

#pragma region : store observable of affected beads
    Ob_relate_old = Ob_m(ind_relate, bond_relate_old);
#pragma endregion

#pragma region : flip bond update
    // [update]
    // remove i-o bond
    mesh[ind_i].nei.erase(mesh[ind_i].nei.begin() + i_nei_o);
    mesh[ind_o].nei.erase(mesh[ind_o].nei.begin() + o_nei_i);
    // one of i_nei_a, i_nei_b changes
    // one of o_nei_a, o_nei_b changes

    // add a-b bond
    if ((a_nei_i - a_nei_o + mesh[ind_a].nei.size()) % mesh[ind_a].nei.size() ==
        1)
    {
        mesh[ind_a].nei.insert(mesh[ind_a].nei.begin() + a_nei_i, ind_b);
        a_nei_b = a_nei_i;
        a_nei_i++;
    }
    else if ((a_nei_o - a_nei_i + mesh[ind_a].nei.size()) %
                 mesh[ind_a].nei.size() ==
             1)
    {
        mesh[ind_a].nei.insert(mesh[ind_a].nei.begin() + a_nei_o, ind_b);
        a_nei_b = a_nei_o;
        a_nei_o++;
    }
    else
    {
        std::cout << "i-o not next in a_nei!"
                  << "a_nei_i,a_nei_o=" << a_nei_i << "," << a_nei_o
                  << " a.nei.size()=" << mesh[ind_a].nei.size() << "\n";
        std::cout << "miao?";
    } // error code 2
    if ((b_nei_i - b_nei_o + mesh[ind_b].nei.size()) % mesh[ind_b].nei.size() ==
        1)
    {
        mesh[ind_b].nei.insert(mesh[ind_b].nei.begin() + b_nei_i, ind_a);
        b_nei_a = b_nei_i;
        b_nei_i++;
    }
    else if ((b_nei_o - b_nei_i + mesh[ind_b].nei.size()) %
                 mesh[ind_b].nei.size() ==
             1)
    {
        mesh[ind_b].nei.insert(mesh[ind_b].nei.begin() + b_nei_o, ind_a);
        b_nei_a = b_nei_o;
        b_nei_o++;
    }
    else
    {
        std::cout << "i-o not next in b_nei!"
                  << "b_nei_i,b_nei_o=" << b_nei_i << "," << b_nei_o
                  << " b.nei.size()=" << mesh[ind_b].nei.size() << "\n";
        std::cout << "huh?";
    } // error code 2
#pragma endregion

#pragma region : get after - update observables
    mesh_bead_info_update(ind_relate);
    Ob_relate_new = Ob_m(ind_relate, bond_relate_new);

#pragma endregion

#pragma region : Metropolis
    if (rand_uni(gen) <=
        std::exp(-beta * (Ob_relate_new.E - Ob_relate_old.E)))
    {
        // [accepted]
        Ob_sys_update(Ob_relate_new, Ob_relate_old);
        // update bulk_bond_list
        delete_bulk_bond_list(ind_i, ind_o);
        std::pair<int, int> bond0, bond1;
        bond0.first = ind_a;
        bond0.second = ind_b;
        bond1.first = ind_b;
        bond1.second = ind_a;
        bulk_bond_list.push_back(bond0);
        bulk_bond_list.push_back(bond1);

        return 1;
    }
    else
    {
        // [rejected]
        // put everything back (so sad T.T)
        for (int i = 0; i < ind_relate.size(); i++)
        {
            mesh[ind_relate[i]] = bead_relate[i];
        }
        return 0;
    }
#pragma endregion
}

int dtmc_lc::edge_metropolis()
{
    std::vector<int> ind_relate;
    int indi_edge_num_old;
    std::vector<vertex> bead_relate;
    std::vector<std::pair<int, int>> bond_relate_old, bond_relate_new;

    // observables
    observable Ob_relate_old, Ob_relate_new;
    // useful variables
    std::pair<int, int> bond;
    std::vector<int> fedge_list;
    int e2update;
    e2update = int(Ne * rand_uni(gen));
    fedge_list = edge_lists[e2update]; // choose a edge to update with equall probability
    //std::cout<<"e2update="<<e2update<<"\n";
    /*
    fedge_list.clear();
    for (int i = 0; i < edge_lists.size(); i++)
    {
        fedge_list.insert(fedge_list.end(), edge_lists[i].begin(),
                          edge_lists[i].end());
    } // flatten edge_lists and put in to one int vector
    */
    if (rand_uni(gen) < 0.5)
    {
        // [shrink] (add a bond:2 bulk bond-1 edge bond)
        // remove ind_i
        /*
                        (j)
        (k)             =
            =         =
            = (i)

        to:
                = =   (j)
        (k)  =        -
            -       -
            - (i)
        */

#pragma region : [shrink] find update and related beads.
        int ve_ind_i, e_ind_i, ind_i;
        ve_ind_i = int(fedge_list.size() * rand_uni(gen));
        ind_i = fedge_list[ve_ind_i];
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
        // check j k connection
        if (list_a_nei_b(mesh[ind_j].nei, ind_k) != -1)
        {
            return 0;
        }
        // check distance, note the bond limit differes for on-edge and in-bulk bond
        /*
        if (distance(ind_j, ind_k) >= l1 || distance(ind_i, ind_j) >= l0 ||
            distance(ind_i, ind_k) >= l0)
        {
            return 0;
        }
        */
        // check distance, no l1 case
        /*
        if (distance(ind_j, ind_k) >= l0)
        {
            return 0;
        }
        */
        if (distance(ind_j, ind_k) >= l0)
        {
            return 0;
        }
        /*if (distance(ind_j, ind_k) >= l1 || distance(ind_i, ind_j) >= l0 ||
            distance(ind_i, ind_k) >= l0) {
            return 0;
        }*/
        // check # nei, can't be greater than 9
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
            1.0 * fedge_list.size() / (fedge_list.size() - 1) *
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
    } // end remove ind_i

    else
    {
        // [extend] (remove a bond: + 1 edge bond - 2 bulk bond)
        // add ind_i
        /*
                = =   (j)
        (k)  =        -
            -       -
            - (i)
        to:
                        (j)
        (k)             =
            =         =
            = (i)
        */
        // [extend] add i (from the bulk) to edge
        std::vector<int> i_nei, j_nei, k_nei, l_nei;
        // std::vector<int> i_e_nei, j_e_nei, k_e_nei;
        int ind_i, e_ind_j, ind_j, ind_k;
        int k_nei_j, k_nei_i, j_nei_k;
        int k_e_nei_j, j_e_nei_k;
        int v2a_ind_i;
        int i_nei_j, i_nei_k;
#pragma region : [extend] find update and related beads
        e_ind_j = int(fedge_list.size() * rand_uni(gen));
        // ind_i = v2add_list[v2a_ind_i];
        ind_j = fedge_list[e_ind_j];
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
        // let's find i
        j_nei = mesh[ind_j].nei;
        k_nei = mesh[ind_k].nei;
        ind_i = -1;
        int j_next, j_previous;
        int ind_l, l_nei_j;
        for (int l = 0; l < mesh[ind_j].nei.size(); l++)
        {
            ind_l = mesh[ind_j].nei[l];
            l_nei = mesh[ind_l].nei;
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

        if(lf!=0){
            // pulling experiment for Ne2
            if(mesh[ind_j].edge_num==0 && mesh[ind_i].R[2]>edge_zlim[0]){
                // bottom edge beads has to be z<=0
                return 0;
            }
            if(mesh[ind_j].edge_num==1 && mesh[ind_i].R[2]<edge_zlim[1]){
                // top edge beads has to be z>=lf
                return 0;
            }
        }
        else if (Epar.g!=0){
            // graviti experiment for edge 0
            if(mesh[ind_j].edge_num==0 && mesh[ind_i].R[2]>0){
                // bottom edge beads has to be z<=0
                return 0;
            }
        }

        // for Ne>1 check if i has neighbour on the other edge
        // same here, this part is unneseccery since not difference were found
        // [Aug7 2021] but, since there is some folding when using simulated tempering start, just add it back
        // check for potential edge bridging bond
        for (int k = 0; k < mesh[ind_i].nei.size(); k++)
        { // check every neighbor of ind_i
            if ((mesh[mesh[ind_i].nei[k]].edge_num != -1) && (mesh[ind_i].nei[k] != ind_j) && (mesh[ind_i].nei[k] != ind_k))
            {
                return 0;
            }
        }

        i_nei = mesh[ind_i].nei;

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
            1.0 * fedge_list.size() / (fedge_list.size() + 1) *
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
    } // end add ind_i
}

/*
int dtmc_lc::swap_metropolis()
{
    // since what's been changed is very simlple. I won't go through the entire process of storaging every related variables, instead I'll just update Ob_sys two un2 variable and the energy E.
    int ind_i, ind_j; // swap i,j cn status
    double Er_old, Er_new;
    double Tun2r_old, Tun2r_new;
    double Tun2pr_old, Tun2pr_new;
#pragma region : find related beads and obserable
    ind_i = rand_pos(gen);
    ind_j = (ind_i + 1 + int((N - 1) * rand_uni(gen))) % N;

    Er_old = 0;
    Tun2r_old = 0;
    Tun2pr_old = 0;
    std::vector<int> inds{ind_i, ind_j};
    for (int k = 0; k < 2; k++)
    {
        if (mesh[inds[k]].is_cnp)
        {
            Tun2pr_old += mesh[inds[k]].un2;
        }
        else
        {
            Tun2r_old += mesh[inds[k]].un2;
        }
    }
    Er_old += Epar.Cn * Tun2r_old + Epar.Cnp * Tun2pr_old;

#pragma endregion
#pragma region : do swap update
    //swap is_cnp parameter
    int is_cnp_cache = mesh[ind_i].is_cnp;
    mesh[ind_i].is_cnp = mesh[ind_j].is_cnp;
    mesh[ind_j].is_cnp = is_cnp_cache;

#pragma endregion
#pragma region : get after update observables
    Er_new = 0;
    Tun2r_new = 0;
    Tun2pr_new = 0;
    for (int k = 0; k < 2; k++)
    {
        if (mesh[inds[k]].is_cnp)
        {
            Tun2pr_new += mesh[inds[k]].un2;
        }
        else
        {
            Tun2r_new += mesh[inds[k]].un2;
        }
    }
    Er_new += Epar.Cn * Tun2r_new + Epar.Cnp * Tun2pr_new;

#pragma endregion
#pragma region : Metropolis
    if (rand_uni(gen) < std::exp(-beta * (Er_new - Er_old)))
    {
        // pass the update
        Ob_sys.Tun2 += Tun2r_new - Tun2r_old;
        Ob_sys.Tun2p += Tun2pr_new - Tun2pr_old;
        Ob_sys.E += Er_new - Er_old;
        return 1;
    }
    else
    {
        // reject, swap back
        is_cnp_cache = mesh[ind_i].is_cnp;
        mesh[ind_i].is_cnp = mesh[ind_j].is_cnp;
        mesh[ind_j].is_cnp = is_cnp_cache;
        return 0;
    }

#pragma endregion
}
*/

/*
int dtmc_lc::hop_metropolis()
{
    // hop update changes the relative density of short rods on the outside vs. inside
    // the density related parameter phi continuously take (-1,1)
    // notice this update is local (for now, without interaction Aug13_2021)
#pragma region : variable declaration
    // vertex info related,
    int index, ind_j;
    double phi_old, phi_new;
    double local_IKphi2_old, local_IKphi2_new;
    double local_I2H2dis_old, local_I2H2dis_new; // only affected variables
    double local_Tphi_nei;
    double Er_old, Er_new;
#pragma endregion

#pragma region : find update, and related beads
    index = rand_pos(gen);
#pragma endregion

#pragma region : store observable of affected beads
    phi_old = mesh[index].phi;
    local_I2H2dis_old = mesh[index].dAn2H[0] * std::pow(mesh[index].dAn2H[1] - mesh[index].phi * Epar.C0, 2);
    local_IKphi2_old = mesh[index].dAK * mesh[index].phi * mesh[index].phi;
    local_Tphi_nei = 0;
    for (int j = 0; j < mesh[index].nei.size(); j++)
    {
        ind_j = mesh[index].nei[j];
        local_Tphi_nei += mesh[ind_j].phi;
    }
    Er_old = 0.5 * Epar.kar * local_I2H2dis_old + Epar.karg * local_IKphi2_old - Epar.J * mesh[index].phi * local_Tphi_nei;

#pragma endregion

#pragma region : hopping MC update proposal
    mesh[index].phi = 2 * rand_uni(gen) - 1; //randomly set within (-1,1)
#pragma endregion

#pragma region : get after - update observables
    local_I2H2dis_new = mesh[index].dAn2H[0] * std::pow(mesh[index].dAn2H[1] - mesh[index].phi * Epar.C0, 2);
    local_IKphi2_new = mesh[index].dAK * mesh[index].phi * mesh[index].phi;
    Er_new = 0.5 * Epar.kar * local_I2H2dis_new + Epar.karg * local_IKphi2_new - Epar.J * mesh[index].phi * local_Tphi_nei;
#pragma endregion

#pragma region : Metropolis
    if (rand_uni(gen) <=
        std::exp(-beta * (Er_new - Er_old)))
    {
        // [accepted]
        Ob_sys.I2H2dis += local_I2H2dis_new - local_I2H2dis_old;
        Ob_sys.IKphi2 += local_IKphi2_new - local_IKphi2_old;
        Ob_sys.E += Er_new - Er_old;
        Ob_sys.Tphi2 += (mesh[index].phi - phi_old) * local_Tphi_nei;
        Ob_sys.Iphi += mesh[index].phi - phi_old;

        return 1;
    }
    else
    {
        // [rejected]
        mesh[index].phi = phi_old;
        return 0;
    }
#pragma endregion
}
*/