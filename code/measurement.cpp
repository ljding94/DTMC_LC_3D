#include "dtmc_lc_3d.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#define PI 3.14159265358979323846

#pragma region : tools
double dtmc_lc::innerproduct(std::vector<double> a, std::vector<double> b)
{
    double result = 0;
    for (int k = 0; k < a.size(); k++)
    {
        result += a[k] * b[k];
    }
    return result;
}
std::vector<double> crossproduct(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> cp{0, 0, 0};
    cp[0] = a[1] * b[2] - a[2] * b[1];
    cp[1] = a[2] * b[0] - a[0] * b[2];
    cp[2] = a[0] * b[1] - a[1] * b[0];
    return cp;
}
double dtmc_lc::distance2(int ind_1, int ind_2)
{
    double result = 0;
    for (int k = 0; k < 3; k++)
    {
        result += std::pow(mesh[ind_2].R[k] - mesh[ind_1].R[k], 2);
    }
    return result;
}
double dtmc_lc::distance(int ind_1, int ind_2)
{
    double result = distance2(ind_1, ind_2);
    result = std::sqrt(result);
    return result;
}
double dtmc_lc::distancefp(int ind_1, std::vector<double> p)
{
    double result = 0;
    for (int k = 0; k < 3; k++)
    {
        result += std::pow(p[k] - mesh[ind_1].R[k], 2);
    }
    result = std::sqrt(result);
    return result;
}

int dtmc_lc::sort_nei(int index)
{
    // nei list can be sorted as long as two edge_nei are next to each other in the nei list
    std::vector<int> index_nei = mesh[index].nei;
    std::vector<int> index_enei = mesh[index].edge_nei;
    if (mesh[index].edge_nei.size() == 0)
    {
        return 0;
    }

    //alternative way involves no while()
    int count = 0;
    for (int count = 0; count < mesh[index].nei.size(); count++)
    {
        mesh[index].nei.push_back(mesh[index].nei[0]);
        mesh[index].nei.erase(mesh[index].nei.begin());
        if (mesh[index].nei[0] == mesh[index].edge_nei[0])
        {
            break;
        }
    }
    /*
    while (mesh[index].nei[0] != mesh[index].edge_nei[0])
    {
        mesh[index].nei.push_back(mesh[index].nei[0]);
        mesh[index].nei.erase(mesh[index].nei.begin());
    }
    */
    if (mesh[index].nei[1] == mesh[index].edge_nei[1])
    {
        mesh[index].nei.push_back(mesh[index].nei[0]);
        mesh[index].nei.erase(mesh[index].nei.begin());
    }
    else if (mesh[index].nei.back() != mesh[index].edge_nei[1])
    {
        std::cout << "wrong, can't sort neighbors!" << index << "em(" << mesh[index].edge_num << " \n ";
        std::cout << "nei size" << mesh[index].nei.size() << " are\n";
        for (int k = 0; k < mesh[index].nei.size(); k++)
        {
            std::cout << "," << mesh[index].nei[k];
            if (mesh[mesh[index].nei[k]].edge_num != -1)
            {
                std::cout << "(e)";
            }
        }
        return 0;
    }
    return 1;
}

int dtmc_lc::list_a_nei_b(std::vector<int> a, int b)
{
    for (int j = 0; j < a.size(); j++)
    {
        if (a[j] == b)
        {
            return j;
        }
    }
    return -1;
}

int dtmc_lc::check_nei_connect()
{
    int j_next;
    for (int i = 0; i < mesh.size(); i++)
    {
        if (mesh[i].edge_nei.size() != 0 && mesh[i].edge_nei.size() != 2)
        {
            std::cout << "mesh[" << i << "].edge_nei.size()=" << mesh[i].edge_nei.size();
        }
    }
    return 0;
}
int dtmc_lc::check_duplication(int ind_i)
{
    for (int j = 0; j < mesh[ind_i].nei.size() - 1; j++)
    {
        for (int k = j + 1; k < mesh[ind_i].nei.size(); k++)
        {
            if (mesh[ind_i].nei[j] == mesh[ind_i].nei[k])
            {
                std::cout << "duplicated nei!\n";
                return 1;
            }
        }
    }
    return 0;
}

#pragma endregion

#pragma region : energy related

void dtmc_lc::mesh_bead_info_update(std::vector<int> ind_relate)
{
    int ind;
    for (int i = 0; i < ind_relate.size(); i++)
    {
        ind = ind_relate[i];
        mesh[ind].n = n_m(ind);
        mesh[ind].dAn2H = dAn2H_m(ind); // signed H depends on n
        mesh[ind].ds = ds_m(ind);
        //mesh[ind].dsk2 = dsk2_m(ind);
        mesh[ind].dAK = dAK_m(ind);
        mesh[ind].un2 = un2_m(ind);
    }
}

observable dtmc_lc::Ob_m(std::vector<int> ind_relate,
                         std::vector<std::pair<int, int>> bond_relate)
{
    //Observable related measurement for local metropolis update usage
    // note: hop_metropolis not using this function [Aug17]
    int ind, ind_i, ind_j;
    observable Ob;
    Ob_init(Ob); // set inital observable value, all to 0

    for (int i = 0; i < ind_relate.size(); i++)
    {
        ind = ind_relate[i];
        // geometric terms
        Ob.I2H2 += mesh[ind].dAn2H[0] * mesh[ind].dAn2H[1] * mesh[ind].dAn2H[1];
        //Ob.I2H2dis += mesh[ind].dAn2H[0] * std::pow(mesh[ind].dAn2H[1] - mesh[ind].phi * Epar.C0, 2);
        Ob.IK += mesh[ind].dAK;
        //Ob.IKphi2 += mesh[ind].dAK * mesh[ind].phi * mesh[ind].phi;
        if (mesh[ind].edge_num != -1)
        {
            Ob.Les[mesh[ind].edge_num] += mesh[ind].ds;
            //Ob.Ik2s[mesh[ind].edge_num] += mesh[ind].dsk2;
            //Ob.Leuns[mesh[ind].edge_num] += mesh[ind].ds * std::sqrt(mesh[ind].un2);
            //Ob.Leuns[mesh[ind].edge_num] += mesh[ind].ds * ut_m(ind);
        }
        // coupling terms
        Ob.Tun2 += mesh[ind].un2;
        //Ob.IKun2 += mesh[ind].dAK * mesh[ind].un2;
        // miscellany terms
        Ob.IdA += mesh[ind].dAn2H[0];
        Ob.I2H += mesh[ind].dAn2H[0] * mesh[ind].dAn2H[1];
    }
    // crystalline terms
    for (int i = 0; i < bond_relate.size(); i++)
    {
        ind_i = bond_relate[i].first;
        ind_j = bond_relate[i].second;
        Ob.Tp2uu += p2uu_m(ind_i, ind_j);
        Ob.Tuuc += uuc_m(ind_i, ind_j);
        // also, Ising-like phi field tems
        //Ob.Tphi2 += mesh[ind_i].phi * mesh[ind_j].phi;
        // use new LC energy that seperate twist from the splay and bend
        //Ob.Tuusb += uusb_m(ind_i, ind_j);
        //Ob.Tuut += uut_m(ind_i, ind_j);
    }
    Ob.Bond_num += bond_relate.size();
    Ob.E = E_m(Ob);
    return Ob;
}
void dtmc_lc::Ob_sys_update(observable Ob_new, observable Ob_old)
{
    Ob_sys.E += Ob_new.E - Ob_old.E;
    Ob_sys.I2H2 += Ob_new.I2H2 - Ob_old.I2H2;
    //Ob_sys.Iphi += Ob_new.Iphi - Ob_old.Iphi;
    //Ob_sys.Tphi2 += Ob_new.Tphi2 - Ob_old.Tphi2;
    //Ob_sys.I2H2dis += Ob_new.I2H2dis - Ob_old.I2H2dis;
    Ob_sys.IK += Ob_new.IK - Ob_old.IK;
    //Ob_sys.IKphi2 += Ob_new.IKphi2 - Ob_old.IKphi2;
    for (int e = 0; e < Ne; e++)
    {
        Ob_sys.Les[e] += Ob_new.Les[e] - Ob_old.Les[e];
        //Ob_sys.Ik2s[e] += Ob_new.Ik2s[e] - Ob_old.Ik2s[e];
        //Ob_sys.Leuns[e] += Ob_new.Leuns[e] - Ob_old.Leuns[e];
    }
    Ob_sys.Tp2uu += Ob_new.Tp2uu - Ob_old.Tp2uu;
    Ob_sys.Tuuc += Ob_new.Tuuc - Ob_old.Tuuc;
    //Ob_sys.Tuusb += Ob_new.Tuusb - Ob_old.Tuusb;
    //Ob_sys.Tuut += Ob_new.Tuut - Ob_old.Tuut;
    Ob_sys.Tun2 += Ob_new.Tun2 - Ob_old.Tun2;
    //Ob_sys.IKun2 += Ob_new.IKun2 - Ob_old.IKun2;
    Ob_sys.IdA += Ob_new.IdA - Ob_old.IdA;
    Ob_sys.I2H += Ob_new.I2H - Ob_old.I2H;

    Ob_sys.Bond_num += Ob_new.Bond_num - Ob_old.Bond_num;
}
double dtmc_lc::E_m(observable Ob)
{
    // Energy measurement for local Ob and global Ob_sys usage
    double E = 0;
    E += 0.5 * Epar.kar * Ob.I2H2;
    //E += -Epar.J * Ob.Tphi2;
    E += Epar.karg * Ob.IK;
    //E += Epar.karg * Ob.IKphi2;
    for (int e = 0; e < Ne; e++)
    {
        E += Epar.lam * Ob.Les[e];
        //E += 0.5 * Epar.B * Ob.Ik2s[e];
        //E += Epar.lamd * Ob.Leuns[e];
    }
    E += -Epar.Kd * (Ob.Tp2uu + Epar.q * Ob.Tuuc);
    //E += Epar.Ksb * Ob.Tuusb + Epar.Kt * Ob.Tuut; // use different moduli for twist and the other two
    E += -0.5 * Epar.Cn * Ob.Tun2;
    //E += Epar.kard * Ob.IKun2;

    //E += -0.5 * (Epar.Cn * Ob.Tun2 + Epar.Cnp * Ob.Tun2p);
    // + kard * Ob.IKun2;00
    return E;
}

double dtmc_lc::ds_m(int index)
{
    if (mesh[index].edge_num == -1)
    {
        return 0;
    }
    double Le_l = 0;
    int nei_ind;
    double distance;
    for (int i = 0; i < mesh[index].edge_nei.size(); i++)
    {
        nei_ind = mesh[index].edge_nei[i];
        distance = 0;
        for (int j = 0; j < 3; j++)
        {
            // Rij[j] = mesh[index].R[j] - mesh[nei_ind].R[j];
            distance += std::pow(mesh[index].R[j] - mesh[nei_ind].R[j], 2);
        }
        distance = std::sqrt(distance);
        Le_l += 0.5 * distance;
    }
    if (Le_l > 2)
    {
        std::cout << "Le_l=" << Le_l << "\n";
    }
    return Le_l;
}

/*
double dtmc_lc::dsk2_m(int index)
{
    if (mesh[index].edge_num == -1)
    {
        return 0;
    }

    //(j)===(i)===(k)

    double cos_jk, theta_jk;
    double l2_ji, l2_ik;
    int ind_j, ind_k;
    std::vector<double> r_ji{0, 0, 0};
    std::vector<double> r_ik{0, 0, 0};
    cos_jk = 0;
    l2_ji = 0;
    l2_ik = 0;
    ind_j = mesh[index].edge_nei[0];
    ind_k = mesh[index].edge_nei[1];
    for (int k = 0; k < 3; k++)
    {
        r_ji[k] = mesh[index].R[k] - mesh[ind_j].R[k];
        l2_ji += r_ji[k] * r_ji[k];
        r_ik[k] = mesh[ind_k].R[k] - mesh[index].R[k];
        l2_ik += r_ik[k] * r_ik[k];
        cos_jk += r_ji[k] * r_ik[k];
    }
    cos_jk = cos_jk / (std::sqrt(l2_ji * l2_ik));
    // make sure cos is not over the bound
    cos_jk = std::min(cos_jk, 1.0);
    cos_jk = std::max(cos_jk, -1.0);

    theta_jk = std::acos(cos_jk);
    return theta_jk * theta_jk / mesh[index].ds;
}
*/

double dtmc_lc::ut_m(int index)
{ // director edge tangent,[don't think I use this anymore]
    if (mesh[index].edge_num == -1)
    {
        return 0;
    }
    double ut = 0;
    std::vector<double> rij = {0, 0, 0};
    for (int k = 0; k < 3; k++)
    {
        rij[k] = mesh[index].R[k] - mesh[mesh[index].edge_nei[0]].R[k];
    }
    ut = std::sqrt(1 - std::pow(innerproduct(mesh[index].u, rij), 2));
    return ut;
}
std::vector<double> dtmc_lc::dAn2H_m(int index)
{
    std::vector<double> dAn2H{0.0, 0.0};
    if (mesh[index].edge_nei.size() != 0)
    {
        return dAn2H;
    }
    double l_ij, l_ij0, l_ij2, l_jj0, l_jj2;
    double dot_ij0j, dot_ij2j, cot0, cot2, theta0, theta2;
    double sigma_ij, sigma_i;
    int ind_j0, ind_j, ind_j2;

    std::vector<double> r_ij{0, 0, 0};
    std::vector<double> r_ij0{0, 0, 0};
    std::vector<double> r_ij2{0, 0, 0};
    std::vector<double> r_jj0{0, 0, 0};
    std::vector<double> r_jj2{0, 0, 0};
    /* vertex configuration:
               - (j2)
             -
        (i)- - - (j)
             -
               - (j0)
    */
    sigma_i = 0;
    std::vector<double> dH_i{0, 0, 0}; // 2H_i, depends on r_ij
    for (int j = 0; j < mesh[index].nei.size(); j++)
    {
        // get nei index
        ind_j = mesh[index].nei[j];
        if (j == 0)
        {
            ind_j0 = mesh[index].nei[mesh[index].nei.size() - 1];
            ind_j2 = mesh[index].nei[j + 1];
        }
        else if (j == mesh[index].nei.size() - 1)
        {
            ind_j0 = mesh[index].nei[j - 1];
            ind_j2 = mesh[index].nei[0];
        }
        else
        {
            ind_j0 = mesh[index].nei[j - 1];
            ind_j2 = mesh[index].nei[j + 1];
        }
        // end set nei_ind
        // triangulation site to site vectors
        for (int k = 0; k < 3; k++)
        {
            r_ij[k] = mesh[index].R[k] - mesh[ind_j].R[k];
            r_ij0[k] = mesh[index].R[k] - mesh[ind_j0].R[k];
            r_ij2[k] = mesh[index].R[k] - mesh[ind_j2].R[k];
            r_jj0[k] = mesh[ind_j].R[k] - mesh[ind_j0].R[k];
            r_jj2[k] = mesh[ind_j].R[k] - mesh[ind_j2].R[k];
        }
        // site-site distance
        l_ij = std::sqrt(innerproduct(r_ij, r_ij));
        l_ij0 = std::sqrt(innerproduct(r_ij0, r_ij0));
        l_ij2 = std::sqrt(innerproduct(r_ij2, r_ij2));
        l_jj0 = std::sqrt(innerproduct(r_jj0, r_jj0));
        l_jj2 = std::sqrt(innerproduct(r_jj2, r_jj2));
        // inner product for cot calculation
        dot_ij0j = innerproduct(r_ij0, r_jj0);
        theta0 = std::acos(dot_ij0j / (l_ij0 * l_jj0));
        dot_ij2j = innerproduct(r_ij2, r_jj2);
        theta2 = std::acos(dot_ij2j / (l_ij2 * l_jj2));
        // get cot
        /*
        cot0 = dot_ij0j / std::sqrt(std::pow(l_ij0 * l_jj0, 2) -
                                    std::pow(dot_ij0j, 2));
        cot2 = dot_ij2j / std::sqrt(std::pow(l_ij2 * l_jj2, 2) -
                                    std::pow(dot_ij2j, 2));
        */
        cot0 = std::cos(theta0) / std::sin(theta0);
        cot2 = std::cos(theta2) / std::sin(theta2);
        sigma_ij = 0.5 * l_ij * (cot0 + cot2);
        sigma_i += 0.25 * sigma_ij * l_ij;
        for (int k = 0; k < 3; k++)
        {
            dH_i[k] += sigma_ij / l_ij * r_ij[k];
        }
    } // end for nei
    for (int k = 0; k < 3; k++)
    {
        dH_i[k] = dH_i[k] / sigma_i;
    }
    dAn2H[0] = sigma_i;
    dAn2H[1] = std::sqrt(innerproduct(dH_i, dH_i));
    // find the sign of H by comparing with local normal
    if (innerproduct(dH_i, mesh[index].n) < 0)
    {
        dAn2H[1] *= -1;
    }
    return dAn2H;
}

double dtmc_lc::dskg_m(int index)
{
    if (mesh[index].edge_nei.size() == 0)
    {
        return 0;
    }
    double kg;
    std::vector<double> r_ij{0, 0, 0};
    std::vector<double> r_ik{0, 0, 0};
    double l_ij, l_ik;
    int ind_j, ind_k;
    double cos_jk;
    /* vertex configuration:
               - (k)
             -
        (i)- - - (j)
    */
    std::vector<int> nei_original;
    // sort neighbors
    nei_original = mesh[index].nei;
    sort_nei(index);
    // get thetai
    // kg = PI - sum(theta_i)
    kg = PI;
    for (int j = 0; j < mesh[index].nei.size() - 1; j++)
    {
        cos_jk = 0;
        l_ij = 0;
        l_ik = 0;
        ind_j = mesh[index].nei[j];
        ind_k = mesh[index].nei[j + 1];
        for (int k = 0; k < 3; k++)
        {
            r_ij[k] = mesh[ind_j].R[k] - mesh[index].R[k];
            l_ij += r_ij[k] * r_ij[k];
            r_ik[k] = mesh[ind_k].R[k] - mesh[index].R[k];
            l_ik += r_ik[k] * r_ik[k];
            cos_jk += r_ij[k] * r_ik[k];
        }
        cos_jk = cos_jk / (std::sqrt(l_ij * l_ik));
        kg -= std::acos(cos_jk);
    }
    mesh[index].nei = nei_original; // put nei sort back;
    return kg;
}

double dtmc_lc::dAK_m(int index)
{
    // account both in-bulk and on-edge beads

    // if (mesh[index].edge_nei.size() != 0) {

    double K;
    std::vector<double> r_ij{0, 0, 0};
    std::vector<double> r_ik{0, 0, 0};
    double l_ij, l_ik;
    int ind_j, ind_k;
    double cos_jk;
    /* vertex configuration:
               - (k)
             -
        (i)- - - (j)
    */
    // get thetai
    // K = 2*PI - sum(theta_i) if in-bulk
    // K = PI - sum(theta_i) if on-edge
    std::vector<int> nei_original;
    // sort neighbors
    nei_original = mesh[index].nei;
    if (mesh[index].edge_nei.size() == 0)
    {
        K = 2 * PI;
    }
    else
    {
        return 0;
        //K = PI;
        //sort_nei(index);
    }

    for (int j = 0; j < mesh[index].nei.size(); j++)
    {
        cos_jk = 0;
        l_ij = 0;
        l_ik = 0;
        ind_j = mesh[index].nei[j];
        if (j == (mesh[index].nei.size() - 1))
        {
            if (mesh[index].edge_nei.size())
            {
                break;
                mesh[index].nei = nei_original;
            }
            ind_k = mesh[index].nei[0];
        }
        else
        {
            ind_k = mesh[index].nei[j + 1];
        }
        for (int k = 0; k < 3; k++)
        {
            r_ij[k] = mesh[ind_j].R[k] - mesh[index].R[k];
            l_ij += r_ij[k] * r_ij[k];
            r_ik[k] = mesh[ind_k].R[k] - mesh[index].R[k];
            l_ik += r_ik[k] * r_ik[k];
            cos_jk += r_ij[k] * r_ik[k];
        }
        cos_jk = cos_jk / (std::sqrt(l_ij * l_ik));
        K -= std::acos(cos_jk);
    }
    return K;
}

std::vector<double> dtmc_lc::n_m(int index)
{
    // measure normal, work for both in-bulk and on-edge beads
    int on_edge = 0; // flag indicating if it's an edge bead
    std::vector<double> n_ind{0, 0, 0};
    double n_ind_abs;
    int j_next, ind_j, ind_k;
    double l_ij, l_ik;
    std::vector<double> r_ij{0, 0, 0};
    std::vector<double> r_ik{0, 0, 0};
    std::vector<double> n_jk{0, 0, 0};
    double n_jk_abs;
    double theta_jk;
    std::vector<int> nei_original;
    nei_original = mesh[index].nei; // thus it has correct size from the beginning
    if (mesh[index].edge_num != -1)
    {
        on_edge = 1;
        sort_nei(index);
    }
    for (int j = 0; j < (mesh[index].nei.size() - on_edge); j++)
    {
        j_next = (j + 1) % mesh[index].nei.size();
        ind_j = mesh[index].nei[j];
        ind_k = mesh[index].nei[j_next];
        for (int k = 0; k < 3; k++)
        {
            r_ij[k] = mesh[ind_j].R[k] - mesh[index].R[k];
            r_ik[k] = mesh[ind_k].R[k] - mesh[index].R[k];
        }
        l_ij = std::sqrt(innerproduct(r_ij, r_ij));
        l_ik = std::sqrt(innerproduct(r_ik, r_ik));
        n_jk = crossproduct(r_ij, r_ik);
        n_jk_abs = std::sqrt(innerproduct(n_jk, n_jk));
        theta_jk = std::asin(n_jk_abs / (l_ij * l_ik));
        for (int k = 0; k < 3; k++)
        {
            n_jk[k] = n_jk[k] / n_jk_abs;
            n_ind[k] += theta_jk * n_jk[k];
        }
    }
    n_ind_abs = std::sqrt(innerproduct(n_ind, n_ind));
    for (int k = 0; k < 3; k++)
    {
        n_ind[k] = n_ind[k] / n_ind_abs;
    }
    // put nei sort back
    mesh[index].nei = nei_original;

    return n_ind;
}

double dtmc_lc::p2uu_m(int ind_i, int ind_j)
{
    // spin spin interaction of i and j
    double p2uu = std::pow(innerproduct(mesh[ind_i].u, mesh[ind_j].u), 2);
    return 1.5 * p2uu - 0.5;
}
double dtmc_lc::uuc_m(int ind_i, int ind_j)
{
    // spin twist interaction of i and j
    double uuc_local, lij;
    std::vector<double> rij{0, 0, 0};
    for (int k = 0; k < mesh[ind_i].R.size(); k++)
    {
        rij[k] = mesh[ind_j].R[k] - mesh[ind_i].R[k];
    }
    lij = std::sqrt(innerproduct(rij, rij));
    uuc_local = innerproduct(crossproduct(mesh[ind_i].u, mesh[ind_j].u), rij);
    uuc_local = uuc_local * innerproduct(mesh[ind_i].u, mesh[ind_j].u);
    uuc_local = uuc_local / lij;
    return uuc_local;
}

double dtmc_lc::uusb_m(int ind_i, int ind_j)
{
    // splay and bend interaction of i and j
    //[(u_i\cross u_j)\cross r_{ij}]^2
    double uusb_local, lij2;
    std::vector<double> uusb_vec{0, 0, 0};
    std::vector<double> rij{0, 0, 0};
    for (int k = 0; k < mesh[ind_i].R.size(); k++)
    {
        rij[k] = mesh[ind_j].R[k] - mesh[ind_i].R[k];
    }
    uusb_vec = crossproduct(crossproduct(mesh[ind_i].u, mesh[ind_j].u), rij);
    uusb_local = innerproduct(uusb_vec, uusb_vec);
    lij2 = innerproduct(rij, rij);
    uusb_local /= lij2;
    return uusb_local;
}
double dtmc_lc::uut_m(int ind_i, int ind_j)
{
    // twist interaction of i and j
    //[(u_i\cross u_j)\cdot lr_{ij}]^2
    double uut_local, lij2;
    std::vector<double> uusb_vec{0, 0, 0};
    std::vector<double> rij{0, 0, 0};
    for (int k = 0; k < mesh[ind_i].R.size(); k++)
    {
        rij[k] = mesh[ind_j].R[k] - mesh[ind_i].R[k];
    }
    uut_local = innerproduct(crossproduct(mesh[ind_i].u, mesh[ind_j].u), rij);
    uut_local = uut_local * uut_local;
    lij2 = innerproduct(rij, rij);
    uut_local /= lij2;
    return uut_local;
}
double dtmc_lc::un2_m(int index)
{
    // spin normal interaction of i
    return std::pow(innerproduct(mesh[index].u, mesh[index].n), 2);
}

#pragma endregion

#pragma region : membrane structure related
std::vector<double> dtmc_lc::Gij_m()
{
    // G_ij = 1/(2N)\sum_n\sum_m{(r_i(n)-r_i(m))(r_j(n)-r_j(m))}
    // G_ij = <r_i r_j> - <r_i><r_j>
    std::vector<double> Gij(9, 0);
    double xc = 0, yc = 0, zc = 0;
    for (int n = 0; n < mesh.size(); n++)
    {
        // G_xx
        Gij[0] += mesh[n].R[0] * mesh[n].R[0];
        // G_xy
        Gij[1] += mesh[n].R[0] * mesh[n].R[1];
        // G_xz
        Gij[2] += mesh[n].R[0] * mesh[n].R[2];
        // G_yx Gij[3] = Gij[1]
        // G_yy
        Gij[4] += mesh[n].R[1] * mesh[n].R[1];
        // G_yz
        Gij[5] += mesh[n].R[1] * mesh[n].R[2];
        // G_zx Gij[6] = Gij[2]
        // G_zy Gij[7] = Gij[5]
        // G_zz
        Gij[8] += mesh[n].R[2] * mesh[n].R[2];
        // center
        xc += mesh[n].R[0];
        yc += mesh[n].R[1];
        zc += mesh[n].R[2];
    }
    xc /= N;
    yc /= N;
    zc /= N;
    Gij[0] = Gij[0] / N - xc * xc;
    Gij[1] = Gij[1] / N - xc * yc;
    Gij[2] = Gij[2] / N - xc * zc;
    Gij[3] = Gij[1];
    Gij[4] = Gij[4] / N - yc * yc;
    Gij[5] = Gij[5] / N - yc * zc;
    Gij[6] = Gij[2];
    Gij[7] = Gij[5];
    Gij[8] = Gij[8] / N - zc * zc;
    return Gij;
}

double dtmc_lc::D_edge_com_m()
{
    // only work for Ne=2 now
    if (Ne != 2)
    {
        return 0;
    }
    // com of 2 edges
    std::vector<double> Rc0{0, 0, 0};
    std::vector<double> Rc1{0, 0, 0};
    double D_edge = 0;
    for (int i = 0; i < edge_lists[0].size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Rc0[j] += mesh[edge_lists[0][i]].R[j];
        }
    }
    for (int i = 0; i < edge_lists[1].size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Rc1[j] += mesh[edge_lists[1][i]].R[j];
        }
    }
    for (int j = 0; j < 3; j++)
    {
        Rc0[j] /= edge_lists[0].size();
        Rc1[j] /= edge_lists[1].size();
        D_edge += std::pow(Rc1[j] - Rc0[j], 2);
    }
    D_edge = std::sqrt(D_edge);
    return D_edge;
}
#pragma endregion
