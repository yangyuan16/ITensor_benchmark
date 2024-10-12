#include "itensor/all.h"
#include <iostream>
#include <chrono>
#include <fstream>
#include <math.h>

using namespace itensor;
std::pair<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>> kagomeAdjacency(int Lx, int Ly);
void Measure_charge_density(int N, SiteSet const& sites, MPS& psi0);
void Measure_densitycor(SiteSet const& sites, MPS& psi0);
void Measure_Particlecor(SiteSet const& sites, MPS& psi0);
void Measure_Szcor(SiteSet const& sites, MPS& psi0);
void Measure_Sccor(int Ly, int N, int Rf_site, SiteSet const& sites, MPS& psi0);
std::pair<MPS, SiteSet> run_dmrg(int Lx, int Ly, Real t1, Real J1, Real t2, Real J2, int numHoles, 
             std::vector<int> const& maxdims, Real cutoff);
MPO constructHamiltonian(SiteSet const& sites, int Lx, int Ly, Real t1, Real J1, Real t2, Real J2);

// 设置系统的参数
int Lx = 5; // x方向正三角个数
int Ly = 6; // y方向三角个数*3
int N = Ly * Lx; //总格点数
Real t1 = 3.0; // 跳跃参数
Real J1 = 1.0; // 交换参数
Real t2 = 3.0 * 0.0; // 跳跃参数
Real J2 = 1.0 * 0.0; // 交换参数
int numHoles = 10; // 空穴数

// 设置DMRG参数
std::vector<int> maxdims = {100, 200, 200, 300, 300, 300, 400, 400, 400, 400};
std::vector<int> continue_sweep = {600, 700, 800, 900};
Real cutoff = 1E-7;

// 1: first run; 2: continue run; 3: just measure
int Switch = 3; 
int Rf_site = 5;

int main()
{
    if(Switch == 1) 
    {
        auto [psi0, sites] = run_dmrg(Lx, Ly, t1, J1, t2, J2, numHoles, maxdims, cutoff);
        writeToFile("final_psi0", psi0);   
        writeToFile("sites_file", sites);
        Measure_charge_density(N, sites, psi0);
        Measure_densitycor(sites, psi0);
	Measure_Particlecor(sites, psi0);
	Measure_Szcor(sites, psi0);
	Measure_Sccor(Ly, N, Rf_site, sites, psi0);
    } 
    else if(Switch == 2)
    {
        auto psi = readFromFile<MPS>("final_psi0");
        auto sites = readFromFile<tJ>("sites_file");
        MPO H = constructHamiltonian(sites, Lx, Ly, t1, J1, t2, J2);
 
        auto sweeps = Sweeps(continue_sweep.size());
        auto maxdimSetter = sweeps.maxdim();
        for(auto const& maxdim : continue_sweep)
        {
            maxdimSetter = maxdim;
        }
        sweeps.cutoff() = cutoff;
        sweeps.niter() = 2;
    
        auto [energy, psi0] = dmrg(H, psi, sweeps, {"Quiet", true, "WriteDim", 100});
        println("Ground State Energy = ", energy);
        writeToFile("final_psi0", psi0);   
        writeToFile("sites_file", sites);
        
        Measure_charge_density(N, sites, psi0);
        Measure_densitycor(sites, psi0);
	Measure_Particlecor(sites, psi0);
	Measure_Szcor(sites, psi0);
	Measure_Sccor(Ly, N, Rf_site, sites, psi0);
    }
    else 
    {
        auto psi0 = readFromFile<MPS>("final_psi0");
        auto sites = readFromFile<tJ>("sites_file");
        Measure_charge_density(N, sites, psi0);
        Measure_densitycor(sites, psi0);
	Measure_Particlecor(sites, psi0);
	Measure_Szcor(sites, psi0);
	Measure_Sccor(Ly, N, Rf_site, sites, psi0);
    }
    return 0;
}

std::pair<MPS, SiteSet> run_dmrg(int Lx, int Ly, Real t1, Real J1, Real t2, Real J2, int numHoles, 
	     std::vector<int> const& maxdims, Real cutoff)
{
    // 创建一个t-J模型的Hilbert空间，启用粒子数守恒和自旋守恒
    auto sites = tJ(N, {"ConserveNf", true, "ConserveSz", true});
 
    // 定义初始状态，确保空穴数固定
    auto state = InitState(sites);
    int holeCount = 0;
    int flag1 = N / numHoles;
    int flag2 = 0;
    for(int i = 1; i <= N; ++i)
    {
        if(i % flag1 == 2 && holeCount < numHoles)
        {
            state.set(i, "0"); // 空穴
            holeCount++;
        }
        else if(flag2 == 0)
        {
            state.set(i, "Up"); // 自旋向上
            flag2 = 1;
	}
        else
        {
            state.set(i, "Dn"); // 自旋向下
            flag2 = 0;
	}
    }

    // 创建初始MPS
    auto psi = MPS(state);

    // 生成格子和哈密顿量
    MPO H = constructHamiltonian(sites, Lx, Ly, t1, J1, t2, J2);    
 
    // 使用DMRG算法计算基态
    auto sweeps = Sweeps(maxdims.size());
    auto maxdimSetter = sweeps.maxdim();
    for(auto const& maxdim : maxdims)
    {
        maxdimSetter = maxdim;
    }
    sweeps.cutoff() = cutoff;
    sweeps.niter() = 2;

    auto [energy, psi0] = dmrg(H, psi, sweeps, {"Quiet", true, "WriteDim", 100});

    println("Ground State Energy = ", energy);

    return {psi0, sites};
}


void Measure_charge_density(int N, SiteSet const& sites, MPS& psi0)
{
    // 电荷密度分布ok
    std::ofstream outfile11("11_charge_density.txt");
    double charge_density[N];
    for(int i = 1; i <= N; ++i)
    {
        psi0.position(i);
        auto n_op = op(sites, "Ntot", i);
        charge_density[i - 1] = elt(dag(prime(psi0.A(i),"Site"))*n_op*psi0.A(i));
        outfile11 << i << " " << charge_density[i - 1] << "\n";
    }
    outfile11.close();
}

void Measure_densitycor(SiteSet const& sites, MPS& psi0)
{
    // 密度密度关联函数ok
    std::ofstream outFile12("12_density_correlation.txt");
    auto DD = correlationMatrix(psi0,sites,"Ntot","Ntot");
    for (size_t i = 0; i < DD.size(); ++i)
    {
        for (size_t j = 0; j < DD[i].size(); ++j)
        {
            outFile12 << i+1 << " " << j+1 << " " << DD[i][j] << " " << "\n"; 
        }
    } 
    outFile12.close();
}

void Measure_Particlecor(SiteSet const& sites, MPS& psi0)
{
    // 单粒子关联函数ok
    std::ofstream outFile13("13_single_particle_correlation.txt");
    auto up = correlationMatrix(psi0,sites,"Cdagup","Cup");
    auto dn = correlationMatrix(psi0,sites,"Cdagdn","Cdn");
    for (size_t i = 0; i < up.size(); ++i)
    {
        for (size_t j = 0; j < up[i].size(); ++j)
        {
            outFile13 << i+1 << " " << j+1 << " " << up[i][j]+dn[i][j] << "\n"; 
        }
    }
 }

void Measure_Szcor(SiteSet const& sites, MPS& psi0)
{
    // 自旋关联函数ok
    std::ofstream outFile21("21_spin_correlation.txt");
    auto czz = correlationMatrix(psi0,sites,"Sz","Sz");
    auto sz = expect(psi0,sites,"Sz"); 
    for (size_t i = 0; i < czz.size(); ++i)
    {
        for (size_t j = 0; j < czz[i].size(); ++j)
        {
            outFile21 << i+1 << " " << j+1 << " " << czz[i][j] << " " << sz[i] << " " << sz[j] << "\n"; 
        }
    } 
    outFile21.close();
   
}

    
void Measure_Sccor(int Ly, int N, int Rf_site, SiteSet const& sites, MPS& psi0)
{
    // 超导关联函数ok
    std::ofstream outFile3("3_Superconductivity_correlation.txt");
    auto sc_cor = AutoMPO(sites);
    for (int i = Rf_site + Ly; i < N; i+=Ly){
        auto sc_cor = AutoMPO(sites);
        sc_cor += 0.5, "Cdagup", Rf_site-1, "Cdagdn", Rf_site, "Cup", i-1, "Cdn", i;
        sc_cor += -0.5, "Cdagup", Rf_site-1, "Cdagdn", Rf_site, "Cdn", i-1, "Cup", i;
        sc_cor += -0.5, "Cdagdn", Rf_site-1, "Cdagup", Rf_site, "Cup", i-1, "Cdn", i;
        sc_cor += 0.5, "Cdagdn", Rf_site-1, "Cdagup", Rf_site, "Cdn", i-1, "Cup", i; 
        auto Sc_cor = toMPO(sc_cor);
        outFile3 << Rf_site-1 << " " << i-1 << " " << inner(psi0, Sc_cor, psi0) << "\n";        
    }
    outFile3.close();
}

MPO constructHamiltonian(SiteSet const& sites, int Lx, int Ly, Real t1, Real J1, Real t2, Real J2)
{
    auto [bonds_nn, bonds_nnn] = kagomeAdjacency(Lx, Ly);
    auto ampo = AutoMPO(sites);
    for(auto b : bonds_nn)
    {
        int i = b.first;
        int j = b.second;
        ampo += -t1, "Cdagup", i, "Cup", j;
        ampo += -t1, "Cdagup", j, "Cup", i;
        ampo += -t1, "Cdagdn", i, "Cdn", j;
        ampo += -t1, "Cdagdn", j, "Cdn", i;
        ampo += 0.5*J1, "S+", i, "S-", j;
        ampo += 0.5*J1, "S-", i, "S+", j;
        ampo += J1, "Sz", i, "Sz", j;
        ampo += -0.25*J1, "Ntot", i, "Ntot", j;
    }
    for(auto b : bonds_nnn)
    {
        int i = b.first;
        int j = b.second;
        ampo += -t2, "Cdagup", i, "Cup", j;
        ampo += -t2, "Cdagup", j, "Cup", i;
        ampo += -t2, "Cdagdn", i, "Cdn", j;
        ampo += -t2, "Cdagdn", j, "Cdn", i;
        ampo += 0.5*J2, "S+", i, "S-", j;
        ampo += 0.5*J2, "S-", i, "S+", j;
        ampo += J2, "Sz", i, "Sz", j;
        ampo += -0.25*J2, "Ntot", i, "Ntot", j;
    }
    
    auto H = toMPO(ampo);
    return H;
}


//------Kagome lattice(Lx,Ly)

  //            1   2   10  11  19  20
  //          9       18      27
  //        7   8   16  17  25  26
  //      6       15      24
  //    4   5   13  14  22  23
  //  3       12      21
  //1   2   10  11  19  20


// 定义 Kagome 格子的邻接表
std::pair<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>> kagomeAdjacency(int Lx, int Ly)
{
    std::vector<std::pair<int, int>> bonds_nn;
    std::vector<std::pair<int, int>> bonds_nnn;

    // 最近邻关系
    //----正三角
    for(int i = 0; i < Lx * Ly/3; ++i)
    {
        bonds_nn.push_back({3 * i + 1, 3 * i + 2});
        bonds_nn.push_back({3 * i + 1, 3 * i + 3});
        bonds_nn.push_back({3 * i + 2, 3 * i + 3});
    }
    //----倒三角1
    for(int i = 0; i < Lx - 1; ++i)
    {   
        for (int j = 0; j < Ly/3 - 1; ++j)
        {
	    bonds_nn.push_back({Ly * i + 3 * j + 5, Ly * i + 3 * j + Ly + 3});
   	    bonds_nn.push_back({Ly * i + 3 * j + 5, Ly * i + 3 * j + Ly + 4});
      	    bonds_nn.push_back({Ly * i + 3 * j + Ly + 3, Ly * i + 3 * j + Ly + 4});
	}
    }
    //----倒三角2
    for(int i = 0; i < Lx - 1; ++i)
    {
        bonds_nn.push_back({Ly * i + 2, Ly * i + Ly + 1});
        bonds_nn.push_back({Ly * i + 2, Ly * i + Ly * 2});
        bonds_nn.push_back({Ly * i + Ly + 1, Ly * i + Ly * 2});
    }
    //----左边线
    for(int i = 0; i < Ly/3 - 1; ++i)
    {
        bonds_nn.push_back({3 * i + 3, 3 * i + 4});
    }
    bonds_nn.push_back({1, Ly});
    

    // 次近邻关系
    //----六边形1
    for(int i = 0; i < Lx - 1; ++i)
    {
        for(int j = 0; j < Ly/3 - 1; ++j)
	{
            bonds_nnn.push_back({Ly * i + 3 * j + 2, Ly * i + 3 * j + 4});
            bonds_nnn.push_back({Ly * i + 3 * j + 3, Ly * i + 3 * j + 5});
            bonds_nnn.push_back({Ly * i + 3 * j + 4, Ly * i + 3 * j + Ly + 3});
            bonds_nnn.push_back({Ly * i + 3 * j + 5, Ly * i + 3 * j + Ly + 1});
            bonds_nnn.push_back({Ly * i + 3 * j + Ly + 3, Ly * i + 3 * j + 2});
            bonds_nnn.push_back({Ly * i + 3 * j + Ly + 1, Ly * i + 3 * j + 3});
 	}
    }
    //----六边形2
    for(int i = 0; i < Lx - 1; ++i)
    {
        bonds_nnn.push_back({Ly * i + 1, Ly * i + Ly - 1});
        bonds_nnn.push_back({Ly * i + 2, Ly * i + Ly});
        bonds_nnn.push_back({Ly * i + Ly * 2, Ly * i + 1});
        bonds_nnn.push_back({Ly * i + Ly * 2 - 2, Ly * i + 2});
        bonds_nnn.push_back({Ly * i + Ly - 1, Ly * i + Ly * 2});
        bonds_nnn.push_back({Ly * i + Ly, Ly * i + Ly * 2 - 2});
    }
    //----右边线
    for(int i = 0; i < Ly/3 - 1; ++i)
    {
	bonds_nnn.push_back({Ly * (Lx - 1) + 3 * i + 2, Ly * (Lx - 1) + 3 * i + 4});
        bonds_nnn.push_back({Ly * (Lx - 1) + 3 * i + 3, Ly * (Lx - 1) + 3 * i + 5});
    }
    bonds_nnn.push_back({Ly * (Lx - 1) + 1, Ly * Lx - 1});
    bonds_nnn.push_back({Ly * (Lx - 1) + 2, Ly * Lx});
    
    return {bonds_nn, bonds_nnn};
}


