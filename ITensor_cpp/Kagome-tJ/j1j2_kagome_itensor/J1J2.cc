#include "itensor/all.h"
#include <iostream>
#include <chrono>
#include <fstream>
#include <math.h>

using namespace itensor;
std::pair<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>> kagomeAdjacency(int Lx, int Ly);
void Measure_Szcor(SiteSet const& sites, MPS& psi0);
std::pair<MPS, SiteSet> run_dmrg(int Lx, int Ly, Real J1, Real J2, 
             std::vector<int> const& maxdims, Real cutoff);
MPO constructHamiltonian(SiteSet const& sites, int Lx, int Ly, Real J1, Real J2);

// 设置系统的参数
int Lx = 12; // x方向正三角个数
int Ly = 9; // y方向三角个数*3
int N = Ly * Lx; //总格点数
Real J1 = 1.0; // 交换参数
Real J2 = 1.0 * 0.4; // 交换参数

// 设置DMRG参数
std::vector<int> maxdims = {1000, 2000, 4000, 6000, 7000, 7000, 7000, 7000, 7000};
std::vector<int> continue_sweep = {600, 700, 800, 900};
Real cutoff = 1E-7;

// 1: first run; 2: continue run; 3: just measure
int Switch = 1; 

int main()
{
    if(Switch == 1) 
    {
        auto [psi0, sites] = run_dmrg(Lx, Ly, J1, J2, maxdims, cutoff);
        writeToFile("final_psi0", psi0);   
        writeToFile("sites_file", sites);
	Measure_Szcor(sites, psi0);
    } 
    else if(Switch == 2)
    {
        auto psi = readFromFile<MPS>("final_psi0");
        auto sites = readFromFile<tJ>("sites_file");
        MPO H = constructHamiltonian(sites, Lx, Ly, J1, J2);
 
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
        
	Measure_Szcor(sites, psi0);
    }
    else 
    {
        auto psi0 = readFromFile<MPS>("final_psi0");
        auto sites = readFromFile<tJ>("sites_file");
	Measure_Szcor(sites, psi0);
    }
    return 0;
}

std::pair<MPS, SiteSet> run_dmrg(int Lx, int Ly, Real J1, Real J2, 
	     std::vector<int> const& maxdims, Real cutoff)
{
    auto sites = SpinHalf(N);
 
    // 定义初始状态，确保空穴数固定
    auto state = InitState(sites);

    for(auto i : range1(N))
    {
        if(i%2 == 1) state.set(i,"Up");
        else         state.set(i,"Dn");
    }
     

    // 创建初始MPS
    auto psi = MPS(state);

    // 生成格子和哈密顿量
    MPO H = constructHamiltonian(sites, Lx, Ly, J1, J2);    
 
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

    
MPO constructHamiltonian(SiteSet const& sites, int Lx, int Ly, Real J1, Real J2)
{
    auto [bonds_nn, bonds_nnn] = kagomeAdjacency(Lx, Ly);
    auto ampo = AutoMPO(sites);
    for(auto b : bonds_nn)
    {
        int i = b.first;
        int j = b.second;
        ampo += 0.5*J1, "S+", i, "S-", j;
        ampo += 0.5*J1, "S-", i, "S+", j;
        ampo += J1, "Sz", i, "Sz", j;
    }
    for(auto b : bonds_nnn)
    {
        int i = b.first;
        int j = b.second;
        ampo += 0.5*J2, "S+", i, "S-", j;
        ampo += 0.5*J2, "S-", i, "S+", j;
        ampo += J2, "Sz", i, "Sz", j;
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


