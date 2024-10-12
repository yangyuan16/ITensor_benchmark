
#include "itensor/all.h"
#include <stdio.h>
#include <iostream>
#include "sample/src/electronk.h"
#include "sample/src/hubbard.h"
#include <vector>

using namespace itensor;

int main(int argc, char* argv[])
{
	// set parameters
	// 2*W*L W=1
	int Nx = 64;// L
	int Ny = 4;// 2 * 2 orbitals 
	double U = 8.0;// 12.0
	if (argc > 3)
		U = std::stof(argv[3]);
	if (argc > 2)
		Ny = std::stoi(argv[2]);
	if (argc > 1)
		Nx = std::stoi(argv[1]);

	string filename;
	std::ofstream fout;

	int N = Nx * Ny;
	//SiteSet sites = ElectronK(N, args);
	auto sites = Electron(N, { "ConserveQNs =", true });

	auto txx = 1.0;
	auto Ez = 0.0;
	auto J = 1.0;
	auto Upri = U - 2*J;
	auto txzpri = 0.07;

	auto tzz = 0.241;
        auto txz = 0.506;
        auto tzz2 = 1.314;
        auto Ex = 1.116;

	string pathname = "Output/3000";

	// set the Hamiltonian MPO
	auto ampo = AutoMPO(sites);
	auto lattice = squareLattice(Nx, Ny, { "YPeriodic = ", false });
	// on-site repulsive
	for (int i = 1; i <= N; ++i)
	{
		ampo += U, "Nupdn", i; //
	}
	// Hund coupling
	for (int i = 1; i <= N; ++i)
	{
		if (i % 2 == 1) {
				ampo += Upri - J, "Nup", i, "Nup", i+1;
				ampo += Upri - J, "Ndn", i, "Ndn", i+1;
				ampo += Upri, "Nup", i, "Ndn", i+1;
				ampo += Upri, "Ndn", i, "Nup", i+1;
		}
	}
	// on-site energy of two orbitals
	for (int i = 1; i <= N; ++i) {
		if (i % 4 == 1 || i % 4 == 0) {
			ampo += Ex, "Cdagup", i, "Cup", i;
			ampo += Ex, "Cdagdn", i, "Cdn", i;
		}
		else {
			ampo += Ez, "Cdagup", i, "Cup", i;
			ampo += Ez, "Cdagdn", i, "Cdn", i;
		}
	}
	// hoping integral of NN
	for (int i = 1; i <= N; ++i) {
		if (i % 2 == 1) {
			if (i + 5 <= N) {
				ampo += -txz, "Cdagup", i, "Cup", i + 5;
				ampo += -txz, "Cdagup", i + 5, "Cup", i;
				ampo += -txz, "Cdagdn", i, "Cdn", i + 5;
				ampo += -txz, "Cdagdn", i + 5, "Cdn", i;
			}			
		}
		else {			
			if (i + 3 <= N) {
				ampo += -txz, "Cdagup", i, "Cup", i + 3;
				ampo += -txz, "Cdagup", i + 3, "Cup", i;
				ampo += -txz, "Cdagdn", i, "Cdn", i + 3;
				ampo += -txz, "Cdagdn", i + 3, "Cdn", i;
			}
		}
		if (i + 4 <= N) {
			if (i % 4 == 1 || i % 4 == 0) {
				ampo += -txx, "Cdagup", i, "Cup", i + 4;
				ampo += -txx, "Cdagup", i + 4, "Cup", i;
				ampo += -txx, "Cdagdn", i, "Cdn", i + 4;
				ampo += -txx, "Cdagdn", i + 4, "Cdn", i;
			}
			if (i % 4 == 2 || i % 4 == 3) {
				ampo += -tzz, "Cdagup", i, "Cup", i + 4;
				ampo += -tzz, "Cdagup", i + 4, "Cup", i;
				ampo += -tzz, "Cdagdn", i, "Cdn", i + 4;
				ampo += -tzz, "Cdagdn", i + 4, "Cdn", i;
			}
		}
		if (i % 4 == 2) {
			ampo += -tzz2, "Cdagup", i, "Cup", i + 1;
			ampo += -tzz2, "Cdagup", i + 1, "Cup", i;
			ampo += -tzz2, "Cdagdn", i, "Cdn", i + 1;
			ampo += -tzz2, "Cdagdn", i + 1, "Cdn", i;
		}
	}
	// hoping integral of NNN
	for (int i = 1; i <= N; ++i) {
		if (i % 4 == 1 && i % 4 == 2) {
			if (i + 6 <= N) {
				ampo += -txzpri, "Cdagup", i, "Cup", i + 6;
				ampo += -txzpri, "Cdagup", i + 6, "Cup", i;
				ampo += -txzpri, "Cdagdn", i, "Cdn", i + 6;
				ampo += -txzpri, "Cdagdn", i + 6, "Cdn", i;
			}			
		}
		else {			
			if (i + 2 <= N) {
				ampo += -txzpri, "Cdagup", i, "Cup", i + 2;
				ampo += -txzpri, "Cdagup", i + 2, "Cup", i;
				ampo += -txzpri, "Cdagdn", i, "Cdn", i + 2;
				ampo += -txzpri, "Cdagdn", i + 2, "Cdn", i;
			}
		}
	}
	auto H = toMPO(ampo);
	//
	// Set the initial wavefunction matrix product state of 3/8 filling
	// Sz = 0
	auto state = InitState(sites);
	for (auto j : range1(N))
	{
		if (j % 8 == 1) {
		 	state.set(j, "Up");
			state.set(j + 1, "Up");
			state.set(j + 2, "Dn");
			state.set(j + 3, "Emp");
			state.set(j + 4, "Dn");
			state.set(j + 5, "Dn");
			state.set(j + 6, "Up");
			state.set(j + 7, "Emp");
		}
	}
	// Set sweeps
	auto sweeps = Sweeps(10);
	sweeps.maxdim() = 100, 500, 1000, 2000,3000;
	sweeps.noise() = 1E-7, 1E-8, 1E-10, 0;
	sweeps.cutoff() = 1E-8;// truncation?

	PrintData(sweeps);

	auto psi0 = randomMPS(state);
	// auto psi0 = MPS(state);
	//auto [energy, psi] = dmrg(H, psi0, sweeps, { "Quiet=",true });
	auto [energy, psi] = dmrg(H, psi0, sweeps, { "Quiet=",true, "WriteDim", 1000, "WriteDir", "Output/PH" });

	PrintData(Nx);
	PrintData(Ny);
	PrintData(U);
	PrintData(tzz2);
	PrintData(tzz);
	PrintData(txz);
	PrintData(txx);
	PrintData(txzpri);
	PrintData(J);
	PrintData(Ex);
	PrintData(totalQN(psi));
	PrintData(maxLinkDim(psi));
	PrintData(energy);
	writeToFile(pathname + "sites_file", sites);
	writeToFile(pathname + "psi_file", psi);

	auto Nt = expect(psi, sites, "Ntot");
	printf(" The value of Nt is in %s Ntot.txt. :\n",pathname);
	std::stringstream ss;
	ss << pathname << "Ntot.txt";
	filename = ss.str();
	fout.precision(12);
	fout.open(filename.c_str(), std::ofstream::app); //start writing file
	for (const auto& el : Nt) fout << el << "\t";
	fout << std::endl;
	fout.close();	


	auto Fzz = correlationMatrix(psi, sites, "Sz", "Sz");
	printf(" The value of Fzz is in %s Fzz.txt.:\n",pathname);
	ss.str("");
	ss << pathname << "Fzz.txt";
	filename = ss.str();
	fout.precision(12);
	fout.open(filename.c_str(), std::ofstream::app); //start writing file
	for (const auto& rowV : Fzz) {
		for (const auto& el : rowV) fout << el << "\t";
		fout << std::endl;
	}
	fout.close();


	// compute the pair correlation of z2 orbital
	//
	int jsta = 18;
	int jend = 64-18;
	ss.str("");
	ss << pathname << "Pzz.txt";
	filename = ss.str();
	fout.precision(12);
	fout.open(filename.c_str(), std::ofstream::app); //start writing file
	auto PCMz = AutoMPO(sites);
	for (int j = jsta; j <= jend; j++){
		int jz1 = 4 * (j - 1) + 2;
	        int jz2 = 4 * (j - 1) + 3;

		for (int r = jsta; r <= jend; r++) {
			int rz1 = 4 * (r - 1) + 2;
			int rz2 = 4 * (r - 1) + 3;
			PCMz = AutoMPO(sites);

			PCMz += 0.5, "Cdagdn", jz2, "Cdagup", jz1, "Cup", rz1, "Cdn", rz2;
			PCMz += -0.5, "Cdagdn", jz2, "Cdagup", jz1, "Cdn", rz1, "Cup", rz2;
			PCMz += -0.5, "Cdagup", jz2, "Cdagdn", jz1, "Cup", rz1, "Cdn", rz2;
			PCMz += 0.5, "Cdagup", jz2, "Cdagdn", jz1, "Cdn", rz1, "Cup", rz2;
			auto Pzz = toMPO(PCMz);
			fout << inner(psi, Pzz, psi) << "\t";
			}
		fout << std::endl;
	}
	fout.close();
        printf(" The value of Pzz is in %s Pzz.txt.:\n",pathname);

	
	ss.str("");
	ss << pathname << "Pxy.txt";
	filename = ss.str();
	fout.precision(12);
	fout.open(filename.c_str(), std::ofstream::app); //start writing file
	auto PCMxy = AutoMPO(sites);
	for (int j = jsta; j <= jend; j++){
		int jxy1 = 4 * (j - 1) + 1;
	        int jxy2 = 4 * (j - 1) + 4;

		for (int r = jsta; r <= jend; r++) {
			int rxy1 = 4 * (r - 1) + 1;
			int rxy2 = 4 * (r - 1) + 4;
			PCMxy = AutoMPO(sites);

			PCMxy += 0.5, "Cdagdn", jxy2, "Cdagup", jxy1, "Cup", rxy1, "Cdn", rxy2;
			PCMxy+= -0.5, "Cdagdn", jxy2, "Cdagup", jxy1, "Cdn", rxy1, "Cup", rxy2;
			PCMxy += -0.5, "Cdagup", jxy2, "Cdagdn", jxy1, "Cup", rxy1, "Cdn", rxy2;
			PCMxy += 0.5, "Cdagup", jxy2, "Cdagdn", jxy1, "Cdn", rxy1, "Cup", rxy2;
			auto Pxy = toMPO(PCMxy);
			fout << log10(fabs(inner(psi, Pxy, psi))) << "\t";
		}
		fout << std::endl;
	}
	fout.close();
	printf(" The value of Pxy is in %s Pxy.txt.:\n",pathname);

	return 0;
}
