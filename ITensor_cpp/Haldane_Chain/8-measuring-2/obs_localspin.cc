#include "obs_localspin.h"
#include "itensor/all.h"
#include <fstream>

using namespace itensor;

void localspin(MPS& psi, const SiteSet& sites, int N){
    println("\nj Sz = ");

    std::ofstream szFile("measure_localspin.dat");

    if (!szFile.is_open()){
        println("Error opening Sz.dat file.");
        return;
    }

    szFile << "j Sz\n";

    for (auto j : range1(N)){
        // re-gauge psi to get ready to measure at position j
        psi.position(j);

        auto ket = psi(j);

        auto bra = dag(prime(ket, "Site"));

        auto Szjop = op(sites, "Sz", j);

        // take an inner product
        auto szj = elt(bra * Szjop * ket);

        szFile << j << " " << szj << "\n";
        //printfln("%d %.12f", j, szj);
    }
    szFile.close();
}
