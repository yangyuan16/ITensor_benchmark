#include "obs_spincorrelation.h"
#include "itensor/all.h"

using namespace itensor;

double spincorrelation_2site(MPS& psi, const SiteSet& sites, int i, int j){

    auto op_i = op(sites, "Sz", i);
    auto op_j = op(sites, "Sz", j);

    // "gauge" the MPS to site i
    psi.position(i);

    // "Create" the bra/dual version of MPS psi
    auto psidag = dag(psi);

    // Prime the link indices to make them distinct from the original ket links
    psidag.prime("Link");

    // index linking i-1 to i:
    auto li_1 = leftLinkIndex(psi, i);

    auto C = prime(psi(i),li_1)*op_i;
    C *= prime(psidag(i), "Site");
    for (int k=i+1; k<j; ++k)
    {
        C *= psi(k);
        C *= psidag(k);   
    }

    // index linking j to j+1;
    auto lj = rightLinkIndex(psi, j);

    C *= prime(psi(j),lj) * op_j;
    C *= prime(psidag(j), "Site");

    double result = elt(C); // or eltC(C) if expecting complex

    return result; 
}

void spincorrelation(MPS& psi, const SiteSet& sites, int istart, int iend){

    std::ofstream File("measure_spincorrelation.dat");
    if (!File.is_open()){
        println("Error opening measure_spincorrelation.dat file.");
        return;
    } 
    for (int it= istart+1; it <= iend; ++it)
    {
        double corre = spincorrelation_2site(psi, sites, istart, it);

        File << istart << " " << it << " " << corre << "\n";
    }
}
