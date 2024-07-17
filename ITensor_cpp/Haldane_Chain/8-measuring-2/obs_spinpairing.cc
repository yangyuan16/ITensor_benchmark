#include "obs_spinpairing.h"
#include "itensor/all.h"

using namespace itensor;

double spincorrelation_4site(MPS& psi, const SiteSet& sites, int i, int j, int m, int n){

    auto op_i = op(sites, "Sz", i);
    auto op_j = op(sites, "Sz", j);
    auto op_m = op(sites, "Sz", m);
    auto op_n = op(sites, "Sz", n);

    // "gauge" the MPS to site i
    psi.position(i);

    // "Create" the bra/dual version of MPS psi
    auto psidag = dag(psi);

    // Prime the link indices to make them distinct from the original ket links
    psidag.prime("Link");

    // index linking i-1 to i;
    auto li_1 = leftLinkIndex(psi, i);

    auto C = prime(psi(i),li_1) * op_i;
    C *= prime(psidag(i),"Site");
    // from i+1 to j-1
    for (int k=i+1; k<j; ++k){
        C *= psi(k);
        C *= psidag(k);
    }
    // add op_j 
    C *= psi(j) * op_j;
    C *= prime(psidag(j),"Site");
    // from j+1 to m-1
    for (int k=j+1; k<m; ++k){
        C *= psi(k);
        C *= psidag(k);
    }
    // add op_m
    C *= psi(m) * op_m;
    C *= prime(psidag(m),"Site");
    // from m+1 to n-1
    for (int k=m+1; k<n; ++k){
        C *= psi(k);
        C *= psidag(k);
    }
    // add op_n
    // index linking n to n+1;
    auto ln = rightLinkIndex(psi, n);

    C *= prime(psi(n), ln) * op_n;
    C *= prime(psidag(n), "Site");

    double result = elt(C); // or eltC(C) if expecting complex 
    

    std::cout << "i " << "j " << "m " << "n " << " " << result << std::endl; 

    return result;

}