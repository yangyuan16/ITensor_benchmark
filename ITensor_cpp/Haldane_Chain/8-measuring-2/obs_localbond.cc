#include "obs_localbond.h"
#include "itensor/all.h"
#include <fstream>

using namespace itensor;

void localbond(MPS& psi, const SiteSet& sites, int N){
    
    Real totalSdS = 0.;
    
    println("\nj S.S = ");

    std::ofstream File("measure_localbond.dat");

    if (!File.is_open()){
        println("Error opening measure_localbond.dat file");
        return;
    }

    for (auto b:range1(N-1)){
        psi.position(b);
        
        auto bondket = psi(b) * psi(b+1);
        auto bondbra = dag(prime(bondket, "Site"));

        auto zzop = op(sites, "Sz", b) * op(sites, "Sz", b+1);
        auto pmop = 0.5 * op(sites, "S+", b) * op(sites, "S-", b+1);
        auto mpop = 0.5 * op(sites, "S-", b) * op(sites, "S+", b+1);

        auto zz = elt(bondbra * zzop * bondket);
        auto pm = elt(bondbra * pmop * bondket);
        auto mp = elt(bondbra * mpop * bondket);

        printfln("%d %.12f", b, zz + pm + mp);

        File << b << " " << (b+1) << " " << (zz + pm + mp) << "\n";

        totalSdS += zz + pm + mp;

    }

    File.close();

    printfln("\nSum of S.S = %.12f", totalSdS);

}
