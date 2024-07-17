
#include "itensor/all.h"
#include "obs_spincorrelation.h"
#include "obs_localspin.h"
#include "obs_localbond.h"
#include "obs_spinpairing.h"
#include "obs_entanglement.h"

using namespace itensor;

int main(){
    
    int N = 40;

    auto sites = SpinOne(N, {"ConserveQNs=", false});

    auto ampo = AutoMPO(sites);
    for (int j = 1; j<N; ++j){
        ampo += 0.5, "S+", j, "S-", j+1;
        ampo += 0.5, "S-", j, "S+", j+1;
        ampo += "Sz", j, "Sz", j+1;
    }
    auto H = toMPO(ampo);

    auto sweeps = Sweeps(5); // number of sweeps is 5
    sweeps.maxdim() = 10, 20, 100, 100, 200;
    sweeps.cutoff() = 1E-10;

    auto psi0 = randomMPS(sites);

    auto [energy, psi] = dmrg(H, psi0, sweeps, {"Quiet=", true});

    // Measuring Sz
    localspin(psi, sites, N); 

    // Measuring SiSi+1
    localbond(psi,sites,N);

    // Measuring Spin Correlation
    spincorrelation(psi, sites, 10, 20);

    // Measuring 4 site spin correlation
    spincorrelation_4site(psi, sites, 3,4, 5,6);
    spincorrelation_4site(psi, sites, 3,4, 6,7);
    spincorrelation_4site(psi, sites, 3,4, 8,9);
    spincorrelation_4site(psi, sites, 3,6, 8,11);
    spincorrelation_4site(psi, sites, 3,6, 9,12);
    spincorrelation_4site(psi, sites, 3,6, 10, 15);
    
    // Measuring entanglement
    entanglement(psi, sites, N);
    
    return 0;

}
