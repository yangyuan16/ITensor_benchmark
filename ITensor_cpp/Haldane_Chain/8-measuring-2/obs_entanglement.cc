#include "obs_entanglement.h"
#include "itensor/all.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

using namespace itensor;

void entanglement(MPS& psi, const SiteSet& sites, int N){
    
    std::string folderName = "entanglement";

    // Create a directory
    if (mkdir(folderName.c_str(),0777)==-1){
        std::cerr << "Error creating directory: " << folderName << std::endl;
    } else {
        std::cout << "Directory created: " << folderName << std::endl;
    }

    // Create entropy file
    std::string entropyfilename = folderName + "/entropy.dat";
    // Open the file
    std::ofstream entropyfile(entropyfilename);
    // Check if the file was opened successfully
    if (!entropyfile){
        std::cerr << "Error creating entropy file: " << entropyfilename << std::endl;
    }

    for (int b=1; b<N; ++b){

        // "Gauge" the MPS to site b
        psi.position(b);

        // SVD this wavefunction to get the spectrum of density-matrix eigenvalues

        auto l = leftLinkIndex(psi,b);
        auto s = siteIndex(psi,b);
        auto [U, S, V] = svd(psi(b),{l,s});
        auto u = commonIndex(U,S);
        
        // Create spectrum file
        std::string spectrumfilename = folderName + "/spectrum" + std::to_string(b) + ".dat";
        // Open the file
        std::ofstream spectrumfile(spectrumfilename);
        //Check if the file was open successfully
        if (!spectrumfile){
            std::cerr << "Error creating spectrum file: " << spectrumfilename << std::endl;
        }

        //Apply von Neumann formula to the square of the singular values
        Real SvN = 0.;
        for (auto n:range1(dim(u)))
        {
            auto Sn = elt(S, n, n);
            auto p = sqr(Sn);
            if (p > 1E-12) SvN += -p*log(p);
            //printfln("Across bond b=%d, the %d-th eigs of entangle. spectrum is: %.10f", b, n, Sn);
            // Write spectrum file 
            spectrumfile << n << " " << Sn << std::endl;
        }
        // Close the spectrum file
        spectrumfile.close();
    printfln("Across bond b=%d, SvN= %.10f", b, SvN); 
    
    //Write entropy into file:
    entropyfile << b << " " << SvN << std::endl;
    }

    //Close entropy file
    entropyfile.close(); 
}