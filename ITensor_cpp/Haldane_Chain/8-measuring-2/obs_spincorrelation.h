#ifndef OBS_SPINCORRELATION_H
#define OBS_SPINCORRELATION_H

#include "itensor/all.h"

double spincorrelation_2site(itensor::MPS& psi, const itensor::SiteSet& sites, int i, int j);
void spincorrelation(itensor::MPS& psi, const itensor::SiteSet& sites, int istart, int iend);
#endif 
