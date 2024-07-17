#ifndef OBS_SPINPAIRING_H
#define OBS_SPINPAIRING_H

#include "itensor/all.h"

double spincorrelation_4site(itensor::MPS& psi, const itensor::SiteSet& sites, int i, int j, int m, int n);

#endif