#ifndef LIVETIME_H
#define LIVETIME_H

#include "TH1.h"
#include "binning.h"

void UpdateZoneLivetime (float Livetime, float Rcutoff, TH1F * esposizionegeo,Binning bins, float prescale=1);

#endif

