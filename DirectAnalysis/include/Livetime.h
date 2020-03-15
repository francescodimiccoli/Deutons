#ifndef LIVETIME_H
#define LIVETIME_H

#include "TH1.h"
#include "binning.h"
#include "Globals.h"

float GetBeta( float R, float mass);
void UpdateZoneLivetime (float Livetime, float Rcutoff, TH1F * esposizionegeo,Binning bins, float prescale=1);

#endif

