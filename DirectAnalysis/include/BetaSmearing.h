#ifndef BETASMEAR_H
#define BETASMEAR_H

#include "BadEventSimulator.h"
#include "Variables.hpp"

float SmearBetaRICH(Variables * vars);

float SmearBeta(Variables * vars);

float GetSmearedBetaTOF  (Variables * vars);   

float GetSmearedBetaRICH  (Variables * vars);   

float GetSmearedRecMassTOF(Variables * vars); 


#endif
