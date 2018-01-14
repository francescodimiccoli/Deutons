#ifndef CUTS_H
#define CUTS_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include "TH1.h"
#include "TH2.h"
#include "Variables.hpp"

using namespace std;

extern std::vector<float> LatEdges;

TH1F* ProjectionXtoTH1F(TH2F* h2, string title, int binmin, int binmax);
TH1F* TH1DtoTH1F(TH1D* hd);

bool IsProtonMC    (Variables * vars);
bool IsDeutonMC    (Variables * vars);
bool IsHeliumMC    (Variables * vars);

bool IsFragmentedPfromHeMC (Variables * vars); 
bool IsFragmentedDfromHeMC (Variables * vars); 
bool IsFragmentedTMC 	   (Variables * vars); 
bool IsPureDMC 		   (Variables * vars); 
bool IsFragmentedPfromDMC  (Variables * vars); 
bool IsFragment		   (Variables * vars);

bool L1LooseCharge1(Variables * vars) ;
bool IsPrimary	   (Variables * vars);
bool IsMC          (Variables * vars); 
bool IsData        (Variables * vars);
bool IsPreselectedHe (Variables * vars);
bool IsPreselectedInner (Variables * vars);
bool IsPreselected (Variables * vars);
bool IsMinimumBias (Variables * vars);

bool IsFromNaF     (Variables * vars);
bool IsFromAgl     (Variables * vars);
bool IsOnlyFromToF (Variables * vars);
bool InnerAndL1Charge2 (Variables * vars);
bool ProtonsMassCut(Variables * vars);
bool DeutonsMassCut(Variables * vars);
bool TemplatesMassCut(Variables * vars);

bool ControlSampleMassCut(Variables * vars);

bool IsNegativeCharged (Variables * vars); 


bool Qualitycut(Variables * vars, float cutvariable, float cutTOF, float cutNaF, float cutAgl);

bool QualChargeCut (Variables * vars);

bool DistanceCut   (Variables * vars);
bool LikelihoodCut (Variables * vars);

bool IsGoodHe      (Variables * vars);

bool IsInLatZone   (Variables * vars, int lat) ;

bool ControlSample (Variables * vars) ;
bool PresControlSample (Variables * vars) ;

bool TofBetaSafetyCut (Variables * vars) ;
bool NafBetaSafetyCut (Variables * vars) ;
bool AglBetaSafetyCut (Variables * vars);


template<typename Out>
void split(const std::string &s, char delim, Out result);
std::vector<std::string> split(const std::string &s, char delim) ;

bool ApplyCuts(std::string cut, Variables * Vars);

int GetLatitude(Variables * vars);

TH1F* ProjectionXtoTH1F(TH2F* h2, string title, int binmin, int binmax);

TH1F* TH1DtoTH1F(TH1D* hd) ;

#endif
