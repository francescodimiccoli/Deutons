#ifndef GLOBALS_H
#define GLOBALS_H


#include "TF1.h"
#include "binning.h"
#include "TH1.h"
//#include "GlobalPaths.h"
#include "particle.h"
#include "RangeMerger.h"
#include <TMVA/Reader.h>
#include <TMVA/Tools.h>
#include "TRandom3.h"


using namespace std;


extern int Ev_Num;
extern int Timebeg;
extern float FRAC;

extern float RcutoffCut;
extern float BetacutoffCut;


extern int nbinsr;
extern int nbinsToF;
extern int nbinsNaF;
extern int nbinsAgl;

extern float ToFsmearSigma;
extern float ToFsmearShift;

extern TRandom3 * Rand;

extern Particle proton;  
extern Particle deuton; 

extern RangeMerger Global;
extern Binning HefragmToF;
extern Binning HefragmNaF;
extern Binning HefragmAgl;

extern Binning DRB;
extern Binning PRB;
extern Binning ForAcceptance;
extern Binning ForEffCorr;
extern Binning ForEffCorr_D;

extern Binning ForCutoff;

extern float con_min;
extern float con_max;

//resolution binning
extern RangeMerger GlobalRig;

extern Binning PResB;

extern TMVA::Reader *readerTOF;
extern TMVA::Reader *readerNaF;
extern TMVA::Reader *readerAgl;

extern TF1 * ResponseTOF ;
extern TF1 * ResponseNaF ;
extern TF1 * ResponseAgl ;

void UpdateProgressBar(int currentevent, int totalentries);
void SetBins();
void SetUpUsualBinning();
void SetUpEffCorrBinning();
void SetUpTOIBinning();
void SetUpRigTOIBinning();

//unfolding binning
extern Binning UnfoldingToF;
extern Binning UnfoldingNaF;
extern Binning UnfoldingAgl;

extern Binning UnfoldingToF_D;
extern Binning UnfoldingNaF_D;
extern Binning UnfoldingAgl_D;



int FindTimeIndex(std::string input);
#endif
