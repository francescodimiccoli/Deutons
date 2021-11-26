#ifndef MACROFORPLOTS_H
#define MACROFORPLOTS_H

#include "TText.h"
#include "TLine.h"
#include "binning.h"
#include "Globals.h"
#include "reweight.h"
#include "Variables.hpp"
#include "filesaver.h"
#include "TCanvas.h"
#include "PlottingFunctions.h"

float GammaFromBeta (float beta);
float EkFromBeta (float beta, float mass) ;


void DrawBins(TVirtualPad * c1, TH2F* histo, std::vector<float> orizontal_set,std::vector<float> vertical_set,int col_oriz, int col_vert);

TH2F * RebinHisto(TH2F * histo,std::vector<float> orizontal_set,std::vector<float> vertical_set,std::string name);

void DrawQvsRBeta(TH2F * h6, TH2F *h9, TH2F * h6_, TH2F * h9_, std::string name, FileSaver finalResults);


TH1D * ReduceOnTimeHisto(TH1D *OnTime);


void DrawCleaning(FileSaver finalHistos,FileSaver finalResults);


void DrawBDT(FileSaver finalHistos,FileSaver finalResults);

void DrawMasses(FileSaver finalHistos,FileSaver finalResults);

TH1D * GetMeans(TH2F * h1);

TH1D * GetSTD(TH2F * h1);

void DrawMassRes(FileSaver finalHistos,FileSaver finalResults);

void DrawBetaSmear(FileSaver finalHistos,FileSaver finalResults);


void DrawBetaRes(FileSaver finalHistos,FileSaver finalResults);


void DrawBetavsRig(FileSaver finalHistos,FileSaver finalResults);



void DrawAcceptanceMatrix(FileSaver finalHistos,FileSaver finalResults);

TH1D * EvalRunningIntegral(TH1D * Distrib);

void DrawVariables_app(FileSaver finalHistos,FileSaver finalResults, std::string app);
 
void DrawVariables(FileSaver finalHistos,FileSaver finalResults);

#endif
