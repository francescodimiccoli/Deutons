#ifndef PLOTTING_H
#define PLOTTING_H

#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "../include/binning.h"
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include <TFile.h>
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "../include/Globals.h"

using namespace std;

void SetCanvas(TCanvas *c1);

TH1F * CreateHisto(std::string name,  Binning Bins,bool IsEkin=false);
TH1D * CreateHistoD(std::string name,  Binning Bins,bool IsEkin=false);



TH1F * ConvertBinnedHisto(TH1F * histo,std::string name,  Binning Bins,bool IsEkin=false);
TH1D * ConvertBinnedHisto(TH1D * histo,std::string name,  Binning Bins,bool IsEkin=false);


TH2F* CreateFrame (TVirtualPad * c,float xmin,float xmax,float ymin, float ymax,std::string Xaxis,std::string Yaxis);


void PlotDistribution(TVirtualPad * c, TH1F * Distribution, std::string Xaxis, std::string Yaxis, int color, std::string options, float ymin=-1,float ymax=-1,float thickline=3,std::string legendname="",bool filled = false,bool dots= false,bool skipleg=false,int rebin=1);

void PlotTH1F(TVirtualPad * c, TH1F * Distribution, std::string Xaxis, std::string Yaxis, int color, std::string options,std::string legendname="", float thickness=3, bool skipleg=false);

void PlotTH2F(TVirtualPad * c, TH2F * Distribution, std::string Xaxis, std::string Yaxis, std::string options);

void PlotGraph(TVirtualPad * c,TGraphErrors * graph,std::string Xaxis, std::string Yaxis, int color,std::string options, float xmin=-1,float xmax=-1,float ymin=-1,float ymax=-1,std::string legendname="",int dotstyle=8,bool skipleg=false) ;


void PlotFunction(TVirtualPad * c, TF1 * Function, std::string Xaxis, std::string Yaxis, int color, std::string options, float xmin=-1,float xmax=-1,float ymin=-1,float ymax=-1,std::string legendname="");

void PlotFunction(TVirtualPad * c, TSpline3 * Function, std::string Xaxis, std::string Yaxis, int color, std::string options, float xmin=-1,float xmax=-1,float ymin=-1,float ymax=-1,std::string legendname="");

void PlotTH1FintoGraph(TVirtualPad * c, Binning bins, TH1F * Values, std::string Xaxis, std::string Yaxis, int color,bool Ekin=false, std::string options="same", float xmin=-1,float xmax=-1,float ymin=-1,float ymax=-1,std::string legendname="",int dotstyle=8,bool skipleg=false, bool cleanhigherrors=false, TH1F* ErrUp=0x0, TH1F* ErrDw =0x0);


void PlotMergedRanges(TVirtualPad * c, TH1F * ValuesTOF, TH1F* ValuesNaF, TH1F* ValuesAgl, std::string Xaxis, std::string Yaxis, int color,bool Ekin=false, std::string options="same", float xmin=-1,float xmax=-1,float ymin=-1,float ymax=-1,std::string legendname="",int dotstyle=8,bool skipleg=false, bool cleanhigherrors=false);


void PlotTH1FRatiointoGraph(TVirtualPad * c, Binning bins, TH1F * Values1, TH1F * Values2, std::string Xaxis, std::string Yaxis, int color,bool Ekin=false, std::string options="same", float xmin=-1,float xmax=-1,float ymin=-1,float ymax=-1,std::string legendname1="",std::string legendname2="");


void PlotRatioWithSplineintoGraph(TVirtualPad * c, Binning bins, TH1F * Values1,TSpline3 * Spline, std::string Xaxis, std::string Yaxis, int color,bool Ekin=false, std::string options="same", float xmin=-1,float xmax=-1,float ymin=-1,float ymax=-1,std::string legendname1="",int dotstyle=8);

void PlotRatioWithSplineiAvg(TVirtualPad * c, Binning bins, TH1F * Values1,TSpline3 * Spline, std::string Xaxis, std::string Yaxis, int color,bool drawError = false, bool Ekin=false, std::string options="same", float xmin=-1,float xmax=-1,float ymin=-1,float ymax=-1,std::string legendname1="",int dotstyle=8);


#endif
