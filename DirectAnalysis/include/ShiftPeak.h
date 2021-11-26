#ifndef SHIFTPEAK_H
#define SHIFTPEAK_H

#include <TSpline.h>
#include <TH1.h>
#include <TF1.h>

TSpline3 * Model_Histo(TH1F * Histo);

TF1 * WeightGausswithSpline(TSpline3 * Spline, float mean,float mu, float sigma);

TH1F * InfinitesimalTermForConvolution(TH1F * Histo, TSpline3 * Spline, float mean, float mu, float sigma,int steps);

void AddInfinitesimalTermForConvolution(TH1F* Term, TH1F * Histo, TSpline3 * Spline, float mean, float mu, float sigma,int steps);


TH1F * ConvolveWithGaus(TH1F * Histo, float mu, float sigma,int steps=500);

TH1F * SimpleShiftHisto(TH1F * Histo, float mu,int steps=500);

#endif
