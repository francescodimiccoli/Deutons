#ifndef UNFOLDING_H
#define UNFOLDING_H

#include <TSpline.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TRandom3.h>
//#include "Globals.h"


using namespace std;

int slicenormalizex(TH2F * h);
double splineIntegral(TSpline3* spline,float xmin,float xmax);

struct UnfoldRes {
	TSpline3 * spline;
	TF1 * funct;
};

class Folder {

        private:
        TH1F * Counts;
        TH2F * Migr_matr;
        TSpline3 * UnfoldedCounts;
	TGraph * Expo;
        int knots;
	TRandom3 * rand;
	float fit_min;
	float fit_max;
	float offset;

        public:
        Folder (TH2 * Migr, TH1* counts,TGraph *expo, float fitmin, float fitmax,int Knots,float Offset) {
                Counts = (TH1F*) counts;
                Migr_matr = (TH2F*) Migr;
              Expo = expo;
		  knots = Knots;
		offset=Offset;
		fit_min=fitmin;
		fit_max=fitmax;
                double x[knots];
                double y[knots];
                for(int i=0;i<knots;i++) {x[i]=0; x[i] = GetXknot(i);  }
                for(int i=0;i<knots;i++) {y[i]=0; y[i] = counts->GetBinContent(2*i+1); }
		rand = new TRandom3(time(0));
                UnfoldedCounts = new TSpline3("UnfoldedCounts",x,y,knots);

        }

        TSpline3* GetSpline() {return UnfoldedCounts;}

        double foldwithmatrix(TSpline3* flux, float x);
	TRandom3* GetRandomGen(){return rand;}
        double Fold(double *x, double *p);
        int GetNknots(){return knots;}

        float GetXknot(int indx);
        void SetKnots(double *y,std::string name);
};

UnfoldRes Unfold(TH2* migr_matr_norm, TH1 * measured_R, TGraph *expo, float fit_min, float fit_max,int Knots, float Offset);

#endif
