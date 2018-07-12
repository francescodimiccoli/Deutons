#include "Livetime.h"

void UpdateZoneLivetime (float Livetime, float Rcutoff, TH1F * esposizionegeo,Binning bins,float prescale){

        for(int i=0;i<esposizionegeo->GetNbinsX();i++)
                        if(bins.RigBins()[i+1]>=1.3*Rcutoff){
                                esposizionegeo -> SetBinContent(i+1, esposizionegeo -> GetBinContent(i+1) + Livetime*prescale) ;
        }
        return;
}


