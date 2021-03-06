#include "Livetime.h"

float GetBeta( float R, float mass) {
	float beta = sqrt(R*R/(R*R + mass*mass));
	return beta;
}

void UpdateZoneLivetime (float Livetime, float Rcutoff, TH1F * esposizionegeo,Binning bins,float prescale){
		float SF;
		if(bins.IsUsingBetaEdges()) SF = BetacutoffCut;
		else SF = RcutoffCut; 

		for(int i=0;i<esposizionegeo->GetNbinsX();i++)
			if(SF*Rcutoff<=bins.RigBins()[i+1]){ 
				//per essere accettato un ev deve essere > di Rlowedge E di R cutoff, quindi vengono acc quelli tra Rcutoff e Rhighedge (è stima per eccesso - il sec viene contato cmq tutto)
				//al contrario è stima per difetto: il sec viene contato SOLO se è tutto osservato
				float partial = 1;	
				if(SF*Rcutoff>bins.RigBins()[i]) partial = (bins.RigBins()[i+1]-SF*Rcutoff)/(bins.RigBins()[i+1]-bins.RigBins()[i]);
				esposizionegeo -> SetBinContent(i+1, esposizionegeo -> GetBinContent(i+1) + Livetime*partial) ;
			}
		
		
	
	return;
}


