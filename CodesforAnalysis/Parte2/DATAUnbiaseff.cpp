#include "PlottingFunctions/DATAUnbiaseff_Plot.h"

using namespace std;

Efficiency * EffUnbiasDATA  = new Efficiency ("EffUnbiasDATA" );
Efficiency * EffUnbiasDATAQ = new Efficiency ("EffUnbiasDATAQ");

void DATAUnbiaseffQ_Fill () {
	if(Tup.Beta<=0||Tup.R<=0||Tup.Beta<=0||Tup.R<1.2*Tup.Rcutoff) return;
	if(!Herejcut) return;
	if(!ProtonsMassWindow) return;
	
	int Kbin=PRB.GetRBin (fabs (Tup.R) );
        if(trgpatt.IsPhysical()) EffUnbiasDATA->beforeR->Fill(Kbin);
        if(trgpatt.IsUnbias()  ) EffUnbiasDATA->beforeR->Fill(Kbin,100);
        if(trgpatt.IsPhysical()) EffUnbiasDATA->afterR->Fill(Kbin);


	if(!(Distcut&&Likcut)) return;
        
        if(trgpatt.IsPhysical()) EffUnbiasDATAQ->beforeR->Fill(Kbin);
        if(trgpatt.IsUnbias()  ) EffUnbiasDATAQ->beforeR->Fill(Kbin,100);
        if(trgpatt.IsPhysical()) EffUnbiasDATAQ->afterR->Fill(Kbin);

        return;
}

void DATAUnbiaseff_Write() {
	EffUnbiasDATA -> Write();
	return;
}

void DATAUnbiaseffQ_Write() {
        EffUnbiasDATAQ -> Write();
        return;
}



void DATAUnbiaseff (string filename) {
	cout<<"********** DATA Tup.Unbias TRIGG. EFFICIENCY ******************************"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
	TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


	Efficiency * EffUnbiasDATA  = new Efficiency (inputHistoFile,"EffUnbiasDATA");
	Efficiency * EffUnbiasDATAQ = new Efficiency (inputHistoFile,"EffUnbiasDATAQ");	

	cout<<"********** DATA Tup.Unbias TRIGG. EFFICIENCY ******************************"<<endl;

	EffUnbiasDATA  -> Eval_Efficiency();
	EffUnbiasDATAQ -> Eval_Efficiency();

	TH1F *EffUnbDATA_R_TH1F = (TH1F *) EffUnbiasDATA -> effR   ->Clone();
	TH1F *EffUnbDATAQ_R_TH1F = (TH1F *) EffUnbiasDATAQ -> effR   ->Clone();

	float EffMean = 0;
	int bin=0;
	for(int n = 0; n< EffUnbDATA_R_TH1F->GetNbinsX();n++){
		if(EffUnbDATA_R_TH1F->GetBinContent(n)>0){
			EffMean+=EffUnbDATA_R_TH1F->GetBinContent(n+1);
			bin++;
			}	
	}
	EffMean=EffMean/(float)bin;	

	float EffMeanQ = 0;
        bin=0;
	for(int n = 0; n< EffUnbDATAQ_R_TH1F->GetNbinsX();n++){
                if(EffUnbDATAQ_R_TH1F->GetBinContent(n)>0){
		EffMeanQ+=EffUnbDATAQ_R_TH1F->GetBinContent(n+1);
        	bin++;
		}
	}
        EffMeanQ=EffMeanQ/(float)bin;
		
	

	TH1F * TriggerGlobalFactor = new TH1F("TriggerGlobalFactor","TriggerGlobalFactor",1,0,1);
	TH1F * TriggerGlobalFactorQ = new TH1F("TriggerGlobalFactorQ","TriggerGlobalFactorQ",1,0,1);

	TriggerGlobalFactor -> SetBinContent(1,EffMean);
	TriggerGlobalFactor -> SetBinError(1,0.001); 

	TriggerGlobalFactorQ -> SetBinContent(1,EffMeanQ);
        TriggerGlobalFactorQ -> SetBinError(1,0.001);

	finalHistos.Add(EffUnbDATA_R_TH1F  );
	finalHistos.Add(TriggerGlobalFactor);
	finalHistos.Add(TriggerGlobalFactorQ);
	finalHistos.writeObjsInFolder("Results");

	cout<<"*** Plotting ...  ****"<<endl;

	DATAUnbiaseff_Plot(EffUnbDATA_R_TH1F , 
			   EffUnbDATAQ_R_TH1F ,
			   TriggerGlobalFactor,
			   TriggerGlobalFactorQ	     
			);


	return;
}

