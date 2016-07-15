#include "PlottingFunctions/DATAUnbiaseff_Plot.h"

using namespace std;

Efficiency * EffUnbiasDATA  = new Efficiency ("EffUnbiasDATA" );
Efficiency * EffUnbiasDATAQ = new Efficiency ("EffUnbiasDATAQ");


void DATAUnbiaseff_Fill () {
	if (!cmask.isPreselected() ||Tup.R_pre<=0||Tup.Beta_pre<=0||Tup.R_pre<1.2*Tup.Rcutoff) return;
	//if (!(Tup.EdepTrack<EdepTrackbeta->Eval (Tup.Beta_pre)+0.2&&Tup.EdepTrack>EdepTrackbeta->Eval (Tup.Beta_pre)-0.2)) return;	
	if(!Herejcut) return;

	int Kbin=PRB.GetRBin (fabs (Tup.R_pre) );
	if(Tup.Unbias==0) EffUnbiasDATA->beforeR->Fill(Kbin);
	if(Tup.Unbias==1) EffUnbiasDATA->beforeR->Fill(Kbin,100);
	if(Tup.Unbias==0) EffUnbiasDATA->afterR->Fill(Kbin);

	return;
}

void DATAUnbiaseffQ_Fill () {
	if(Tup.Beta<=0||Tup.R<=0||Tup.Beta<=0||Tup.R<1.2*Tup.Rcutoff) return;
	if(!Herejcut) return;
	if(!(Tup.Dist5D_P<6&&Likcut)) return;

	int Kbin=PRB.GetRBin (fabs (Tup.R) );
        if(Tup.Unbias==0) EffUnbiasDATAQ->beforeR->Fill(Kbin);
        if(Tup.Unbias==1) EffUnbiasDATAQ->beforeR->Fill(Kbin,100);
        if(Tup.Unbias==0) EffUnbiasDATAQ->afterR->Fill(Kbin);

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
	TriggerGlobalFactor -> SetBinError(1,0.01); 

	TriggerGlobalFactorQ -> SetBinContent(1,EffMeanQ);
        TriggerGlobalFactorQ -> SetBinError(1,0.01);

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

