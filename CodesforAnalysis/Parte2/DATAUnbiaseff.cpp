#include "PlottingFunctions/DATAUnbiaseff_Plot.h"

using namespace std;

Efficiency * EffUnbiasDATA = new Efficiency ("EffUnbiasDATA");

void DATAUnbiaseff_Fill () {
	if (!cmask.isPreselected() ||Tup.R_pre<=0||Tup.R_pre<1.2*Tup.Rcutoff) return;
	if (!(Tup.EdepTrack<EdepTrackbeta->Eval (Tup.Beta_pre)+0.2&&Tup.EdepTrack>EdepTrackbeta->Eval (Tup.Beta_pre)-0.2)) return;	

	int Kbin=PRB.GetRBin (fabs (Tup.R_pre) );
	if(Tup.Unbias==0) EffUnbiasDATA->beforeR->Fill(Kbin);
	if(Tup.Unbias==1) EffUnbiasDATA->beforeR->Fill(Kbin,100);
	if(Tup.Unbias==0) EffUnbiasDATA->afterR->Fill(Kbin);

	return;
}


void DATAUnbiaseff_Write() {
	EffUnbiasDATA -> Write();
	return;
}


void DATAUnbiaseff (string filename) {
	cout<<"********** DATA Tup.Unbias TRIGG. EFFICIENCY ******************************"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
	TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


	Efficiency * EffUnbiasDATA = new Efficiency (inputHistoFile,"EffUnbiasDATA");

	cout<<"********** DATA Tup.Unbias TRIGG. EFFICIENCY ******************************"<<endl;

	EffUnbiasDATA -> Eval_Efficiency();

	TH1F *EffUnbDATA_R_TH1F = (TH1F *) EffUnbiasDATA -> effR   ->Clone();


	float EffMean = 0;
	for(int n = 0; n< EffUnbDATA_R_TH1F->GetNbinsX();n++){
		EffMean+=EffUnbDATA_R_TH1F->GetBinContent(n+1);
	}
	EffMean=EffMean/(float)EffUnbDATA_R_TH1F->GetNbinsX();	

	TH1F * TriggerGlobalFactor = new TH1F("TriggerGlobalFactor","TriggerGlobalFactor",1,0,1);
	TriggerGlobalFactor -> SetBinContent(1,EffMean);
	TriggerGlobalFactor -> SetBinError(1,0.01); 

	finalHistos.Add(EffUnbDATA_R_TH1F  );
	finalHistos.Add(TriggerGlobalFactor);
	finalHistos.writeObjsInFolder("Results");

	cout<<"*** Plotting ...  ****"<<endl;

	DATAUnbiaseff_Plot(EffUnbDATA_R_TH1F , 
			TriggerGlobalFactor     

			);


	return;
}

