#include "PlottingFunctions/Hecut_Plot.h"


TH2F * Hecut_D=new TH2F("Hecut_D","Hecut_D",1000,0,40,1000,0,40);
TH2F * HecutMC_P=new TH2F("HecutMC_P","HecutMC_P",1000,0,40,1000,0,40);
TH2F * HecutMC_He=new TH2F("HecutMC_He","HecutMC_He",1000,0,40,1000,0,40);

Efficiency * HecutMCP = new Efficiency("HecutMCP");
Efficiency * HecutMCHe = new Efficiency("HecutMCHe");

Efficiency * HeEff = new Efficiency("HeEff");

Efficiency * DataHeCont = new Efficiency("DataHeCont");



void HecutMC_Fill() {

	 if(!trgpatt.IsPhysical()) return;
         if(Tup.Beta<=0||Tup.R<=0) return;


	float EdepTOFud=(Tup.EdepTOFU+Tup.EdepTOFD)/2;
	int Kbin=PRB.GetRBin(Tup.R);
	
	if(Massa_gen<1) {
		HecutMC_P->Fill( fabs(EdepTOFbeta  ->Eval(Tup.Beta)-EdepTOFud) / (pow(EdepTOFbeta  ->Eval(Tup.Beta),2)*etofu ->Eval(Tup.Beta)) ,
				fabs(EdepTrackbeta->Eval(Tup.Beta)-Tup.EdepTrack) / (pow(EdepTrackbeta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta)) );
		HecutMCP->beforeR->Fill(Kbin);
		if(Herejcut) HecutMCP->afterR->Fill(Kbin);

	}
	if(Massa_gen>2) {
		Kbin=ToFDB.GetBin(RUsed);
		if(Likcut&&!(Distcut)) HeEff->beforeTOF->Fill(Kbin);
		if(Distcut && Likcut)    HeEff->afterTOF->Fill(Kbin);	

		if(cmask.isFromNaF()) {
			Kbin=NaFDB.GetBin(RUsed);
			if(Likcut&&!(Distcut)) HeEff->beforeNaF->Fill(Kbin);
			if(Distcut && Likcut)    HeEff->afterNaF->Fill(Kbin);		
		}
		if(cmask.isFromAgl()) {
			Kbin=AglDB.GetBin(RUsed);
			if(Likcut&&!(Distcut))   HeEff->beforeAgl->Fill(Kbin);
			if(Distcut && Likcut)    HeEff->afterAgl->Fill(Kbin);
		}
	}
}

void HecutD_Fill() {

	if(!trgpatt.IsPhysical()) return;
        if(Tup.Beta<=0||Tup.R<=0) return;
	if(!Tup.R>1.2*Tup.Rcutoff) return;

	float EdepTOFud=(Tup.EdepTOFU+Tup.EdepTOFD)/2;
	Hecut_D->Fill(fabs(EdepTOFbeta->Eval(Tup.Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Tup.Beta),2)*etofu->Eval(Tup.Beta)),fabs(EdepTrackbeta->Eval(Tup.Beta)-Tup.EdepTrack)/(pow(EdepTrackbeta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta)));
	int Kbin=0;
	
	Kbin=ToFDB.GetBin(RUsed);
	if(Likcut&&!(Distcut)) DataHeCont->beforeTOF->Fill(Kbin);
	if(Distcut && Likcut)    DataHeCont->afterTOF ->Fill(Kbin);	

	if(cmask.isFromNaF()) {	
		Kbin=NaFDB.GetBin(RUsed);
		if(Likcut&&!(Distcut)) DataHeCont->beforeNaF->Fill(Kbin);
		if(Distcut && Likcut)    DataHeCont->afterNaF ->Fill(Kbin);
	}
	if(cmask.isFromAgl()) {
		Kbin=AglDB.GetBin(RUsed);
		if(Likcut&&!(Distcut)) DataHeCont->beforeAgl->Fill(Kbin);
		if(Distcut && Likcut)    DataHeCont->afterAgl ->Fill(Kbin);
	}

}


void HecutMC_Write() {
	HecutMC_P->Write();
	HecutMC_He->Write();
	Hecut_D->Write();
	HecutMCP->Write();
	HecutMCHe->Write();
	HeEff -> Write();
	DataHeCont-> Write();	
}


TH1F * Eval_Contamination(TH1F * Q1Data, TH1F * Q2Data, TH1F * Hefragmeff){
		TH1F * Contamination = (TH1F *) Q2Data -> Clone();
		Contamination -> Multiply (Hefragmeff);
		Contamination -> Divide (Q1Data );
		return Contamination;
}


void Hecut(string filename) {
	cout<<"*************** He control sample cut Efficiency on P*******************"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
	TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


	TH2F* HecutMC_P =(TH2F*)inputHistoFile->Get("HecutMC_P");
	TH2F* HecutMC_He=(TH2F*)inputHistoFile->Get("HecutMC_He");
	TH2F* Hecut_D   =(TH2F*)inputHistoFile->Get("Hecut_D");
	
	Efficiency * HecutMCP = new Efficiency(inputHistoFile,"HecutMCP");
	Efficiency * HecutMCHe = new Efficiency(inputHistoFile,"HecutMCHe");

	Efficiency * HeEff      = new Efficiency(inputHistoFile,"HeEff");
	Efficiency * DataHeCont = new Efficiency(inputHistoFile,"DataHeCont");



	cout<<"*************** He control sample cut Efficiency on P*******************"<<endl;

	HecutMCP->Eval_Efficiency();
	HecutMCHe->Eval_Efficiency();
	
	HeEff->Eval_Efficiency();
	DataHeCont->Eval_Efficiency();

	TH1F * HecutMCP_TH1F = 	(TH1F *)HecutMCP ->effR->Clone();
	TH1F * HecutMCHe_TH1F=  (TH1F *)HecutMCHe->effR->Clone();

	cout<<"*************** He Contamination ******************"<<endl;

	TH1F * ContaminationTOF = Eval_Contamination((TH1F *)DataHeCont->afterTOF,(TH1F *)DataHeCont->beforeTOF,(TH1F *)HeEff->effTOF);
	TH1F * ContaminationNaF = Eval_Contamination((TH1F *)DataHeCont->afterNaF,(TH1F *)DataHeCont->beforeNaF,(TH1F *)HeEff->effNaF);
	TH1F * ContaminationAgl = Eval_Contamination((TH1F *)DataHeCont->afterAgl,(TH1F *)DataHeCont->beforeAgl,(TH1F *)HeEff->effAgl);

	ContaminationTOF -> SetName("ContaminationTOF");
        ContaminationNaF -> SetName("ContaminationNaF");
        ContaminationAgl -> SetName("ContaminationAgl");
	
	finalHistos.Add(ContaminationTOF);
	finalHistos.Add(ContaminationNaF);
	finalHistos.Add(ContaminationAgl);
   	finalHistos.writeObjsInFolder("Results");

	cout<<"*** Plotting ...  ****"<<endl;

	Hecut_Plot(

        HecutMC_He        ,
        Hecut_D    ,
        HecutMCP_TH1F,
        HecutMCHe_TH1F,
        HecutMC_P,
	(TH1F *)HeEff->effTOF,
	(TH1F *)HeEff->effNaF,
	(TH1F *)HeEff->effAgl,
	ContaminationTOF,	
        ContaminationNaF,
	ContaminationAgl);


	return; 
}
