#include "PlottingFunctions/Hecut_Plot.h"


TH2F * Hecut_D=new TH2F("Hecut_D","Hecut_D",1000,0,40,1000,0,40);
TH2F * HecutMC_P=new TH2F("HecutMC_P","HecutMC_P",1000,0,40,1000,0,40);
TH2F * HecutMC_He=new TH2F("HecutMC_He","HecutMC_He",1000,0,40,1000,0,40);

TH2F * L1TOF_MC = new TH2F("L1TOF_MC","L1TOF_MC",100,0,1.5,nbinsToF,0,nbinsToF);
TH2F * L1TOF_DATA = new TH2F("L1TOF_DATA","L1TOF_DATA",100,0,1.5,nbinsNaF,0,nbinsNaF);
TH2F * L1TOF_DATAcutoff = new TH2F("L1TOF_DATAcutoff","L1TOF_DATAcutoff",100,0,1.5,nbinsAgl,0,nbinsAgl);

TH2F * L1NaF_MC = new TH2F("L1NaF_MC","L1NaF_MC",100,0,1.5,nbinsToF,0,nbinsToF);
TH2F * L1NaF_DATA = new TH2F("L1NaF_DATA","L1NaF_DATA",100,0,1.5,nbinsNaF,0,nbinsNaF);
TH2F * L1NaF_DATAcutoff = new TH2F("L1NaF_DATAcutoff","L1NaF_DATAcutoff",100,0,1.5,nbinsAgl,0,nbinsAgl);

TH2F * L1Agl_MC = new TH2F("L1Agl_MC","L1Agl_MC",100,0,1.5,nbinsToF,0,nbinsToF);
TH2F * L1Agl_DATA = new TH2F("L1Agl_DATA","L1Agl_DATA",100,0,1.5,nbinsNaF,0,nbinsNaF);
TH2F * L1Agl_DATAcutoff = new TH2F("L1Agl_DATAcutoff","L1Agl_DATAcutoff",100,0,1.5,nbinsAgl,0,nbinsAgl);


TH1F * L1TOFs_MC	 = new TH1F("L1TOFs_MC","L1TOFs_MC",1000,-20,40);
TH1F * L1TOFs_DATA 	 = new TH1F("L1TOFs_DATA","L1TOFs_DATA",1000,-20,40);
TH1F * L1TOFs_DATAcutoff = new TH1F("L1TOFs_DATAcutoff","L1TOFs_DATAcutoff",1000,-20,40);

TH1F * L1NaFs_MC 	 = new TH1F("L1NaFs_MC","L1NaFs_MC",1000,-20,40);
TH1F * L1NaFs_DATA 	 = new TH1F("L1NaFs_DATA","L1NaFs_DATA",1000,-20,40);
TH1F * L1NaFs_DATAcutoff = new TH1F("L1NaFs_DATAcutoff","L1NaFs_DATAcutoff",1000,-20,40);

TH1F * L1Agls_MC 	 = new TH1F("L1Agls_MC","L1Agls_MC",1000,-20,40);
TH1F * L1Agls_DATA 	 = new TH1F("L1Agls_DATA","L1Agls_DATA",1000,-20,40);
TH1F * L1Agls_DATAcutoff = new TH1F("L1Agls_DATAcutoff","L1Agls_DATAcutoff",1000,-20,40);


Efficiency * HecutMCP = new Efficiency("HecutMCP");
Efficiency * HecutMCHe = new Efficiency("HecutMCHe");

Efficiency * HeEff = new Efficiency("HeEff");
Efficiency * DataHeCont = new Efficiency("DataHeCont");

Efficiency * HeL1Cont = new Efficiency("HeL1Cont");

Efficiency * HeTRDCont_MC = new Efficiency("HeTRDCont_MC");
Efficiency * HeTRDCont_DATA = new Efficiency("HeTRDCont_DATA");


TH2F * L1RvsBetaTest = new TH2F("L1RvsBetaTest","L1RvsBetaTest",1000,0,20,100,0,1.3);


void HecutMC_Fill() {

	 if(!trgpatt.IsPhysical()) return;
         if(Tup.Beta<=0||Tup.R<=0) return;

	float EdepTOFud=(Tup.EdepTOFU+Tup.EdepTOFD)/2;

	int Kbin=0;
	float sigma=0;
	
	bool HeL1sample = false;
        if(Tup.EdepL1>0.015 && Massa_gen>2 ) HeL1sample = true;

	if(HeL1sample)
		{
			Kbin=ToFDB.GetBin(RUsed);

			if(cmask.isOnlyFromToF()){	
				L1TOF_MC ->Fill(Tup.EdepL1,Kbin);
				sigma = (pow(EdepL1beta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta));
				L1TOFs_MC ->Fill((Tup.EdepL1-EdepL1beta->Eval(Tup.Beta))/sigma );
				}
			Kbin=NaFDB.GetBin(RUsed);

			if(cmask.isFromNaF()){
				L1NaF_MC ->Fill(Tup.EdepL1,Kbin);
				sigma = (pow(EdepL1beta->Eval(Tup.BetaRICH),2)*etrack->Eval(Tup.BetaRICH));
				L1NaFs_MC ->Fill((Tup.EdepL1-EdepL1beta->Eval(Tup.BetaRICH))/sigma );
				}
			Kbin=AglDB.GetBin(RUsed);

			if(cmask.isFromAgl()){
                                L1Agl_MC ->Fill(Tup.EdepL1,Kbin);
				sigma = (pow(EdepL1beta->Eval(Tup.BetaRICH),2)*etrack->Eval(Tup.BetaRICH));
				L1Agls_MC ->Fill((Tup.EdepL1-EdepL1beta->Eval(Tup.BetaRICH))/sigma );
				}
		}

	Kbin=PRB.GetRBin(Tup.R);
	
	if(Massa_gen<1) {
		HecutMC_P->Fill( fabs(EdepTOFbeta  ->Eval(Tup.Beta)-EdepTOFud) / (pow(EdepTOFbeta  ->Eval(Tup.Beta),2)*etofu ->Eval(Tup.Beta)) ,
				fabs(EdepTrackbeta->Eval(Tup.Beta)-Tup.EdepTrack) / (pow(EdepTrackbeta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta)) );
		HecutMCP->beforeR->Fill(Kbin);
		if(Herejcut) HecutMCP->afterR->Fill(Kbin);

	}
	
	
	//total He->P 
	if(Massa_gen>2){	
			Kbin=ToFDB.GetBin(RUsed);
			if(Likcut&&!(Distcut)) HeEff->beforeTOF->Fill(Kbin);
			if(Likcut&&Distcut) HeEff->afterTOF->Fill(Kbin);	
			
			Kbin=NaFDB.GetBin(RUsed);
			if(Likcut&&!(Distcut)) HeEff->beforeNaF->Fill(Kbin);
			if(cmask.isFromNaF()) if(Likcut&&Distcut)    HeEff->afterNaF->Fill(Kbin);		
	
			Kbin=AglDB.GetBin(RUsed);
			if(Likcut&&!(Distcut)) HeEff->beforeAgl->Fill(Kbin);
			if(cmask.isFromAgl()) if(Likcut&&Distcut)    HeEff->afterAgl->Fill(Kbin);
	}
	
	if(HeL1sample) {	
		//He->P in TRD
		if(cmask.isOnlyFromToF()){ 
			Kbin=ToFDB.GetBin(RUsed);
			if(IsHeL1 ) {
				HeTRDCont_MC ->beforeTOF->Fill(Kbin);
				if(Distcut && Likcut) HeTRDCont_MC ->afterTOF->Fill(Kbin);
			}
		}
		if(cmask.isFromNaF()) {
			Kbin=NaFDB.GetBin(RUsed);
			if(IsHeL1 ) {
				HeTRDCont_MC ->beforeNaF->Fill(Kbin);
				if(Distcut && Likcut) HeTRDCont_MC ->afterNaF->Fill(Kbin);
			}
		}		

		if(cmask.isFromAgl()) {
			Kbin=AglDB.GetBin(RUsed);
			if(IsHeL1) {
				HeTRDCont_MC ->beforeAgl->Fill(Kbin);
				if(Distcut && Likcut) HeTRDCont_MC ->afterAgl->Fill(Kbin);
			}
		}


	}
}

void HecutD_Fill() {

	if(!trgpatt.IsPhysical()) return;
        if(Tup.Beta<=0||Tup.R<=0) return;
	if(!(Tup.R>1.3*Tup.Rcutoff)) return;

	float EdepTOFud=(Tup.EdepTOFU+Tup.EdepTOFD)/2;
	Hecut_D->Fill(fabs(EdepTOFbeta->Eval(Tup.Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Tup.Beta),2)*etofu->Eval(Tup.Beta)),fabs(EdepTrackbeta->Eval(Tup.Beta)-Tup.EdepTrack)/(pow(EdepTrackbeta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta)));
	

	int Kbin=0;
	float sigma=0;


	if(Tup.R<=1.3*Tup.Rcutoff) cout<<"ecco"<<endl;

	bool HeL1sample = false;
	if(Tup.EdepL1>0.015 && Tup.R>1.2*Tup.Rcutoff ) HeL1sample = true;



	if(HeL1sample)
                {
   		       Kbin=ToFDB.GetBin(RUsed);
 
                       if(cmask.isOnlyFromToF()) {                 
				L1TOF_DATA  ->Fill(Tup.EdepL1,Kbin);
				sigma = (pow(EdepL1beta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta));
                        	L1TOFs_DATA ->Fill((Tup.EdepL1-EdepL1beta->Eval(Tup.Beta))/sigma );
				if(Tup.Rcutoff > 1.3*RBeta->Eval(Tup.Beta) ){
					L1TOF_DATAcutoff->Fill(Tup.EdepL1,Kbin);
					L1TOFs_DATAcutoff ->Fill((Tup.EdepL1-EdepL1beta->Eval(Tup.Beta))/sigma );
					
					}
				}
			
			Kbin=NaFDB.GetBin(RUsed);
			
			if(cmask.isFromNaF()){
				L1NaF_DATA ->Fill(Tup.EdepL1,Kbin);
                        	sigma = (pow(EdepL1beta->Eval(Tup.BetaRICH),2)*etrack->Eval(Tup.BetaRICH));
				L1NaFs_DATA->Fill((Tup.EdepL1-EdepL1beta->Eval(Tup.BetaRICH))/sigma );
				if(Tup.Rcutoff > 1.3*RBeta->Eval(Tup.BetaRICH) ){
					L1NaF_DATAcutoff ->Fill(Tup.EdepL1,Kbin);
					L1NaFs_DATAcutoff->Fill((Tup.EdepL1-EdepL1beta->Eval(Tup.BetaRICH))/sigma );
					}
				}
			
			Kbin=AglDB.GetBin(RUsed);

			if(cmask.isFromAgl()){
				L1Agl_DATA ->Fill(Tup.EdepL1,Kbin);
				sigma = (pow(EdepL1beta->Eval(Tup.BetaRICH),2)*etrack->Eval(Tup.BetaRICH));
				L1Agls_DATA->Fill((Tup.EdepL1-EdepL1beta->Eval(Tup.BetaRICH))/sigma );
				if(Tup.Rcutoff > 1.3*RBeta->Eval(Tup.BetaRICH)){
					L1Agl_DATAcutoff->Fill(Tup.EdepL1,Kbin);
					L1Agls_DATAcutoff->Fill((Tup.EdepL1-EdepL1beta->Eval(Tup.BetaRICH))/sigma );
					}
				}
			if(cmask.isOnlyFromToF()) L1RvsBetaTest->Fill(Tup.R,Tup.Beta);
			if(cmask.isFromNaF()||cmask.isFromAgl()) L1RvsBetaTest->Fill(Tup.R,Tup.BetaRICH);			



                }



		Kbin=ToFDB.GetBin(RUsed);
		if(Likcut&&!(Distcut))   DataHeCont->beforeTOF->Fill(Kbin);
		if(Distcut && Likcut)    DataHeCont->afterTOF ->Fill(Kbin);	

		Kbin=NaFDB.GetBin(RUsed);
		if(Likcut&&!(Distcut)) DataHeCont->beforeNaF->Fill(Kbin);
		if(cmask.isFromNaF()) if(Distcut && Likcut)    DataHeCont->afterNaF ->Fill(Kbin);
	
		Kbin=AglDB.GetBin(RUsed);
		if(Likcut&&!(Distcut)) DataHeCont->beforeAgl->Fill(Kbin);
		if(cmask.isFromAgl()) if(Distcut && Likcut)    DataHeCont->afterAgl ->Fill(Kbin);


	//He->P in TRD
	if(HeL1sample) {
		if(cmask.isOnlyFromToF()){ 
			Kbin=ToFDB.GetBin(RUsed);
			if(Tup.Rcutoff > 1.3*RBeta->Eval(Tup.Beta) &&  IsHeL1 ) {
				HeTRDCont_DATA ->beforeTOF->Fill(Kbin);
				if(Distcut && Likcut) HeTRDCont_DATA ->afterTOF->Fill(Kbin);
			}
		}
		if(cmask.isFromNaF()) {
			Kbin=NaFDB.GetBin(RUsed);
			if(Tup.Rcutoff > 1.3*RBeta->Eval(Tup.BetaRICH) && IsHeL1 ) {
				HeTRDCont_DATA ->beforeNaF->Fill(Kbin);
				if(Distcut && Likcut) HeTRDCont_DATA ->afterNaF->Fill(Kbin);
			}
		}               

		if(cmask.isFromAgl()) {
			Kbin=AglDB.GetBin(RUsed);
			if(Tup.Rcutoff > 1.3*RBeta->Eval(Tup.BetaRICH) && IsHeL1 ) {
				HeTRDCont_DATA ->beforeAgl->Fill(Kbin);
				if(Distcut && Likcut) HeTRDCont_DATA ->afterAgl->Fill(Kbin);
			}
		}

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
	

	HeTRDCont_MC ->Write();
	HeTRDCont_DATA ->Write();


	L1TOF_MC->Write();
	L1TOF_DATA->Write();
	L1TOF_DATAcutoff->Write();
	L1NaF_MC->Write();
        L1NaF_DATA->Write();
        L1NaF_DATAcutoff->Write();
	L1Agl_MC->Write();
	L1Agl_DATA->Write();
        L1Agl_DATAcutoff->Write();


	L1TOFs_MC->Write();
	L1TOFs_DATA->Write();
	L1TOFs_DATAcutoff->Write();
	L1NaFs_MC->Write();
        L1NaFs_DATA->Write();
        L1NaFs_DATAcutoff->Write();
	L1Agls_MC->Write();
	L1Agls_DATA->Write();
        L1Agls_DATAcutoff->Write();
	
	L1RvsBetaTest->Write();


}


TH1F * Eval_Contamination(TH1F * Q1Data, TH1F * Q2Data, TH1F * HefragmL1, TH1F * HefragmTRD){
		TH1F * Contamination = (TH1F *) Q2Data -> Clone();
		Contamination -> Multiply ( HefragmL1 );
		/*for (int i = 0; i<Contamination ->GetNbinsX();i++){
			float Heliums = Contamination->GetBinContent(i+1)/(1 - HefragmL1->GetBinContent(i+1) - HefragmTRD->GetBinContent(i+1));
			Contamination->SetBinContent(i+1,Heliums*(HefragmL1->GetBinContent(i+1) + HefragmTRD->GetBinContent(i+1)));
		}*/
		Contamination -> Divide (Q1Data );
		return Contamination;
}


void Hecut(string filename) {
	cout<<"*************** He control sample cut Efficiency on P*******************"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
	TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


	TF1 * f2=new TF1("f2","ROOT::Math::crystalball_function(x, 2, 1, 1, 0)",-5,5);


	TH2F* HecutMC_P =(TH2F*)inputHistoFile->Get("HecutMC_P");
	TH2F* HecutMC_He=(TH2F*)inputHistoFile->Get("HecutMC_He");
	TH2F* Hecut_D   =(TH2F*)inputHistoFile->Get("Hecut_D");

	TH1F * L1TOFs_MC	 =(TH1F*)inputHistoFile->Get("L1TOFs_MC"	);	
	TH1F * L1TOFs_DATA 	 =(TH1F*)inputHistoFile->Get("L1TOFs_DATA" 	);
        TH1F * L1TOFs_DATAcutoff =(TH1F*)inputHistoFile->Get("L1TOFs_DATAcutoff");
                                                                               
        TH1F * L1NaFs_MC 	 =(TH1F*)inputHistoFile->Get("L1NaFs_MC" 	);	
        TH1F * L1NaFs_DATA 	 =(TH1F*)inputHistoFile->Get("L1NaFs_DATA" 	);
        TH1F * L1NaFs_DATAcutoff =(TH1F*)inputHistoFile->Get("L1NaFs_DATAcutoff");
                                                                               
        TH1F * L1Agls_MC 	 =(TH1F*)inputHistoFile->Get("L1Agls_MC" 	);
        TH1F * L1Agls_DATA 	 =(TH1F*)inputHistoFile->Get("L1Agls_DATA" 	);
	TH1F * L1Agls_DATAcutoff =(TH1F*)inputHistoFile->Get("L1Agls_DATAcutoff");

	TH1F * L1TOF_MC	 =(TH1F*)inputHistoFile->Get("L1TOF_MC"	);	
	TH1F * L1TOF_DATA 	 =(TH1F*)inputHistoFile->Get("L1TOF_DATA" 	);
        TH1F * L1TOF_DATAcutoff =(TH1F*)inputHistoFile->Get("L1TOF_DATAcutoff");
                                                                               
        TH1F * L1NaF_MC 	 =(TH1F*)inputHistoFile->Get("L1NaF_MC" 	);	
        TH1F * L1NaF_DATA 	 =(TH1F*)inputHistoFile->Get("L1NaF_DATA" 	);
        TH1F * L1NaF_DATAcutoff =(TH1F*)inputHistoFile->Get("L1NaF_DATAcutoff");
                                                                               
        TH1F * L1Agl_MC 	 =(TH1F*)inputHistoFile->Get("L1Agl_MC" 	);
        TH1F * L1Agl_DATA 	 =(TH1F*)inputHistoFile->Get("L1Agl_DATA" 	);
	TH1F * L1Agl_DATAcutoff =(TH1F*)inputHistoFile->Get("L1Agl_DATAcutoff");
	
	
	Efficiency * HecutMCP = new Efficiency(inputHistoFile,"HecutMCP");
	Efficiency * HecutMCHe = new Efficiency(inputHistoFile,"HecutMCHe");

	Efficiency * HeEff      = new Efficiency(inputHistoFile,"HeEff");
	Efficiency * DataHeCont = new Efficiency(inputHistoFile,"DataHeCont");

	Efficiency * HeTRDCont_MC   = new Efficiency(inputHistoFile,"HeTRDCont_MC");
	Efficiency * HeTRDCont_DATA = new Efficiency(inputHistoFile,"HeTRDCont_DATA");




	cout<<"*************** He control sample cut Efficiency on P*******************"<<endl;

	HecutMCP ->Eval_Efficiency();
	HecutMCHe->Eval_Efficiency();
	
	HeTRDCont_MC   ->Eval_Efficiency();
        HeTRDCont_DATA ->Eval_Efficiency();

	HeEff->Eval_Efficiency();
	DataHeCont->Eval_Efficiency();

	TH1F * HecutMCP_TH1F = 	(TH1F *)HecutMCP ->effR->Clone();
	TH1F * HecutMCHe_TH1F=  (TH1F *)HecutMCHe->effR->Clone();

	cout<<"*************** He Contamination ******************"<<endl;


	//L1 fragm eff from MC
	TH1F * efffragmL1_TOF = (TH1F *)HeEff->effTOF->Clone();
        TH1F * efffragmL1_NaF = (TH1F *)HeEff->effNaF->Clone();
	TH1F * efffragmL1_Agl = (TH1F *)HeEff->effAgl->Clone();

/*	efffragmL1_TOF -> Add ((TH1F *)HeTRDCont_MC -> effTOF->Clone(),-1);
	efffragmL1_NaF -> Add ((TH1F *)HeTRDCont_MC -> effNaF->Clone(),-1);
        efffragmL1_Agl -> Add ((TH1F *)HeTRDCont_MC -> effAgl->Clone(),-1);
*/
	TH1F * efffragmL1cut_TOF= new TH1F("efffragmL1cut_TOF","efffragmL1cut_TOF",nbinsToF,0,nbinsToF); 
        TH1F * efffragmL1cut_NaF= new TH1F("efffragmL1cut_NaF","efffragmL1cut_NaF",nbinsNaF,0,nbinsNaF);  
        TH1F * efffragmL1cut_Agl= new TH1F("efffragmL1cut_Agl","efffragmL1cut_Agl",nbinsAgl,0,nbinsAgl);  

	//TRD fragm. eff from data
	TH1F * efffragmTRD_TOF = (TH1F *)HeTRDCont_DATA   ->effTOF -> Clone();
	TH1F * efffragmTRD_NaF = (TH1F *)HeTRDCont_DATA   ->effNaF -> Clone();
	TH1F * efffragmTRD_Agl = (TH1F *)HeTRDCont_DATA   ->effAgl -> Clone();

	TH1F * ContaminationTOF = Eval_Contamination((TH1F *)DataHeCont->afterTOF,(TH1F *)DataHeCont->beforeTOF,efffragmL1_TOF ,efffragmTRD_TOF);
	TH1F * ContaminationNaF = Eval_Contamination((TH1F *)DataHeCont->afterNaF,(TH1F *)DataHeCont->beforeNaF,efffragmL1_NaF ,efffragmTRD_NaF);
	TH1F * ContaminationAgl = Eval_Contamination((TH1F *)DataHeCont->afterAgl,(TH1F *)DataHeCont->beforeAgl,efffragmL1_Agl ,efffragmTRD_Agl);

	TH1F * ContaminationTOF_cut = Eval_Contamination((TH1F *)DataHeCont->afterTOF,(TH1F *)DataHeCont->beforeTOF,efffragmL1cut_TOF ,efffragmTRD_TOF);
	TH1F * ContaminationNaF_cut = Eval_Contamination((TH1F *)DataHeCont->afterNaF,(TH1F *)DataHeCont->beforeNaF,efffragmL1cut_NaF ,efffragmTRD_NaF);
	TH1F * ContaminationAgl_cut = Eval_Contamination((TH1F *)DataHeCont->afterAgl,(TH1F *)DataHeCont->beforeAgl,efffragmL1cut_Agl ,efffragmTRD_Agl);

	ContaminationTOF -> SetName("ContaminationTOF");
        ContaminationNaF -> SetName("ContaminationNaF");
        ContaminationAgl -> SetName("ContaminationAgl");

	ContaminationTOF_cut -> SetName("ContaminationTOF_cut");
        ContaminationNaF_cut -> SetName("ContaminationNaF_cut");
        ContaminationAgl_cut -> SetName("ContaminationAgl_cut");
	
	finalHistos.Add(ContaminationTOF);
	finalHistos.Add(ContaminationNaF);
	finalHistos.Add(ContaminationAgl);
   	finalHistos.Add(ContaminationTOF_cut);
	finalHistos.Add(ContaminationNaF_cut);
	finalHistos.Add(ContaminationAgl_cut);
   	
	finalHistos.writeObjsInFolder("Results");

	cout<<"*** Plotting ...  ****"<<endl;

	Hecut_Plot(

	L1TOF_DATA 	, 
        L1TOF_DATAcutoff, 
        L1NaF_DATA 	 ,
        L1NaF_DATAcutoff ,
        L1Agl_DATA 	 ,
        L1Agl_DATAcutoff ,

	L1TOFs_DATA 	, 
        L1TOFs_DATAcutoff, 
        L1NaFs_DATA 	 ,
        L1NaFs_DATAcutoff ,
        L1Agls_DATA 	 ,
        L1Agls_DATAcutoff ,
        

	HecutMC_He        ,
        Hecut_D    ,
        HecutMCP_TH1F,
        HecutMCHe_TH1F,
        HecutMC_P,
	
	(TH1F *)HeTRDCont_MC   ->effTOF,
        (TH1F *)HeTRDCont_MC   ->effNaF,
	(TH1F *)HeTRDCont_MC   ->effAgl,
	(TH1F *)HeTRDCont_DATA ->effTOF,	
	(TH1F *)HeTRDCont_DATA ->effNaF,
        (TH1F *)HeTRDCont_DATA ->effAgl,
	
	efffragmL1_TOF ,
	efffragmL1_NaF ,
	efffragmL1_Agl ,
	ContaminationTOF,	
        ContaminationNaF,
	ContaminationAgl);


	return; 
}
