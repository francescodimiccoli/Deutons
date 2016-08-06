#include "PlottingFunctions/Hecut_Plot.h"


TH2F * Hecut_D=new TH2F("Hecut_D","Hecut_D",1000,0,40,1000,0,40);
TH2F * HecutMC_P=new TH2F("HecutMC_P","HecutMC_P",1000,0,40,1000,0,40);
TH2F * HecutMC_He=new TH2F("HecutMC_He","HecutMC_He",1000,0,40,1000,0,40);

Efficiency * HecutMCP = new Efficiency("HecutMCP");
Efficiency * HecutMCHe = new Efficiency("HecutMCHe");

TemplateFIT * HeliumContaminationTOF = new TemplateFIT("HeliumContaminationTOF",1,0,60);
TemplateFIT * HeliumContaminationNaF = new TemplateFIT("HeliumContaminationNaF",1,0,60);
TemplateFIT * HeliumContaminationAgl = new TemplateFIT("HeliumContaminationAgl",1,0,60);



void HecutMC_Fill() {

	 if(!trgpatt.IsPhysical()) return;
         if(Tup.Beta<=0||Tup.R<=0) return;


	float EdepTOFud=(Tup.EdepTOFU+Tup.EdepTOFD)/2;
	int Kbin=PRB.GetRBin(Tup.R);
	float dummy=0;
	if(Massa_gen<1) {
		HecutMC_P->Fill( fabs(EdepTOFbeta  ->Eval(Tup.Beta)-EdepTOFud) / (pow(EdepTOFbeta  ->Eval(Tup.Beta),2)*etofu ->Eval(Tup.Beta)) ,
				fabs(EdepTrackbeta->Eval(Tup.Beta)-Tup.EdepTrack) / (pow(EdepTrackbeta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta)) );
		HecutMCP->beforeR->Fill(Kbin);
		if(Herejcut) HecutMCP->afterR->Fill(Kbin);

		if(Betastrongcut&&Likcut){		
			if(cmask.isOnlyFromToF()&&Tup.R<3)  ((TH2*)HeliumContaminationTOF -> TemplateP)->Fill(Tup.Dist5D_P,dummy,Tup.mcweight);
			if(cmask.isFromNaF()&&Tup.R<6)		((TH2*)  HeliumContaminationNaF -> TemplateP) -> Fill(Tup.Dist5D_P,dummy,Tup.mcweight);	
			if(cmask.isFromAgl()&&Tup.R<14)  	((TH2*)  HeliumContaminationAgl -> TemplateP )-> Fill(Tup.Dist5D_P,dummy,Tup.mcweight);
		}	
	}
	if(Massa_gen>1&&Massa_gen<2){
		if(Betastrongcut&&Likcut){
			if(cmask.isOnlyFromToF()&&Tup.R<3)   ((TH2*) HeliumContaminationTOF -> TemplateD) -> Fill(Tup.Dist5D_P,dummy,Tup.mcweight);
			if(cmask.isFromNaF()&&Tup.R<6)		((TH2*)  HeliumContaminationNaF -> TemplateD )-> Fill(Tup.Dist5D_P,dummy,Tup.mcweight);
			if(cmask.isFromAgl()&&Tup.R<14) 	((TH2*)	  HeliumContaminationAgl -> TemplateD )-> Fill(Tup.Dist5D_P,dummy,Tup.mcweight);
		}
		
	}
	if(Massa_gen>2) {
		HecutMC_He->Fill(fabs(EdepTOFbeta  ->Eval(Tup.Beta)-EdepTOFud)/(pow(EdepTOFbeta  ->Eval(Tup.Beta),2)*etofu ->Eval(Tup.Beta)),
				fabs(EdepTrackbeta->Eval(Tup.Beta)-Tup.EdepTrack)/(pow(EdepTrackbeta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta)));
		HecutMCHe->beforeR->Fill(Kbin);
		if(Herejcut) HecutMCHe->afterR->Fill(Kbin);

		if(Betastrongcut&&Likcut){
			if(cmask.isOnlyFromToF()&&Tup.R<3)   ((TH2*) HeliumContaminationTOF -> TemplateHe) -> Fill(Tup.Dist5D_P,dummy,Tup.mcweight);
			if(cmask.isFromNaF()&&Tup.R<6)		((TH2*)  HeliumContaminationNaF -> TemplateHe )-> Fill(Tup.Dist5D_P,dummy,Tup.mcweight);
			if(cmask.isFromAgl()&&Tup.R<14) 	((TH2*)	  HeliumContaminationAgl -> TemplateHe) -> Fill(Tup.Dist5D_P,dummy,Tup.mcweight);
		}
	}
}

void HecutD_Fill() {

	if(!trgpatt.IsPhysical()) return;
        if(Tup.Beta<=0||Tup.R<=0) return;
	if(!Tup.R>1.2*Tup.Rcutoff) return;

	float EdepTOFud=(Tup.EdepTOFU+Tup.EdepTOFD)/2;
	Hecut_D->Fill(fabs(EdepTOFbeta->Eval(Tup.Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Tup.Beta),2)*etofu->Eval(Tup.Beta)),fabs(EdepTrackbeta->Eval(Tup.Beta)-Tup.EdepTrack)/(pow(EdepTrackbeta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta)));
	float dummy=0;
	if(Betastrongcut&&Likcut){
		if(cmask.isOnlyFromToF()&&Tup.R<3)  HeliumContaminationTOF -> DATA -> Fill(Tup.Dist5D_P,dummy);
		if(cmask.isFromNaF()&&Tup.R<6)		HeliumContaminationNaF -> DATA -> Fill(Tup.Dist5D_P,dummy);
		if(cmask.isFromAgl()&&Tup.R<14)  	HeliumContaminationAgl -> DATA -> Fill(Tup.Dist5D_P,dummy);	
	}
}


void HecutMC_Write() {
	HecutMC_P->Write();
	HecutMC_He->Write();
	Hecut_D->Write();
	HecutMCP->Write();
	HecutMCHe->Write();
	HeliumContaminationTOF -> Write();
	HeliumContaminationNaF -> Write();	
	HeliumContaminationAgl -> Write();
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

	TemplateFIT * HeliumContaminationTOF = new TemplateFIT(inputHistoFile,"HeliumContaminationTOF","HeliumContaminationTOF");
	TemplateFIT * HeliumContaminationNaF = new TemplateFIT(inputHistoFile,"HeliumContaminationNaF","HeliumContaminationNaF");
	TemplateFIT * HeliumContaminationAgl = new TemplateFIT(inputHistoFile,"HeliumContaminationAgl","HeliumContaminationAgl");

	cout<<"*************** He control sample cut Efficiency on P*******************"<<endl;

	HecutMCP->Eval_Efficiency();
	HecutMCHe->Eval_Efficiency();

	TH1F * HecutMCP_TH1F = 	(TH1F *)HecutMCP ->effR->Clone();
	TH1F * HecutMCHe_TH1F=  (TH1F *)HecutMCHe->effR->Clone();

	cout<<"*************** He Contamination ******************"<<endl;

	 HeliumContaminationTOF -> SetFitConstraints(0.0,1,0.0,1,0.00,1);
         HeliumContaminationNaF -> SetFitConstraints(0.0,1,0.0,1,0.00,1);
         HeliumContaminationAgl -> SetFitConstraints(0.0,1,0.0,1,0.00,1);
	
	HeliumContaminationTOF-> TemplateFits();
	HeliumContaminationNaF-> TemplateFits();
	HeliumContaminationAgl-> TemplateFits();

	float HeCountsTOF= HeliumContaminationTOF->GetResult_He(0)->Integral(0,20);//HeliumContaminationTOF->GetResult_He(0)->FindBin(4));
	float PCountsTOF = HeliumContaminationTOF->GetResult_P(0)->Integral(0,20);//HeliumContaminationTOF->GetResult_P(0)->FindBin(4)); 
	float DCountsTOF = HeliumContaminationTOF->GetResult_D(0)->Integral(0,20);
	float HeCont_TOF = HeCountsTOF/(PCountsTOF+DCountsTOF);

	float HeCountsNaF= HeliumContaminationNaF->GetResult_He(0)->Integral(0,20);//HeliumContaminationNaF->GetResult_He(0)->FindBin(4));
        float PCountsNaF = HeliumContaminationNaF->GetResult_P(0)->Integral(0,20);//HeliumContaminationNaF->GetResult_P(0)->FindBin(4));
        float DCountsNaF = HeliumContaminationNaF->GetResult_D(0)->Integral(0,20);
	float HeCont_NaF = HeCountsNaF/(PCountsNaF+DCountsNaF);
	
	float HeCountsAgl= HeliumContaminationAgl->GetResult_He(0)->Integral(0,20);//HeliumContaminationAgl->GetResult_He(0)->FindBin(4));
        float PCountsAgl = HeliumContaminationAgl->GetResult_P(0)->Integral(0,20);//HeliumContaminationAgl->GetResult_P(0)->FindBin(4));
	float DCountsAgl = HeliumContaminationNaF->GetResult_D(0)->Integral(0,20);
	float HeCont_Agl = HeCountsAgl/(PCountsAgl+DCountsAgl);	


	finalHistos.Add(HecutMC_He	);
	finalHistos.Add(Hecut_D   	);
	finalHistos.Add(HecutMCP_TH1F 	);
	finalHistos.Add(HecutMCHe_TH1F	);
	finalHistos.Add(HecutMC_P       );
	finalHistos.writeObjsInFolder("Results");

	cout<<"*** Plotting ...  ****"<<endl;

	Hecut_Plot(HecutMC_He	 ,	
		   Hecut_D	 ,   	     
                   HecutMCP_TH1F ,	
	           HecutMCHe_TH1F,	
                   HecutMC_P     ,
		  HeliumContaminationTOF,	
		  HeliumContaminationNaF,     
                  HeliumContaminationAgl,
		  HeCont_TOF,HeCont_NaF,HeCont_Agl
		  );    
	return; 
}
