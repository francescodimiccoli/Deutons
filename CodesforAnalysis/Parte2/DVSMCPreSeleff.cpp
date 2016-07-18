#include "PlottingFunctions/DVSMCPreSeleff_Plot.h"


using namespace std;

DatavsMC * PreSel_DvsMC_P = new DatavsMC("PreSel_DvsMC_P",11,3);

void DVSMCPreSeleff_D_Fill(int zona){

	//cuts
	if(!trgpatt.IsPhysical()||Tup.R_pre<=0||Tup.R_pre<1.2*Tup.Rcutoff||!ProtonsMassWindow) return;
	if(!((Tup.R_pre>Rcut[zona]&&zona<10)||(zona==10)))  return;
	if(!Herejcut) return;
	if(!(Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+1&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-1)) return;
	//
	int Kbin;
	for(int S=0;S<3;S++){
		//R bins
		Kbin = PRB.GetRBin(RUsed);
		if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_P -> DataEff -> beforeR) -> Fill(Kbin,zona,S);	
		if(cmask.passed(S))	      ((TH3*)PreSel_DvsMC_P -> DataEff -> afterR ) -> Fill(Kbin,zona,S);     

		//Beta bins
		//ToF
		Kbin=ToFDB.GetBin(RUsed);	
		if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_P -> DataEff -> beforeTOF) -> Fill(Kbin,zona,S);
                if(cmask.passed(S))       ((TH3*)PreSel_DvsMC_P -> DataEff -> afterTOF ) -> Fill(Kbin,zona,S);
		//NaF
		if(cmask.isFromNaF()) {	
			Kbin=NaFDB.GetBin(RUsed);
			if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_P -> DataEff -> beforeNaF) -> Fill(Kbin,zona,S);
                	if(cmask.passed(S))       ((TH3*)PreSel_DvsMC_P -> DataEff -> afterNaF ) -> Fill(Kbin,zona,S);
		}
		//Agl
		if(cmask.isFromAgl()) {
			Kbin=AglDB.GetBin(RUsed);
			if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_P -> DataEff -> beforeAgl) -> Fill(Kbin,zona,S);
                	if(cmask.passed(S))       ((TH3*)PreSel_DvsMC_P -> DataEff -> afterAgl ) -> Fill(Kbin,zona,S);
		}
	}
	return;
	
}
void DVSMCPreSeleff_Fill(){

	if(!trgpatt.IsPhysical()||Tup.R_pre<=0||!ProtonsMassWindow) return;
	if(!Herejcut) return;
	if(!(Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+1&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-1)) return;
	//cuts
	//
	int Kbin;
	for(int S=0;S<3;S++){
		if(Massa_gen<1) {
			//R bins
			Kbin = PRB.GetRBin(RUsed);	
			 if(cmask.notPassed(S))		  ((TH2*)PreSel_DvsMC_P -> MCEff -> beforeR) -> Fill(Kbin,S,Tup.mcweight);
			 if(cmask.passed(S))                  ((TH2*)PreSel_DvsMC_P -> MCEff -> afterR ) -> Fill(Kbin,S,Tup.mcweight);
			//Beta bins

			//ToF
			Kbin=ToFDB.GetBin(RUsed);	
			if(cmask.notPassed(S)) ((TH2*)PreSel_DvsMC_P -> MCEff -> beforeTOF) -> Fill(Kbin,S,Tup.mcweight);
			if(cmask.passed(S))   ((TH2*)PreSel_DvsMC_P -> MCEff -> afterTOF ) -> Fill(Kbin,S,Tup.mcweight);

			//NaF
			if(cmask.isFromNaF()) {	
				Kbin=NaFDB.GetBin(RUsed);	
				if(cmask.notPassed(S)) ((TH2*)PreSel_DvsMC_P -> MCEff -> beforeNaF) -> Fill(Kbin,S,Tup.mcweight);
				if(cmask.passed(S))       ((TH2*)PreSel_DvsMC_P -> MCEff -> afterNaF ) -> Fill(Kbin,S,Tup.mcweight);
			}
			//Agl
			if(cmask.isFromAgl()) {	
				Kbin=AglDB.GetBin(RUsed);
				if(cmask.notPassed(S)) ((TH2*)PreSel_DvsMC_P -> MCEff -> beforeAgl) -> Fill(Kbin,S,Tup.mcweight);
				if(cmask.passed(S))       ((TH2*)PreSel_DvsMC_P -> MCEff -> afterAgl ) -> Fill(Kbin,S,Tup.mcweight);
			}

		}

	}
	return;
}


void DVSMCPreSeleff_Write(){

	PreSel_DvsMC_P -> Write();
	return;
}


void DVSMCPreSeleff(string filename){

	cout<<"******* Data vs MC:  PRESELECTIONS ********"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");	

	DatavsMC * PreSel_DvsMC_P = new DatavsMC(inputHistoFile,"PreSel_DvsMC_P");

	LATcorr * LATpreSelDATA = new LATcorr(inputHistoFile,"LATpreSelDATA"      ,"Results");


	cout<<"******* Data vs MC:  PRESELECTIONS ********"<<endl;

	PreSel_DvsMC_P -> Assign_LatCorr( LATpreSelDATA   ->  LATcorrR_fit , 
					  LATpreSelDATA   ->  LATcorrR_fit ,
					  LATpreSelDATA   ->  LATcorrR_fit ,
					  LATpreSelDATA   ->  LATcorrR_fit );


	PreSel_DvsMC_P ->Eval_DandMC_Eff();  

	PreSel_DvsMC_P ->Eval_Corrections();


	TH2F* PreSel_Correction_R   =(TH2F*) PreSel_DvsMC_P -> GetCorrection_R()  ;
	TH2F* PreSel_Correction_TOF =(TH2F*) PreSel_DvsMC_P -> GetCorrection_TOF();
	TH2F* PreSel_Correction_NaF =(TH2F*) PreSel_DvsMC_P -> GetCorrection_NaF();
	TH2F* PreSel_Correction_Agl =(TH2F*) PreSel_DvsMC_P -> GetCorrection_Agl();

	TH2F* EffData_R   =(TH2F*) PreSel_DvsMC_P -> DataEff_corr -> effR -> Clone();
	TH2F* EffMC_R     =(TH2F*) PreSel_DvsMC_P -> MCEff 	  -> effR -> Clone();


	PreSel_Correction_R    -> SetName("PreSel_DvsMC_P_CorrectionR"  );
	PreSel_Correction_TOF  -> SetName("PreSel_DvsMC_P_CorrectionTOF");
	PreSel_Correction_NaF  -> SetName("PreSel_DvsMC_P_CorrectionNaF");
	PreSel_Correction_Agl  -> SetName("PreSel_DvsMC_P_CorrectionAgl");

	
	finalHistos.Add(PreSel_Correction_R   );
        finalHistos.Add(PreSel_Correction_TOF );
        finalHistos.Add( PreSel_Correction_NaF );
        finalHistos.Add( PreSel_Correction_Agl );

	finalHistos.writeObjsInFolder("Results");

        cout<<"*** Plotting ...  ****"<<endl;


	DVSMCPreSeleff_Plot(PreSel_Correction_R  ,
	                    PreSel_Correction_TOF,
                            PreSel_Correction_NaF,
                            PreSel_Correction_Agl,
				EffData_R, 
                                EffMC_R   
	);

	return;
}


