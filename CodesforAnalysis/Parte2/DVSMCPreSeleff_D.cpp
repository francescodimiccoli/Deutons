#include "PlottingFunctions/DVSMCPreSeleffD.h"


using namespace std;

DatavsMC * PreSel_DvsMC_D = new DatavsMC("PreSel_DvsMC_D",11,3,6);

void DVSMCPreSeleffD_D_Fill(int zona){


	if(Tup.R_pre<=0||Tup.R_pre<1.2*Tup.Rcutoff||Tup.Beta_pre>protons->Eval(Tup.R_pre)+0.1||Tup.Beta_pre<protons->Eval(Tup.R_pre)-0.1) return;
	if(!((Tup.R_pre>Rcut[zona]&&zona<10)||(zona==10)))  return;
	if(!Herejcut) return;
	if(!(Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+1&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-1)) return;
	//
	int Kbin;
	for(int S=0;S<3;S++){
		//Beta bins
		//ToF
		Kbin=ToFDB.GetRBin(RUsed);	
		if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_D -> DataEff -> beforeTOF) -> Fill(Kbin,zona,S);
		if(cmask.passed(S))       ((TH3*)PreSel_DvsMC_D -> DataEff -> afterTOF ) -> Fill(Kbin,zona,S);
		//NaF
		if(cmask.isFromNaF()) {	
			Kbin=NaFDB.GetRBin(RUsed);
			if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_D -> DataEff -> beforeNaF) -> Fill(Kbin,zona,S);
			if(cmask.passed(S))       ((TH3*)PreSel_DvsMC_D -> DataEff -> afterNaF ) -> Fill(Kbin,zona,S);
		}
		//Agl
		if(cmask.isFromAgl()) {
			Kbin=AglDB.GetRBin(RUsed);
			if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_D -> DataEff -> beforeAgl) -> Fill(Kbin,zona,S);
			if(cmask.passed(S))       ((TH3*)PreSel_DvsMC_D -> DataEff -> afterAgl ) -> Fill(Kbin,zona,S);
		}
	}
	return;

}
void DVSMCPreSeleffD_Fill(){
	//cuts
	if(Tup.Beta_pre<=0||Tup.R_pre<=0||Tup.Beta_pre>protons->Eval(Tup.R_pre)+0.1||Tup.Beta_pre<protons->Eval(Tup.R_pre)-0.1) return;
	if(!Herejcut) return;
	if(!(Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+1&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-1)) return;
	//
	int Kbin;
	for(int S=0;S<3;S++){
		if(Massa_gen>1&&Massa_gen<2) {
			//Beta bins
			//ToF
			Kbin=ToFDB.GetRBin(RUsed);	
			if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_D -> MCEff -> beforeTOF) -> Fill(Kbin,ReturnMCGenType(),S);
			if(cmask.passed(S))       ((TH3*)PreSel_DvsMC_D -> MCEff -> afterTOF ) -> Fill(Kbin,ReturnMCGenType(),S);

			//NaF
			if(cmask.isFromNaF()) {	
				Kbin=NaFDB.GetRBin(RUsed);	
				if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_D -> MCEff -> beforeNaF) -> Fill(Kbin,ReturnMCGenType(),S);
				if(cmask.passed(S))       ((TH3*)PreSel_DvsMC_D -> MCEff -> afterNaF ) -> Fill(Kbin,ReturnMCGenType(),S);
			}
			//Agl
			if(cmask.isFromAgl()) {	
				Kbin=AglDB.GetRBin(RUsed);
				if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_D -> MCEff -> beforeAgl) -> Fill(Kbin,ReturnMCGenType(),S);
				if(cmask.passed(S))       ((TH3*)PreSel_DvsMC_D -> MCEff -> afterAgl ) -> Fill(Kbin,ReturnMCGenType(),S);
			}

		}

	}
	return;
}


void DVSMCPreSeleffD_Write(){

	PreSel_DvsMC_D -> Write();
	return;
}


void DVSMCPreSeleffD(string filename){

	cout<<"******* Data vs MC:  PRESELECTIONS (D) ********"<<endl;
	        cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


	DatavsMC * PreSel_DvsMC_D = new DatavsMC(inputHistoFile,"PreSel_DvsMC_D",6);

	LATcorr * LATpreSelDATA = new LATcorr(inputHistoFile,"LATpreSelDATA"      ,"Results");


	cout<<"******* Data vs MC:  PRESELECTIONS (D) ********"<<endl;

	PreSel_DvsMC_D -> Assign_LatCorr( LATpreSelDATA   ->  LATcorrR_fit , 
			LATpreSelDATA   ->  LATcorrR_fit ,
			LATpreSelDATA   ->  LATcorrR_fit ,
			LATpreSelDATA   ->  LATcorrR_fit );


	PreSel_DvsMC_D ->Eval_DandMC_Eff();  

	PreSel_DvsMC_D ->Eval_Corrections();


	TH3F* PreSelD_Correction_R   =(TH3F*) PreSel_DvsMC_D -> GetCorrection_R()  ;
	TH3F* PreSelD_Correction_TOF =(TH3F*) PreSel_DvsMC_D -> GetCorrection_TOF();
	TH3F* PreSelD_Correction_NaF =(TH3F*) PreSel_DvsMC_D -> GetCorrection_NaF();
	TH3F* PreSelD_Correction_Agl =(TH3F*) PreSel_DvsMC_D -> GetCorrection_Agl();



	PreSelD_Correction_R    -> SetName("PreSel_DvsMC_D_CorrectionR"  );
	PreSelD_Correction_TOF  -> SetName("PreSel_DvsMC_D_CorrectionTOF");
	PreSelD_Correction_NaF  -> SetName("PreSel_DvsMC_D_CorrectionNaF");
	PreSelD_Correction_Agl  -> SetName("PreSel_DvsMC_D_CorrectionAgl");

	
	finalHistos.Add(PreSelD_Correction_R   );
       finalHistos.Add( PreSelD_Correction_TOF );
       finalHistos.Add( PreSelD_Correction_NaF );
      finalHistos.Add(  PreSelD_Correction_Agl );
	finalHistos.writeObjsInFolder("Results");

        cout<<"*** Plotting ...  ****"<<endl;
	DVSMCPreSeleffD_Plot(	PreSelD_Correction_R  ,
	                     	PreSelD_Correction_TOF,
                             	PreSelD_Correction_NaF,
                             	PreSelD_Correction_Agl

	);

	return;
}

