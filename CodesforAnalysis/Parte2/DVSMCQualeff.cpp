#include "PlottingFunctions/DVSMCQualeff_Plot.h"

using namespace std;

DatavsMC * Dist_DvsMC_P = new DatavsMC("Dist_DvsMC_P",11);
DatavsMC * Lik_DvsMC_P  = new DatavsMC("Lik_DvsMC_P" ,11);

void DVSMCQualeff2_D_Fill(int zona){

	//cuts
	if(Tup.Beta<=0||Tup.R<=0||Tup.R<1.2*Tup.Rcutoff||Tup.Beta>protons->Eval(Tup.R)+0.1||Tup.Beta<protons->Eval(Tup.R)-0.1) return;
	if(!(Tup.R>1.2*Tup.Rcutoff)) return;
	if(!Herejcut) return;
	//
	int Kbin;
	
	//R bins
	Kbin = PRB.GetRBin(Tup.R);
	Dist_DvsMC_P -> DataEff -> beforeR -> Fill(Kbin,zona);	
	if(Tup.Dist5D_P<6) Lik_DvsMC_P  -> DataEff -> beforeR -> Fill(Kbin,zona);
	
	if(Tup.Dist5D_P<6){
		Dist_DvsMC_P -> DataEff -> afterR -> Fill(Kbin,zona);     
        	if(Likcut) Lik_DvsMC_P  -> DataEff -> afterR -> Fill(Kbin,zona);
	}

	//Beta bins
	//ToF
	Kbin=ToFDB.GetRBin(RUsed);	
	Dist_DvsMC_P -> DataEff -> beforeTOF -> Fill(Kbin,zona);
	if(Distcut) Lik_DvsMC_P  -> DataEff -> beforeTOF -> Fill(Kbin,zona);

	if(Distcut){
		Dist_DvsMC_P -> DataEff -> afterTOF -> Fill(Kbin,zona);
		if(Likcut) Lik_DvsMC_P  -> DataEff -> afterTOF -> Fill(Kbin,zona);
	}
	//NaF
	if(cmask.isFromNaF()) {	
		Kbin=NaFDB.GetRBin(RUsed);
		Dist_DvsMC_P -> DataEff -> beforeNaF -> Fill(Kbin,zona);
		if(Distcut) Lik_DvsMC_P  -> DataEff -> beforeNaF -> Fill(Kbin,zona);

		if(Distcut){
			Dist_DvsMC_P -> DataEff -> afterNaF -> Fill(Kbin,zona);
			if(Likcut) Lik_DvsMC_P  -> DataEff -> afterNaF -> Fill(Kbin,zona);
		}
	}
	//Agl
	if(cmask.isFromAgl()) {
		Kbin=AglDB.GetRBin(RUsed);
		Dist_DvsMC_P -> DataEff -> beforeAgl -> Fill(Kbin,zona);
		if(Distcut) Lik_DvsMC_P  -> DataEff -> beforeAgl -> Fill(Kbin,zona);

		if(Distcut){
			Dist_DvsMC_P -> DataEff -> afterAgl -> Fill(Kbin,zona);
			if(Likcut) Lik_DvsMC_P  -> DataEff -> afterAgl -> Fill(Kbin,zona);
		}
	}
	return;

}

void DVSMCQualeff2_Fill(){

	//cuts
	if(Tup.Beta<=0||Tup.R<=0||Tup.Beta>protons->Eval(Tup.R)+0.1||Tup.Beta<protons->Eval(Tup.R)-0.1) return;
        if(!Herejcut) return;
	//
	int Kbin;

	//R bins
	Kbin = PRB.GetRBin(Tup.R);

	if(Massa_gen<1) {
		//R bins
		Kbin = PRB.GetRBin(Tup.R);	
		Dist_DvsMC_P -> MCEff -> beforeR -> Fill(Kbin);
		if(Tup.Dist5D_P<6) Lik_DvsMC_P  -> MCEff -> beforeR -> Fill(Kbin);

		if(Tup.Dist5D_P<6){
			Dist_DvsMC_P -> MCEff -> afterR -> Fill(Kbin);
			if(Likcut) Lik_DvsMC_P  -> MCEff -> afterR -> Fill(Kbin);
		}
		//Beta bins

		//ToF
		Kbin=ToFDB.GetRBin(RUsed);	
		Dist_DvsMC_P -> MCEff -> beforeTOF -> Fill(Kbin);
		if(Distcut) Lik_DvsMC_P  -> MCEff -> beforeTOF -> Fill(Kbin);

		if(Distcut){
			Dist_DvsMC_P -> MCEff -> afterTOF -> Fill(Kbin);
			if(Likcut) Lik_DvsMC_P  -> MCEff -> afterTOF -> Fill(Kbin);
		}
		//NaF
		if(cmask.isFromNaF()) {	
			Kbin=NaFDB.GetRBin(RUsed);	
			Dist_DvsMC_P -> MCEff -> beforeNaF -> Fill(Kbin);
			if(Distcut) Lik_DvsMC_P  -> MCEff -> beforeNaF -> Fill(Kbin);

			if(Distcut){
				Dist_DvsMC_P -> MCEff -> afterNaF -> Fill(Kbin);
				if(Likcut) Lik_DvsMC_P  -> MCEff -> afterNaF -> Fill(Kbin);
			}

		}
		//Agl
		if(cmask.isFromAgl()) {	
			Kbin=AglDB.GetRBin(RUsed);
			Dist_DvsMC_P -> MCEff -> beforeAgl -> Fill(Kbin);
			if(Distcut) Lik_DvsMC_P  -> MCEff -> beforeAgl -> Fill(Kbin);

			if(Distcut){
				Dist_DvsMC_P -> MCEff -> afterAgl -> Fill(Kbin);
				if(Likcut) Lik_DvsMC_P  -> MCEff -> afterAgl -> Fill(Kbin);
			}

		}

	}                        
}

void DVSMCQualeff2_Write(){

	Dist_DvsMC_P -> Write();
	Lik_DvsMC_P  -> Write();

	return;
}


void DVSMCQualeff2(string filename){

	cout<<"******* Data vs MC: QUALITY SEL ********"<<endl;

        cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


	DatavsMC * Dist_DvsMC_P = new DatavsMC(inputHistoFile,"Dist_DvsMC_P");
	DatavsMC * Lik_DvsMC_P  = new DatavsMC(inputHistoFile,"Lik_DvsMC_P" );

	LATcorr * LATLikelihoodDATA_TOF = new LATcorr(inputHistoFile,"LATLikDATA_TOF"   	 ,"Results");
	LATcorr * LATDistanceDATA_TOF   = new LATcorr(inputHistoFile,"LATDistDATA_TOF" 	 ,"Results");

	LATcorr * LATLikelihoodDATA_NaF = new LATcorr(inputHistoFile,"LATLikDATA_NaF"  	 ,"Results");
	LATcorr * LATDistanceDATA_NaF   = new LATcorr(inputHistoFile,"LATDistDATA_NaF" 	 ,"Results");

	LATcorr * LATLikelihoodDATA_Agl = new LATcorr(inputHistoFile,"LATLikDATA_Agl"  	 ,"Results");
	LATcorr * LATDistanceDATA_Agl   = new LATcorr(inputHistoFile,"LATDistDATA_Agl" 	 ,"Results");




	cout<<"******* Data vs MC: QUALITY SEL ********"<<endl;

	Dist_DvsMC_P -> Assign_LatCorr( LATDistanceDATA_TOF   ->  LATcorrR_fit , 
			LATDistanceDATA_TOF   ->  LATcorrR_fit ,
			LATDistanceDATA_NaF   ->  LATcorrR_fit ,
			LATDistanceDATA_Agl   ->  LATcorrR_fit );

	Lik_DvsMC_P  ->	Assign_LatCorr( LATLikelihoodDATA_TOF ->  LATcorrR_fit , 	
			LATLikelihoodDATA_TOF ->  LATcorrR_fit ,
			LATLikelihoodDATA_NaF ->  LATcorrR_fit ,
			LATLikelihoodDATA_Agl ->  LATcorrR_fit );



	Dist_DvsMC_P ->Eval_DandMC_Eff();  
	Lik_DvsMC_P  ->Eval_DandMC_Eff();

	Dist_DvsMC_P ->Eval_Corrections();
	Lik_DvsMC_P  ->Eval_Corrections();


	TH1F* DistP_Correction_R   =(TH1F*) Dist_DvsMC_P -> GetCorrection_R()  ;
	TH1F* DistP_Correction_TOF =(TH1F*) Dist_DvsMC_P -> GetCorrection_TOF();
	TH1F* DistP_Correction_NaF =(TH1F*) Dist_DvsMC_P -> GetCorrection_NaF();
	TH1F* DistP_Correction_Agl =(TH1F*) Dist_DvsMC_P -> GetCorrection_Agl();

	TH1F* LikP_Correction_R    =(TH1F*) Lik_DvsMC_P -> GetCorrection_R()  ;
	TH1F* LikP_Correction_TOF  =(TH1F*) Lik_DvsMC_P -> GetCorrection_TOF();
	TH1F* LikP_Correction_NaF  =(TH1F*) Lik_DvsMC_P -> GetCorrection_NaF();
	TH1F* LikP_Correction_Agl  =(TH1F*) Lik_DvsMC_P -> GetCorrection_Agl();



	DistP_Correction_R    -> SetName("Dist_DvsMC_P_CorrectionR"  );
	DistP_Correction_TOF  -> SetName("Dist_DvsMC_P_CorrectionTOF");
	DistP_Correction_NaF  -> SetName("Dist_DvsMC_P_CorrectionNaF");
	DistP_Correction_Agl  -> SetName("Dist_DvsMC_P_CorrectionAgl");

	LikP_Correction_R     -> SetName("Lik_DvsMC_P_CorrectionR"    );
	LikP_Correction_TOF   -> SetName("Lik_DvsMC_P_CorrectionTOF"  );
	LikP_Correction_NaF   -> SetName("Lik_DvsMC_P_CorrectionNaF"  );
	LikP_Correction_Agl   -> SetName("Lik_DvsMC_P_CorrectionAgl"  );



	finalHistos.Add(DistP_Correction_R  );
        finalHistos.Add(DistP_Correction_TOF);
       finalHistos.Add( DistP_Correction_NaF);
       finalHistos.Add( DistP_Correction_Agl);
                            
       finalHistos.Add( LikP_Correction_R   );
       finalHistos.Add( LikP_Correction_TOF );
       finalHistos.Add( LikP_Correction_NaF );
       finalHistos.Add( LikP_Correction_Agl );
	finalHistos.writeObjsInFolder("Results");

        cout<<"*** Plotting ...  ****"<<endl;

	DVSMCQualeff_Plot(	DistP_Correction_R  	,
                                DistP_Correction_TOF,
                                DistP_Correction_NaF,
                                DistP_Correction_Agl,
                                    
                                LikP_Correction_R   ,
                                LikP_Correction_TOF ,
	                        LikP_Correction_NaF ,
                                LikP_Correction_Agl 
	);
	

	return;


}

