#include "PlottingFunctions/DVSMCQualeffD_Plot.h"


using namespace std;

DatavsMC * Dist_DvsMC_D = new DatavsMC("Dist_DvsMC_D",11,1,6);
DatavsMC * Lik_DvsMC_D  = new DatavsMC("Lik_DvsMC_D" ,11,1,6);

void DVSMCQualeffD_D_Fill(int zona){

	
	if(Tup.R<1.2*Tup.Rcutoff||Tup.Beta>deutons->Eval(Tup.R)+0.1||Tup.Beta<deutons->Eval(Tup.R)-0.1) return;
	if(!((Tup.R>Rcut[zona]&&zona<10)||(zona==10)))  return;
	if(!Herejcut) return;
	if(!Betastrongcut) return;
	//
	int Kbin;
	
	//Beta bins
	//ToF
	Kbin=ToFDB.GetRBin(RUsed);	
	Dist_DvsMC_D -> DataEff -> beforeTOF -> Fill(Kbin,zona);
	if(Distcut) Lik_DvsMC_D  -> DataEff -> beforeTOF -> Fill(Kbin,zona);

	if(Distcut){
		Dist_DvsMC_D -> DataEff -> afterTOF -> Fill(Kbin,zona);
		if(Likcut) Lik_DvsMC_D  -> DataEff -> afterTOF -> Fill(Kbin,zona);
	}
	//NaF
	if(cmask.isFromNaF()) {	
		Kbin=NaFDB.GetRBin(RUsed);
		Dist_DvsMC_D -> DataEff -> beforeNaF -> Fill(Kbin,zona);
		if(Distcut) Lik_DvsMC_D  -> DataEff -> beforeNaF -> Fill(Kbin,zona);

		if(Distcut){
			Dist_DvsMC_D -> DataEff -> afterNaF -> Fill(Kbin,zona);
			if(Likcut) Lik_DvsMC_D  -> DataEff -> afterNaF -> Fill(Kbin,zona);
		}
	}
	//Agl
	if(cmask.isFromAgl()) {
		Kbin=AglDB.GetRBin(RUsed);
		Dist_DvsMC_D -> DataEff -> beforeAgl -> Fill(Kbin,zona);
		if(Distcut) Lik_DvsMC_D  -> DataEff -> beforeAgl -> Fill(Kbin,zona);

		if(Distcut){
			Dist_DvsMC_D -> DataEff -> afterAgl -> Fill(Kbin,zona);
			if(Likcut) Lik_DvsMC_D  -> DataEff -> afterAgl -> Fill(Kbin,zona);
		}
	}
	return;

}

void DVSMCQualeffD_Fill(){


	//cuts
	if(!Betastrongcut) return;
	if(!(Massa_gen>1&&Massa_gen<2)) return;
	
	int Kbin;


		//Beta bins

		//ToF
		Kbin=ToFDB.GetRBin(RUsed);	
		Dist_DvsMC_D -> MCEff -> beforeTOF -> Fill(Kbin,ReturnMCGenType());
		if(Distcut) Lik_DvsMC_D  -> MCEff -> beforeTOF -> Fill(Kbin,ReturnMCGenType());

		if(Distcut){
			Dist_DvsMC_D -> MCEff -> afterTOF -> Fill(Kbin,ReturnMCGenType());
			if(Likcut) Lik_DvsMC_D  -> MCEff -> afterTOF -> Fill(Kbin,ReturnMCGenType());
		}
		//NaF
		if(cmask.isFromNaF()) {	
			Kbin=NaFDB.GetRBin(RUsed);	
			Dist_DvsMC_D -> MCEff -> beforeNaF -> Fill(Kbin,ReturnMCGenType());
			if(Distcut) Lik_DvsMC_D  -> MCEff -> beforeNaF -> Fill(Kbin,ReturnMCGenType());

			if(Distcut){
				Dist_DvsMC_D -> MCEff -> afterNaF -> Fill(Kbin,ReturnMCGenType());
				if(Likcut) Lik_DvsMC_D  -> MCEff -> afterNaF -> Fill(Kbin,ReturnMCGenType());
			}

		}
		//Agl
		if(cmask.isFromAgl()) {	
			Kbin=AglDB.GetRBin(RUsed);
			Dist_DvsMC_D -> MCEff -> beforeAgl -> Fill(Kbin,ReturnMCGenType());
			if(Distcut) Lik_DvsMC_D  -> MCEff -> beforeAgl -> Fill(Kbin,ReturnMCGenType());

			if(Distcut){
				Dist_DvsMC_D -> MCEff -> afterAgl -> Fill(Kbin,ReturnMCGenType());
				if(Likcut) Lik_DvsMC_D  -> MCEff -> afterAgl -> Fill(Kbin,ReturnMCGenType());
			}

		}

                  
}

void DVSMCQualeffD_Write(){

	Dist_DvsMC_D -> Write();
	Lik_DvsMC_D  -> Write();

	return;
}


void DVSMCQualeffD(string filename){

	cout<<"******* Data vs MC: QUALITY SEL (D)********"<<endl;

	 cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


	DatavsMC * Dist_DvsMC_D = new DatavsMC(inputHistoFile,"Dist_DvsMC_D",6);
	DatavsMC * Lik_DvsMC_D  = new DatavsMC(inputHistoFile,"Lik_DvsMC_D" ,6);

	LATcorr * LATLikelihoodDATA_TOF = new LATcorr(inputHistoFile,"LATLikDATA_TOF"   	 ,"Results");
	LATcorr * LATDistanceDATA_TOF   = new LATcorr(inputHistoFile,"LATDistDATA_TOF" 	 ,"Results");

	LATcorr * LATLikelihoodDATA_NaF = new LATcorr(inputHistoFile,"LATLikDATA_NaF"  	 ,"Results");
	LATcorr * LATDistanceDATA_NaF   = new LATcorr(inputHistoFile,"LATDistDATA_NaF" 	 ,"Results");

	LATcorr * LATLikelihoodDATA_Agl = new LATcorr(inputHistoFile,"LATLikDATA_Agl"  	 ,"Results");
	LATcorr * LATDistanceDATA_Agl   = new LATcorr(inputHistoFile,"LATDistDATA_Agl" 	 ,"Results");




	cout<<"******* Data vs MC: QUALITY SEL (D)********"<<endl;

	Dist_DvsMC_D -> Assign_LatCorr( LATDistanceDATA_TOF   ->  LATcorrR_fit , 
					LATDistanceDATA_TOF   ->  LATcorrR_fit ,
					LATDistanceDATA_NaF   ->  LATcorrR_fit ,
					LATDistanceDATA_Agl   ->  LATcorrR_fit );

	Lik_DvsMC_D  ->	Assign_LatCorr( LATLikelihoodDATA_TOF ->  LATcorrR_fit , 	
					LATLikelihoodDATA_TOF ->  LATcorrR_fit ,
					LATLikelihoodDATA_NaF ->  LATcorrR_fit ,
					LATLikelihoodDATA_Agl ->  LATcorrR_fit );



	Dist_DvsMC_D ->Eval_DandMC_Eff();  
	Lik_DvsMC_D  ->Eval_DandMC_Eff();

	Dist_DvsMC_D ->Eval_Corrections();
	Lik_DvsMC_D  ->Eval_Corrections();


	TH2F* DistD_Correction_R   =(TH2F*) Dist_DvsMC_D -> GetCorrection_R()  ;
	TH2F* DistD_Correction_TOF =(TH2F*) Dist_DvsMC_D -> GetCorrection_TOF();
	TH2F* DistD_Correction_NaF =(TH2F*) Dist_DvsMC_D -> GetCorrection_NaF();
	TH2F* DistD_Correction_Agl =(TH2F*) Dist_DvsMC_D -> GetCorrection_Agl();

	TH2F* LikD_Correction_R    =(TH2F*) Lik_DvsMC_D -> GetCorrection_R()  ;
	TH2F* LikD_Correction_TOF  =(TH2F*) Lik_DvsMC_D -> GetCorrection_TOF();
	TH2F* LikD_Correction_NaF  =(TH2F*) Lik_DvsMC_D -> GetCorrection_NaF();
	TH2F* LikD_Correction_Agl  =(TH2F*) Lik_DvsMC_D -> GetCorrection_Agl();


	DistD_Correction_R    -> SetName("Dist_DvsMC_D_CorrectionR"  );
	DistD_Correction_TOF  -> SetName("Dist_DvsMC_D_CorrectionTOF");
	DistD_Correction_NaF  -> SetName("Dist_DvsMC_D_CorrectionNaF");
	DistD_Correction_Agl  -> SetName("Dist_DvsMC_D_CorrectionAgl");

	LikD_Correction_R    -> SetName("Lik_DvsMC_D_CorrectionR"  );
	LikD_Correction_TOF  -> SetName("Lik_DvsMC_D_CorrectionTOF");
	LikD_Correction_NaF  -> SetName("Lik_DvsMC_D_CorrectionNaF");
	LikD_Correction_Agl  -> SetName("Lik_DvsMC_D_CorrectionAgl");


	finalHistos.Add(DistD_Correction_R  );
       finalHistos.Add( DistD_Correction_TOF);
       finalHistos.Add( DistD_Correction_NaF);
	finalHistos.Add(DistD_Correction_Agl);	
                            
       finalHistos.Add( LikD_Correction_R   );
       finalHistos.Add( LikD_Correction_TOF );
       finalHistos.Add( LikD_Correction_NaF );
	finalHistos.Add(LikD_Correction_Agl );
	finalHistos.writeObjsInFolder("Results");

        cout<<"*** Plotting ...  ****"<<endl;

	DVSMCQualeffD_Plot (DistD_Correction_R  ,
                            DistD_Correction_TOF,
                            DistD_Correction_NaF,
                            DistD_Correction_Agl,
                                
                            LikD_Correction_R   ,
                            LikD_Correction_TOF ,
                            LikD_Correction_NaF ,
                            LikD_Correction_Agl 

	);

	return;
}

