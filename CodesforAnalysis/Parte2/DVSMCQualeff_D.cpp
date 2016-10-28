#include "PlottingFunctions/DVSMCQualeffD_Plot.h"


using namespace std;

DatavsMC * Dist_DvsMC_D = new DatavsMC("Dist_DvsMC_D",11,1,6);
DatavsMC * Lik_DvsMC_D  = new DatavsMC("Lik_DvsMC_D" ,11,1,6);

void DVSMCQualeffD_D_Fill(int zona){

	if(Tup.Beta<=0||Tup.R<=0||Tup.R<1.2*Tup.Rcutoff) return;
        if(!trgpatt.IsPhysical()) return;
        if(!Herejcut) return;
	if(!Betastrongcut) return;
	if(fabs(Tup.EdepTOFU-Tup.EdepTOFD)>1) return;



	//
	int Kbin;
	float mass = 0;

	//Beta bins
	Kbin=ToFDB.GetBin(RUsed);
	mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));

	float Lik = -log(1-Tup.LDiscriminant);
	
	if(mass>2&&mass<4)	{
	
		if(Lik > 1.7 && (Tup.Dist5D_P<3.5||Tup.Dist5D<3.5) ) {
			Dist_DvsMC_D -> DataEff[0] -> beforeTOF -> Fill(Kbin,zona);
			if(Distcut) Dist_DvsMC_D -> DataEff[0] -> afterTOF -> Fill(Kbin,zona);
		}

		if(Tup.Dist5D<2){
			Lik_DvsMC_D  -> DataEff[0] -> beforeTOF -> Fill(Kbin,zona);
			if(Likcut) Lik_DvsMC_D  -> DataEff[0] -> afterTOF -> Fill(Kbin,zona);
		}

	}
	//NaF
	if(cmask.isFromNaF() ) {	
		Kbin=NaFDB.GetBin(RUsed);
		mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
		if(mass>2&&mass<4 ) {

			if(Lik>2.85 && (Tup.Dist5D<2.5||Tup.Dist5D<2.5) ) {
				Dist_DvsMC_D -> DataEff[0] -> beforeNaF -> Fill(Kbin,zona);
				if(Distcut) Dist_DvsMC_D -> DataEff[0] -> afterNaF -> Fill(Kbin,zona);
			}

			if(Tup.Dist5D<1.5 && Lik > 2 ){
				Lik_DvsMC_D  -> DataEff[0] -> beforeNaF -> Fill(Kbin,zona);
				if(Lik > 2.4) Lik_DvsMC_D  -> DataEff[0] -> afterNaF -> Fill(Kbin,zona);
			}

		}
	}
	//Agl
	if(cmask.isFromAgl()) {
		Kbin=AglDB.GetBin(RUsed);
		mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
		if(mass>2&&mass<4) {

			if(Lik>3 && (Tup.Dist5D<4||Tup.Dist5D<4)  ) {
				Dist_DvsMC_D -> DataEff[0] -> beforeAgl -> Fill(Kbin,zona);
				if(Distcut) Dist_DvsMC_D -> DataEff[0] -> afterAgl -> Fill(Kbin,zona);
			}

			if(Tup.Dist5D<2 && Lik > 2){
				Lik_DvsMC_D  -> DataEff[0] -> beforeAgl -> Fill(Kbin,zona);
				if(Lik > 2.4) Lik_DvsMC_D  -> DataEff[0] -> afterAgl -> Fill(Kbin,zona);
			}

		}	
	}


	return;

}


void DVSMCQualeffD_Fill(){


	//cuts
	if(Tup.Beta<=0||Tup.R<=0) return;
	if(!trgpatt.IsPhysical()) return;
	if(!Herejcut) return;
	if(!Betastrongcut) return;
	if(fabs(Tup.EdepTOFU-Tup.EdepTOFD)>1) return;
 		

	if(!(Massa_gen>1&&Massa_gen<2)) return;

	int Kbin;
	float mass=0;	

	float Lik = -log(1-Tup.LDiscriminant);	
	//Beta bins

	//ToF
	Kbin=ToFDB.GetBin(RUsed);	
	mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));

	if(mass>2&&mass<4) {

		if(Lik> 1.7 && (Tup.Dist5D_P<3.5||Tup.Dist5D<3.5)) {
			((TH2*)Dist_DvsMC_D -> MCEff[0] -> beforeTOF) -> Fill(Kbin,ReturnMCGenType(),Tup.mcweight);
			if(Distcut)  ((TH2*)Dist_DvsMC_D -> MCEff[0] -> afterTOF) -> Fill(Kbin,ReturnMCGenType(),Tup.mcweight);
		}

		if(Tup.Dist5D<2){
			((TH2*)Lik_DvsMC_D  -> MCEff[0] -> beforeTOF) -> Fill(Kbin,ReturnMCGenType(),Tup.mcweight);
			if(Likcut)  ((TH2*)Lik_DvsMC_D  -> MCEff[0] -> afterTOF) -> Fill(Kbin,ReturnMCGenType(),Tup.mcweight);
		}

	}
	//NaF
	if(cmask.isFromNaF()) {	
		Kbin=NaFDB.GetBin(RUsed);	
		mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
		if(mass>2&&mass<4) {

			if(Lik>2.85  && (Tup.Dist5D_P<2.5||Tup.Dist5D<2.5)) {
				((TH2*)Dist_DvsMC_D -> MCEff[0] -> beforeNaF) -> Fill(Kbin,ReturnMCGenType(),Tup.mcweight);
				if(Distcut)  ((TH2*)Dist_DvsMC_D -> MCEff[0] -> afterNaF) -> Fill(Kbin,ReturnMCGenType(),Tup.mcweight);
			}

			if(Tup.Dist5D<1.5 && Lik>2 ){
				((TH2*)Lik_DvsMC_D  -> MCEff[0] -> beforeNaF) -> Fill(Kbin,ReturnMCGenType(),Tup.mcweight);
				if(Lik > 2.4)  ((TH2*)Lik_DvsMC_D  -> MCEff[0] -> afterNaF) -> Fill(Kbin,ReturnMCGenType(),Tup.mcweight);
			}
		}
	}


	//Agl
	if(cmask.isFromAgl()) {
		Kbin=AglDB.GetBin(RUsed);
		mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
		if(mass>2&&mass<4) {

			if(Lik>3 && (Tup.Dist5D_P<4||Tup.Dist5D<4)) {
				((TH2*)Dist_DvsMC_D -> MCEff[0] -> beforeAgl) -> Fill(Kbin,ReturnMCGenType(),Tup.mcweight);
				if(Distcut)  ((TH2*)Dist_DvsMC_D -> MCEff[0] -> afterAgl) -> Fill(Kbin,ReturnMCGenType(),Tup.mcweight);
			}

			if(Tup.Dist5D<2 &&Lik>2 ){
				((TH2*)Lik_DvsMC_D  -> MCEff[0] -> beforeAgl) -> Fill(Kbin,ReturnMCGenType(),Tup.mcweight);
				if(Lik>2.4)  ((TH2*)Lik_DvsMC_D  -> MCEff[0] -> afterAgl) -> Fill(Kbin,ReturnMCGenType(),Tup.mcweight);
			}

		}
	}

	return;
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

