#include "PlottingFunctions/DVSMCQualeff_Plot.h"

using namespace std;

DatavsMC * Dist_DvsMC_P = new DatavsMC("Dist_DvsMC_P",11,1,1,20);
DatavsMC * Lik_DvsMC_P  = new DatavsMC("Lik_DvsMC_P" ,11,1,1,20);

void DVSMCQualeff2_D_Fill(int zona){

	//cuts
	if(Tup.Beta<=0||Tup.R<=0||Tup.R<1.2*Tup.Rcutoff) return;
	if(!trgpatt.IsPhysical()) return;
	if(!Herejcut) return;
	//if(Tup.EdepL1<=0.01) return;

	bool Likcutloose = Qualitycut(log(1-Tup.LDiscriminant),-0.3,-1.3,-1.3); 

	int Kbin;
	float mass = 0;

	//R bins

	mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));


	Kbin = PRB.GetRBin(Tup.R);

	for(int s=0;s<20;s++) {
		if(Qualitycut(log(1-Tup.LDiscriminant),-0.1+s*0.02,-0.9+s*0.02,-0.9+s*0.02)) {
			Dist_DvsMC_P -> DataEff[s] -> beforeR -> Fill(Kbin,zona);
			if(Distcut) Dist_DvsMC_P -> DataEff[s] -> afterR -> Fill(Kbin,zona);
		}
		float cut = 2+0.1*2;
		if(Tup.Dist5D<cut||Tup.Dist5D_P<cut){
			Lik_DvsMC_P  -> DataEff[s] -> beforeR -> Fill(Kbin,zona);
			if(Likcut) Lik_DvsMC_P  -> DataEff[s] -> afterR -> Fill(Kbin,zona);
		}
	}

	//Beta bins
	Kbin=ToFPB.GetBin(RUsed);
	mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));	

	for(int s=0;s<20;s++) {
		if(Qualitycut(log(1-Tup.LDiscriminant),-0.1+s*0.02,-0.9+s*0.02,-0.9+s*0.02)) {	
			Dist_DvsMC_P -> DataEff[s] -> beforeTOF -> Fill(Kbin,zona);
			if(Distcut) Dist_DvsMC_P -> DataEff[s] -> afterTOF -> Fill(Kbin,zona);
		}
		float cut = 2+0.1*2;
		if(Tup.Dist5D<cut||Tup.Dist5D_P<cut){
			Lik_DvsMC_P  -> DataEff[s] -> beforeTOF -> Fill(Kbin,zona);
			if(Likcut) Lik_DvsMC_P  -> DataEff[s] -> afterTOF -> Fill(Kbin,zona);
		}
	}

	//NaF
	if(cmask.isFromNaF()) {	
		Kbin=NaFPB.GetBin(RUsed);
		mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
		for(int s=0;s<20;s++) {

			if(Qualitycut(log(1-Tup.LDiscriminant),-0.1+s*0.02,-0.9+s*0.02,-0.9+s*0.02)) {
				Dist_DvsMC_P -> DataEff[s] -> beforeNaF -> Fill(Kbin,zona);
				if(Distcut) Dist_DvsMC_P -> DataEff[s] -> afterNaF -> Fill(Kbin,zona);
			}
			float cut = 2+0.1*2;
			if(Tup.Dist5D<cut||Tup.Dist5D_P<cut){

				Lik_DvsMC_P  -> DataEff[s] -> beforeNaF -> Fill(Kbin,zona);
				if(Likcut) Lik_DvsMC_P  -> DataEff[s] -> afterNaF -> Fill(Kbin,zona);
			}

		}
	}
	//Agl
	if(cmask.isFromAgl()) {
		Kbin=AglPB.GetBin(RUsed);
		mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
		for(int s=0;s<20;s++) {

			if(Qualitycut(log(1-Tup.LDiscriminant),-0.1+s*0.02,-0.9+s*0.02,-0.9+s*0.02)) {
				Dist_DvsMC_P -> DataEff[s] -> beforeAgl -> Fill(Kbin,zona);
				if(Distcut) Dist_DvsMC_P -> DataEff[s] -> afterAgl -> Fill(Kbin,zona);
			}

			float cut = 2+0.1*2;
			if(Tup.Dist5D<cut||Tup.Dist5D_P<cut){
				Lik_DvsMC_P  -> DataEff[s] -> beforeAgl -> Fill(Kbin,zona);
				if(Likcut) Lik_DvsMC_P  -> DataEff[s] -> afterAgl -> Fill(Kbin,zona);
			}
		}

	}
	return;

}

void DVSMCQualeff2_Fill(){

	//cuts
	if(Tup.Beta<=0||Tup.R<=0) return;
	if(!Herejcut) return;
	if(!trgpatt.IsPhysical()) return;       
	//if(Tup.EdepL1<=0.01) return;	

	//
	int Kbin;
	float mass=0;

	//R bins
	Kbin = PRB.GetRBin(Tup.R);

	bool Likcutloose = Qualitycut(log(1-Tup.LDiscriminant),-0.3,-1.3,-1.3);


	if(Massa_gen<1) {
		//R bins
		Kbin = PRB.GetRBin(Tup.R);
		for(int s=0;s<20;s++) {

			if(Qualitycut(log(1-Tup.LDiscriminant),-0.1+s*0.02,-0.9+s*0.02,-0.9+s*0.02)) { 
				Dist_DvsMC_P -> MCEff[s] -> beforeR -> Fill(Kbin,Tup.mcweight);
				if(Distcut) Dist_DvsMC_P -> MCEff[s]-> afterR -> Fill(Kbin,Tup.mcweight);
			}

			float cut = 2+0.1*2;
			if(Tup.Dist5D<cut||Tup.Dist5D_P<cut){
				Lik_DvsMC_P  -> MCEff[s] -> beforeR -> Fill(Kbin,Tup.mcweight);
				if(Likcut) Lik_DvsMC_P  -> MCEff[s] -> afterR -> Fill(Kbin,Tup.mcweight);
			}
		}
		//Beta bins

		//ToF
		Kbin=ToFPB.GetBin(RUsed);	
		mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));
		for(int s=0;s<20;s++) {

			if(Qualitycut(log(1-Tup.LDiscriminant),-0.1+s*0.02,-0.9+s*0.02,-0.9+s*0.02)) {                
				Dist_DvsMC_P -> MCEff[s] -> beforeTOF -> Fill(Kbin,Tup.mcweight);
				if(Distcut) Dist_DvsMC_P -> MCEff[s] -> afterTOF -> Fill(Kbin,Tup.mcweight);
			}

			float cut = 2+0.1*2;
			if(Tup.Dist5D<cut||Tup.Dist5D_P<cut){
				Lik_DvsMC_P  -> MCEff[s] -> beforeTOF -> Fill(Kbin,Tup.mcweight);
				if(Likcut) Lik_DvsMC_P  -> MCEff[s] -> afterTOF -> Fill(Kbin,Tup.mcweight);
			}
		}
		//NaF
		if(cmask.isFromNaF()) {	
			Kbin=NaFPB.GetBin(RUsed);	
			mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
			for(int s=0;s<20;s++) {

				if(Qualitycut(log(1-Tup.LDiscriminant),-0.1+s*0.02,-0.9+s*0.02,-0.9+s*0.02)) {
					Dist_DvsMC_P -> MCEff[s] -> beforeNaF -> Fill(Kbin,Tup.mcweight);
					if(Distcut) Dist_DvsMC_P -> MCEff[s] -> afterNaF -> Fill(Kbin,Tup.mcweight);
				}
				float cut = 2+0.1*2;
				if(Tup.Dist5D<cut||Tup.Dist5D_P<cut){
					Lik_DvsMC_P  -> MCEff[s] -> beforeNaF -> Fill(Kbin,Tup.mcweight);
					if(Likcut) Lik_DvsMC_P  -> MCEff[s] -> afterNaF -> Fill(Kbin,Tup.mcweight);
				}
			}
		}
		//Agl
		if(cmask.isFromAgl()) {	
			Kbin=AglPB.GetBin(RUsed);
			mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
			for(int s=0;s<20;s++) {

				if(Qualitycut(log(1-Tup.LDiscriminant),-0.1+s*0.02,-0.9+s*0.02,-0.9+s*0.02)) {
					Dist_DvsMC_P -> MCEff[s] -> beforeAgl -> Fill(Kbin,Tup.mcweight);
					if(Distcut) Dist_DvsMC_P -> MCEff[s] -> afterAgl -> Fill(Kbin,Tup.mcweight);
				}
				float cut = 2+0.1*2;
				if(Tup.Dist5D<cut||Tup.Dist5D_P<cut){
					Lik_DvsMC_P  -> MCEff[s] -> beforeAgl -> Fill(Kbin,Tup.mcweight);
					if(Likcut) Lik_DvsMC_P  -> MCEff[s] -> afterAgl -> Fill(Kbin,Tup.mcweight);
				}
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


	DatavsMC * Dist_DvsMC_P = new DatavsMC(inputHistoFile,"Dist_DvsMC_P",1,20);
	DatavsMC * Lik_DvsMC_P  = new DatavsMC(inputHistoFile,"Lik_DvsMC_P" ,1,20);

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


	TH1F* DistP_Correction_R   =(TH1F*) Dist_DvsMC_P -> GetCorrection_R(10)  ;
	TH1F* DistP_Correction_TOF =(TH1F*) Dist_DvsMC_P -> GetCorrection_TOF(10);
	TH1F* DistP_Correction_NaF =(TH1F*) Dist_DvsMC_P -> GetCorrection_NaF(10);
	TH1F* DistP_Correction_Agl =(TH1F*) Dist_DvsMC_P -> GetCorrection_Agl(10);

	TH1F* LikP_Correction_R    =(TH1F*) Lik_DvsMC_P -> GetCorrection_R(10)  ;
	TH1F* LikP_Correction_TOF  =(TH1F*) Lik_DvsMC_P -> GetCorrection_TOF(10);
	TH1F* LikP_Correction_NaF  =(TH1F*) Lik_DvsMC_P -> GetCorrection_NaF(10);
	TH1F* LikP_Correction_Agl  =(TH1F*) Lik_DvsMC_P -> GetCorrection_Agl(10);


	Dist_DvsMC_P ->Eval_FittedCorrections();
	Lik_DvsMC_P  ->Eval_FittedCorrections();


	TH1F* DistP_CorrectionFit_R   =(TH1F*) Dist_DvsMC_P -> GetCorrection_R(10)  ;
	TH1F* DistP_CorrectionFit_TOF =(TH1F*) Dist_DvsMC_P -> GetCorrection_TOF(10);
	TH1F* DistP_CorrectionFit_NaF =(TH1F*) Dist_DvsMC_P -> GetCorrection_NaF(10);
	TH1F* DistP_CorrectionFit_Agl =(TH1F*) Dist_DvsMC_P -> GetCorrection_Agl(10);

	TH1F* LikP_CorrectionFit_R    =(TH1F*) Lik_DvsMC_P -> GetCorrection_R(10)  ;
	TH1F* LikP_CorrectionFit_TOF  =(TH1F*) Lik_DvsMC_P -> GetCorrection_TOF(10);
	TH1F* LikP_CorrectionFit_NaF  =(TH1F*) Lik_DvsMC_P -> GetCorrection_NaF(10);
	TH1F* LikP_CorrectionFit_Agl  =(TH1F*) Lik_DvsMC_P -> GetCorrection_Agl(10);

	DistP_CorrectionFit_R    -> SetName("Dist_DvsMC_P_CorrectionR"  );
	DistP_CorrectionFit_TOF  -> SetName("Dist_DvsMC_P_CorrectionTOF");
	DistP_CorrectionFit_NaF  -> SetName("Dist_DvsMC_P_CorrectionNaF");
	DistP_CorrectionFit_Agl  -> SetName("Dist_DvsMC_P_CorrectionAgl");

	LikP_CorrectionFit_R     -> SetName("Lik_DvsMC_P_CorrectionR"    );
	LikP_CorrectionFit_TOF   -> SetName("Lik_DvsMC_P_CorrectionTOF"  );
	LikP_CorrectionFit_NaF   -> SetName("Lik_DvsMC_P_CorrectionNaF"  );
	LikP_CorrectionFit_Agl   -> SetName("Lik_DvsMC_P_CorrectionAgl"  );



	finalHistos.Add(DistP_CorrectionFit_R  );
	finalHistos.Add(DistP_CorrectionFit_TOF);
	finalHistos.Add( DistP_CorrectionFit_NaF);
	finalHistos.Add( DistP_CorrectionFit_Agl);

	finalHistos.Add( LikP_CorrectionFit_R   );
	finalHistos.Add( LikP_CorrectionFit_TOF );
	finalHistos.Add( LikP_CorrectionFit_NaF );
	finalHistos.Add( LikP_CorrectionFit_Agl );
	finalHistos.writeObjsInFolder("Results");

	cout<<"*** Plotting ...  ****"<<endl;

	DVSMCQualeff_Plot(	
			DistP_CorrectionFit_R  ,
			DistP_CorrectionFit_TOF,
			DistP_CorrectionFit_NaF,
			DistP_CorrectionFit_Agl,

			LikP_CorrectionFit_R   ,
			LikP_CorrectionFit_TOF ,
			LikP_CorrectionFit_NaF ,
			LikP_CorrectionFit_Agl, 

			DistP_Correction_R  ,
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

