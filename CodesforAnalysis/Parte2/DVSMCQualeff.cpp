#include "PlottingFunctions/DVSMCQualeff_Plot.h"

using namespace std;

DatavsMC * Dist_DvsMC_P = new DatavsMC("Dist_DvsMC_P",11,1,20);
DatavsMC * Lik_DvsMC_P  = new DatavsMC("Lik_DvsMC_P" ,11,1,20);

void DVSMCQualeff2_D_Fill(int zona){

	//cuts
	if(Tup.Beta<=0||Tup.R<=0||Tup.R<SF*Tup.Rcutoff) return;
	if(!trgpatt.IsPhysical()) return;
	if(!Herejcut) return;
	//if(Tup.EdepL1<=0.01) return;
	if(!ProtonsMassWindow) return;

	int Kbin;
	float mass = 0;

	//R bins

	mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));


	Kbin = PRB.GetRBin(Tup.R);

	for(int s=0;s<20;s++) {
		if(Qualitycut(log(1-Tup.LDiscriminant),-0.5-0.05*s,-1.5-0.05*s,-1.8-0.05*s)) {
			Dist_DvsMC_P -> DataEff[s] -> beforeR -> Fill(Kbin,zona);
			if(Distcut) Dist_DvsMC_P -> DataEff[s] -> afterR -> Fill(Kbin,zona);
		}
		float cut = 43 -s;
		if(Tup.Dist5D<cut||Tup.Dist5D_P<cut){
			Lik_DvsMC_P  -> DataEff[s] -> beforeR -> Fill(Kbin,zona);
			if(Likcut) Lik_DvsMC_P  -> DataEff[s] -> afterR -> Fill(Kbin,zona);
		}
	}

	//Beta bins
	Kbin=ToFPB.GetBin(RUsed);
	mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));	

	for(int s=0;s<20;s++) {
		if(Qualitycut(log(1-Tup.LDiscriminant),-0.5-0.05*s,-1.5-0.05*s,-1.8-0.05*s)) {	
			Dist_DvsMC_P -> DataEff[s] -> beforeTOF -> Fill(Kbin,zona);
			if(Distcut) Dist_DvsMC_P -> DataEff[s] -> afterTOF -> Fill(Kbin,zona);
		}
		float cut = 43 -s;
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

			if(Qualitycut(log(1-Tup.LDiscriminant),-0.5-0.05*s,-1.5-0.05*s,-1.8-0.05*s)) {
				Dist_DvsMC_P -> DataEff[s] -> beforeNaF -> Fill(Kbin,zona);
				if(Distcut) Dist_DvsMC_P -> DataEff[s] -> afterNaF -> Fill(Kbin,zona);
			}
			float cut = 43 -s;
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

			if(Qualitycut(log(1-Tup.LDiscriminant),-0.5-0.05*s,-1.5-0.05*s,-1.8-0.05*s)) {
				Dist_DvsMC_P -> DataEff[s] -> beforeAgl -> Fill(Kbin,zona);
				if(Distcut) Dist_DvsMC_P -> DataEff[s] -> afterAgl -> Fill(Kbin,zona);
			}

			float cut = 43 -s;
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
	//if(!trgpatt.IsPhysical()) return;       
	//if(Tup.EdepL1<=0.01) return;	
	if(!ProtonsMassWindow) return;

	//
	int Kbin;
	float mass=0;

	//R bins
	Kbin = PRB.GetRBin(Tup.R);



	if(Massa_gen<1) {
		//R bins
		Kbin = PRB.GetRBin(Tup.R);
		for(int s=0;s<20;s++) {

			if(Qualitycut(log(1-Tup.LDiscriminant),-0.5-0.05*s,-1.5-0.05*s,-1.8-0.05*s)) { 
				Dist_DvsMC_P -> MCEff[s] -> beforeR -> Fill(Kbin,Tup.mcweight);
				if(Distcut) Dist_DvsMC_P -> MCEff[s]-> afterR -> Fill(Kbin,Tup.mcweight);
			}

			float cut = 43 -s;
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

			if(Qualitycut(log(1-Tup.LDiscriminant),-0.5-0.05*s,-1.5-0.05*s,-1.8-0.05*s)) {                
				Dist_DvsMC_P -> MCEff[s] -> beforeTOF -> Fill(Kbin,Tup.mcweight);
				if(Distcut) Dist_DvsMC_P -> MCEff[s] -> afterTOF -> Fill(Kbin,Tup.mcweight);
			}

			float cut = 43 -s;
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

				if(Qualitycut(log(1-Tup.LDiscriminant),-0.5-0.05*s,-1.5-0.05*s,-1.8-0.05*s)) {
					Dist_DvsMC_P -> MCEff[s] -> beforeNaF -> Fill(Kbin,Tup.mcweight);
					if(Distcut) Dist_DvsMC_P -> MCEff[s] -> afterNaF -> Fill(Kbin,Tup.mcweight);
				}
				float cut = 43 -s;
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

				if(Qualitycut(log(1-Tup.LDiscriminant),-0.5-0.05*s,-1.5-0.05*s,-1.8-0.05*s)) {
					Dist_DvsMC_P -> MCEff[s] -> beforeAgl -> Fill(Kbin,Tup.mcweight);
					if(Distcut) Dist_DvsMC_P -> MCEff[s] -> afterAgl -> Fill(Kbin,Tup.mcweight);
				}
				float cut = 43 -s;
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

	LATcorr * LATLikelihoodDATA_TOF = new LATcorr(inputHistoFile,"LATLikDATA_TOF"  	 ,"Results");
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


	TH1F* DistP_MCEff_R   =(TH1F*) Dist_DvsMC_P -> GetMCEff_R(10)  ;
	TH1F* DistP_MCEff_TOF =(TH1F*) Dist_DvsMC_P -> GetMCEff_TOF(10);
	TH1F* DistP_MCEff_NaF =(TH1F*) Dist_DvsMC_P -> GetMCEff_NaF(10);
	TH1F* DistP_MCEff_Agl =(TH1F*) Dist_DvsMC_P -> GetMCEff_Agl(10);

	TH1F* DistP_DataEff_R   =(TH1F*) Dist_DvsMC_P -> GetDataEff_R(10)  ;
	TH1F* DistP_DataEff_TOF =(TH1F*) Dist_DvsMC_P -> GetDataEff_TOF(10);
	TH1F* DistP_DataEff_NaF =(TH1F*) Dist_DvsMC_P -> GetDataEff_NaF(10);
	TH1F* DistP_DataEff_Agl =(TH1F*) Dist_DvsMC_P -> GetDataEff_Agl(10);

	TH1F* LikP_MCEff_R   =(TH1F*) Lik_DvsMC_P -> GetMCEff_R(10)  ;
	TH1F* LikP_MCEff_TOF =(TH1F*) Lik_DvsMC_P -> GetMCEff_TOF(10);
	TH1F* LikP_MCEff_NaF =(TH1F*) Lik_DvsMC_P -> GetMCEff_NaF(10);
	TH1F* LikP_MCEff_Agl =(TH1F*) Lik_DvsMC_P -> GetMCEff_Agl(10);

	TH1F* LikP_DataEff_R   =(TH1F*) Lik_DvsMC_P -> GetDataEff_R(10)  ;
	TH1F* LikP_DataEff_TOF =(TH1F*) Lik_DvsMC_P -> GetDataEff_TOF(10);
	TH1F* LikP_DataEff_NaF =(TH1F*) Lik_DvsMC_P -> GetDataEff_NaF(10);
	TH1F* LikP_DataEff_Agl =(TH1F*) Lik_DvsMC_P -> GetDataEff_Agl(10);

	TH1F* DistP_Correction_R   =(TH1F*) Dist_DvsMC_P -> GetCorrection_R(10)  ;
	TH1F* DistP_Correction_TOF =(TH1F*) Dist_DvsMC_P -> GetCorrection_TOF(10);
	TH1F* DistP_Correction_NaF =(TH1F*) Dist_DvsMC_P -> GetCorrection_NaF(10);
	TH1F* DistP_Correction_Agl =(TH1F*) Dist_DvsMC_P -> GetCorrection_Agl(10);

	TH1F* LikP_Correction_R    =(TH1F*) Lik_DvsMC_P -> GetCorrection_R(10)  ;
	TH1F* LikP_Correction_TOF  =(TH1F*) Lik_DvsMC_P -> GetCorrection_TOF(10);
	TH1F* LikP_Correction_NaF  =(TH1F*) Lik_DvsMC_P -> GetCorrection_NaF(10);
	TH1F* LikP_Correction_Agl  =(TH1F*) Lik_DvsMC_P -> GetCorrection_Agl(10);


	cout<<"*************SYST ERR**************"<<endl;
        Dist_DvsMC_P ->Initialize_SystError();
        Lik_DvsMC_P  ->Initialize_SystError();

	Dist_DvsMC_P ->Eval_SystError();
        Lik_DvsMC_P  ->Eval_SystError();

        TH2F* SystDist_R   =(TH2F*) Dist_DvsMC_P -> GetSystPlot_R()  ;
        TH2F* SystDist_TOF =(TH2F*) Dist_DvsMC_P -> GetSystPlot_TOF();
        TH2F* SystDist_NaF =(TH2F*) Dist_DvsMC_P -> GetSystPlot_NaF();
        TH2F* SystDist_Agl =(TH2F*) Dist_DvsMC_P -> GetSystPlot_Agl();

        TH2F* SystLik_R   =(TH2F*) Lik_DvsMC_P -> GetSystPlot_R()  ;
        TH2F* SystLik_TOF =(TH2F*) Lik_DvsMC_P -> GetSystPlot_TOF();
        TH2F* SystLik_NaF =(TH2F*) Lik_DvsMC_P -> GetSystPlot_NaF();
        TH2F* SystLik_Agl =(TH2F*) Lik_DvsMC_P -> GetSystPlot_Agl();



	cout<<"************* FIT **************"<<endl;
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

        finalHistos.Add(SystDist_R  );
        finalHistos.Add(SystDist_TOF);
        finalHistos.Add(SystDist_NaF);
        finalHistos.Add(SystDist_Agl);

        finalHistos.Add(SystLik_R   );
        finalHistos.Add(SystLik_TOF );
        finalHistos.Add(SystLik_NaF );
        finalHistos.Add(SystLik_Agl );


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
			LikP_Correction_Agl , 
			

			DistP_MCEff_R  ,
                        DistP_MCEff_TOF,
                        DistP_MCEff_NaF,
                        DistP_MCEff_Agl,
                                                 
                        LikP_MCEff_R   ,
                        LikP_MCEff_TOF ,
                        LikP_MCEff_NaF ,
                        LikP_MCEff_Agl, 

			DistP_DataEff_R  ,
			DistP_DataEff_TOF,
                        DistP_DataEff_NaF,
                        DistP_DataEff_Agl,
                                                
                        LikP_DataEff_R   ,
                        LikP_DataEff_TOF ,
                        LikP_DataEff_NaF ,
                        LikP_DataEff_Agl 


);


	return;


}

