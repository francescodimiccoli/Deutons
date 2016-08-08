#include "PlottingFunctions/DVSMCRICHeff_Plot.h"


using namespace std;

DatavsMC * RICH_DvsMC_P = new DatavsMC("RICH_DvsMC_P",11);
DatavsMC * RICH_DvsMC_D = new DatavsMC("RICH_DvsMC_D",11,1,6);


void DVSMCRICHeff_D_Fill(int zona){

	//cuts
	
	if(Tup.Beta<=0||Tup.R<=0||Tup.R<1.2*Tup.Rcutoff) return;
        if(!trgpatt.IsPhysical()) return;
        if(!Herejcut) return;
        if(!ProtonsMassThres) return;	//
	//if(!(Distcut&&Likcut)) return;
	float mass = 0;	
 
	int Kbin;
	
	//Beta bins
	//NaF
	Kbin=NaFPB.GetBin(RUsed);
	RICH_DvsMC_P -> DataEff -> beforeNaF -> Fill(Kbin,zona);	
	RICH_DvsMC_D -> DataEff -> beforeNaF -> Fill(Kbin,zona);
	mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
	if(cmask.isFromNaF()/*&&mass>0&&mass<3*/) {
		RICH_DvsMC_P -> DataEff -> afterNaF -> Fill(Kbin,zona);
		RICH_DvsMC_D -> DataEff -> afterNaF -> Fill(Kbin,zona);
	}
	//Agl
	Kbin=AglPB.GetBin(RUsed);
        RICH_DvsMC_P -> DataEff -> beforeAgl -> Fill(Kbin,zona);	
	RICH_DvsMC_D -> DataEff -> beforeAgl -> Fill(Kbin,zona);
	if(cmask.isFromAgl()/*&&mass>0&&mass<3*/) { 
		RICH_DvsMC_P -> DataEff -> afterAgl -> Fill(Kbin,zona); 
		RICH_DvsMC_D -> DataEff -> afterAgl -> Fill(Kbin,zona); 
	}
	return;
}

void DVSMCRICHeff_Fill(){

	if(Tup.Beta<=0||Tup.R<=0) return;
        if(!trgpatt.IsPhysical()) return;
        if(!Herejcut) return;
        if(!ProtonsMassThres) return;	//
	//if(!(Distcut&&Likcut)) return;        
	//
	int Kbin;
	float mass=0;

	if(Massa_gen<1) {
		//Beta bins
		//NaF
		Kbin=NaFPB.GetBin(RUsed);
		RICH_DvsMC_P -> MCEff -> beforeNaF -> Fill(Kbin,Tup.mcweight);
		for(int mc_type=0;mc_type<6;mc_type++) ((TH2*)RICH_DvsMC_D -> MCEff -> beforeNaF) -> Fill(Kbin,mc_type,Tup.mcweight);
		mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
		if(cmask.isFromNaF()/*&&mass>0&&mass<3*/){	
			RICH_DvsMC_P -> MCEff -> afterNaF -> Fill(Kbin,Tup.mcweight);
			for(int mc_type=0;mc_type<6;mc_type++) ((TH2*)RICH_DvsMC_D -> MCEff -> afterNaF) -> Fill(Kbin,mc_type,Tup.mcweight);
		}
		//Agl
    		Kbin=AglPB.GetBin(RUsed);
		RICH_DvsMC_P -> MCEff -> beforeAgl -> Fill(Kbin,Tup.mcweight);
		for(int mc_type=0;mc_type<6;mc_type++) ((TH2*)RICH_DvsMC_D -> MCEff -> beforeAgl) -> Fill(Kbin,mc_type,Tup.mcweight);	

		if(cmask.isFromAgl()/*&&mass>0&&mass<3*/) {
			RICH_DvsMC_P -> MCEff -> afterAgl -> Fill(Kbin,Tup.mcweight); 	
			for(int mc_type=0;mc_type<6;mc_type++) 	((TH2*)RICH_DvsMC_D -> MCEff -> afterAgl) -> Fill(Kbin,mc_type,Tup.mcweight);
		}
	}
}                        


void DVSMCRICHeff_Write(){

	RICH_DvsMC_P -> Write();
	RICH_DvsMC_D -> Write();

	return;
}


void DVSMCRICHeff(string filename){

	cout<<"******* Data vs MC: RICH ********"<<endl;
	
	cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

	DatavsMC * RICH_DvsMC_P = new DatavsMC(inputHistoFile,"RICH_DvsMC_P");
	DatavsMC * RICH_DvsMC_D = new DatavsMC(inputHistoFile,"RICH_DvsMC_D",6);


	LATcorr * LATrichDATA_TOF       = new LATcorr(inputHistoFile,"LATrichDATA_Agl" 	 ,"Results");
	LATcorr * LATrichDATA_NaF       = new LATcorr(inputHistoFile,"LATrichDATA_NaF" 	 ,"Results");
	LATcorr * LATrichDATA_Agl       = new LATcorr(inputHistoFile,"LATrichDATA_Agl" 	 ,"Results");	


	cout<<"******* Data vs MC: RICH ********"<<endl;

	RICH_DvsMC_P -> Assign_LatCorr( LATrichDATA_TOF   ->  LATcorrR_fit , 
			LATrichDATA_TOF   ->  LATcorrR_fit ,
			LATrichDATA_NaF   ->  LATcorrR_fit ,
			LATrichDATA_Agl   ->  LATcorrR_fit );

	RICH_DvsMC_D -> Assign_LatCorr( LATrichDATA_TOF   ->  LATcorrR_fit , 
			LATrichDATA_TOF   ->  LATcorrR_fit ,
			LATrichDATA_NaF   ->  LATcorrR_fit ,
			LATrichDATA_Agl   ->  LATcorrR_fit );

	RICH_DvsMC_P ->Eval_DandMC_Eff();  
	RICH_DvsMC_D ->Eval_DandMC_Eff();

	RICH_DvsMC_P ->Eval_Corrections();
	RICH_DvsMC_D ->Eval_Corrections();


	TH1F* RICH_Correction_P_NaF =(TH1F*) RICH_DvsMC_P -> GetCorrection_NaF();
	TH1F* RICH_Correction_P_Agl =(TH1F*) RICH_DvsMC_P -> GetCorrection_Agl();

	TH2F* RICH_Correction_D_NaF =(TH2F*) RICH_DvsMC_D -> GetCorrection_NaF();
	TH2F* RICH_Correction_D_Agl =(TH2F*) RICH_DvsMC_D -> GetCorrection_Agl();


	RICH_Correction_P_NaF  -> SetName("RICH_DvsMC_P_CorrectionNaF");
	RICH_Correction_P_Agl  -> SetName("RICH_DvsMC_P_CorrectionAgl");

	RICH_Correction_D_NaF  -> SetName("RICH_DvsMC_D_CorrectionNaF");
	RICH_Correction_D_Agl  -> SetName("RICH_DvsMC_D_CorrectionAgl");

	finalHistos.Add(RICH_Correction_P_NaF);
       finalHistos.Add( RICH_Correction_P_Agl);
                             
       finalHistos.Add( RICH_Correction_D_NaF);
       finalHistos.Add( RICH_Correction_D_Agl);

	finalHistos.writeObjsInFolder("Results");

        cout<<"*** Plotting ...  ****"<<endl;

	DVSMCRICHeff_Plot (RICH_Correction_P_NaF,
                           RICH_Correction_P_Agl,
                                
                           RICH_Correction_D_NaF,
                           RICH_Correction_D_Agl
	);
	
	return;
}

