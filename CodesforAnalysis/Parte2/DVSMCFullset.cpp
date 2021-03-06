#include "PlottingFunctions/DVSMCFullset.h"



void DVSMCFullset(string filename){

        cout<<"******* Data vs MC: Efficiency Corrections Calculation********"<<endl;

        cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");



        DatavsMC * Dist_DvsMC_P = new DatavsMC(inputHistoFile,"Dist_DvsMC_P",1,20);
        DatavsMC * Lik_DvsMC_P  = new DatavsMC(inputHistoFile,"Lik_DvsMC_P" ,1,20);

        LATcorr * LATLikelihoodDATA_TOF = new LATcorr(inputHistoFile,"LATLikDATA_TOF"    ,"Results");
        LATcorr * LATDistanceDATA_TOF   = new LATcorr(inputHistoFile,"LATDistDATA_TOF"   ,"Results");

        LATcorr * LATLikelihoodDATA_NaF = new LATcorr(inputHistoFile,"LATLikDATA_NaF"    ,"Results");
        LATcorr * LATDistanceDATA_NaF   = new LATcorr(inputHistoFile,"LATDistDATA_NaF"   ,"Results");

        LATcorr * LATLikelihoodDATA_Agl = new LATcorr(inputHistoFile,"LATLikDATA_Agl"    ,"Results");
        LATcorr * LATDistanceDATA_Agl   = new LATcorr(inputHistoFile,"LATDistDATA_Agl"   ,"Results");


	DatavsMC * PreSel_DvsMC_P[3];

        PreSel_DvsMC_P[0]= new DatavsMC(inputHistoFile,"PreSel_DvsMC_P0_",1,3);
        PreSel_DvsMC_P[1]= new DatavsMC(inputHistoFile,"PreSel_DvsMC_P1_",1,3);
        PreSel_DvsMC_P[2]= new DatavsMC(inputHistoFile,"PreSel_DvsMC_P2_",1,3);

	LATcorr * LATpreSelDATA = new LATcorr(inputHistoFile,"LATpreSelDATA"      ,"Results");
		
	cout<<"****PRESELECTIONS****"<<endl;

	for(int sel=0;sel<3;sel++){
                   PreSel_DvsMC_P[sel] -> Assign_LatCorr( LATpreSelDATA   ->  LATcorrR_fit ,
                                   LATpreSelDATA   ->  LATcorrR_fit ,
                                   LATpreSelDATA   ->  LATcorrR_fit ,
                                   LATpreSelDATA   ->  LATcorrR_fit );

                   PreSel_DvsMC_P[sel] ->Initialize_SystError();
		   PreSel_DvsMC_P[sel] ->Eval_DandMC_Eff(sel);
                   PreSel_DvsMC_P[sel] ->Eval_Corrections();
			
	}



        cout<<"****QUALITY SEL *****"<<endl;

        Dist_DvsMC_P -> Assign_LatCorr( LATDistanceDATA_TOF   ->  LATcorrR_fit ,
                        LATDistanceDATA_TOF   ->  LATcorrR_fit ,
                        LATDistanceDATA_NaF   ->  LATcorrR_fit ,
                        LATDistanceDATA_Agl   ->  LATcorrR_fit );

        Lik_DvsMC_P  -> Assign_LatCorr( LATLikelihoodDATA_TOF ->  LATcorrR_fit ,
                        LATLikelihoodDATA_TOF ->  LATcorrR_fit ,
                        LATLikelihoodDATA_NaF ->  LATcorrR_fit ,
                        LATLikelihoodDATA_Agl ->  LATcorrR_fit );


	Dist_DvsMC_P ->Initialize_SystError();
        Lik_DvsMC_P  ->Initialize_SystError();

        Dist_DvsMC_P ->Eval_DandMC_Eff();
        Lik_DvsMC_P  ->Eval_DandMC_Eff();

        Dist_DvsMC_P ->Eval_Corrections();
        Lik_DvsMC_P  ->Eval_Corrections();



	cout<<"******* Data vs MC: FULL SET ********"<<endl;
	
	DatavsMC * Presel_DvsMC_P = new DatavsMC(PreSel_DvsMC_P[0],"Presel_DvsMC_P");	
	
	for(int sel=1;sel<3;sel++) Presel_DvsMC_P -> ComposeCorrection(PreSel_DvsMC_P[sel],1,1);
	


	DatavsMC * Fullset_DvsMC_P = new DatavsMC(Dist_DvsMC_P,"Fullset_DvsMC_P");

	Fullset_DvsMC_P -> ComposeCorrection(Lik_DvsMC_P,10,10);
		
	for(int sel=0;sel<3;sel++) Fullset_DvsMC_P -> ComposeCorrection(PreSel_DvsMC_P[sel],10,1);


	TH1F* Fullset_Correction_R    =(TH1F*) Fullset_DvsMC_P -> GetCorrection_R(10)  ;
	TH1F* Fullset_Correction_TOF  =(TH1F*) Fullset_DvsMC_P -> GetCorrection_TOF(10);
	TH1F* Fullset_Correction_NaF  =(TH1F*) Fullset_DvsMC_P -> GetCorrection_NaF(10);
	TH1F* Fullset_Correction_Agl  =(TH1F*) Fullset_DvsMC_P -> GetCorrection_Agl(10);

	cout<<"*************SYST ERR**************"<<endl;

	Presel_DvsMC_P -> Eval_SystError();
	Fullset_DvsMC_P -> Eval_SystError();

	TH2F* SystFullset_R   =(TH2F*) Fullset_DvsMC_P -> GetSystPlot_R()  ;
        TH2F* SystFullset_TOF =(TH2F*) Fullset_DvsMC_P -> GetSystPlot_TOF();
        TH2F* SystFullset_NaF =(TH2F*) Fullset_DvsMC_P -> GetSystPlot_NaF();
        TH2F* SystFullset_Agl =(TH2F*) Fullset_DvsMC_P -> GetSystPlot_Agl();


	cout<<"************* FIT **************"<<endl;

	Presel_DvsMC_P -> Eval_FittedCorrections();	
	
	TH1F* Presel_CorrectionFit_R    =(TH1F*) Presel_DvsMC_P -> GetCorrection_R(1)  ;
	TH1F* Presel_CorrectionFit_TOF  =(TH1F*) Presel_DvsMC_P -> GetCorrection_TOF(1);
	TH1F* Presel_CorrectionFit_NaF  =(TH1F*) Presel_DvsMC_P -> GetCorrection_NaF(1);
	TH1F* Presel_CorrectionFit_Agl  =(TH1F*) Presel_DvsMC_P -> GetCorrection_Agl(1);

	TH1F* Presel_StatFit_R    =(TH1F*) Presel_DvsMC_P -> GetStatErr_R()  ; 
	TH1F* Presel_StatFit_TOF  =(TH1F*) Presel_DvsMC_P -> GetStatErr_TOF(); 
	TH1F* Presel_StatFit_NaF  =(TH1F*) Presel_DvsMC_P -> GetStatErr_NaF(); 
	TH1F* Presel_StatFit_Agl  =(TH1F*) Presel_DvsMC_P -> GetStatErr_Agl(); 

	TH1F* Presel_SystFit_R    =(TH1F*) Presel_DvsMC_P -> GetSystErr_R()  ; 
	TH1F* Presel_SystFit_TOF  =(TH1F*) Presel_DvsMC_P -> GetSystErr_TOF(); 
	TH1F* Presel_SystFit_NaF  =(TH1F*) Presel_DvsMC_P -> GetSystErr_NaF(); 
	TH1F* Presel_SystFit_Agl  =(TH1F*) Presel_DvsMC_P -> GetSystErr_Agl(); 

	Presel_CorrectionFit_R    -> SetName("Presel_DvsMC_P_CorrectionR"  );
	Presel_CorrectionFit_TOF  -> SetName("Presel_DvsMC_P_CorrectionTOF");
	Presel_CorrectionFit_NaF  -> SetName("Presel_DvsMC_P_CorrectionNaF");
	Presel_CorrectionFit_Agl  -> SetName("Presel_DvsMC_P_CorrectionAgl");

	Fullset_DvsMC_P -> Eval_FittedCorrections();

	TH1F* Fullset_CorrectionFit_R    =(TH1F*) Fullset_DvsMC_P -> GetCorrection_R(10)  ;
	TH1F* Fullset_CorrectionFit_TOF  =(TH1F*) Fullset_DvsMC_P -> GetCorrection_TOF(10);
	TH1F* Fullset_CorrectionFit_NaF  =(TH1F*) Fullset_DvsMC_P -> GetCorrection_NaF(10);
	TH1F* Fullset_CorrectionFit_Agl  =(TH1F*) Fullset_DvsMC_P -> GetCorrection_Agl(10);

	TH1F* Fullset_StatFit_R    =(TH1F*) Fullset_DvsMC_P -> GetStatErr_R()  ; 
	TH1F* Fullset_StatFit_TOF  =(TH1F*) Fullset_DvsMC_P -> GetStatErr_TOF(); 
	TH1F* Fullset_StatFit_NaF  =(TH1F*) Fullset_DvsMC_P -> GetStatErr_NaF(); 
	TH1F* Fullset_StatFit_Agl  =(TH1F*) Fullset_DvsMC_P -> GetStatErr_Agl(); 

	TH1F* Fullset_SystFit_R    =(TH1F*) Fullset_DvsMC_P -> GetSystErr_R()  ; 
	TH1F* Fullset_SystFit_TOF  =(TH1F*) Fullset_DvsMC_P -> GetSystErr_TOF(); 
	TH1F* Fullset_SystFit_NaF  =(TH1F*) Fullset_DvsMC_P -> GetSystErr_NaF(); 
	TH1F* Fullset_SystFit_Agl  =(TH1F*) Fullset_DvsMC_P -> GetSystErr_Agl(); 

	Fullset_CorrectionFit_R    -> SetName("Fullset_DvsMC_P_CorrectionR"  );
	Fullset_CorrectionFit_TOF  -> SetName("Fullset_DvsMC_P_CorrectionTOF");
	Fullset_CorrectionFit_NaF  -> SetName("Fullset_DvsMC_P_CorrectionNaF");
	Fullset_CorrectionFit_Agl  -> SetName("Fullset_DvsMC_P_CorrectionAgl");

	TH2F * SystPlot_R   = (TH2F *)Fullset_DvsMC_P -> GetSystPlot_R();
	TH2F * SystPlot_TOF = (TH2F *)Fullset_DvsMC_P -> GetSystPlot_TOF();
	TH2F * SystPlot_NaF = (TH2F *)Fullset_DvsMC_P -> GetSystPlot_NaF();
	TH2F * SystPlot_Agl = (TH2F *)Fullset_DvsMC_P -> GetSystPlot_Agl();
	
	

	finalHistos.Add(Presel_CorrectionFit_R   );
	finalHistos.Add(Presel_CorrectionFit_TOF );
	finalHistos.Add(Presel_CorrectionFit_NaF );
	finalHistos.Add(Presel_CorrectionFit_Agl );

	finalHistos.Add(Presel_StatFit_R   );
	finalHistos.Add(Presel_StatFit_TOF );
	finalHistos.Add(Presel_StatFit_NaF );
	finalHistos.Add(Presel_StatFit_Agl );
	
	finalHistos.Add(Presel_SystFit_R   );
	finalHistos.Add(Presel_SystFit_TOF );
	finalHistos.Add(Presel_SystFit_NaF );
	finalHistos.Add(Presel_SystFit_Agl );
	
	finalHistos.Add(Fullset_CorrectionFit_R   );
	finalHistos.Add(Fullset_CorrectionFit_TOF );
	finalHistos.Add(Fullset_CorrectionFit_NaF );
	finalHistos.Add(Fullset_CorrectionFit_Agl );

	finalHistos.Add(Fullset_StatFit_R   );
	finalHistos.Add(Fullset_StatFit_TOF );
	finalHistos.Add(Fullset_StatFit_NaF );
	finalHistos.Add(Fullset_StatFit_Agl );
	
	finalHistos.Add(Fullset_SystFit_R   );
	finalHistos.Add(Fullset_SystFit_TOF );
	finalHistos.Add(Fullset_SystFit_NaF );
	finalHistos.Add(Fullset_SystFit_Agl );
	
	finalHistos.Add(SystFullset_R   );
        finalHistos.Add(SystFullset_TOF );
        finalHistos.Add(SystFullset_NaF );
        finalHistos.Add(SystFullset_Agl );



	finalHistos.writeObjsInFolder("Results");

	cout<<"*** Plotting ...  ****"<<endl;

        DVSMCFullset_Plot(
		SystPlot_R,  	
                SystPlot_TOF,
                SystPlot_NaF,
		SystPlot_Agl	);
	



	return;
}
