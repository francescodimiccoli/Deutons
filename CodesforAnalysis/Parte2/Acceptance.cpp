#include "PlottingFunctions/Acceptance_Plot.h"

using namespace std;


void Acceptance(string filename){
	cout<<"****************** ACCEPTANCE CALCULATION ******************"<<endl;
	
	 cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


	ACCEPTANCE * AcceptanceP = new ACCEPTANCE (inputHistoFile,"Results","EffpreselMCP_Now","EffFullsetMCP","TOTLATCorr","CorrezioneLATp",1);
	ACCEPTANCE * AcceptanceD = new ACCEPTANCE (inputHistoFile,"Results","EffpreselMCD_Now","EffFullsetMCD","TOTLATCorr","CorrezioneLATd",6);
	
	ACCEPTANCE * AcceptancePreP = new ACCEPTANCE (inputHistoFile,"Results","EffpreselMCP_Now","EffCSCFullsetMCP","PreLATCorr","CorrezioneLATPrep",1);
	ACCEPTANCE * AcceptancePreD = new ACCEPTANCE (inputHistoFile,"Results","EffpreselMCD_Now","EffpreselMCD","PreLATCorr","CorrezioneLATPred",6);
       
	TH1F * TrackerGlobalFactor =  (TH1F *) inputHistoFile -> Get("Results/TrackerGlobalFactor");	
	TH1F * TriggerGlobalFactor =  (TH1F *) inputHistoFile -> Get("Results/TriggerGlobalFactor"); 
	TH1F * TriggerGlobalFactorQ = (TH1F *) inputHistoFile -> Get("Results/TriggerGlobalFactorQ");
	
	


	cout<<"****************** ACCEPTANCE CALCULATION ******************"<<endl;


	enum {protons, deutons};

 	AcceptanceP -> Set_MC_Par  (0.0308232619, 0.5, 100); 
 	AcceptanceP -> Set_Binning (false);

	AcceptanceD -> Set_MC_Par  (0.0242236931, 1, 20); 
	AcceptanceD -> Set_Binning (true);
	
	AcceptanceP-> Eval_Gen_Acceptance(1);
	AcceptanceP-> Eval_MC_Acceptance();
	AcceptanceP-> Eval_Geomag_Acceptance(1);
	AcceptanceP-> Eval_Corrected_Acceptance(1);

	AcceptanceD-> Eval_Gen_Acceptance(6);
	AcceptanceD-> Eval_MC_Acceptance();
	AcceptanceD-> Eval_Geomag_Acceptance(6);
	AcceptanceD-> Eval_Corrected_Acceptance(6);
	

	AcceptancePreP -> Set_MC_Par  (0.0308232619, 0.5, 100); 
 	AcceptancePreP -> Set_Binning (protons);

	AcceptancePreD -> Set_MC_Par  (0.0242236931, 0.5, 20); 
	AcceptancePreD -> Set_Binning (deutons);

	AcceptancePreP-> Eval_Gen_Acceptance(1);
	AcceptancePreP-> Eval_MC_Acceptance();
	AcceptancePreP-> Eval_Geomag_Acceptance(1);
	AcceptancePreP-> Eval_Corrected_Acceptance(1);

	AcceptancePreD-> Eval_Gen_Acceptance(6);
	AcceptancePreD-> Eval_MC_Acceptance();
	AcceptancePreD-> Eval_Geomag_Acceptance(6);	
	AcceptancePreD-> Eval_Corrected_Acceptance(6);	

	cout<<"****** DVSMC APPLICATION *********"<<endl;

	TH1F* PreSel_Correction_R[3];  
        TH1F* PreSel_Correction_TOF[3];
        TH1F* PreSel_Correction_NaF[3];
	TH1F* PreSel_Correction_Agl[3];


	TH1F* PreselP_Correction_R   =(TH1F*) inputHistoFile -> Get ("Results/Presel_DvsMC_P_CorrectionR"             );	

	TH1F* FullsetP_Correction_R   =(TH1F*) inputHistoFile -> Get ("Results/Fullset_DvsMC_P_CorrectionR"  		); 
	TH1F* FullsetP_Correction_TOF =(TH1F*) inputHistoFile -> Get ("Results/Fullset_DvsMC_P_CorrectionTOF"		);
	TH1F* FullsetP_Correction_NaF =(TH1F*) inputHistoFile -> Get ("Results/Fullset_DvsMC_P_CorrectionNaF"		);
	TH1F* FullsetP_Correction_Agl =(TH1F*) inputHistoFile -> Get ("Results/Fullset_DvsMC_P_CorrectionAgl"		);
                                                                                      
	TH1F* RICH_Correction_P_NaF =(TH1F*) inputHistoFile -> Get ("Results/RICH_DvsMC_P_CorrectionNaF"		);
	TH1F* RICH_Correction_P_Agl =(TH1F*) inputHistoFile -> Get ("Results/RICH_DvsMC_P_CorrectionAgl"		);
	
	TH2F* RICH_Correction_D_NaF =(TH2F*) inputHistoFile -> Get ("Results/RICH_DvsMC_D_CorrectionNaF"         );
        TH2F* RICH_Correction_D_Agl =(TH2F*) inputHistoFile -> Get ("Results/RICH_DvsMC_D_CorrectionAgl"         );


	//global factors
	
	AcceptanceP    -> ApplyGlobalFactor(TriggerGlobalFactorQ -> GetBinContent(1), TriggerGlobalFactorQ -> GetBinError(1));
	AcceptanceD    -> ApplyGlobalFactor(TriggerGlobalFactorQ -> GetBinContent(1), TriggerGlobalFactorQ -> GetBinError(1));
	AcceptancePreP -> ApplyGlobalFactor(TriggerGlobalFactorQ -> GetBinContent(1), TriggerGlobalFactorQ  -> GetBinError(1));
 
	AcceptanceP    -> ApplyGlobalFactor(TrackerGlobalFactor -> GetBinContent(1), TrackerGlobalFactor -> GetBinError(1));
	AcceptanceD    -> ApplyGlobalFactor(TrackerGlobalFactor -> GetBinContent(1), TrackerGlobalFactor -> GetBinError(1));
	AcceptancePreP -> ApplyGlobalFactor(TrackerGlobalFactor -> GetBinContent(1), TrackerGlobalFactor -> GetBinError(1));


	//pre
	AcceptancePreP -> Apply_DvsMCcorrection_R  (PreselP_Correction_R);	
	
	//qual
	AcceptanceP -> Apply_DvsMCcorrection_R  (FullsetP_Correction_R);
	AcceptanceP -> Apply_DvsMCcorrection_TOF(FullsetP_Correction_TOF);
        AcceptanceP -> Apply_DvsMCcorrection_NaF(FullsetP_Correction_NaF);
        AcceptanceP -> Apply_DvsMCcorrection_Agl(FullsetP_Correction_Agl);
/*	
	AcceptanceD -> Apply_DvsMCcorrection_R  (FullsetP_Correction_R,6);
	AcceptanceD -> Apply_DvsMCcorrection_TOF(FullsetP_Correction_TOF,6);
        AcceptanceD -> Apply_DvsMCcorrection_NaF(FullsetP_Correction_NaF,6);
        AcceptanceD -> Apply_DvsMCcorrection_Agl(FullsetP_Correction_Agl,6);
*/		

	//rich
	AcceptanceP -> Apply_DvsMCcorrection_NaF(RICH_Correction_P_NaF);
	AcceptanceP -> Apply_DvsMCcorrection_Agl(RICH_Correction_P_Agl);
	

	AcceptanceD -> Apply_DvsMCcorrection_NaF(RICH_Correction_D_NaF,6);
        AcceptanceD -> Apply_DvsMCcorrection_Agl(RICH_Correction_D_Agl,6);	
	
	//Protons
	AcceptanceP   ->Gen_Acceptance_R  ->SetName("Gen_AcceptanceP_R"  ); 
	AcceptanceP   ->Gen_Acceptance_TOF->SetName("Gen_AcceptanceP_TOF");
	AcceptanceP   ->Gen_Acceptance_NaF->SetName("Gen_AcceptanceP_NaF"); 
	AcceptanceP   ->Gen_Acceptance_Agl->SetName("Gen_AcceptanceP_Agl");

	AcceptancePreP->MCAcceptance_R  ->SetName("MC_AcceptancePreP_R"  );
	AcceptancePreD->MCAcceptance_R  ->SetName("MC_AcceptancePreD_R"  );
	AcceptancePreP->CorrectedAcceptance_R  ->SetName("Corr_AcceptancePreP_R"  );
	AcceptancePreD->CorrectedAcceptance_R  ->SetName("Corr_AcceptancePreD_R"  );

	AcceptanceP   ->MCAcceptance_R  ->SetName("MC_AcceptanceP_R"  );       
	AcceptanceP   ->MCAcceptance_TOF->SetName("MC_AcceptanceP_TOF");
	AcceptanceP   ->MCAcceptance_NaF->SetName("MC_AcceptanceP_NaF");
	AcceptanceP   ->MCAcceptance_Agl->SetName("MC_AcceptanceP_Agl");
	
	AcceptanceP   ->Geomag_Acceptance_R  ->SetName("Geomag_AcceptanceP_R"  );
	AcceptanceP   ->Geomag_Acceptance_TOF->SetName("Geomag_AcceptanceP_TOF");
	AcceptanceP   ->Geomag_Acceptance_NaF->SetName("Geomag_AcceptanceP_NaF");
	AcceptanceP   ->Geomag_Acceptance_Agl->SetName("Geomag_AcceptanceP_Agl");

	AcceptanceP   ->CorrectedAcceptance_R    ->SetName("Corr_AcceptanceP_R"  );
	AcceptanceP   ->CorrectedAcceptance_TOF  ->SetName("Corr_AcceptanceP_TOF");	
	AcceptanceP   ->CorrectedAcceptance_NaF  ->SetName("Corr_AcceptanceP_NaF");
	AcceptanceP   ->CorrectedAcceptance_Agl  ->SetName("Corr_AcceptanceP_Agl");

	
	//deutons
	AcceptanceD   ->Gen_Acceptance_R  ->SetName("Gen_AcceptanceD_R"  ); 
	AcceptanceD   ->Gen_Acceptance_TOF->SetName("Gen_AcceptanceD_TOF");
	AcceptanceD   ->Gen_Acceptance_NaF->SetName("Gen_AcceptanceD_NaF");      	
	AcceptanceD   ->Gen_Acceptance_Agl->SetName("Gen_AcceptanceD_Agl");

	AcceptanceD   ->MCAcceptance_R  ->SetName("MC_AcceptanceD_R"  );       
	AcceptanceD   ->MCAcceptance_TOF->SetName("MC_AcceptanceD_TOF");
	AcceptanceD   ->MCAcceptance_NaF->SetName("MC_AcceptanceD_NaF");
	AcceptanceD   ->MCAcceptance_Agl->SetName("MC_AcceptanceD_Agl");

	AcceptanceD   ->Geomag_Acceptance_R  ->SetName("Geomag_AcceptanceD_R"  );
	AcceptanceD   ->Geomag_Acceptance_TOF->SetName("Geomag_AcceptanceD_TOF");
	AcceptanceD   ->Geomag_Acceptance_NaF->SetName("Geomag_AcceptanceD_NaF");
	AcceptanceD   ->Geomag_Acceptance_Agl->SetName("Geomag_AcceptanceD_Agl");	

	AcceptanceD   ->CorrectedAcceptance_R    ->SetName("Corr_AcceptanceD_R"  );
	AcceptanceD   ->CorrectedAcceptance_TOF  ->SetName("Corr_AcceptanceD_TOF");
	AcceptanceD   ->CorrectedAcceptance_NaF  ->SetName("Corr_AcceptanceD_NaF");
	AcceptanceD   ->CorrectedAcceptance_Agl  ->SetName("Corr_AcceptanceD_Agl");

	
	cout<<"*************** ERRORS EXTRACTION *******************"<<endl;

	TH1F * LatCorrErrR_P   = (TH1F*)SetErrors((TH1*)AcceptanceP->LATcorrW_R);
	TH1F * LatCorrErrTOF_P = (TH1F*)SetErrors((TH1*)AcceptanceP->LATcorrW_TOF);
	TH1F * LatCorrErrNaF_P = (TH1F*)SetErrors((TH1*)AcceptanceP->LATcorrW_NaF);
	TH1F * LatCorrErrAgl_P = (TH1F*)SetErrors((TH1*)AcceptanceP->LATcorrW_Agl);
	
	TH1F * EffStatErrR_P   = (TH1F*)inputHistoFile -> Get ("Results/EffFullsetMCP_statR"  		);;
	TH1F * EffStatErrTOF_P = (TH1F*)inputHistoFile -> Get ("Results/EffFullsetMCP_statTOF"		);;
	TH1F * EffStatErrNaF_P = (TH1F*)inputHistoFile -> Get ("Results/EffFullsetMCP_statNaF"		);;
	TH1F * EffStatErrAgl_P = (TH1F*)inputHistoFile -> Get ("Results/EffFullsetMCP_statAgl"		);;

	TH1F * EffSystErrR_P   = (TH1F*)inputHistoFile -> Get ("Results/EffFullsetMCP_systR"  		);
	TH1F * EffSystErrTOF_P = (TH1F*)inputHistoFile -> Get ("Results/EffFullsetMCP_systTOF"		);
	TH1F * EffSystErrNaF_P = (TH1F*)inputHistoFile -> Get ("Results/EffFullsetMCP_systNaF"		);
	TH1F * EffSystErrAgl_P = (TH1F*)inputHistoFile -> Get ("Results/EffFullsetMCP_systAgl"		);

	TH1F * CorrStatErrR_P  = (TH1F*)inputHistoFile -> Get ("Results/Fullset_DvsMC_P_statR"  		);;
	TH1F * CorrStatErrTOF_P= (TH1F*)inputHistoFile -> Get ("Results/Fullset_DvsMC_P_statTOF"		);;
	TH1F * CorrStatErrNaF_P= (TH1F*)inputHistoFile -> Get ("Results/Fullset_DvsMC_P_statNaF"		);;
	TH1F * CorrStatErrAgl_P= (TH1F*)inputHistoFile -> Get ("Results/Fullset_DvsMC_P_statAgl"		);;

	TH1F * CorrSystErrR_P  = (TH1F*)inputHistoFile -> Get ("Results/Fullset_DvsMC_P_systR"  		);;
	TH1F * CorrSystErrTOF_P= (TH1F*)inputHistoFile -> Get ("Results/Fullset_DvsMC_P_systTOF"		);;
	TH1F * CorrSystErrNaF_P= (TH1F*)inputHistoFile -> Get ("Results/Fullset_DvsMC_P_systNaF"		);;
	TH1F * CorrSystErrAgl_P= (TH1F*)inputHistoFile -> Get ("Results/Fullset_DvsMC_P_systAgl"		);;

	TH1F * RICHStatErrNaF_P= (TH1F*)inputHistoFile -> Get ("Results/RICH_DvsMC_P_statNaF"		);;
	TH1F * RICHStatErrAgl_P= (TH1F*)inputHistoFile -> Get ("Results/RICH_DvsMC_P_statAgl"		);;

	TH1F * RICHSystErrNaF_P= (TH1F*)inputHistoFile -> Get ("Results/RICH_DvsMC_P_systNaF"		);;
	TH1F * RICHSystErrAgl_P= (TH1F*)inputHistoFile -> Get ("Results/RICH_DvsMC_P_systAgl"		);;



	TH2F * EffStatErrR_D2   = (TH2F*)inputHistoFile -> Get ("Results/EffFullsetMCD_statR"  		);;
	TH2F * EffStatErrTOF_D2 = (TH2F*)inputHistoFile -> Get ("Results/EffFullsetMCD_statTOF"		);;
	TH2F * EffStatErrNaF_D2 = (TH2F*)inputHistoFile -> Get ("Results/EffFullsetMCD_statNaF"		);;
	TH2F * EffStatErrAgl_D2 = (TH2F*)inputHistoFile -> Get ("Results/EffFullsetMCD_statAgl"		);;

	TH2F * EffSystErrR_D2   = (TH2F*)inputHistoFile -> Get ("Results/EffFullsetMCD_systR"  		);
	TH2F * EffSystErrTOF_D2 = (TH2F*)inputHistoFile -> Get ("Results/EffFullsetMCD_systTOF"		);
	TH2F * EffSystErrNaF_D2 = (TH2F*)inputHistoFile -> Get ("Results/EffFullsetMCD_systNaF"		);
	TH2F * EffSystErrAgl_D2 = (TH2F*)inputHistoFile -> Get ("Results/EffFullsetMCD_systAgl"		);




	TH1F * EffStatErrR_D   =ProjectionXtoTH1F(EffStatErrR_D2,  "EffFullsetMCD_statR",1,1);
        TH1F * EffStatErrTOF_D =ProjectionXtoTH1F(EffStatErrTOF_D2,"EffFullsetMCD_statTOF",1,1);
        TH1F * EffStatErrNaF_D =ProjectionXtoTH1F(EffStatErrNaF_D2,"EffFullsetMCD_statNaF",1,1);
        TH1F * EffStatErrAgl_D =ProjectionXtoTH1F(EffStatErrAgl_D2,"EffFullsetMCD_statAgl",1,1);
                                 
        TH1F * EffSystErrR_D   =ProjectionXtoTH1F(EffSystErrR_D2,  "EffFullsetMCD_systR",1,1);
        TH1F * EffSystErrTOF_D =ProjectionXtoTH1F(EffSystErrTOF_D2,"EffFullsetMCD_systTOF",1,1);
        TH1F * EffSystErrNaF_D =ProjectionXtoTH1F(EffSystErrNaF_D2,"EffFullsetMCD_systNaF",1,1);
        TH1F * EffSystErrAgl_D =ProjectionXtoTH1F(EffSystErrAgl_D2,"EffFullsetMCD_systAgl",1,1);




	cout<<RICHStatErrNaF_P<<endl;
        cout<<RICHStatErrAgl_P<<endl;
                        
        cout<<RICHSystErrNaF_P<<endl;
        cout<<RICHSystErrAgl_P<<endl;

	finalHistos.Add(AcceptanceP   ->	Gen_Acceptance_R   	);
	finalHistos.Add(AcceptanceP   ->	Gen_Acceptance_TOF	);
	finalHistos.Add(AcceptanceP   ->	Gen_Acceptance_NaF 	);
	finalHistos.Add(AcceptanceP   ->	Gen_Acceptance_Agl	);
	finalHistos.Add(AcceptancePreP->	MCAcceptance_R  	);
	finalHistos.Add(AcceptancePreD->	MCAcceptance_R  	);
	finalHistos.Add(AcceptancePreP->	CorrectedAcceptance_R 	);
	finalHistos.Add(AcceptancePreD->	CorrectedAcceptance_R 	);
	finalHistos.Add(AcceptanceP   ->	MCAcceptance_R        	); 
	finalHistos.Add(AcceptanceP   ->	MCAcceptance_TOF	);
	finalHistos.Add(AcceptanceP   ->	MCAcceptance_NaF	);
	finalHistos.Add(AcceptanceP   ->	MCAcceptance_Agl	);
	finalHistos.Add(AcceptanceP   ->	Geomag_Acceptance_R 	);
	finalHistos.Add(AcceptanceP   ->	Geomag_Acceptance_TOF	);
	finalHistos.Add(AcceptanceP   ->	Geomag_Acceptance_NaF	);
	finalHistos.Add(AcceptanceP   ->	Geomag_Acceptance_Agl	);
	finalHistos.Add(AcceptanceP   ->	CorrectedAcceptance_R   );
	finalHistos.Add(AcceptanceP   ->	CorrectedAcceptance_TOF	);
	finalHistos.Add(AcceptanceP   ->	CorrectedAcceptance_NaF	);
	finalHistos.Add(AcceptanceP   ->	CorrectedAcceptance_Agl	);
	finalHistos.Add(AcceptanceD   ->	Gen_Acceptance_R  	);
	finalHistos.Add(AcceptanceD   ->	Gen_Acceptance_TOF	);
	finalHistos.Add(AcceptanceD   ->	Gen_Acceptance_NaF     	);	
	finalHistos.Add(AcceptanceD   ->	Gen_Acceptance_Agl	);
	finalHistos.Add(AcceptanceD   ->	MCAcceptance_R         	);
	finalHistos.Add(AcceptanceD   ->	MCAcceptance_TOF	);
	finalHistos.Add(AcceptanceD   ->	MCAcceptance_NaF	);
	finalHistos.Add(AcceptanceD   ->	MCAcceptance_Agl	);
	finalHistos.Add(AcceptanceD   ->	Geomag_Acceptance_R  	);
	finalHistos.Add(AcceptanceD   ->	Geomag_Acceptance_TOF	);
	finalHistos.Add(AcceptanceD   ->	Geomag_Acceptance_NaF	);
	finalHistos.Add(AcceptanceD   ->	Geomag_Acceptance_Agl	);
	finalHistos.Add(AcceptanceD   ->	CorrectedAcceptance_R  	);
	finalHistos.Add(AcceptanceD   ->	CorrectedAcceptance_TOF	);
	finalHistos.Add(AcceptanceD   ->	CorrectedAcceptance_NaF	);
	finalHistos.Add(AcceptanceD   ->	CorrectedAcceptance_Agl	);


	finalHistos.writeObjsInFolder("Results");
        cout<<"*** Plotting ...  ****"<<endl;

	Acceptance_Plot (AcceptanceP   ->Gen_Acceptance_R   	,
                         AcceptanceP   ->Gen_Acceptance_TOF	,
                         AcceptanceP   ->Gen_Acceptance_NaF 	,
                         AcceptanceP   ->Gen_Acceptance_Agl	,
                         AcceptancePreP->MCAcceptance_R  	,
                         AcceptancePreP->Geomag_Acceptance_R 	,
                         AcceptanceP   ->MCAcceptance_R        	,
                         AcceptanceP   ->MCAcceptance_TOF	,
                         AcceptanceP   ->MCAcceptance_NaF	,
                         AcceptanceP   ->MCAcceptance_Agl	,
                         AcceptanceP   ->Geomag_Acceptance_R 	,
                         AcceptanceP   ->Geomag_Acceptance_TOF	,
                         AcceptanceP   ->Geomag_Acceptance_NaF	,
                         AcceptanceP   ->Geomag_Acceptance_Agl	,
                         AcceptanceP   ->CorrectedAcceptance_R   ,
                         AcceptanceP   ->CorrectedAcceptance_TOF,	
                         AcceptanceP   ->CorrectedAcceptance_NaF,	
                         AcceptanceP   ->CorrectedAcceptance_Agl,	
                         AcceptanceD   ->Gen_Acceptance_R  	,
                         AcceptanceD   ->Gen_Acceptance_TOF	,
                         AcceptanceD   ->Gen_Acceptance_NaF     ,	
                         AcceptanceD   ->Gen_Acceptance_Agl	,
                         AcceptanceD   ->MCAcceptance_R         ,	
                         AcceptanceD   ->MCAcceptance_TOF	,
                         AcceptanceD   ->MCAcceptance_NaF	,
                         AcceptanceD   ->MCAcceptance_Agl	,
                         AcceptanceD   ->Geomag_Acceptance_R  	,
                         AcceptanceD   ->Geomag_Acceptance_TOF	,
                         AcceptanceD   ->Geomag_Acceptance_NaF	,
                         AcceptanceD   ->Geomag_Acceptance_Agl	,
                         AcceptanceD   ->CorrectedAcceptance_R  ,	
                         AcceptanceD   ->CorrectedAcceptance_TOF,	
                         AcceptanceD   ->CorrectedAcceptance_NaF,	
	                 AcceptanceD   ->CorrectedAcceptance_Agl,
			 LatCorrErrR_P   		,
                         LatCorrErrTOF_P ,
                         LatCorrErrNaF_P, 
                         LatCorrErrAgl_P ,
                                         
                         EffStatErrR_P   ,
                         EffStatErrTOF_P ,
                         EffStatErrNaF_P ,
                         EffStatErrAgl_P ,
                                         
                         EffSystErrR_P   ,
                         EffSystErrTOF_P ,
                         EffSystErrNaF_P ,
                         EffSystErrAgl_P ,
                                         
                         CorrStatErrR_P  ,
                         CorrStatErrTOF_P,
                         CorrStatErrNaF_P,
                         CorrStatErrAgl_P,
                                         
                         CorrSystErrR_P  ,
                         CorrSystErrTOF_P,
                         CorrSystErrNaF_P,
                         CorrSystErrAgl_P,
			
			RICHStatErrNaF_P,
                        RICHStatErrAgl_P,
                                        
                        RICHSystErrNaF_P,
                        RICHSystErrAgl_P,


			EffStatErrR_D  , 
                        EffStatErrTOF_D,
                        EffStatErrNaF_D, 
                        EffStatErrAgl_D, 
                                        
                        EffSystErrR_D  , 
                        EffSystErrTOF_D, 
                        EffSystErrNaF_D, 
                        EffSystErrAgl_D 


);	
	return;
}


