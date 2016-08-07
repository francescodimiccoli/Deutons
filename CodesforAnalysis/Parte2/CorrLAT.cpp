#include "PlottingFunctions/CorrLAT_Plot.h"

TH1 * Weighted_CorrLAT(TH2F * esposizionegeo, TH1 * LATcorr);

void CorrLAT(string filename) {
	cout<<"*************** He control sample cut Efficiency on P*******************"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
	TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

	TH2F * esposizionegeo_R    = (TH2F*)inputHistoFile->Get( 	"esposizionegeo_R"       );
	TH2F * esposizionepgeoTOF  = (TH2F*)inputHistoFile->Get(	"esposizionepgeoTOF"	);
	TH2F * esposizionepgeoNaF  = (TH2F*)inputHistoFile->Get(	"esposizionepgeoNaF"	);
	TH2F * esposizionepgeoAgl  = (TH2F*)inputHistoFile->Get(	"esposizionepgeoAgl"	);
	TH2F * esposizionedgeoTOF  = (TH2F*)inputHistoFile->Get(	"esposizionedgeoTOF"	);
	TH2F * esposizionedgeoNaF  = (TH2F*)inputHistoFile->Get(	"esposizionedgeoNaF"	);
	TH2F * esposizionedgeoAgl  = (TH2F*)inputHistoFile->Get(	"esposizionedgeoAgl"	);

	LATcorr * LATpreSelDATA      	= new LATcorr(inputHistoFile,"LATpreSelDATA"  	 ,"Results");

	LATcorr * LATLikelihoodDATA_TOF = new LATcorr(inputHistoFile,"LATLikDATA_TOF"   	 ,"Results");
	LATcorr * LATDistanceDATA_TOF   = new LATcorr(inputHistoFile,"LATDistDATA_TOF" 	 ,"Results");

	LATcorr * LATLikelihoodDATA_NaF = new LATcorr(inputHistoFile,"LATLikDATA_NaF"  	 ,"Results");
	LATcorr * LATDistanceDATA_NaF   = new LATcorr(inputHistoFile,"LATDistDATA_NaF" 	 ,"Results");

	LATcorr * LATLikelihoodDATA_Agl = new LATcorr(inputHistoFile,"LATLikDATA_Agl"  	 ,"Results");
	LATcorr * LATDistanceDATA_Agl   = new LATcorr(inputHistoFile,"LATDistDATA_Agl" 	 ,"Results");

	LATcorr * LATrichDATA_NaF       = new LATcorr(inputHistoFile,"LATrichDATA_NaF" 	 ,"Results");
	LATcorr * LATrichDATA_Agl       = new LATcorr(inputHistoFile,"LATrichDATA_Agl" 	 ,"Results");	

	cout<<"******* TOTAL LAT. CORRECTION *************"<<endl;

	//product of LAT corr. of 3 preselections
	TH1F * PreLATCorr = (TH1F *)(( (TH2F *)LATpreSelDATA ->  LATcorrR_fit ) -> ProjectionX("",1,1)) -> Clone();
	PreLATCorr -> Multiply( (TH1F *)( ( (TH2F *)LATpreSelDATA ->  LATcorrR_fit ) -> ProjectionX("",2,2))->Clone()	);
	PreLATCorr -> Multiply( (TH1F *)( ( (TH2F *)LATpreSelDATA ->  LATcorrR_fit ) -> ProjectionX("",3,3))->Clone()	);

	TH1F *  PreLATCorrR   = (TH1F *) PreLATCorr -> Clone();
	TH1F *  PreLATCorrTOF = (TH1F *) PreLATCorr -> Clone();   
	TH1F *  PreLATCorrNaF = (TH1F *) PreLATCorr -> Clone();
	TH1F *  PreLATCorrAgl = (TH1F *) PreLATCorr -> Clone();

	TH1F *  TOTLATCorrR   = (TH1F *) PreLATCorr -> Clone();
	TH1F *  TOTLATCorrTOF = (TH1F *) PreLATCorr -> Clone();
	TH1F *  TOTLATCorrNaF = (TH1F *) PreLATCorr -> Clone();
	TH1F *  TOTLATCorrAgl = (TH1F *) PreLATCorr -> Clone();

	TOTLATCorrR    -> Multiply ( ( (TH1F *)LATLikelihoodDATA_TOF ->  LATcorrR_fit) );
	TOTLATCorrR    -> Multiply ( ( (TH1F *)LATDistanceDATA_TOF   ->  LATcorrR_fit) );

	TOTLATCorrTOF  -> Multiply ( ( (TH1F *)LATLikelihoodDATA_TOF ->  LATcorrR_fit) );
	TOTLATCorrTOF  -> Multiply ( ( (TH1F *)LATDistanceDATA_TOF   ->  LATcorrR_fit) );

	TOTLATCorrNaF  -> Multiply ( ( (TH1F *)LATLikelihoodDATA_NaF ->  LATcorrR_fit) );
	TOTLATCorrNaF  -> Multiply ( ( (TH1F *)LATDistanceDATA_NaF   ->  LATcorrR_fit) );
//	TOTLATCorrNaF  -> Multiply ( ( (TH1F *)LATrichDATA_NaF       ->  LATcorrR_fit) );   

	TOTLATCorrAgl  -> Multiply ( ( (TH1F *)LATLikelihoodDATA_Agl ->  LATcorrR_fit) );
	TOTLATCorrAgl  -> Multiply ( ( (TH1F *)LATDistanceDATA_Agl   ->  LATcorrR_fit) );
//	TOTLATCorrAgl  -> Multiply ( ( (TH1F *)LATrichDATA_Agl       ->  LATcorrR_fit) );

	//Only pres.
	TH1F * CorrezioneLATpre_pR = (TH1F *) Weighted_CorrLAT ( esposizionegeo_R , PreLATCorr);

	TH1F * CorrezioneLATpre_dR = (TH1F *) Weighted_CorrLAT ( esposizionegeo_R , PreLATCorr);

	//Full set
	TH1F * CorrezioneLAT_pR   = (TH1F *) Weighted_CorrLAT ( esposizionegeo_R   , TOTLATCorrTOF 	);
	TH1F * CorrezioneLAT_pTOF = (TH1F *) Weighted_CorrLAT ( esposizionepgeoTOF , TOTLATCorrTOF 	);
	TH1F * CorrezioneLAT_pNaF = (TH1F *) Weighted_CorrLAT ( esposizionepgeoNaF , TOTLATCorrNaF 	);
	TH1F * CorrezioneLAT_pAgl = (TH1F *) Weighted_CorrLAT ( esposizionepgeoAgl , TOTLATCorrAgl 	);

	TH1F * CorrezioneLAT_dR   = (TH1F *) Weighted_CorrLAT ( esposizionegeo_R   , TOTLATCorrTOF 	);
	TH1F * CorrezioneLAT_dTOF = (TH1F *) Weighted_CorrLAT ( esposizionedgeoTOF , TOTLATCorrTOF	);
	TH1F * CorrezioneLAT_dNaF = (TH1F *) Weighted_CorrLAT ( esposizionedgeoNaF , TOTLATCorrNaF	);
	TH1F * CorrezioneLAT_dAgl = (TH1F *) Weighted_CorrLAT ( esposizionedgeoAgl , TOTLATCorrAgl	);



	PreLATCorrR         -> SetName(     "PreLATCorr_LATcorrR_fit" 	);
	PreLATCorrTOF       -> SetName(     "PreLATCorr_LATcorrTOF_fit");
	PreLATCorrNaF       -> SetName(     "PreLATCorr_LATcorrNaF_fit");
	PreLATCorrAgl       -> SetName(     "PreLATCorr_LATcorrAgl_fit");

	TOTLATCorrR  	-> SetName(     "TOTLATCorr_LATcorrR_fit"  );
	TOTLATCorrTOF 	 -> SetName(     "TOTLATCorr_LATcorrTOF_fit");
	TOTLATCorrNaF        -> SetName(     "TOTLATCorr_LATcorrNaF_fit");
	TOTLATCorrAgl        -> SetName(     "TOTLATCorr_LATcorrAgl_fit");

	CorrezioneLATpre_pR -> SetName(  "CorrezioneLATPrep_R"	);
	CorrezioneLATpre_dR -> SetName(  "CorrezioneLATPred_R"	);

	CorrezioneLAT_pR    -> SetName(  "CorrezioneLATp_R"  	);
	CorrezioneLAT_pTOF  -> SetName(  "CorrezioneLATp_TOF"	);
	CorrezioneLAT_pNaF  -> SetName(  "CorrezioneLATp_NaF"       );
	CorrezioneLAT_pAgl  -> SetName(  "CorrezioneLATp_Agl"	);

	CorrezioneLAT_dR     -> SetName(  "CorrezioneLATd_R"  	);
	CorrezioneLAT_dTOF   -> SetName(  "CorrezioneLATd_TOF"	);
	CorrezioneLAT_dNaF   -> SetName(  "CorrezioneLATd_NaF"	);
	CorrezioneLAT_dAgl   -> SetName(  "CorrezioneLATd_Agl"       );


	finalHistos.Add(PreLATCorrR        );
	finalHistos.Add(PreLATCorrTOF      );
	finalHistos.Add(PreLATCorrNaF      );
	finalHistos.Add(PreLATCorrAgl      );

	finalHistos.Add(TOTLATCorrR        );
	finalHistos.Add(TOTLATCorrTOF      );
	finalHistos.Add(TOTLATCorrNaF      );
	finalHistos.Add(TOTLATCorrAgl      );

	finalHistos.Add(CorrezioneLATpre_pR);
	finalHistos.Add(CorrezioneLATpre_dR);

	finalHistos.Add(CorrezioneLAT_pR   );
	finalHistos.Add(CorrezioneLAT_pTOF );
	finalHistos.Add(CorrezioneLAT_pNaF );
	finalHistos.Add(CorrezioneLAT_pAgl );

	finalHistos.Add(CorrezioneLAT_dR   );
	finalHistos.Add(CorrezioneLAT_dTOF );
	finalHistos.Add(CorrezioneLAT_dNaF );
	finalHistos.Add(CorrezioneLAT_dAgl );



	finalHistos.writeObjsInFolder("Results");

	cout<<"*** Plotting ...  ****"<<endl;





	CorrLAT_Plot(TOTLATCorrTOF,
			PreLATCorr,
			CorrezioneLATpre_pR,
			CorrezioneLAT_pR,
			CorrezioneLAT_pTOF,
			CorrezioneLAT_dTOF,
			CorrezioneLAT_pNaF,
			CorrezioneLAT_dNaF,
			CorrezioneLAT_pAgl,
			CorrezioneLAT_dAgl			
		    );

	return;

}



TH1 * Weighted_CorrLAT(TH2F * esposizionegeo, TH1 * LATcorr) {
	TH2F * temp    = (TH2F *)esposizionegeo -> Clone();
	TH2F * temp_err= (TH2F *)esposizionegeo -> Clone(); 
	for(int m=0; m<11; m++) {
		for(int i=0; i< temp -> GetNbinsX(); i++) {
			temp    ->SetBinContent(i+1,m,esposizionegeo->GetBinContent(i+1,m)*LATcorr -> GetBinContent(m+1));
			temp_err->SetBinContent(i+1,m,pow(esposizionegeo->GetBinContent(i+1,m)*LATcorr -> GetBinError(m+1),2));
		}
	}
	//summ all over latitudes	
	TH1F * temp2     = ProjectionXtoTH1F(temp     , "",0,10);
	TH1F * temp2_err = ProjectionXtoTH1F(temp_err , "",0,10);

	for(int i=0; i< temp2_err -> GetNbinsX(); i++) temp2_err->SetBinContent(i+1,pow(temp2_err->GetBinContent(i+1),0.5));
	//	

	//Divide by total exposure time
	TH1F * Exptime =ProjectionXtoTH1F( esposizionegeo , "",0,10);
	temp2     -> Divide ( Exptime );
	temp2_err -> Divide ( Exptime );
	//

	//setting errors
	for(int i=0; i< temp2 -> GetNbinsX(); i++)	temp2 -> SetBinError(i+1,temp2_err -> GetBinContent(i+1));

	return (TH1 *)temp2 -> Clone();
}
