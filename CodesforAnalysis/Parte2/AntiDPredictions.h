#include "PlottingFunctions/AntiDPrediction.h"
#include "AntiPFlux.h"

using namespace std;





void AntiDEfficiencies(string filename){
	cout<<"******************** ANTI-D EFFICIENCIES ************************"<<endl;
        cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

	Efficiency * EffpreselMCP = new Efficiency(inputHistoFile, "EffpreselMCP");
	Efficiency * EffpreselMCD = new Efficiency(inputHistoFile, "EffpreselMCD");

	Efficiency * EffLikMCP  = new Efficiency (inputHistoFile,"EffLikMCP");
	Efficiency * EffLikMCD  = new Efficiency (inputHistoFile,"EffLikMCD");
        
	OptimizationCut * DiscriminantCutTOF = new   OptimizationCut(inputHistoFile,"D_DiscrCUT_TOF","Results" );
        OptimizationCut * DiscriminantCutNaF = new   OptimizationCut(inputHistoFile,"D_DiscrCUT_NaF","Results" );
        OptimizationCut * DiscriminantCutAgl = new   OptimizationCut(inputHistoFile,"D_DiscrCUT_Agl","Results" );

        cout<<"******************** ANTI-D EFFICIENCIES ************************"<<endl;

	EffpreselMCP -> Eval_Efficiency();
	EffpreselMCD -> Eval_Efficiency();

	EffLikMCP   -> Eval_Efficiency();
	EffLikMCD   -> Eval_Efficiency();

	// preselection efficiencies
		
	TH1F * EffPreMCP_TH1F    =  (TH1F *)EffpreselMCP->effTOF->Clone();
	TH1F * EffPreMCPNaF_TH1F =  (TH1F *)EffpreselMCP->effNaF->Clone();
	TH1F * EffPreMCPAgl_TH1F =  (TH1F *)EffpreselMCP->effAgl->Clone();
	
	TH2F * EffPreMCD_TH2F    =  (TH2F *)EffpreselMCD->effTOF->Clone();
	TH2F * EffPreMCDNaF_TH2F =  (TH2F *)EffpreselMCD->effNaF->Clone();
	TH2F * EffPreMCDAgl_TH2F =  (TH2F *)EffpreselMCD->effAgl->Clone();

	
	//Likelihood efficiencies
	
	TH1F * EffFullSetP_TOF 	=(TH1F *)EffLikMCP ->effTOF->Clone();
   	TH2F * EffFullSetD_T_TOF=(TH2F *)EffLikMCD ->effTOF->Clone();
   	TH1F * EffFullSetP_NaF 	=(TH1F *)EffLikMCP ->effNaF->Clone();
   	TH2F * EffFullSetD_T_NaF=(TH2F *)EffLikMCD ->effNaF->Clone();
   	TH1F * EffFullSetP_Agl 	=(TH1F *)EffLikMCP ->effAgl->Clone();
   	TH2F * EffFullSetD_T_Agl=(TH2F *)EffLikMCD ->effAgl->Clone();

	EffFullSetP_TOF  -> Multiply( EffPreMCP_TH1F    ); 
	EffFullSetD_T_TOF-> Multiply( EffPreMCD_TH2F	); 
        EffFullSetP_NaF  -> Multiply( EffPreMCPNaF_TH1F ); 
        EffFullSetD_T_NaF-> Multiply( EffPreMCDNaF_TH2F ); 
        EffFullSetP_Agl  -> Multiply( EffPreMCPAgl_TH1F ); 
        EffFullSetD_T_Agl-> Multiply( EffPreMCDAgl_TH2F ); 
	
	// select only one MC_Type for deuterons	
	TH1F * EffFullSetD_TOF =  ProjectionXtoTH1F( (TH2F *)EffFullSetD_T_TOF,"EffFullSetD_TOF",2,3);
	TH1F * EffFullSetD_NaF =  ProjectionXtoTH1F( (TH2F *)EffFullSetD_T_NaF,"EffFullSetD_NaF",2,3);
	TH1F * EffFullSetD_Agl =  ProjectionXtoTH1F( (TH2F *)EffFullSetD_T_Agl,"EffFullSetD_Agl",2,3);
	
	// Final Distance cut
	DiscriminantCutTOF->Eval_Efficiencies();	
        DiscriminantCutNaF->Eval_Efficiencies();
        DiscriminantCutAgl->Eval_Efficiencies();

	TH1F * EffCutTOF_D = DiscriminantCutTOF-> GetEfficiency_D();
	TH1F * EffCutNaF_D = DiscriminantCutNaF-> GetEfficiency_D();
	TH1F * EffCutAgl_D = DiscriminantCutAgl-> GetEfficiency_D();

	TH1F * EffCutTOF_P = DiscriminantCutTOF-> GetEfficiency_P();
	TH1F * EffCutNaF_P = DiscriminantCutNaF-> GetEfficiency_P();
	TH1F * EffCutAgl_P = DiscriminantCutAgl-> GetEfficiency_P();

	EffFullSetP_TOF-> Multiply(	EffCutTOF_P);
	EffFullSetD_TOF-> Multiply(	EffCutTOF_D);
        EffFullSetP_NaF-> Multiply(	EffCutNaF_P);
        EffFullSetD_NaF-> Multiply(	EffCutNaF_D);
        EffFullSetP_Agl-> Multiply(	EffCutAgl_P);
        EffFullSetD_Agl-> Multiply(	EffCutAgl_D);

	

	EffFullSetP_TOF -> SetName ("AntiP_EffTOF"); 
        EffFullSetD_TOF -> SetName ("AntiD_EffTOF"); 
        EffFullSetP_NaF -> SetName ("AntiP_EffNaF"); 
        EffFullSetD_NaF -> SetName ("AntiD_EffNaF"); 
        EffFullSetP_Agl -> SetName ("AntiP_EffAgl"); 
        EffFullSetD_Agl -> SetName ("AntiD_EffAgl"); 

	finalHistos.Add(EffFullSetP_TOF);
	finalHistos.Add(EffFullSetD_TOF);
	finalHistos.Add(EffFullSetP_NaF);
	finalHistos.Add(EffFullSetD_NaF);
	finalHistos.Add(EffFullSetP_Agl);
	finalHistos.Add(EffFullSetD_Agl);

	finalHistos.writeObjsInFolder("Results");		

	return;
}

void AntiDpredictions(string filename){

	cout<<"******************** ANTI-D PREDICTIONS ************************"<<endl;
        cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

	ACCEPTANCE * AcceptanceAntiP = new ACCEPTANCE (inputHistoFile,"Results","EffpreselMCP","AntiP","TOTLATCorr","CorrezioneLATp",1);	
	ACCEPTANCE * AcceptanceAntiD = new ACCEPTANCE (inputHistoFile,"Results","EffpreselMCP","AntiD","TOTLATCorr","CorrezioneLATp",1);

	Tempi = (TH1F *)inputHistoFile->Get("Tempi");	

	TH2F * esposizionegeo_R    = (TH2F*)inputHistoFile->Get(      "esposizionegeo_R"      );
	TH2F * esposizionedgeoTOF  = (TH2F*)inputHistoFile->Get(      "esposizionedgeoTOF"    );
	TH2F * esposizionedgeoNaF  = (TH2F*)inputHistoFile->Get(      "esposizionedgeoNaF"    );
	TH2F * esposizionedgeoAgl  = (TH2F*)inputHistoFile->Get(      "esposizionedgeoAgl"    );
	TH2F * esposizionepgeoTOF  = (TH2F*)inputHistoFile->Get(      "esposizionepgeoTOF"    );
	TH2F * esposizionepgeoNaF  = (TH2F*)inputHistoFile->Get(      "esposizionepgeoNaF"    );
	TH2F * esposizionepgeoAgl  = (TH2F*)inputHistoFile->Get(      "esposizionepgeoAgl"    );



	cout<<"******************** ANTI-D PREDICTIONS ************************"<<endl;

	AcceptanceAntiD -> Set_MC_Par  (0.0242236931, 0.5, 20); 
	AcceptanceAntiD -> Set_Binning (deutons);
	
	AcceptanceAntiP -> Set_MC_Par  (0.0308232619, 0.5, 100); 
	AcceptanceAntiP -> Set_Binning (deutons);
	
	AcceptanceAntiP-> Eval_Gen_Acceptance(1);
	AcceptanceAntiP-> Eval_MC_Acceptance();

	AcceptanceAntiD-> Eval_Gen_Acceptance(1);
	AcceptanceAntiD-> Eval_MC_Acceptance();

	//Anti - D Acceptance
	
	AcceptanceAntiD   ->MCAcceptance_TOF->SetName("AcceptanceAntiD_TOF");
	AcceptanceAntiD   ->MCAcceptance_NaF->SetName("AcceptanceAntiD_NaF"); 
	AcceptanceAntiD   ->MCAcceptance_Agl->SetName("AcceptanceAntiD_Agl");


	TH1F * Acceptance_TOF = (TH1F*)AcceptanceAntiD   ->MCAcceptance_TOF -> Clone();
	TH1F * Acceptance_NaF = (TH1F*)AcceptanceAntiD   ->MCAcceptance_NaF -> Clone();
	TH1F * Acceptance_Agl = (TH1F*)AcceptanceAntiD   ->MCAcceptance_Agl -> Clone();		
	

	
	//Anti - P rejection power
	//
	TH1F * AcceptanceP_TOF = (TH1F*)AcceptanceAntiP   ->MCAcceptance_TOF -> Clone();
	TH1F * AcceptanceP_NaF = (TH1F*)AcceptanceAntiP   ->MCAcceptance_NaF -> Clone();
	TH1F * AcceptanceP_Agl = (TH1F*)AcceptanceAntiP   ->MCAcceptance_Agl -> Clone();		
	

	
	//Anti-D thresholds

	TH1F * ExposureD_TOF =  ProjectionXtoTH1F( (TH2F *)esposizionedgeoTOF,"ExposureD_TOF",0,11);
	TH1F * ExposureD_NaF =  ProjectionXtoTH1F( (TH2F *)esposizionedgeoNaF,"ExposureD_TOF",0,11);
	TH1F * ExposureD_Agl =  ProjectionXtoTH1F( (TH2F *)esposizionedgeoAgl,"ExposureD_TOF",0,11);

	 ExposureD_TOF ->Scale(1/Tempi->Integral());
	 ExposureD_NaF ->Scale(1/Tempi->Integral());
	 ExposureD_Agl ->Scale(1/Tempi->Integral());

	 ExposureD_TOF ->Scale(10*3.15*1e7);
	 ExposureD_NaF ->Scale(10*3.15*1e7);
	 ExposureD_Agl ->Scale(10*3.15*1e7);



	float TOF_Threshold = 0;
	for (int i=0;i<nbinsToF;i++) 
		if(Acceptance_TOF->GetBinContent(i+1)>0)  
			TOF_Threshold+=((ToFPB.EkPerMassBin(i+1)-ToFPB.EkPerMassBin(i))*Acceptance_TOF->GetBinContent(i+1)*0.8*ExposureD_TOF->GetBinContent(i+1));
	TOF_Threshold=1/TOF_Threshold;
	
	
	float NaF_Threshold = 0;
	for (int i=0;i<nbinsNaF;i++) 
			NaF_Threshold+=((NaFPB.EkPerMassBin(i+1)-NaFPB.EkPerMassBin(i))*Acceptance_NaF->GetBinContent(i+1)*0.8*ExposureD_TOF->GetBinContent(i+1));
	NaF_Threshold=1/NaF_Threshold;
	
	float Agl_Threshold = 0;
        for (int i=0;i<nbinsAgl;i++)
                	Agl_Threshold+=((AglPB.EkPerMassBin(i+1)-AglPB.EkPerMassBin(i))*Acceptance_Agl->GetBinContent(i+1)*0.8*ExposureD_TOF->GetBinContent(i+1));
	Agl_Threshold=1/Agl_Threshold;


	//Anti-P Contamination

	TH1F * ExposureP_TOF =  ProjectionXtoTH1F( (TH2F *)esposizionepgeoTOF,"ExposureP_TOF",0,11);
	TH1F * ExposureP_NaF =  ProjectionXtoTH1F( (TH2F *)esposizionepgeoNaF,"ExposureP_TOF",0,11);
	TH1F * ExposureP_Agl =  ProjectionXtoTH1F( (TH2F *)esposizionepgeoAgl,"ExposureP_TOF",0,11);

	 ExposureP_TOF ->Scale(1/Tempi->Integral());
	 ExposureP_NaF ->Scale(1/Tempi->Integral());
	 ExposureP_Agl ->Scale(1/Tempi->Integral());

	 ExposureP_TOF ->Scale(10*3.15*1e7);
	 ExposureP_NaF ->Scale(10*3.15*1e7);
	 ExposureP_Agl ->Scale(10*3.15*1e7);

	
	TSpline3 * AntiP_Flux  = new TSpline3("AntiP_Flux",AntiP_FluxX,AntiP_FluxY,24,0,120);	
	
	float TOF_Pbar_Threshold = 0;
	for (int i=0;i<nbinsToF;i++)
                if(Acceptance_TOF->GetBinContent(i+1)>0) 
			TOF_Pbar_Threshold +=  AntiP_Flux->Eval(ToFPB.EkPerMassBin(i))*AcceptanceP_TOF->GetBinContent(i+1)*((ToFPB.EkPerMassBin(i+1)-ToFPB.EkPerMassBin(i))*0.8*ExposureD_TOF->GetBinContent(i+1)); 	
			TOF_Pbar_Threshold*=TOF_Threshold;


	float NaF_Pbar_Threshold = 0;
	for (int i=0;i<nbinsNaF;i++)
                if(Acceptance_NaF->GetBinContent(i+1)>0) 
			NaF_Pbar_Threshold +=  AntiP_Flux->Eval(NaFPB.EkPerMassBin(i))*AcceptanceP_NaF->GetBinContent(i+1)*((ToFPB.EkPerMassBin(i+1)-ToFPB.EkPerMassBin(i))*0.8*ExposureD_NaF->GetBinContent(i+1)); 	
			NaF_Pbar_Threshold*=NaF_Threshold;


	float Agl_Pbar_Threshold = 0;
	for (int i=0;i<nbinsAgl;i++)
                if(Acceptance_Agl->GetBinContent(i+1)>0) 
			Agl_Pbar_Threshold +=  AntiP_Flux->Eval(AglPB.EkPerMassBin(i))*AcceptanceP_Agl->GetBinContent(i+1)*((ToFPB.EkPerMassBin(i+1)-ToFPB.EkPerMassBin(i))*0.8*ExposureD_Agl->GetBinContent(i+1)); 	
			Agl_Pbar_Threshold*=Agl_Threshold;





	cout<<"*** Plotting ...  ****"<<endl;
	AntiDpredictions_Plot(  Acceptance_TOF,
			  	Acceptance_NaF,  
				Acceptance_Agl,
			AcceptanceP_TOF ,	
			AcceptanceP_NaF ,
			AcceptanceP_Agl ,
			TOF_Threshold,
			NaF_Threshold,
			Agl_Threshold,
			TOF_Pbar_Threshold,
			NaF_Pbar_Threshold,
			Agl_Pbar_Threshold
			);


	return;




}








