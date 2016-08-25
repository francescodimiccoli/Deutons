#include "PlottingFunctions/DVSMCTrackeff_Plot.h"

using namespace std;

TH2F * ECALvsR_D=new TH2F("ECALvsR_D","ECALvsR_D",1000,0,100,1000,0,100);
TH2F * ECALvsR_MC=new TH2F("ECALvsR_MC","ECALvsR_MC",1000,0,100,1000,0,100);

Efficiency * TrackerEfficiencyMCP = new Efficiency("TrackerEfficiencyMCP");
Efficiency * TrackerEfficiencyD   = new Efficiency("TrackerEfficiencyC"  );

void DVSMCTrackeff_D_Fill(){
	
	//cuts
	if(!trgpatt.IsPhysical()) return;
	if(!(Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+1&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-1)) return;

	if(cmask.isPreselected()&&Tup.EdepECAL>1)
		ECALvsR_D->Fill(Tup.R_pre,Tup.EdepECAL);
        //R bins
	int Kbin=PRB.GetRBin (20) ;
	if(((int) Tup.Cutmask&3 ) == 3   )             TrackerEfficiencyD -> beforeR -> Fill(Kbin);
	if(((int) Tup.Cutmask&11) == 11  ) 	       TrackerEfficiencyD -> afterR  -> Fill(Kbin); 	

	return;
}


void DVSMCTrackeff_Fill(){
        //cuts
	if(!trgpatt.IsPhysical()||Tup.Beta_pre<=0) return;
	if(!(Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+1&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-1)) return;
	if(!(Massa_gen<1&&Massa_gen>0.5)) return;
	if(cmask.isPreselected()&&Tup.EdepECAL>1)
		ECALvsR_MC->Fill(Tup.R_pre,Tup.EdepECAL);	
        //R bins
        int Kbin=PRB.GetRBin (20) ;
	if(((int) Tup.Cutmask&3)  == 3   ) 	      TrackerEfficiencyMCP -> beforeR -> Fill(Kbin,Tup.mcweight);
        if(((int) Tup.Cutmask&11) == 11  )	      TrackerEfficiencyMCP -> afterR  -> Fill(Kbin,Tup.mcweight);
        
	return;

}


void DVSMCTrackeff_Write(){
        ECALvsR_D->Write(); 
        ECALvsR_MC->Write();
        TrackerEfficiencyMCP->Write();
	TrackerEfficiencyD  ->Write();
	return;
}


void Set_GlobalDatavsMCCorr(TH1F * Correction, TH1F * DataEff, TH1F * MCEff){
	for(int i =0; i< Correction -> GetNbinsX(); i++) {
		Correction -> SetBinContent(i,DataEff -> GetBinContent(PRB.GetRBin(20)+1)/(float)MCEff -> GetBinContent(PRB.GetRBin(20)+1) );
		Correction -> SetBinError(i,DataEff -> GetBinError(PRB.GetRBin(20)+1));
	}
	return;
}

void DVSMCTrackeff(string filename){

	cout<<"*************** Tracker Eff: Data vs MC ***************"<<endl;
	 TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


	Efficiency * TrackerEfficiencyMCP = new Efficiency(inputHistoFile,"TrackerEfficiencyMCP");
	Efficiency * TrackerEfficiencyD   = new Efficiency(inputHistoFile,"TrackerEfficiencyC"  );

	TH2F * ECALvsR_D =(TH2F*) inputHistoFile->Get("ECALvsR_D");
        TH2F * ECALvsR_MC =(TH2F*) inputHistoFile->Get("ECALvsR_MC");

	cout<<"*************** Tracker Eff: Data vs MC ***************"<<endl;	
 	TrackerEfficiencyMCP  -> Eval_Efficiency();
   	TrackerEfficiencyD    -> Eval_Efficiency();

	TH1F * TrackerEfficiencyMC    = (TH1F *)TrackerEfficiencyMCP  -> effR    ->    Clone();
	TH1F * TrackerEfficiencyData  = (TH1F *)TrackerEfficiencyD    -> effR    ->    Clone();		

	TH1F * TrackerGlobalFactor = new TH1F ("TrackerGlobalFactor","TrackerGlobalFactor",1,0,1);
	TrackerGlobalFactor -> SetBinContent(1,TrackerEfficiencyData -> GetBinContent(PRB.GetRBin(20)+1)/(float)TrackerEfficiencyMC -> GetBinContent(PRB.GetRBin(20)+1));
	TrackerGlobalFactor -> SetBinError(1,TrackerEfficiencyData -> GetBinError(PRB.GetRBin(20)+1));		
	
	TrackerEfficiencyData -> SetName("TrackerEfficiencyData");	


	finalHistos.Add(TrackerEfficiencyData);
	finalHistos.Add(TrackerEfficiencyMC);
	finalHistos.Add(TrackerGlobalFactor  );
	finalHistos.writeObjsInFolder("Results");

	cout<<"*** Plotting ...  ****"<<endl;

	DVSMCTrackeff_Plot(TrackerEfficiencyData,
			TrackerEfficiencyMC,
			TrackerGlobalFactor,
			ECALvsR_D, 
			ECALvsR_MC

			);

	return;
}

