#include "PlottingFunctions/MCpreeff_Plot.h"

using namespace std;



Efficiency * EffpreselMCP = new Efficiency("EffpreselMCP");
Efficiency * EffpreselMCD = new Efficiency("EffpreselMCD", 6);

Efficiency * EffpreselMCP_L1 = new Efficiency("EffpreselMCP_L1");
Efficiency * EffpreselMCD_L1 = new Efficiency("EffpreselMCD_L1", 6);

Efficiency * EffpreselMCP_Now = new Efficiency("EffpreselMCP_Now");	//withouth reweighting, for acceptance
Efficiency * EffpreselMCD_Now = new Efficiency("EffpreselMCD_Now", 6);





Efficiency * EffCascadeMCP[6];
Efficiency * EffCascadeRICHMCP[4];

Efficiency * EffPreMCP = new Efficiency("EffPreMCP");
Efficiency * EffAglMCP = new Efficiency("EffAglMCP");
Efficiency * EffNaFMCP = new Efficiency("EffNaFMCP");

Efficiency * RInner = new Efficiency("RInner");
Efficiency * R_L1   = new Efficiency("R_L1");

void InizializeEff(){
	for(int i=0;i<6;i++) EffCascadeMCP[i]= new Efficiency(("EffCascadeMCP"+to_string(i)));
	for(int i=0;i<4;i++) EffCascadeRICHMCP[i]= new Efficiency(("EffCascadeRICHMCP"+to_string(i)));
}

void MCpreseff_Fill() {

	int kbin;
	if(Massa_gen<1&&Massa_gen>0.5) {
		//R bins
		EffpreselMCP->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
		EffpreselMCP_L1->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
		EffpreselMCP_Now->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),1);



		if(cmask.isPreselected()&&Tup.Beta_pre>0&&Tup.R_pre>0) {
			EffpreselMCP->afterR->Fill(PRB.GetRBin(RUsed),Tup.mcweight);
			if(Tup.R_L1>0) EffpreselMCP_L1->afterR->Fill(PRB.GetRBin(RUsed),Tup.mcweight);
		}

		// Beta bins
		EffpreselMCP->beforeTOF->Fill(ToFPB.GetBin(Tup.Momento_gen),Tup.mcweight);	
		EffpreselMCP->beforeNaF->Fill(NaFPB.GetBin(Tup.Momento_gen),Tup.mcweight);	
		EffpreselMCP->beforeAgl->Fill(AglPB.GetBin(Tup.Momento_gen),Tup.mcweight);	

		EffpreselMCP_L1->beforeTOF->Fill(ToFPB.GetBin(Tup.Momento_gen),Tup.mcweight);	
		EffpreselMCP_L1->beforeNaF->Fill(NaFPB.GetBin(Tup.Momento_gen),Tup.mcweight);	
		EffpreselMCP_L1->beforeAgl->Fill(AglPB.GetBin(Tup.Momento_gen),Tup.mcweight);	

		EffpreselMCP_Now->beforeTOF->Fill(ToFPB.GetBin(Tup.Momento_gen),1);	
		EffpreselMCP_Now->beforeNaF->Fill(NaFPB.GetBin(Tup.Momento_gen),1);	
		EffpreselMCP_Now->beforeAgl->Fill(AglPB.GetBin(Tup.Momento_gen),1);	




		if(cmask.isPreselected() && Tup.Beta_pre>0 && Tup.R_pre>0)
		{
			kbin=ToFPB.GetBin(RUsed);
			EffpreselMCP->afterTOF->Fill(kbin,Tup.mcweight); 
			if(Tup.R_L1>0) EffpreselMCP->afterTOF->Fill(kbin,Tup.mcweight);
			
			if(cmask.isFromNaF()) {
				kbin=NaFPB.GetBin(RUsed);	
				EffpreselMCP->afterNaF->Fill(kbin,Tup.mcweight); 
				if(Tup.R_L1>0) EffpreselMCP_L1->afterNaF->Fill(kbin,Tup.mcweight);
				
			} 
			if(cmask.isFromAgl()) {
				kbin=AglPB.GetBin(RUsed);   
				EffpreselMCP->afterAgl->Fill(kbin,Tup.mcweight);
				if(Tup.R_L1>0) EffpreselMCP_L1->afterAgl->Fill(kbin,Tup.mcweight);
			}  
		}

	}
	//Cascade Pres.
	if(Massa_gen<1&&Massa_gen>0.5&&trgpatt.IsPhysical()){
		EffCascadeMCP[0]->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
		EffCascadeMCP[1]->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
		EffCascadeMCP[2]->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
		EffCascadeMCP[3]->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
		EffCascadeMCP[4]->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
		EffCascadeMCP[5]->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
		
		if(cmask.isMinimumBiasToF3or4Layers()){
			EffCascadeMCP[0]->afterR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
			if(cmask.isMinimumBiasTracker()){
				EffCascadeMCP[1]->afterR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
				if(cmask.isGoldenToF3or4Layers()){
					EffCascadeMCP[2]->afterR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
					if(cmask.isGoldenTracker()){
						EffCascadeMCP[3]->afterR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
						if(cmask.hasSingleTrTrack()){
							EffCascadeMCP[4]->afterR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
							if(Tup.R_L1>0) EffCascadeMCP[5]->afterR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
						}	
					}
				}
			}
		}	
	}
	//Cascade RICH
	if(Massa_gen<1&&Massa_gen>0.5&&trgpatt.IsPhysical()&&!(cmask.isFromNaF())){
		EffCascadeRICHMCP[0]->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
		EffCascadeRICHMCP[1]->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
		EffCascadeRICHMCP[2]->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
		EffCascadeRICHMCP[3]->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
		
		if(cmask.hasRICHRing()){
			EffCascadeRICHMCP[0]->afterR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
			if((((int)Tup.Cutmask>>11)&4)==0){
				EffCascadeRICHMCP[1]->afterR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
				if((((int)Tup.Cutmask>>11)&24)==0){
					EffCascadeRICHMCP[2]->afterR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
					if((((int)Tup.Cutmask>>11)&480)==0){
						EffCascadeRICHMCP[3]->afterR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
					}
				}
			}
		}
	}

	//RICH eff.
	if(Massa_gen<1&&Massa_gen>0.5&&trgpatt.IsPhysical()){
		EffPreMCP->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
                EffNaFMCP->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
                EffAglMCP->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);	
		if(cmask.isPreselected()&&Tup.Beta_pre>0&&Tup.R_pre>0)	{
			EffPreMCP->afterR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
                	if(cmask.isFromNaF()) EffNaFMCP->afterR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
                	if(cmask.isFromAgl()) EffAglMCP->afterR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);	
		}
	}

	//Inner vs L1
	if(Massa_gen<1&&Massa_gen>0.5&&trgpatt.IsPhysical()){
		RInner->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
                R_L1->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
		if(Tup.R_pre>0) RInner->afterR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
		if(Tup.R_L1>0) R_L1->afterR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
	}


	if(Massa_gen>1&&Massa_gen<2) {
		// R bins      
		((TH2*)EffpreselMCD->beforeR) ->Fill( DRB.GetRBin(Tup.Momento_gen),ReturnMCGenType(),Tup.mcweight);
		((TH2*)EffpreselMCD_L1->beforeR) ->Fill( DRB.GetRBin(Tup.Momento_gen),ReturnMCGenType(),Tup.mcweight);
		((TH2*)EffpreselMCD_Now->beforeR) ->Fill( DRB.GetRBin(Tup.Momento_gen),ReturnMCGenType(),1);
	


	
		if(cmask.isPreselected()&&Tup.Beta_pre>0&&Tup.R_pre>0){
			((TH2*) EffpreselMCD->afterR) ->Fill( DRB.GetRBin(RUsed),ReturnMCGenType(),Tup.mcweight);
			if(Tup.R_L1>0) ((TH2*) EffpreselMCD_L1->afterR) ->Fill( DRB.GetRBin(RUsed),ReturnMCGenType(),Tup.mcweight);	
		}
		// Beta bins

		((TH2*)EffpreselMCD->beforeTOF)->Fill( ToFDB.GetBin(Tup.Momento_gen),ReturnMCGenType(),Tup.mcweight); 
		((TH2*)EffpreselMCD->beforeNaF)->Fill( NaFDB.GetBin(Tup.Momento_gen),ReturnMCGenType(),Tup.mcweight); 
		((TH2*)EffpreselMCD->beforeAgl)->Fill( AglDB.GetBin(Tup.Momento_gen),ReturnMCGenType(),Tup.mcweight); 

		((TH2*)EffpreselMCD_L1->beforeTOF)->Fill( ToFDB.GetBin(Tup.Momento_gen),ReturnMCGenType(),Tup.mcweight); 
		((TH2*)EffpreselMCD_L1->beforeNaF)->Fill( NaFDB.GetBin(Tup.Momento_gen),ReturnMCGenType(),Tup.mcweight); 
		((TH2*)EffpreselMCD_L1->beforeAgl)->Fill( AglDB.GetBin(Tup.Momento_gen),ReturnMCGenType(),Tup.mcweight); 

		((TH2*)EffpreselMCD_Now->beforeTOF)->Fill( ToFDB.GetBin(Tup.Momento_gen),ReturnMCGenType(),1); 
		((TH2*)EffpreselMCD_Now->beforeNaF)->Fill( NaFDB.GetBin(Tup.Momento_gen),ReturnMCGenType(),1); 
		((TH2*)EffpreselMCD_Now->beforeAgl)->Fill( AglDB.GetBin(Tup.Momento_gen),ReturnMCGenType(),1); 





		if(cmask.isPreselected() && Tup.Beta_pre>0 && Tup.R_pre>0)
		{
			kbin=ToFDB.GetBin(RUsed);
			((TH2*)EffpreselMCD->afterTOF)->Fill(kbin,ReturnMCGenType(),Tup.mcweight);
			if(Tup.R_L1>0) ((TH2*) EffpreselMCD_L1->afterTOF)->Fill(kbin,ReturnMCGenType(),Tup.mcweight);

			if(cmask.isFromNaF()) {
				kbin=NaFDB.GetBin(RUsed);
				((TH2*)EffpreselMCD->afterNaF)->Fill(kbin,ReturnMCGenType(),Tup.mcweight);
				if(Tup.R_L1>0) ((TH2*) EffpreselMCD_L1->afterNaF)->Fill(kbin,ReturnMCGenType(),Tup.mcweight);
			}
			if(cmask.isFromAgl()) {
				kbin=AglDB.GetBin(RUsed);
				((TH2*)EffpreselMCD->afterAgl)->Fill(kbin,ReturnMCGenType(),Tup.mcweight);
				if(Tup.R_L1>0) ((TH2*) EffpreselMCD_L1->afterAgl)->Fill(kbin,ReturnMCGenType(),Tup.mcweight);
			}
		}

	}
	return;
}


void MCpreeff_Write() {
	EffpreselMCP->Write();
	EffpreselMCD->Write();
	EffpreselMCP_L1->Write();
	EffpreselMCD_L1->Write();
	EffpreselMCP_Now->Write();
	EffpreselMCD_Now->Write();
	
	
	for(int i=0;i<6;i++) EffCascadeMCP[i]->Write();
	for(int i=0;i<4;i++) EffCascadeRICHMCP[i]->Write();
	
	EffPreMCP ->Write();
        EffAglMCP ->Write();
	EffNaFMCP ->Write(); 
	RInner 	->Write();
        R_L1    ->Write(); 

	return;
}



void MCpreeff(string filename) {

	cout<<"**** MC PRESELECTIONS EFFICIENCY (FULL SET) ****"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
	TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

	Efficiency * EffpreselMCP = new Efficiency(inputHistoFile, "EffpreselMCP");
	Efficiency * EffpreselMCD = new Efficiency(inputHistoFile, "EffpreselMCD");

	Efficiency * EffpreselMCP_L1 = new Efficiency(inputHistoFile, "EffpreselMCP_L1");
	Efficiency * EffpreselMCD_L1 = new Efficiency(inputHistoFile, "EffpreselMCD_L1");

	
	Efficiency * EffCascadeMCP[6];
	Efficiency * EffCascadeRICHMCP[4];

	for(int i=0;i<6;i++) EffCascadeMCP[i]=     new Efficiency(inputHistoFile,("EffCascadeMCP"+to_string(i)));
        for(int i=0;i<4;i++) EffCascadeRICHMCP[i]= new Efficiency(inputHistoFile,("EffCascadeRICHMCP"+to_string(i)));

	Efficiency * EffPreMCP = new Efficiency(inputHistoFile,"EffPreMCP");
	Efficiency * EffAglMCP = new Efficiency(inputHistoFile,"EffAglMCP");
	Efficiency * EffNaFMCP = new Efficiency(inputHistoFile,"EffNaFMCP");

	Efficiency *  RInner 	= new Efficiency(inputHistoFile,"RInner");
	Efficiency *  R_L1      = new Efficiency(inputHistoFile,"R_L1");
	

	string tagli[10]= {"Trigger","3of4 TOF","TRD Segments","Rigidity exists","Chi^2 R","Matching TOF","Matching TRD","In TRD Accept.","1 Particle","1 Tr. Track"};
	string nome;

	cout<<"**** MC PRESELECTIONS EFFICIENCY (FULL SET) ****"<<endl;

	EffpreselMCP -> Eval_Efficiency();
	EffpreselMCD -> Eval_Efficiency();

	EffpreselMCP_L1 -> Eval_Efficiency();
	EffpreselMCD_L1 -> Eval_Efficiency();

	for(int i=0;i<6;i++) EffCascadeMCP[i]    -> Eval_Efficiency();  
        for(int i=0;i<4;i++) EffCascadeRICHMCP[i]-> Eval_Efficiency();

	EffPreMCP  -> Eval_Efficiency();
        EffAglMCP  -> Eval_Efficiency();
        EffNaFMCP  -> Eval_Efficiency();

	RInner -> Eval_Efficiency();
        R_L1   -> Eval_Efficiency();

	TH1F * EffPreMCP_R_TH1F  =  (TH1F *)EffpreselMCP->effR	->Clone();
	TH1F * EffPreMCP_TH1F    =  (TH1F *)EffpreselMCP->effTOF->Clone();
	TH1F * EffPreMCPNaF_TH1F =  (TH1F *)EffpreselMCP->effNaF->Clone();
	TH1F * EffPreMCPAgl_TH1F =  (TH1F *)EffpreselMCP->effAgl->Clone();
	TH2F * EffPreMCD_R_TH2F  =  (TH2F *)EffpreselMCD->effR  ->Clone();
	TH2F * EffPreMCD_TH2F    =  (TH2F *)EffpreselMCD->effTOF->Clone();
	TH2F * EffPreMCDNaF_TH2F =  (TH2F *)EffpreselMCD->effNaF->Clone();
	TH2F * EffPreMCDAgl_TH2F =  (TH2F *)EffpreselMCD->effAgl->Clone();

	TH1F * EffCascade[6];
	TH1F * EffCascadeRICH[4];

	for(int i=0;i<6;i++) EffCascade[i] = (TH1F *)EffCascadeMCP[i]->effR->Clone();
	for(int i=0;i<4;i++) EffCascadeRICH[i] = (TH1F *)EffCascadeRICHMCP[i]->effR->Clone();

	TH1F * Effpre =  (TH1F *)EffPreMCP->effR->Clone();
	TH1F * Effagl =  (TH1F *)EffAglMCP->effR->Clone();
	TH1F * Effnaf =  (TH1F *)EffNaFMCP->effR->Clone();

	TH1F * EffRInner =  (TH1F *)RInner->effR->Clone();
        TH1F * EffR_L1   =  (TH1F *)R_L1  ->effR->Clone();
	
	finalHistos.Add(EffPreMCP_R_TH1F  	);
	finalHistos.Add(EffPreMCP_TH1F       );
	finalHistos.Add(EffPreMCPNaF_TH1F    );
	finalHistos.Add(EffPreMCPAgl_TH1F 	);
	finalHistos.Add(EffPreMCD_R_TH2F   	);
	finalHistos.Add(EffPreMCD_TH2F    	);
	finalHistos.Add(EffPreMCDNaF_TH2F 	);
	finalHistos.Add(EffPreMCDAgl_TH2F 	);
	finalHistos.writeObjsInFolder("Results");


	cout<<"*** Plotting ...  ****"<<endl;

	MCpreeff_Plot(

			EffPreMCP_R_TH1F , 
			EffPreMCP_TH1F   , 
			EffPreMCPNaF_TH1F, 
			EffPreMCPAgl_TH1F, 
			EffPreMCD_R_TH2F , 
			EffPreMCD_TH2F   , 
			EffPreMCDNaF_TH2F, 
			EffPreMCDAgl_TH2F, 
			EffCascade,     			
	                EffCascadeRICH,
			Effpre, 		
                        Effagl,
                        Effnaf, 
			EffRInner, 
	                EffR_L1   
		     );	


	return;

	
}


