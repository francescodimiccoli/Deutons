#include "PlottingFunctions/MCpreeff_Plot.h"

using namespace std;



Efficiency * EffpreselMCP = new Efficiency("EffpreselMCP");
Efficiency * EffpreselMCD = new Efficiency("EffpreselMCD", 6);


void MCpreseff_Fill() {

	int kbin;
	if(Massa_gen<1&&Massa_gen>0.5) {
		//R bins
		EffpreselMCP->beforeR->Fill(PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
		if(cmask.isPreselected()&&Tup.Beta_pre>0&&Tup.R_pre>0) 
			EffpreselMCP->afterR->Fill(PRB.GetRBin(RUsed),Tup.mcweight);

		// Beta bins
		EffpreselMCP->beforeTOF->Fill(ToFPB.GetBin(Tup.Momento_gen),Tup.mcweight);	
		EffpreselMCP->beforeNaF->Fill(NaFPB.GetBin(Tup.Momento_gen),Tup.mcweight);	
		EffpreselMCP->beforeAgl->Fill(AglPB.GetBin(Tup.Momento_gen),Tup.mcweight);	

		if(cmask.isPreselected() && Tup.Beta_pre>0 && Tup.R_pre>0)
		{
			kbin=ToFPB.GetBin(RUsed);
			EffpreselMCP->afterTOF->Fill(kbin,Tup.mcweight); 
			if(cmask.isFromNaF()) {
				kbin=NaFPB.GetBin(RUsed);	
				EffpreselMCP->afterNaF->Fill(kbin,Tup.mcweight); 
			} 
			if(cmask.isFromAgl()) {
				kbin=AglPB.GetBin(RUsed);   
				EffpreselMCP->afterAgl->Fill(kbin,Tup.mcweight);
			}  
		}

	}

	if(Massa_gen>1&&Massa_gen<2) {
		// R bins      
		((TH2*)EffpreselMCD->beforeR) ->Fill( PRB.GetRBin(Tup.Momento_gen),ReturnMCGenType(),Tup.mcweight);
		if(cmask.isPreselected()&&Tup.Beta_pre>0&&Tup.R_pre>0)
			((TH2*) EffpreselMCD->beforeR) ->Fill( PRB.GetRBin(RUsed),ReturnMCGenType(),Tup.mcweight);

		// Beta bins

		((TH2*)EffpreselMCD->beforeTOF)->Fill( ToFDB.GetBin(Tup.Momento_gen),ReturnMCGenType(),Tup.mcweight); 
		((TH2*)EffpreselMCD->beforeNaF)->Fill( NaFDB.GetBin(Tup.Momento_gen),ReturnMCGenType(),Tup.mcweight); 
		((TH2*)EffpreselMCD->beforeAgl)->Fill( AglDB.GetBin(Tup.Momento_gen),ReturnMCGenType(),Tup.mcweight); 

		if(cmask.isPreselected() && Tup.Beta_pre>0 && Tup.R_pre>0)
		{
			kbin=ToFDB.GetBin(RUsed);
			((TH2*)EffpreselMCD->afterTOF)->Fill(kbin,ReturnMCGenType(),Tup.mcweight);
			if(cmask.isFromNaF()) {
				kbin=NaFDB.GetBin(RUsed);
				((TH2*)EffpreselMCD->afterNaF)->Fill(kbin,ReturnMCGenType(),Tup.mcweight);
			}
			if(cmask.isFromAgl()) {
				kbin=AglDB.GetBin(RUsed);
				((TH2*)EffpreselMCD->afterAgl)->Fill(kbin,ReturnMCGenType(),Tup.mcweight);
			}
		}

	}
	return;
}


void MCpreeff_Write() {
	EffpreselMCP->Write();
	EffpreselMCD->Write();
	return;
}



void MCpreeff(string filename) {

	cout<<"**** MC PRESELECTIONS EFFICIENCY (FULL SET) ****"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
	TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

	Efficiency * EffpreselMCP = new Efficiency(inputHistoFile, "EffpreselMCP");
	Efficiency * EffpreselMCD = new Efficiency(inputHistoFile, "EffpreselMCD");

	string tagli[10]= {"Trigger","3of4 TOF","TRD Segments","Rigidity exists","Chi^2 R","Matching TOF","Matching TRD","In TRD Accept.","1 Particle","1 Tr. Track"};
	string nome;

	cout<<"**** MC PRESELECTIONS EFFICIENCY (FULL SET) ****"<<endl;

	EffpreselMCP -> Eval_Efficiency();
	EffpreselMCD -> Eval_Efficiency();

	TH1F * EffPreMCP_R_TH1F  =  (TH1F *)EffpreselMCP->effR	->Clone();
	TH1F * EffPreMCP_TH1F    =  (TH1F *)EffpreselMCP->effTOF->Clone();
	TH1F * EffPreMCPNaF_TH1F =  (TH1F *)EffpreselMCP->effNaF->Clone();
	TH1F * EffPreMCPAgl_TH1F =  (TH1F *)EffpreselMCP->effAgl->Clone();
	TH2F * EffPreMCD_R_TH2F  =  (TH2F *)EffpreselMCD->effR  ->Clone();
	TH2F * EffPreMCD_TH2F    =  (TH2F *)EffpreselMCD->effTOF->Clone();
	TH2F * EffPreMCDNaF_TH2F =  (TH2F *)EffpreselMCD->effNaF->Clone();
	TH2F * EffPreMCDAgl_TH2F =  (TH2F *)EffpreselMCD->effAgl->Clone();

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
			EffPreMCDAgl_TH2F 

		     );	


	return;

	
}


