#include "PlottingFunctions/MCFullseteff_Plot.h"

using namespace std;




void MCFullseteff(string filename){

	cout<<"**** MC FULL SET  EFFICIENCY  ****"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

	Efficiency * EffpreselMCP = new Efficiency(inputHistoFile, "EffpreselMCP");
	Efficiency * EffpreselMCD = new Efficiency(inputHistoFile, "EffpreselMCD");

	Efficiency * EffLikMCP = new Efficiency (inputHistoFile,"EffLikMCP");
	Efficiency * EffLikMCD = new Efficiency (inputHistoFile,"EffLikMCD");
	
	Efficiency * EffDistMCP = new Efficiency (inputHistoFile,"EffDistMCP");
	Efficiency * EffDistMCD = new Efficiency (inputHistoFile,"EffDistMCD");
	
	cout<<"**** MC FULL SET  EFFICIENCY  ****"<<endl;
	
	EffpreselMCP -> Eval_Efficiency();
	EffpreselMCD -> Eval_Efficiency();

	EffLikMCP   -> Eval_Efficiency();
	EffLikMCD   -> Eval_Efficiency();

	EffDistMCP   -> Eval_Efficiency();
	EffDistMCD   -> Eval_Efficiency();

	TH1F * EffPreMCP_R_TH1F  =  (TH1F *)EffpreselMCP->effR	->Clone();  
	TH1F * EffPreMCP_TH1F    =  (TH1F *)EffpreselMCP->effTOF->Clone();
	TH1F * EffPreMCPNaF_TH1F =  (TH1F *)EffpreselMCP->effNaF->Clone();
	TH1F * EffPreMCPAgl_TH1F =  (TH1F *)EffpreselMCP->effAgl->Clone();
	TH2F * EffPreMCD_R_TH2F  =  (TH2F *)EffpreselMCD->effR  ->Clone();
	TH2F * EffPreMCD_TH2F    =  (TH2F *)EffpreselMCD->effTOF->Clone();
	TH2F * EffPreMCDNaF_TH2F =  (TH2F *)EffpreselMCD->effNaF->Clone();
	TH2F * EffPreMCDAgl_TH2F =  (TH2F *)EffpreselMCD->effAgl->Clone();


	TH1F * EffFullsetMCP_R_TH1F  =  (TH1F *)EffDistMCP->effR  ->Clone();
	TH1F * EffFullsetMCP_TH1F    =  (TH1F *)EffDistMCP->effTOF->Clone();
	TH1F * EffFullsetMCPNaF_TH1F =  (TH1F *)EffDistMCP->effNaF->Clone();
	TH1F * EffFullsetMCPAgl_TH1F =  (TH1F *)EffDistMCP->effAgl->Clone();
	TH2F * EffFullsetMCD_R_TH2F  =  (TH2F *)EffDistMCD->effR  ->Clone();
	TH2F * EffFullsetMCD_TH2F    =  (TH2F *)EffDistMCD->effTOF->Clone();
	TH2F * EffFullsetMCDNaF_TH2F =  (TH2F *)EffDistMCD->effNaF->Clone();
	TH2F * EffFullsetMCDAgl_TH2F =  (TH2F *)EffDistMCD->effAgl->Clone();

	EffFullsetMCP_R_TH1F -> Multiply( (TH1F *)EffLikMCP->effR  ->Clone()	); 
	EffFullsetMCP_TH1F   -> Multiply( (TH1F *)EffLikMCP->effTOF->Clone()	); 
	EffFullsetMCPNaF_TH1F-> Multiply( (TH1F *)EffLikMCP->effNaF->Clone()	); 
	EffFullsetMCPAgl_TH1F-> Multiply( (TH1F *)EffLikMCP->effAgl->Clone()	); 
	EffFullsetMCD_R_TH2F -> Multiply( (TH2F *)EffLikMCD->effR  ->Clone()	); 
	EffFullsetMCD_TH2F   -> Multiply( (TH2F *)EffLikMCD->effTOF->Clone()	); 
	EffFullsetMCDNaF_TH2F-> Multiply( (TH2F *)EffLikMCD->effNaF->Clone()	); 
	EffFullsetMCDAgl_TH2F-> Multiply( (TH2F *)EffLikMCD->effAgl->Clone()	); 


	EffFullsetMCP_R_TH1F -> Multiply( EffPreMCP_R_TH1F  	); 
	EffFullsetMCP_TH1F   -> Multiply( EffPreMCP_TH1F    	); 
	EffFullsetMCPNaF_TH1F-> Multiply( EffPreMCPNaF_TH1F 	); 
	EffFullsetMCPAgl_TH1F-> Multiply( EffPreMCPAgl_TH1F 	); 
	EffFullsetMCD_R_TH2F -> Multiply( EffPreMCD_R_TH2F  	); 
	EffFullsetMCD_TH2F   -> Multiply( EffPreMCD_TH2F    	); 
	EffFullsetMCDNaF_TH2F-> Multiply( EffPreMCDNaF_TH2F 	); 
	EffFullsetMCDAgl_TH2F-> Multiply( EffPreMCDAgl_TH2F 	); 


	EffFullsetMCP_R_TH1F	->SetName("EffFullsetMCP_EffR");
	EffFullsetMCP_TH1F 	->SetName("EffFullsetMCP_EffTOF");
	EffFullsetMCPNaF_TH1F	->SetName("EffFullsetMCP_EffNaF");
	EffFullsetMCPAgl_TH1F 	->SetName("EffFullsetMCP_EffAgl");
	EffFullsetMCD_R_TH2F	->SetName("EffFullsetMCD_EffR");
	EffFullsetMCD_TH2F	->SetName("EffFullsetMCD_EffTOF");
	EffFullsetMCDNaF_TH2F 	->SetName("EffFullsetMCD_EffNaF");
	EffFullsetMCDAgl_TH2F	->SetName("EffFullsetMCD_EffAgl");


	finalHistos.Add(EffFullsetMCP_R_TH1F );
	finalHistos.Add(EffFullsetMCP_TH1F   );
	finalHistos.Add(EffFullsetMCPNaF_TH1F);
	finalHistos.Add(EffFullsetMCPAgl_TH1F);	
	finalHistos.Add(EffFullsetMCD_R_TH2F );
	finalHistos.Add(EffFullsetMCD_TH2F   );
	finalHistos.Add(EffFullsetMCDNaF_TH2F);
	finalHistos.Add(EffFullsetMCDAgl_TH2F);

	finalHistos.writeObjsInFolder("Results");

        cout<<"*** Plotting ...  ****"<<endl;

	MCFullseteff_Plot(EffFullsetMCP_R_TH1F, 
                          EffFullsetMCP_TH1F   ,
                          EffFullsetMCPNaF_TH1F,
                          EffFullsetMCPAgl_TH1F,
                          EffFullsetMCD_R_TH2F ,
                          EffFullsetMCD_TH2F   ,
                          EffFullsetMCDNaF_TH2F,
                          EffFullsetMCDAgl_TH2F);


	return;
}


