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


	Efficiency * EffFullsetMCP = new Efficiency(EffpreselMCP,"EffFullsetMCP"); 
	Efficiency * EffFullsetMCD = new Efficiency(EffpreselMCD,"EffFullsetMCD"); 

	EffFullsetMCP->Compose_Efficiency(EffLikMCP); 
        EffFullsetMCD->Compose_Efficiency(EffLikMCD);  

	EffFullsetMCP->Compose_Efficiency(EffDistMCP); 
        EffFullsetMCD->Compose_Efficiency(EffDistMCD);  

	EffFullsetMCP->Eval_FittedEfficiency();
	EffFullsetMCD->Eval_FittedEfficiency();

	finalHistos.Add(EffFullsetMCP->effR_fit  );
        finalHistos.Add(EffFullsetMCP->effTOF_fit);
        finalHistos.Add(EffFullsetMCP->effNaF_fit);
        finalHistos.Add(EffFullsetMCP->effAgl_fit);
        finalHistos.Add(EffFullsetMCD->effR_fit  );
        finalHistos.Add(EffFullsetMCD->effTOF_fit);
        finalHistos.Add(EffFullsetMCD->effNaF_fit);
        finalHistos.Add(EffFullsetMCD->effAgl_fit);
	
	finalHistos.Add(EffFullsetMCP->err_systR  );
	finalHistos.Add(EffFullsetMCP->err_systTOF);
	finalHistos.Add(EffFullsetMCP->err_systNaF);
	finalHistos.Add(EffFullsetMCP->err_systAgl);
	finalHistos.Add(EffFullsetMCP->err_statR  );
	finalHistos.Add(EffFullsetMCP->err_statTOF);
	finalHistos.Add(EffFullsetMCP->err_statNaF);
	finalHistos.Add(EffFullsetMCP->err_statAgl);


	finalHistos.Add(EffFullsetMCP->effR  );
	finalHistos.Add(EffFullsetMCP->effTOF);
	finalHistos.Add(EffFullsetMCP->effNaF);
	finalHistos.Add(EffFullsetMCP->effAgl);	
	finalHistos.Add(EffFullsetMCD->effR  );
	finalHistos.Add(EffFullsetMCD->effTOF);
	finalHistos.Add(EffFullsetMCD->effNaF);
	finalHistos.Add(EffFullsetMCD->effAgl);
	
	finalHistos.Add(EffFullsetMCD->err_systR  );
	finalHistos.Add(EffFullsetMCD->err_systTOF);
	finalHistos.Add(EffFullsetMCD->err_systNaF);
	finalHistos.Add(EffFullsetMCD->err_systAgl);
	finalHistos.Add(EffFullsetMCD->err_statR  );
	finalHistos.Add(EffFullsetMCD->err_statTOF);
	finalHistos.Add(EffFullsetMCD->err_statNaF);
	finalHistos.Add(EffFullsetMCD->err_statAgl);



	


	finalHistos.writeObjsInFolder("Results");

        cout<<"*** Plotting ...  ****"<<endl;

	MCFullseteff_Plot(

			  EffFullsetMCP->effR_fit  , 
                          EffFullsetMCP->effTOF_fit,
                          EffFullsetMCP->effNaF_fit,
                          EffFullsetMCP->effAgl_fit,
                          EffFullsetMCD->effR_fit  ,
                          EffFullsetMCD->effTOF_fit,
                          EffFullsetMCD->effNaF_fit,
                          EffFullsetMCD->effAgl_fit,


			  EffFullsetMCP->effR  , 
                          EffFullsetMCP->effTOF,
                          EffFullsetMCP->effNaF,
                          EffFullsetMCP->effAgl,
                          EffFullsetMCD->effR  ,
                          EffFullsetMCD->effTOF,
                          EffFullsetMCD->effNaF,
                          EffFullsetMCD->effAgl);


	return;
}


