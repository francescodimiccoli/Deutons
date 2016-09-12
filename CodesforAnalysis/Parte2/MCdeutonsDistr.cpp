#include "PlottingFunctions/MCdeutonsDistr_Plot.h"



TH2F * MassDMCTOF = new TH2F("MassDMCTOF","MassDMCTOF",50,0,3,6,0,6);
TH2F * MassDMCNaF = new TH2F("MassDMCNaF","MassDMCNaF",50,0,3,6,0,6);
TH2F * MassDMCAgl = new TH2F("MassDMCAgl","MassDMCAgl",50,0,3,6,0,6);

TH2F * LikDMCTOF  = new TH2F( "LikDMCTOF", "LikDMCTOF",50,0,3,6,0,6);
TH2F * LikDMCNaF  = new TH2F( "LikDMCNaF", "LikDMCNaF",50,0,4,6,0,6);
TH2F * LikDMCAgl  = new TH2F( "LikDMCAgl", "LikDMCAgl",50,0,5,6,0,6);
                                                     
TH2F * DistDMCTOF = new TH2F("DistDMCTOF","DistDMCTOF",50,0,6,6,0,6);
TH2F * DistDMCNaF = new TH2F("DistDMCNaF","DistDMCNaF",50,0,6,6,0,6);
TH2F * DistDMCAgl = new TH2F("DistDMCAgl","DistDMCAgl",50,0,6,6,0,6);


void MCdeutonsDistr_Fill(){
	float mass=0;
	//cuts
	if(Tup.Beta<=0||Tup.R<=0) return;
	if(!Betastrongcut) return;
	if(!(Massa_gen<2 && Massa_gen>1)) return;

	if(cmask.isOnlyFromToF()){
		mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));	
		MassDMCTOF -> Fill (mass,ReturnMCGenType());	
		LikDMCTOF -> Fill (-log(1-Tup.LDiscriminant),ReturnMCGenType());
		DistDMCTOF-> Fill (Tup.Dist5D,ReturnMCGenType());
	}	
 	if(cmask.isFromNaF()){
		mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));	
		MassDMCNaF -> Fill (mass,ReturnMCGenType());	
		LikDMCNaF -> Fill (-log(1-Tup.LDiscriminant),ReturnMCGenType());
                DistDMCNaF-> Fill (Tup.Dist5D,ReturnMCGenType());

	}			
	if(cmask.isFromAgl()){
		mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));	
		MassDMCAgl -> Fill (mass,ReturnMCGenType());	
		LikDMCAgl -> Fill (-log(1-Tup.LDiscriminant),ReturnMCGenType());
                DistDMCAgl-> Fill (Tup.Dist5D,ReturnMCGenType());

	}

	return;
}


void MCdeutonsDistr_Write(){

	MassDMCTOF ->Write();
	MassDMCNaF ->Write();
	MassDMCAgl ->Write();

	LikDMCTOF  ->Write();
	LikDMCNaF  ->Write();
	LikDMCAgl  ->Write();
                  
	DistDMCTOF ->Write();
	DistDMCNaF ->Write();
	DistDMCAgl ->Write();

}



void MCdeutonsDistr(string filename){


	cout<<"******* Deutons MC distributions comparisons********"<<endl;

         cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

	TH2F * MassDMCTOF  =  (TH2F*) inputHistoFile ->Get ("MassDMCTOF");
        TH2F * MassDMCNaF  =  (TH2F*) inputHistoFile ->Get ("MassDMCNaF");
        TH2F * MassDMCAgl  =  (TH2F*) inputHistoFile ->Get ("MassDMCAgl");
                                                                        
        TH2F * LikDMCTOF   =  (TH2F*) inputHistoFile ->Get ("LikDMCTOF");
        TH2F * LikDMCNaF   =  (TH2F*) inputHistoFile ->Get ("LikDMCNaF");
        TH2F * LikDMCAgl   =  (TH2F*) inputHistoFile ->Get ("LikDMCAgl");
                                                                       
        TH2F * DistDMCTOF  =  (TH2F*) inputHistoFile ->Get ("DistDMCTOF");
        TH2F * DistDMCNaF  =  (TH2F*) inputHistoFile ->Get ("DistDMCNaF");
        TH2F * DistDMCAgl  =  (TH2F*) inputHistoFile ->Get ("DistDMCAgl");
	

	cout<<"******* Deutons MC distributions comparisons********"<<endl;   

	TH1F * SliceMassDMCTOF[6];
	TH1F * SliceMassDMCNaF[6];
	TH1F * SliceMassDMCAgl[6];

	TH1F * SliceLikDMCTOF[6];
	TH1F * SliceLikDMCNaF[6];
	TH1F * SliceLikDMCAgl[6];

	TH1F * SliceDistDMCTOF[6];
	TH1F * SliceDistDMCNaF[6];
	TH1F * SliceDistDMCAgl[6];

	for (int mc_type=0;mc_type<6;mc_type++){
			
		SliceMassDMCTOF[mc_type] = ProjectionXtoTH1F(MassDMCTOF , ("SliceMassDMCTOF"+to_string(mc_type)).c_str(),mc_type+1,mc_type+1);	
                SliceMassDMCNaF[mc_type] = ProjectionXtoTH1F(MassDMCNaF , ("SliceMassDMCNaF"+to_string(mc_type)).c_str(),mc_type+1,mc_type+1);
                SliceMassDMCAgl[mc_type] = ProjectionXtoTH1F(MassDMCAgl , ("SliceMassDMCAgl"+to_string(mc_type)).c_str(),mc_type+1,mc_type+1);
	                                                                             
		SliceLikDMCTOF[mc_type]  = ProjectionXtoTH1F(LikDMCTOF , ("SliceLikDMCTOF" +to_string(mc_type)).c_str(),mc_type+1,mc_type+1);
                SliceLikDMCNaF[mc_type]  = ProjectionXtoTH1F(LikDMCNaF , ("SliceLikDMCNaF" +to_string(mc_type)).c_str(),mc_type+1,mc_type+1);
                SliceLikDMCAgl[mc_type]  = ProjectionXtoTH1F(LikDMCAgl , ("SliceLikDMCAgl" +to_string(mc_type)).c_str(),mc_type+1,mc_type+1); 
                                                                                     
		SliceDistDMCTOF[mc_type] = ProjectionXtoTH1F(DistDMCTOF , ("SliceDistDMCTOF"+to_string(mc_type)).c_str(),mc_type+1,mc_type+1);
                SliceDistDMCNaF[mc_type] = ProjectionXtoTH1F(DistDMCNaF , ("SliceDistDMCNaF"+to_string(mc_type)).c_str(),mc_type+1,mc_type+1);
                SliceDistDMCAgl[mc_type] = ProjectionXtoTH1F(DistDMCAgl , ("SliceDistDMCAgl"+to_string(mc_type)).c_str(),mc_type+1,mc_type+1);

	}
	


	for (int mc_type=0;mc_type<6;mc_type++){
	
		SliceMassDMCTOF[mc_type] ->Scale(1/ SliceMassDMCTOF[mc_type]->GetEntries());
                SliceMassDMCNaF[mc_type] ->Scale(1/ SliceMassDMCNaF[mc_type]->GetEntries());
                SliceMassDMCAgl[mc_type] ->Scale(1/ SliceMassDMCAgl[mc_type]->GetEntries());
                                                                
	        SliceLikDMCTOF[mc_type]  ->Scale(1/ SliceLikDMCTOF[mc_type] ->GetEntries());
                SliceLikDMCNaF[mc_type]  ->Scale(1/ SliceLikDMCNaF[mc_type] ->GetEntries());
                SliceLikDMCAgl[mc_type]  ->Scale(1/ SliceLikDMCAgl[mc_type] ->GetEntries());
                                                                
                SliceDistDMCTOF[mc_type] ->Scale(1/ SliceDistDMCTOF[mc_type]->GetEntries());
                SliceDistDMCNaF[mc_type] ->Scale(1/ SliceDistDMCNaF[mc_type]->GetEntries());
                SliceDistDMCAgl[mc_type] ->Scale(1/ SliceDistDMCAgl[mc_type]->GetEntries());

	}



	TH1F * MeanMassDMCTOF = (TH1F*) SliceMassDMCTOF[0] ->Clone();
        TH1F * MeanMassDMCNaF = (TH1F*) SliceMassDMCNaF[0] ->Clone();
        TH1F * MeanMassDMCAgl = (TH1F*) SliceMassDMCAgl[0] ->Clone();
                                                         
        TH1F * MeanLikDMCTOF  = (TH1F*) SliceLikDMCTOF[0] ->Clone();
        TH1F * MeanLikDMCNaF  = (TH1F*) SliceLikDMCNaF[0] ->Clone();
        TH1F * MeanLikDMCAgl  = (TH1F*) SliceLikDMCAgl[0] ->Clone();
                                                         
        TH1F * MeanDistDMCTOF = (TH1F*) SliceDistDMCTOF[0] ->Clone();
        TH1F * MeanDistDMCNaF = (TH1F*) SliceDistDMCNaF[0] ->Clone();
        TH1F * MeanDistDMCAgl = (TH1F*) SliceDistDMCAgl[0] ->Clone();

	for (int mc_type=1;mc_type<6;mc_type++){
	
			MeanMassDMCTOF ->Add(SliceMassDMCTOF[mc_type] 	);
	                MeanMassDMCNaF ->Add(SliceMassDMCNaF[mc_type] 	);
                        MeanMassDMCAgl ->Add(SliceMassDMCAgl[mc_type] 	);
                                                                      
                        MeanLikDMCTOF  ->Add(SliceLikDMCTOF[mc_type]  	);
                        MeanLikDMCNaF  ->Add(SliceLikDMCNaF[mc_type]  	);
                        MeanLikDMCAgl  ->Add(SliceLikDMCAgl[mc_type]  	);
                                                                      
                        MeanDistDMCTOF ->Add(SliceDistDMCTOF[mc_type] 	);
	                MeanDistDMCNaF ->Add(SliceDistDMCNaF[mc_type] 	);
                        MeanDistDMCAgl ->Add(SliceDistDMCAgl[mc_type] 	);
	}


		MeanMassDMCTOF -> Scale(1/6.);
                MeanMassDMCNaF ->Scale(1/6.);
                MeanMassDMCAgl ->Scale(1/6.);
                                 
                MeanLikDMCTOF  ->Scale(1/6.);
                MeanLikDMCNaF  ->Scale(1/6.);
                MeanLikDMCAgl  ->Scale(1/6.);
                                 
                MeanDistDMCTOF ->Scale(1/6.);
                MeanDistDMCNaF ->Scale(1/6.);
                MeanDistDMCAgl ->Scale(1/6.);

	for (int mc_type=0;mc_type<6;mc_type++){
	
		SliceMassDMCTOF[mc_type] -> Divide(MeanMassDMCTOF );
                SliceMassDMCNaF[mc_type] -> Divide(MeanMassDMCNaF );
                SliceMassDMCAgl[mc_type] -> Divide(MeanMassDMCAgl );
                                                                  
                SliceLikDMCTOF[mc_type]  -> Divide(MeanLikDMCTOF  );
                SliceLikDMCNaF[mc_type]  -> Divide(MeanLikDMCNaF  );
                SliceLikDMCAgl[mc_type]  -> Divide(MeanLikDMCAgl  );
                                                                  
                SliceDistDMCTOF[mc_type] -> Divide(MeanDistDMCTOF );
                SliceDistDMCNaF[mc_type] -> Divide(MeanDistDMCNaF );
                SliceDistDMCAgl[mc_type] -> Divide(MeanDistDMCAgl );

	}



	MCdeutonsDistr_Plot(	

			SliceMassDMCTOF,	
                        SliceMassDMCNaF,
                        SliceMassDMCAgl,
                                       
                        SliceLikDMCTOF,
                        SliceLikDMCNaF,
                        SliceLikDMCAgl,
                                       
                        SliceDistDMCTOF,
                        SliceDistDMCNaF,
                        SliceDistDMCAgl);


}


