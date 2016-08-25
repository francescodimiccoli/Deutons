#include "PlottingFunctions/AntiDCutOptimization_Plot.h"


OptimizationCut * DiscriminantCutTOF = new   OptimizationCut("D_DiscrCUT_TOF" , nbinsToF,-1.2,1.2);
OptimizationCut * DiscriminantCutNaF = new   OptimizationCut("D_DiscrCUT_NaF" , nbinsNaF,-1.2,1.2);
OptimizationCut * DiscriminantCutAgl = new   OptimizationCut("D_DiscrCUT_Agl" , nbinsAgl,-1.2,1.2);


void AntiDCutOptimization_Fill(){

	//cuts
	if(Tup.Beta<=0||Tup.R<=0) return;
	if(!(Likcut)) return;
	//      	
	      	
	float Distance_Discr = 0;
	int Kbin;
	//TOF
	if(!cmask.isFromNaF()&&!cmask.isFromAgl()){
		Kbin=ToFDB.GetBin(RUsed);
		if(Tup.Beta<0.9){
		Distance_Discr = ((Tup.Dist5D_P-Tup.Dist5D)/(Tup.Dist5D_P+Tup.Dist5D));
		if(Massa_gen<1&&Massa_gen>0.5) ((TH2*)DiscriminantCutTOF -> Distrib_P) -> Fill(Distance_Discr,Kbin,Tup.mcweight);	
		if(Massa_gen<2&&Massa_gen>1.5) ((TH2*)DiscriminantCutTOF -> Distrib_D) -> Fill(Distance_Discr,Kbin,Tup.mcweight);
	}}

	//NaF
	if(cmask.isFromNaF()) {
		Kbin=NaFDB.GetBin(RUsed);
		Distance_Discr = ((Tup.Dist5D_P-Tup.Dist5D)/(Tup.Dist5D_P+Tup.Dist5D));
		if(Massa_gen<1&&Massa_gen>0.5) ((TH2*)DiscriminantCutNaF -> Distrib_P) -> Fill(Distance_Discr,Kbin,Tup.mcweight);
		if(Massa_gen<2&&Massa_gen>1.5) ((TH2*)DiscriminantCutNaF -> Distrib_D) -> Fill(Distance_Discr,Kbin,Tup.mcweight);
	}

	//Agl
	if(cmask.isFromAgl()) {
		Kbin=AglDB.GetBin(RUsed);
		Distance_Discr = ((Tup.Dist5D_P-Tup.Dist5D)/(Tup.Dist5D_P+Tup.Dist5D));
		if(Massa_gen<1&&Massa_gen>0.5) ((TH2*)DiscriminantCutAgl -> Distrib_P) -> Fill(Distance_Discr,Kbin,Tup.mcweight);
		if(Massa_gen<2&&Massa_gen>1.5) ((TH2*)DiscriminantCutAgl -> Distrib_D) -> Fill(Distance_Discr,Kbin,Tup.mcweight);
	}

}




void  AntiDCutOptimization_Write()
{

	DiscriminantCutTOF->Write();
	DiscriminantCutNaF->Write();
	DiscriminantCutAgl->Write();

	return;
}


void AntiDCutOptimization(string filename){

	cout<<"******************** ANTI-D CUT OPTIMIZATION ************************"<<endl;
	cout<<"*** Reading  P1 file ****"<<endl;
	TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

	OptimizationCut * DiscriminantCutTOF = new   OptimizationCut(inputHistoFile,"D_DiscrCUT_TOF" );
	OptimizationCut * DiscriminantCutNaF = new   OptimizationCut(inputHistoFile,"D_DiscrCUT_NaF" );
	OptimizationCut * DiscriminantCutAgl = new   OptimizationCut(inputHistoFile,"D_DiscrCUT_Agl" );

	cout<<"******************** ANTI-D CUT OPTIMIZATION ************************"<<endl;
	
	DiscriminantCutTOF -> NormalizeDistributions();
        DiscriminantCutNaF -> NormalizeDistributions();
        DiscriminantCutAgl -> NormalizeDistributions();
	
	DiscriminantCutTOF->Optimization(); 
        DiscriminantCutNaF->Optimization();   
        DiscriminantCutAgl->Optimization();   

	TH1F * OptimizedCutsTOF = DiscriminantCutTOF-> Getcuts();
	TH1F * OptimizedCutsNaF = DiscriminantCutNaF-> Getcuts();
	TH1F * OptimizedCutsAgl = DiscriminantCutAgl-> Getcuts();


	finalHistos.Add(	OptimizedCutsTOF);
	finalHistos.Add(        OptimizedCutsNaF);
	finalHistos.Add(        OptimizedCutsAgl);

	finalHistos.writeObjsInFolder("Results");	


	cout<<"*** Plotting ...  ****"<<endl;

	AntiDCutOptimization_Plot(DiscriminantCutTOF,
                                  DiscriminantCutNaF,
                                  DiscriminantCutAgl
					);

}

