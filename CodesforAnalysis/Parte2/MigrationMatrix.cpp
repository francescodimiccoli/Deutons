#include "PlottingFunctions/MigrationMatrix_Plot.h"

using namespace std;


TH2F * MigrMatrix = new TH2F("MigrMatrix","MigrMatrix",nbinsr,0,nbinsr,nbinsr,0,nbinsr);

void MigrationMatrix_Fill(){
	if(!(trgpatt.IsUnbias()&&cmask.isPreselected()&&Tup.Beta_pre>0&&Tup.R_pre>0))	return;
	if(Massa_gen<1&&Massa_gen>0.5) 
		MigrMatrix->Fill(PRB.GetRBin(fabs(Tup.R_pre)) , PRB.GetRBin(Tup.Momento_gen),Tup.mcweight);
}

void MigrationMatrix_Write(){
        MigrMatrix->Write();
}


void MigrationMatrix(string filename){
	cout<<"****** Migration Matrix **********"<<endl;
	
	cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

	
	TH2F * MigrMatrix= (TH2F*) inputHistoFile->Get("MigrMatrix");
	
	cout<<"****** Migration Matrix **********"<<endl;
	float norm[nbinsr]={0};
	for(int M=0;M<nbinsr;M++)for(int l=0;l<nbinsr;l++) norm[M]+=MigrMatrix->GetBinContent(l+1,M+1);
	for(int M=0;M<nbinsr;M++)for(int l=0;l<nbinsr;l++) MigrMatrix->SetBinContent(l+1,M+1,MigrMatrix->GetBinContent(l+1,M+1)/norm[M]);
	
	
	cout<<"*** Plotting ...  ****"<<endl;

	MigrationMatrix_Plot( MigrMatrix);

	return;
}

