using namespace std;


TH2F * MigrMatrix = new TH2F("MigrMatrix","MigrMatrix",nbinsr,0,nbinsr,nbinsr,0,nbinsr);

void MigrationMatrix_Fill(){
	if(!(Tup.Unbias==0&&cmask.isPreselected()&&Tup.Beta_pre>0&&Tup.R_pre>0))	return;
	if(Massa_gen<1&&Massa_gen>0.5) 
		MigrMatrix->Fill(RB.GetRBin(fabs(Tup.R_pre)) , RB.GetRBin(Tup.Momento_gen));
}

void MigrationMatrix_Write(){
        MigrMatrix->Write();
}


void MigrationMatrix(TFile * file1){
	TH2F * MigrMatrix= (TH2F*) file1->Get("MigrMatrix");
	
	cout<<"****** Migration Matrix **********"<<endl;
	float norm[nbinsr]={0};
	for(int M=0;M<nbinsr;M++)for(int l=0;l<nbinsr;l++) norm[M]+=MigrMatrix->GetBinContent(l+1,M+1);
	for(int M=0;M<nbinsr;M++)for(int l=0;l<nbinsr;l++) MigrMatrix->SetBinContent(l+1,M+1,MigrMatrix->GetBinContent(l+1,M+1)/norm[M]);
	
	TCanvas * c27 = new TCanvas("Rigidity Migration matrix");
	c27->cd();
	gPad->SetLogz();
	MigrMatrix->GetXaxis()->SetTitle("n.bin (R meas)");
	MigrMatrix->GetYaxis()->SetTitle("n.bin (R gen)");
	MigrMatrix->Draw("col");

	cout<<"*** Updating Results file ***"<<endl;
        fileFinalPlots->cd("MC Results");
	c27 -> Write();
        fileFinalPlots->Write();
	
}
