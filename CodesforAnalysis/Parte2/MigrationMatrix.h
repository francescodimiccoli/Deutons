using namespace std;


TH2F * MigrMatrix = new TH2F("MigrMatrix","MigrMatrix",nbinsr,0,nbinsr,nbinsr,0,nbinsr);

void MigrationMatrix_Fill(TNtuple *ntupla, int l){
	int k = ntupla->GetEvent(l);
	if(!(Unbias==0&&((int)Cutmask&187)==187&&Beta_pre>0&&R_pre>0))	return;
	if(Massa_gen<1&&Massa_gen>0.5) 
		MigrMatrix->Fill(GetArrayBin(fabs(R_pre),       bin, nbinsr),
		                 GetArrayBin(fabs(Momento_gen), bin, nbinsr) );
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
        string nomefile="./Final_plots/"+mese+".root";
        TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
        f_out->cd("MC Results");
	c27 -> Write();
        f_out->Write();
        f_out->Close();
	
}
