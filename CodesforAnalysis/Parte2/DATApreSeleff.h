using namespace std;


LATcorr *EffpreSelDATA = new LATcorr("EffpreSelDATA",3);
 

void DATApreSeleff_Fill(TNtuple *ntupla, int l,int zona){
	int k = ntupla->GetEvent(l);
	if(Unbias!=0||Beta_pre<=0||R_pre<=0||Beta_pre>protons->Eval(R_pre)+0.1||Beta_pre<protons->Eval(R_pre)-0.1) return;
	if(!(EdepL1>0&&EdepL1<EdepL1beta->Eval(Beta)+0.1&&EdepL1>EdepL1beta->Eval(Beta)-0.1)) return;	
	for(int S=0;S<3;S++){
			for(int M=0;M<43;M++) 
				if(fabs(R_pre)<bin[M+1]&&fabs(R_pre)>bin[M]&&R_pre>Rcut[zona]) {
					if(((int)Cutmask&notpassed[S])==notpassed[S]) ((TH3 *)EffpreSelDATA->beforeR)->Fill(M,zona,S);
					if(((int)Cutmask&passed[S])==passed[S]) ((TH3 *)EffpreSelDATA->afterR)->Fill(M,zona,S);
				}
			}	
	return;
}


void DATApreSeleff_Write(){
        EffpreSelDATA->Write();
        return;
}

TF1 * CorrLAT_pre[3];
TH2F *CorrLATpre_spl[3];
TGraphErrors *CorrLATpre_Spl[3];


void DATApreSeleff(TFile * file1){

	LATcorr * EffpreSelDATA = new LATcorr(file1,"EffpreSelDATA",3);	
	
	string tagli[3]={"Matching TOF","Chi^2 R","1 Tr. Track"};
        string nome;

	cout<<"********************** DATA PRESELECTIONS EFFICIENCIES ****************************************"<<endl;

	EffpreSelDATA -> Eval_Efficiency();

	TH3F *EffpreSelDATA_R = (TH3F *) EffpreSelDATA     -> effR -> Clone();	
	
	cout<<"********************** LAT. Eff. CORRECTION **************************************************"<<endl;


	cout<<"*** Updating P1 file ****"<<endl;
        string nomefile=percorso + "/Risultati/risultati/"+mese+"_"+frac+"_P1.root";
        file1 =TFile::Open(nomefile.c_str(),"UPDATE");
        if(!file1){
                nomefile=percorso + "/Risultati/"+mese+"/"+mese+"_"+frac+"_P1.root";
                file1 =TFile::Open(nomefile.c_str(),"UPDATE");
        }

        file1->cd("Results");
	EffpreSelDATA_R -> Write();
	file1->Write();
	file1->Close();

	TGraphErrors * Eff_preSelLAT[3][11];
	TCanvas *c14[4];
	for(int S=0;S<3;S++){
		nome="Latitude Efficiency: "+tagli[S];
		c14[S]=new TCanvas(nome.c_str());
		c14[S]->Divide(2,1);
		c14[S]->cd(1);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		for(int l=1;l<11;l++){
			Eff_preSelLAT[S][l]=new TGraphErrors();
			for(int i=1;i<43;i++){
				Eff_preSelLAT[S][l]->SetPoint(i,R_cent[i],EffpreSelDATA_R->GetBinContent(i+1,l+1,S+1));
				Eff_preSelLAT[S][l]->SetPointError(i,0,EffpreSelDATA_R->GetBinError(i+1,l+1,S+1));
				}
			}

		Eff_preSelLAT[S][10]->SetMarkerColor(1);
		Eff_preSelLAT[S][10]->SetLineColor(1);
		Eff_preSelLAT[S][10]->SetMarkerStyle(8);
		Eff_preSelLAT[S][10]->GetXaxis()->SetTitle("R [GV]");
		Eff_preSelLAT[S][10]->GetYaxis()->SetTitle("Efficiency");
		Eff_preSelLAT[S][10]->GetYaxis()->SetRangeUser(0.1,1);
		Eff_preSelLAT[S][10]->Draw("AP");
		for(int l=0;l<10;l++){
			Eff_preSelLAT[S][l]->SetMarkerColor(l);
			Eff_preSelLAT[S][l]->SetLineColor(l);
			Eff_preSelLAT[S][l]->SetMarkerStyle(8);
			Eff_preSelLAT[S][l]->Draw("Psame");
		}
	}

/*

	TGraphErrors *CorrLATpre[3];
	for(int S=0;S<3;S++){
	c14[S]->cd(2);
	gPad->SetGridy();
	gPad->SetGridx();
	nome="Latitude Efficiency: "+tagli[S];
	nome=nome+": Fit";
	CorrLAT_pre[S]=new TF1(nome.c_str(),"pol3");
	CorrLATpre[S]=new TGraphErrors();
	CorrLATpre[S]->SetTitle("Latitude Efficiency Corr.");
	CorrLATpre[S]->GetXaxis()->SetTitle("Latitude");
	CorrLATpre[S]->GetYaxis()->SetTitle("Eff. Corr. Factor");
	CorrLATpre[S]->SetMarkerStyle(8);
	for(int i=1;i<11;i++) {
			CorrLATpre[S]->SetPoint(i-1,geomagC[i],Tau_D[S][i]);
			CorrLATpre[S]->SetPointError(i-1,0,Tau_D_err[S][i]);	
		}
	CorrLATpre[S]->Fit(nome.c_str());
	CorrLATpre[S]->Draw("AP");	

	nome="CorrLATpre_spl"+tagli[S];
        CorrLATpre_spl[S]=new TH2F(nome.c_str(),nome.c_str(),11,0,11,2,0,2);
	for(int i=0;i<11;i++) CorrLATpre_spl[S]->SetBinContent(i+1,1,CorrLAT_pre[S]->Eval(geomagC[i]));
	for(int i=0;i<11;i++) CorrLATpre_spl[S]->SetBinContent(i+1,2,Tau_D_err[S][i]);
	CorrLATpre_Spl[S]=new TGraphErrors(tagli[S].c_str());
	CorrLATpre_Spl[S]->SetName(tagli[S].c_str());
	int j=0;
	for(int i=1;i<11;i++) { CorrLATpre_Spl[S]->SetPoint(j,geomagC[i],CorrLAT_pre[S]->Eval(geomagC[i]));
				CorrLATpre_Spl[S]->SetPointError(j,0,Tau_D_err[S][i]);
				j++;
				}
	CorrLATpre_Spl[S]->SetLineColor(2);
	CorrLATpre_Spl[S]->SetMarkerColor(2);
	CorrLATpre_Spl[S]->SetFillColor(2);
	CorrLATpre_Spl[S]->SetFillStyle(3001);
	CorrLATpre_Spl[S]->SetLineWidth(2);
	CorrLATpre_Spl[S]->Draw("PLsame");	

	}
*/
	cout<<"*** Updating Results file ***"<<endl;
        nomefile=percorso + "/CodesforAnalysis/Final_plots/"+mese+".root";
        TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
        f_out->mkdir("DATA-driven Results/Latitude effect/\"Clean-event\" Selections");
        f_out->cd("DATA-driven Results/Latitude effect/\"Clean-event\" Selections");
        for(int S=0;S<3;S++) c14[S]->Write();
        f_out->Write();
        f_out->Close();

}

