using namespace std;

TCanvas * c22 = new TCanvas("Protons Acceptance");

TH2F * AcceptPzone = new TH2F("AcceptPzone","AcceptPzone",43,0,43,11,0,11);
TH2F * AcceptPzone_pre = new TH2F("AcceptPzone_pre","AcceptPzone_pre",43,0,43,11,0,11);

void AcceptanceP(TFile * file1){

	float AcceptgeoP[43]={0};
	float AcceptSelMCP[43]={0};
	float AcceptpreMCP[43]={0};
	float AcceptSelP[43][11]={0};
	float AcceptpreP[43][11]={0};	

	float eventiprova=0;
	for(int i=1;i<43;i++) eventiprova+=EffpreselMCP1_R->GetBinContent(i+1);
	float triggerbin=(pow(0.0308232619,-1))*eventiprova/43;
	float triggertot=(pow(0.0308232619,-1))*eventiprova;
	for(int i=0;i<43;i++) 
		AcceptgeoP[i]=EffTriggerMCP_R_TH1F->GetBinContent(i+1)*EffTOF_MCP_R_TH1F->GetBinContent(i+1)*EffTriggMCP1_R->GetBinContent(i+1)/triggerbin*47.78;
	for(int j=0;j<11;j++) 
		for(int i=0;i<43;i++)	{
			//AcceptSelP[i][j]=AcceptgeoP[i]*EffTrackerMCP_R_TH1F->GetBinContent(i+1);
			//for(int S=0;S<3;S++) AcceptSelP[i][j]*=EffPreSelMCP_R_TH2F->GetBinContent(i+1,S+1);
			AcceptpreMCP[i]=EffPreMCP_R_TH1F->GetBinContent(i+1)*EffpreselMCP1_R->GetBinContent(i+1)/triggerbin*47.78;
			//AcceptSelP[i][j]=AcceptSelP[i][j]*EffUnbDATA_R_TH1F->GetBinContent(i+1);
			AcceptSelP[i][j]=AcceptpreMCP[i];
			AcceptpreP[i][j]=AcceptSelP[i][j];
			AcceptSelP[i][j]*=EffMCDistP_TH1F->GetBinContent(i+1);
			AcceptSelMCP[i]=AcceptSelP[i][j];
			// Data-driven corrections
			//latitude
			for(int S=0;S<3;S++) AcceptSelP[i][j]/=CorrLAT_pre[S]->Eval(geomagC[j]);
			for(int S=0;S<3;S++) AcceptpreP[i][j]/=CorrLAT_pre[S]->Eval(geomagC[j]);
			AcceptSelP[i][j]/=CorrLAT_Lik->Eval(geomagC[j])*CorrLAT_Dist->Eval(geomagC[j]);
			//DvsMC
			/*for(int S=0;S<3;S++) AcceptSelP[i][j]*=PreDVSMC_P[S]->GetBinContent(i+1,1);
			for(int S=0;S<3;S++) AcceptpreP[i][j]*=PreDVSMC_P[S]->GetBinContent(i+1,1);		
			AcceptSelP[i][j]*=DistDVSMC_P->Eval(encinprot[i])*LikDVSMC_P->Eval(encinprot[i]);*/	
			}
	//c22->Divide(2,1);
	c22->cd();
	gPad->SetLogx();
        gPad->SetLogy();
	gPad->SetGridx();
        gPad->SetGridy();
	TGraphErrors * AccgeoP= new TGraphErrors();
	TGraphErrors * AccPreMCP= new TGraphErrors();
	TGraphErrors * AccSelMCP= new TGraphErrors();
	TGraphErrors * AccSelP[11];
	TGraphErrors * AccpreP[11];
	int p=0;
	for(int i=0;i<43;i++) if(AcceptgeoP[i]>0) {AccgeoP->SetPoint(p,encinprot[i],AcceptgeoP[i]);p++;}	
	p=0;
        for(int i=0;i<43;i++) if(AcceptpreMCP[i]>0) {AccPreMCP->SetPoint(p,encinprot[i],AcceptpreMCP[i]);p++;}
	p=0;
	for(int i=0;i<43;i++) if(AcceptSelMCP[i]>0) {AccSelMCP->SetPoint(p,encinprot[i],AcceptSelMCP[i]);p++;}
	
	for(int j=0;j<11;j++) for(int i=0;i<43;i++)   AcceptPzone->SetBinContent(i+1,j+1,AcceptSelP[i][j]);
	for(int j=0;j<11;j++) for(int i=0;i<43;i++) AcceptPzone_pre->SetBinContent(i+1,j+1,AcceptpreP[i][j]);

	for(int j=0;j<11;j++) {
		AccSelP[j]=new TGraphErrors();
		p=0;
		for(int i=0;i<43;i++) if(AcceptSelP[i][j]>0) {AccSelP[j]->SetPoint(p,encinprot[i],AcceptSelP[i][j]);p++;}
		AccSelP[j]->SetMarkerStyle(8);
        	AccSelP[j]->SetMarkerColor(j-1);
        	AccSelP[j]->SetLineColor(j-1);
        	AccSelP[j]->SetLineWidth(2);
	}
	for(int j=0;j<11;j++) {
                AccpreP[j]=new TGraphErrors();
                p=0;
		for(int i=0;i<43;i++) if(AcceptpreP[i][j]>0) {AccpreP[j]->SetPoint(p,encinprot[i],AcceptpreP[i][j]);p++;}
                AccpreP[j]->SetMarkerStyle(8);
                AccpreP[j]->SetMarkerColor(j-1);
                AccpreP[j]->SetLineColor(j-1);
                AccpreP[j]->SetLineWidth(2);
        }

	AccgeoP->SetMarkerStyle(8);
	AccgeoP->SetMarkerColor(2);
	AccgeoP->SetLineColor(2);
	AccgeoP->SetLineWidth(4);
	AccPreMCP->SetMarkerStyle(8);
        AccPreMCP->SetMarkerColor(1);
        AccPreMCP->SetLineColor(2);
        AccPreMCP->SetLineWidth(4);
	AccSelMCP->SetMarkerStyle(8);
        AccSelMCP->SetMarkerColor(1);
        AccSelMCP->SetLineColor(2);
        AccSelMCP->SetLineWidth(4);
	AccgeoP->SetTitle("Protons Acceptance");
        AccgeoP->GetXaxis()->SetTitle("R [GV]");
        AccgeoP->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
        AccgeoP->GetXaxis()->SetTitleSize(0.045);
        AccgeoP->GetYaxis()->SetTitleSize(0.045);
	AccgeoP->GetYaxis()->SetRangeUser(1e-2,1.3);
	AccgeoP->Draw("AC");
	for(int j=0;j<11;j++) AccSelP[j]->Draw("PCsame");
	for(int j=0;j<11;j++) AccpreP[j]->Draw("PCsame");
	AccPreMCP->Draw("Csame");
	AccSelMCP->Draw("Csame");
}

