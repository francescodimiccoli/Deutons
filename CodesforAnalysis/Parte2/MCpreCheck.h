using namespace std;


TH1F * EffpreCheckP1_R=new TH1F("EffpreCheckP1_R","EffpreselMCP1_R",43,0,43);
TH1F * EffpreCheckP2_R=new TH1F("EffpreCheckP2_R","EffpreselMCP2_R",43,0,43);


void MCpreCheck_Fill(TNtuple *ntupla, int l){
        int k = ntupla->GetEvent(l);
        if(Massa_gen<1&&Massa_gen>0.5) {
                for(int M=0;M<43;M++) if(fabs(Momento_gen)<bin[M+1]&&fabs(Momento_gen)>bin[M])
                        EffpreCheckP1_R->Fill(M);
                if(Unbias==0&&((int)Cutmask&11)==11&&Beta_pre>0&&R_pre>0)
                        for(int M=0;M<43;M++) if(fabs(R_pre)<bin[M+1]&&fabs(R_pre)>bin[M])
                                //if(EdepTOFU<EdepTOFbeta->Eval(Beta_pre)+1&&EdepTOFU>EdepTOFbeta->Eval(Beta_pre)-1) 
                                  //          if(Beta_pre<protons->Eval(R_pre)+0.1&&Beta_pre>protons->Eval(R_pre)-0.1)                    
						EffpreCheckP2_R->Fill(M);
	}
} 

void MCpreCheck_Copy(TFile * file){
        EffpreCheckP1_R= (TH1F*) file->Get("EffpreCheckP1_R");
        EffpreCheckP2_R= (TH1F*) file->Get("EffpreCheckP2_R");
}

void MCpreCheck_Write(){
        EffpreCheckP1_R->Write();
        EffpreCheckP2_R->Write();
}


TCanvas *c10 =new TCanvas("Cascade Pres. Eff.");
TH1F *EffPreCheckP_R_TH1F = new TH1F("EffPreCheckP_R_TH1F","EffPreCheckP_R_TH1F",43,0,43);
                                                               
void MCpreCheck(){

	cout<<"**** MC P Eff. CHECK ****"<<endl;
	float EffpreCheckP_R[43]={0};
        for(int i=1;i<43;i++) EffpreCheckP_R[i]=EffpreCheckP2_R->GetBinContent(i+1)/(float)EffpreCheckP1_R->GetBinContent(i+1);
	for(int i=1;i<43;i++) EffPreCheckP_R_TH1F->SetBinContent(i+1,EffpreCheckP_R[i]);
	
	float EffPreCompP[43][9]={{0}};
	for(int i=1;i<43;i++) EffPreCompP[i][0]=EffTriggerMCP_R_TH1F->GetBinContent(i+1);
	for(int i=1;i<43;i++) EffPreCompP[i][1]=EffPreCompP[i][0]*EffTOF_MCP_R_TH1F->GetBinContent(i+1);
	for(int i=1;i<43;i++) EffPreCompP[i][2]=EffPreCompP[i][1]*EffTrackerMCP_R_TH1F->GetBinContent(i+1);
	for(int i=1;i<43;i++) EffPreCompP[i][3]=EffPreCompP[i][2]*EffQTOFerMCP_R_TH1F->GetBinContent(i+1);
	for(int i=1;i<43;i++)EffPreCompP[i][4]=EffPreCompP[i][3]*EffMassMCP_R_TH1F->GetBinContent(i+1);	
	for(int i=1;i<43;i++) EffPreCompP[i][5]=EffPreCompP[i][4]*EffPreSelMCP_R_TH2F->GetBinContent(i+1,3);
	for(int i=1;i<43;i++) EffPreCompP[i][6]=EffPreCompP[i][5]*EffPreSelMCP_R_TH2F->GetBinContent(i+1,2);
	for(int i=1;i<43;i++) EffPreCompP[i][7]=EffPreCompP[i][6]*EffPreSelMCP_R_TH2F->GetBinContent(i+1,1);
	for(int i=1;i<43;i++)EffPreCompP[i][8]=EffPreCompP[i][7]*EffUnbMCP_R_TH1F->GetBinContent(i+1);	

	c10->cd();
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        string MCLegend[9]={"3/4 TOF Activity","Beta from 3/4 cluster","Track Exists","TOF Charge","Mass cut","Only 1 Track","ChisquareY<5","TOF-Track Matching","Unbias==0"};
        TGraph * EffPreCompMCP_R[9]; 
	for(int h=0;h<9;h++) EffPreCompMCP_R[h]= new TGraph();
        TGraph * EffPreFullMCP_R = new TGraph();
	for(int h=0;h<9;h++){
	EffPreCompMCP_R[h]->SetTitle(MCLegend[h].c_str());
        for(int i=0;i<43;i++) EffPreCompMCP_R[h]->SetPoint(i,R_cent[i],EffPreCompP[i][h]);
        }
	for(int i=0;i<43;i++) EffPreFullMCP_R->SetPoint(i,R_cent[i],EffPreCheckP_R_TH1F->GetBinContent(i+1));
	//TGraph * EffTrackerMCD_R[6];
      	for(int h=0;h<9;h++){
	EffPreCompMCP_R[h]->SetMarkerColor(h);
	EffPreFullMCP_R->SetMarkerColor(2);
        EffPreCompMCP_R[h]->SetMarkerStyle(8);
        EffPreFullMCP_R->SetMarkerStyle(4);
	EffPreCompMCP_R[h]->SetLineColor(h);
        EffPreCompMCP_R[h]->SetLineWidth(2);
	EffPreFullMCP_R->SetLineColor(2);
        EffPreFullMCP_R->SetLineWidth(2);
        EffPreFullMCP_R->GetXaxis()->SetTitle("R [GV]");
        EffPreFullMCP_R->GetYaxis()->SetTitle("Pres. Efficiency");
	EffPreFullMCP_R->GetYaxis()->SetRangeUser(0,1.1);
        EffPreFullMCP_R->GetXaxis()->SetTitleSize(0.045);
        EffPreFullMCP_R->GetYaxis()->SetTitleSize(0.045);
        }
	{
                EffPreFullMCP_R->Draw("ACP");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffPreFullMCP_R,"Full Set Cuts", "ep");
		for(int h=0;h<9;h++){
			EffPreCompMCP_R[h]->Draw("CPsame");
			leg->AddEntry(EffPreCompMCP_R[h],MCLegend[h].c_str(), "ep");
        	}
		leg->Draw("same");
	}

}
