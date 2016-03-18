using namespace std;

class Efficiency
{
public:
    TH1 * beforeR,   * afterR;
    TH1 * beforeTOF, * afterTOF;
    TH1 * beforeNaF, * afterNaF;
    TH1 * beforeAgl, * afterAgl;

    // Creation constructors:
    Efficiency(std::string basename){
        beforeTOF = new TH1F((basename + "1"   ).c_str(),(basename + "1"   ).c_str(),18,0,18);
        afterTOF  = new TH1F((basename + "2"   ).c_str(),(basename + "2"   ).c_str(),18,0,18);
        beforeNaF = new TH1F((basename + "1NaF").c_str(),(basename + "1NaF").c_str(),18,0,18);
        afterNaF  = new TH1F((basename + "2NaF").c_str(),(basename + "2NaF").c_str(),18,0,18);
        beforeAgl = new TH1F((basename + "1Agl").c_str(),(basename + "1Agl").c_str(),18,0,18);
        afterAgl  = new TH1F((basename + "2Agl").c_str(),(basename + "2Agl").c_str(),18,0,18);
        beforeR   = new TH1F((basename + "1_R" ).c_str(),(basename + "1_R" ).c_str(),43,0,43);
        afterR    = new TH1F((basename + "2_R" ).c_str(),(basename + "2_R" ).c_str(),43,0,43);
    }

    Efficiency(std::string basename, int n){
        beforeTOF = new TH2F((basename + "1"   ).c_str(),(basename + "1"   ).c_str(),18,0,18, n, 0 ,n);
        afterTOF  = new TH2F((basename + "2"   ).c_str(),(basename + "2"   ).c_str(),18,0,18, n, 0 ,n);
        beforeNaF = new TH2F((basename + "1NaF").c_str(),(basename + "1NaF").c_str(),18,0,18, n, 0 ,n);
        afterNaF  = new TH2F((basename + "2NaF").c_str(),(basename + "2NaF").c_str(),18,0,18, n, 0 ,n);
        beforeAgl = new TH2F((basename + "1Agl").c_str(),(basename + "1Agl").c_str(),18,0,18, n, 0 ,n);
        afterAgl  = new TH2F((basename + "2Agl").c_str(),(basename + "2Agl").c_str(),18,0,18, n, 0 ,n);
        beforeR   = new TH2F((basename + "1_R" ).c_str(),(basename + "1_R" ).c_str(),43,0,43, n, 0 ,n);
        afterR    = new TH2F((basename + "2_R" ).c_str(),(basename + "2_R" ).c_str(),43,0,43, n, 0 ,n);
    }

    Efficiency(std::string basename, int n, int m){
        beforeTOF = new TH3F((basename + "1"   ).c_str(),(basename + "1"   ).c_str(),18,0,18, n, 0 ,n, m, 0, m);
        afterTOF  = new TH3F((basename + "2"   ).c_str(),(basename + "2"   ).c_str(),18,0,18, n, 0 ,n, m, 0, m);
        beforeNaF = new TH3F((basename + "1NaF").c_str(),(basename + "1NaF").c_str(),18,0,18, n, 0 ,n, m, 0, m);
        afterNaF  = new TH3F((basename + "2NaF").c_str(),(basename + "2NaF").c_str(),18,0,18, n, 0 ,n, m, 0, m);
        beforeAgl = new TH3F((basename + "1Agl").c_str(),(basename + "1Agl").c_str(),18,0,18, n, 0 ,n, m, 0, m);
        afterAgl  = new TH3F((basename + "2Agl").c_str(),(basename + "2Agl").c_str(),18,0,18, n, 0 ,n, m, 0, m);
        beforeR   = new TH3F((basename + "1_R" ).c_str(),(basename + "1_R" ).c_str(),43,0,43, n, 0 ,n, m, 0, m);
        afterR    = new TH3F((basename + "2_R" ).c_str(),(basename + "2_R" ).c_str(),43,0,43, n, 0 ,n, m, 0, m);
    }

    // Reading constructors
    
    Efficiency(TFile * file, std::string basename){
        beforeTOF = (TH1F *)file->Get((basename + "1"   ).c_str());
        afterTOF  = (TH1F *)file->Get((basename + "2"   ).c_str());
        beforeNaF = (TH1F *)file->Get((basename + "1NaF").c_str());
        afterNaF  = (TH1F *)file->Get((basename + "2NaF").c_str());
        beforeAgl = (TH1F *)file->Get((basename + "1Agl").c_str());
        afterAgl  = (TH1F *)file->Get((basename + "2Agl").c_str());
        beforeR   = (TH1F *)file->Get((basename + "1_R" ).c_str());
        afterR    = (TH1F *)file->Get((basename + "2_R" ).c_str());
    }

    Efficiency(TFile * file, std::string basename, int i){
        beforeTOF = (TH2F *)file->Get((basename + "1"   ).c_str());
        afterTOF  = (TH2F *)file->Get((basename + "2"   ).c_str());
        beforeNaF = (TH2F *)file->Get((basename + "1NaF").c_str());
        afterNaF  = (TH2F *)file->Get((basename + "2NaF").c_str());
        beforeAgl = (TH2F *)file->Get((basename + "1Agl").c_str());
        afterAgl  = (TH2F *)file->Get((basename + "2Agl").c_str());
        beforeR   = (TH2F *)file->Get((basename + "1_R" ).c_str());
        afterR    = (TH2F *)file->Get((basename + "2_R" ).c_str());
    }

    Efficiency(TFile * file, std::string basename, int i, int j){
        beforeTOF = (TH3F *)file->Get((basename + "1"   ).c_str());
        afterTOF  = (TH3F *)file->Get((basename + "2"   ).c_str());
        beforeNaF = (TH3F *)file->Get((basename + "1NaF").c_str());
        afterNaF  = (TH3F *)file->Get((basename + "2NaF").c_str());
        beforeAgl = (TH3F *)file->Get((basename + "1Agl").c_str());
        afterAgl  = (TH3F *)file->Get((basename + "2Agl").c_str());
        beforeR   = (TH3F *)file->Get((basename + "1_R" ).c_str());
        afterR    = (TH3F *)file->Get((basename + "2_R" ).c_str());
    }

    void Write();

};


void Efficiency::Write()
{
    beforeTOF->Write();
    beforeNaF->Write();
    beforeAgl->Write();
    afterTOF->Write();
    afterNaF->Write();
    afterAgl->Write();
}



Efficiency * EffpreselMCP = new Efficiency("EfficiencyP");
Efficiency * EffpreselMCD = new Efficiency("EfficiencyD", 6);


void MCpreseff_Fill(TNtuple *ntupla, int l){
	int k = ntupla->GetEvent(l);
	
	if(Massa_gen<1&&Massa_gen>0.5) {
        //R bins
		for(int M=0;M<43;M++) 
        {
            if(fabs(Momento_gen)<bin[M+1]&&fabs(Momento_gen)>bin[M]) EffpreselMCP->beforeR->Fill(M);
			if(fabs(R_pre)<bin[M+1]&&fabs(R_pre)>bin[M]) 	
            {
                if(Unbias==0&&((int)Cutmask&187)==187&&Beta_pre>0&&R_pre>0) EffpreselMCP->afterR->Fill(M);
            }
        }
        
        // Beta bins
		for(int m=0;m<18;m++) { 
            if(Var3>BetaP[m]   &&Var3<=BetaP[m+1])    EffpreselMCP->beforeTOF->Fill(m);
            if(Var3>BetaNaFP[m]&&Var3<=BetaNaFP[m+1]) EffpreselMCP->beforeNaF->Fill(m);
            if(Var3>BetaAglP[m]&&Var3<=BetaAglP[m+1]) EffpreselMCP->beforeAgl->Fill(m);
            
            if(Unbias==0&&((int)Cutmask&187)==187&&Beta_pre>0&&R_pre>0)
            {
                if(Var>BetaP[m]&&Var<=BetaP[m+1])EffpreselMCP->afterTOF->Fill(m);	
                if((((int)Cutmask)>>11)==512&&Var2>BetaNaFP[m]&&Var2<=BetaNaFP[m+1]) EffpreselMCP->afterNaF->Fill(m);	
                if((((int)Cutmask)>>11)==0&&Var2>BetaAglP[m]&&Var2<=BetaAglP[m+1]) EffpreselMCP->afterAgl->Fill(m);
            }
        }
	}				 

	if(Massa_gen>1&&Massa_gen<2) {
        // R bins
		for(int M=0;M<43;M++) 
        {
            if(fabs(Momento_gen)<bin[M+1]&&fabs(Momento_gen)>bin[M]) 
                EffpreselMCD->beforeR->Fill(M,(int)(10000*Massa_gen-18570));
            
		    if(((int)Cutmask&187)==187&&Beta_pre>0&&Unbias==0&&R_pre>0) 
            {
			    if(fabs(R_pre)<bin[M+1]&&fabs(R_pre)>bin[M])	
				    EffpreselMCD->afterR->Fill(M,(int)(10000*Massa_gen-18570));
            }
        }

        // Beta bins
		for(int m=0;m<18;m++) 
        {
            if(Var3>BetaD[m]&&Var3<=BetaD[m+1])
                EffpreselMCD->beforeTOF->Fill(m,(int)(10000*Massa_gen-18570));
            if(Var3>BetaNaFD[m]&&Var3<=BetaNaFD[m+1])
			    EffpreselMCD->beforeNaF->Fill(m,(int)(10000*Massa_gen-18570));
            if(Var3>BetaAglD[m]&&Var3<=BetaAglD[m+1])
			    EffpreselMCD->beforeAgl->Fill(m,(int)(10000*Massa_gen-18570));
        
			if(((int)Cutmask&187)==187&&Beta_pre>0&&Unbias==0&&R_pre>0) 
            {
                if(Var>BetaD[m]&&Var<=BetaD[m+1])	
                    EffpreselMCD->afterTOF->Fill(m,(int)(10000*Massa_gen-18570));
                if((((int)Cutmask)>>11)==512&&Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1])
                    EffpreselMCD->afterNaF->Fill(m,(int)(10000*Massa_gen-18570));
                if((((int)Cutmask)>>11)==0&&Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1])
                    EffpreselMCD->afterAgl->Fill(m,(int)(10000*Massa_gen-18570));
            }
	    }
    }
	return;
}


void MCpreeff_Write(){
    EffpreselMCP->Write();
    EffpreselMCD->Write();
    return;
}



void MCpreeff(TFile * file1){
	Efficiency * EffpreselMCP = new Efficiency(file1, "EffpreselMCP");
	Efficiency * EffpreselMCD = new Efficiency(file1, "EffpreselMCD", 6);

	string numero[18]={"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","18"};
	string tagli[10]={"Trigger","3of4 TOF","TRD Segments","Rigidity exists","Chi^2 R","Matching TOF","Matching TRD","In TRD Accept.","1 Particle","1 Tr. Track"};
	string nome;

	cout<<"**** MC PRESELECTIONS EFFICIENCY (FULL SET) ****"<<endl;

	TH1F * EffPreMCP_R_TH1F  = (TH1F *)EffpreselMCP->afterR  ->Clone(); 
	TH1F * EffPreMCP_TH1F    = (TH1F *)EffpreselMCP->afterTOF->Clone(); 
	TH1F * EffPreMCPNaF_TH1F = (TH1F *)EffpreselMCP->afterNaF->Clone(); 
	TH1F * EffPreMCPAgl_TH1F = (TH1F *)EffpreselMCP->afterAgl->Clone(); 
	TH2F * EffPreMCD_R_TH2F  = (TH2F *)EffpreselMCD->afterR  ->Clone(); 
	TH2F * EffPreMCD_TH2F    = (TH2F *)EffpreselMCD->afterTOF->Clone(); 
	TH2F * EffPreMCDNaF_TH2F = (TH2F *)EffpreselMCD->afterNaF->Clone(); 
	TH2F * EffPreMCDAgl_TH2F = (TH2F *)EffpreselMCD->afterAgl->Clone(); 

	EffPreMCP_R_TH1F  -> Divide( EffpreselMCP->beforeR   ); 
	EffPreMCP_TH1F    -> Divide( EffpreselMCP->beforeTOF ); 
	EffPreMCPNaF_TH1F -> Divide( EffpreselMCP->beforeNaF ); 
	EffPreMCPAgl_TH1F -> Divide( EffpreselMCP->beforeAgl ); 
	EffPreMCD_R_TH2F  -> Divide( EffpreselMCD->beforeR   ); 
	EffPreMCD_TH2F    -> Divide( EffpreselMCD->beforeTOF ); 
	EffPreMCDNaF_TH2F -> Divide( EffpreselMCD->beforeNaF ); 
	EffPreMCDAgl_TH2F -> Divide( EffpreselMCD->beforeAgl ); 

	cout<<"*** Updating P1 file ****"<<endl;
	string nomefile=percorso + "/Risultati/risultati/"+mese+"_"+frac+"_P1.root";
	file1 =TFile::Open(nomefile.c_str(),"UPDATE");
	if(!file1){
		nomefile=percorso + "/Risultati/"+mese+"/"+mese+"_"+frac+"_P1.root";
		file1 =TFile::Open(nomefile.c_str(),"UPDATE");
	}

	TCanvas *c4=new TCanvas("Preselections Efficiency (R bins)");
	TCanvas *c4_bis=new TCanvas("Preselections Efficiency (Beta bins)");	
	c4_bis->Divide(3,1);

	c4->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	string MCLegend[7]={"protons.B800","d.pl1.0_520_GG_Blic","d.pl1.0_520_GG_BlicDPMJet","d.pl1.0_520_GG_QMD","d.pl1.0_520_Shen_Blic","d.pl1.0_520_Shen_BlicDPMJet","d.pl1.0_520_Shen_QMD"};
	TGraph * EffPreMCP_R = new TGraph();
	EffPreMCP_R->SetTitle(MCLegend[0].c_str());
	for(int i=0;i<43;i++) EffPreMCP_R->SetPoint(i,R_cent[i],EffPreMCP_R_TH1F->GetBinContent(i+1));
	TGraph * EffPreMCD_R[6];
	EffPreMCP_R->SetMarkerColor(2);
	EffPreMCP_R->SetMarkerStyle(8);
	EffPreMCP_R->SetLineColor(2);
	EffPreMCP_R->SetLineWidth(2);
	EffPreMCP_R->SetTitle("Preselections Efficiency MC (R bins)");
	EffPreMCP_R->GetXaxis()->SetTitle("R [GV]");
	EffPreMCP_R->GetYaxis()->SetTitle("Pres. Efficiency");
	EffPreMCP_R->GetXaxis()->SetTitleSize(0.045);
	EffPreMCP_R->GetYaxis()->SetTitleSize(0.045);
	{
		EffPreMCP_R->Draw("ACP");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(EffPreMCP_R,MCLegend[0].c_str(), "ep");

		for(int h=0;h<6;h++){
			EffPreMCD_R[h]= new TGraph();
			EffPreMCD_R[h]->SetTitle(MCLegend[h+1].c_str());
			for(int i=1;i<43;i++) EffPreMCD_R[h]->SetPoint(i,R_cent[i],EffPreMCD_R_TH2F->GetBinContent(i+1,h+1));
			leg->AddEntry(EffPreMCD_R[h],MCLegend[h+1].c_str(), "ep");
			EffPreMCD_R[h]->SetMarkerColor(4);
			EffPreMCD_R[h]->SetMarkerStyle(h+3);
			EffPreMCD_R[h]->SetMarkerSize(2);
			EffPreMCD_R[h]->SetLineColor(4);
			EffPreMCD_R[h]->SetLineWidth(2);
			EffPreMCD_R[h]->Draw("Psame");
			leg->Draw();
		}
	}

	c4_bis->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraph * EffPreMCP = new TGraph();
	for(int i=0;i<18;i++) EffPreMCP->SetPoint(i,Ekincent[i],EffPreMCP_TH1F->GetBinContent(i+1));
	TGraph * EffPreMCD[6];
	EffPreMCP->SetMarkerColor(2);
	EffPreMCP->SetMarkerStyle(8);
	EffPreMCP->SetLineColor(2);
	EffPreMCP->SetLineWidth(2);
	EffPreMCP->SetTitle("Preselections Efficiency MC (Beta bins)");
	EffPreMCP->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	EffPreMCP->GetYaxis()->SetTitle("Pres. Efficiency");
	EffPreMCP->GetXaxis()->SetTitleSize(0.045);
	EffPreMCP->GetYaxis()->SetTitleSize(0.045);
	{
		EffPreMCP->Draw("ACP");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(EffPreMCP,MCLegend[0].c_str(), "ep");

		for(int h=0;h<6;h++){
			EffPreMCD[h]= new TGraph();
			for(int i=0;i<18;i++) EffPreMCD[h]->SetPoint(i,Ekincent[i], EffPreMCD_TH2F->GetBinContent(i+1,h+1));
			EffPreMCD[h]->SetMarkerColor(4);
			EffPreMCD[h]->SetMarkerStyle(h+3);
			leg->AddEntry(EffPreMCD[h],MCLegend[h+1].c_str(), "ep");
			EffPreMCD[h]->SetMarkerSize(2);
			EffPreMCD[h]->SetLineColor(4);
			EffPreMCD[h]->SetLineWidth(2);
			EffPreMCD[h]->Draw("Psame");
			leg->Draw();
		}
	}
	
	c4_bis->cd(2);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraph * EffPreMCPNaF = new TGraph();
        for(int i=0;i<18;i++) EffPreMCPNaF->SetPoint(i,EkincentNaF[i],EffPreMCPNaF_TH1F->GetBinContent(i+1));
        TGraph * EffPreMCDNaF[6];
        EffPreMCPNaF->SetMarkerColor(2);
        EffPreMCPNaF->SetMarkerStyle(8);
        EffPreMCPNaF->SetLineColor(2);
        EffPreMCPNaF->SetLineWidth(2);
        EffPreMCPNaF->SetTitle("Preselections Efficiency MC (Beta bins NaF)");
        EffPreMCPNaF->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        EffPreMCPNaF->GetYaxis()->SetTitle("Pres. Efficiency");
        EffPreMCPNaF->GetXaxis()->SetTitleSize(0.045);
        EffPreMCPNaF->GetYaxis()->SetTitleSize(0.045);
        {
                EffPreMCPNaF->Draw("ACP");
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffPreMCPNaF,MCLegend[0].c_str(), "ep");

                for(int h=0;h<6;h++){
                        EffPreMCDNaF[h]= new TGraph();
                        for(int i=0;i<18;i++) EffPreMCDNaF[h]->SetPoint(i,EkincentNaF[i], 
                            EffPreMCDNaF_TH2F->GetBinContent(i+1,h+1));
                        EffPreMCDNaF[h]->SetMarkerColor(4);
                        EffPreMCDNaF[h]->SetMarkerStyle(h+3);
                        leg->AddEntry(EffPreMCDNaF[h],MCLegend[h+1].c_str(), "ep");
                        EffPreMCDNaF[h]->SetMarkerSize(2);
                        EffPreMCDNaF[h]->SetLineColor(4);
                        EffPreMCDNaF[h]->SetLineWidth(2);
                        EffPreMCDNaF[h]->Draw("Psame");
                        leg->Draw();
                }
        }

	c4_bis->cd(3);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraph * EffPreMCPAgl = new TGraph();
        for(int i=0;i<18;i++) EffPreMCPAgl->SetPoint(i,EkincentAgl[i],EffPreMCPAgl_TH1F->GetBinContent(i+1));
        TGraph * EffPreMCDAgl[6];
        EffPreMCPAgl->SetMarkerColor(2);
        EffPreMCPAgl->SetMarkerStyle(8);
        EffPreMCPAgl->SetLineColor(2);
        EffPreMCPAgl->SetLineWidth(2);
        EffPreMCPAgl->SetTitle("Preselections Efficiency MC (Beta bins Agl)");
        EffPreMCPAgl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        EffPreMCPAgl->GetYaxis()->SetTitle("Pres. Efficiency");
        EffPreMCPAgl->GetXaxis()->SetTitleSize(0.045);
        EffPreMCPAgl->GetYaxis()->SetTitleSize(0.045);
        {
                EffPreMCPAgl->Draw("ACP");
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffPreMCPAgl,MCLegend[0].c_str(), "ep");

                for(int h=0;h<6;h++){
                        EffPreMCDAgl[h]= new TGraph();
                        for(int i=0;i<18;i++) EffPreMCDAgl[h]->SetPoint(i,EkincentAgl[i], EffPreMCDAgl_TH2F->GetBinContent(i+1,h+1));
                        EffPreMCDAgl[h]->SetMarkerColor(4);
                        EffPreMCDAgl[h]->SetMarkerStyle(h+3);
                        leg->AddEntry(EffPreMCDAgl[h],MCLegend[h+1].c_str(), "ep");
                        EffPreMCDAgl[h]->SetMarkerSize(2);
                        EffPreMCDAgl[h]->SetLineColor(4);
                        EffPreMCDAgl[h]->SetLineWidth(2);
                        EffPreMCDAgl[h]->Draw("Psame");
                        leg->Draw();
                }
        }
	
	file1->mkdir("Result");
	file1->cd("Result");
	EffPreMCP_R_TH1F->Write();
	EffPreMCP_TH1F ->Write();
	EffPreMCPNaF_TH1F->Write();
	EffPreMCPAgl_TH1F ->Write();
	EffPreMCD_R_TH2F->Write();
	EffPreMCD_TH2F->Write();
	EffPreMCDNaF_TH2F ->Write();
	EffPreMCDAgl_TH2F->Write();
	file1->Close();


	cout<<"*** Updating Results file ***"<<endl;
	nomefile=percorso + "/CodesforAnalysis/Final_plots/"+mese+".root";
        TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
	f_out->mkdir("MC Results");
	f_out->mkdir("MC Results/Preselections");
	f_out->cd("MC Results/Preselections");
	c4->Write();
	c4_bis->Write();
	f_out->Close();

}
