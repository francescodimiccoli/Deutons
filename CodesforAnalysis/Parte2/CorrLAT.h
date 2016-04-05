TH1 * ExposureTime(TH2 * esposizionegeo);

TH1 * Weighted_CorrLAT(TH2 * esposizionegeo, TH1 * LATcorr);

void CorrLAT(){
	string nomefile=percorso + "/Risultati/risultati/"+mese+"_"+frac+"_P1.root";
        TFile * file1 = TFile::Open(nomefile.c_str(),"READ");
        if(!file1){
                nomefile=percorso + "/Risultati/"+mese+"/"+mese+"_"+frac+"_P1.root";
                file1 =TFile::Open(nomefile.c_str(),"READ");
        }

	TH2F * esposizionegeo_R    = (TH2F*)file1->Get( "esposizionegeo"       );
	TH2F * esposizionepgeoTOF  = (TH2F*)file1->Get(	"esposizionepgeo"	);
        TH2F * esposizionepgeoNaF  = (TH2F*)file1->Get(	"esposizionepgeoNaF"	);
        TH2F * esposizionepgeoAgl  = (TH2F*)file1->Get(	"esposizionepgeoAgl"	);
        TH2F * esposizionedgeoTOF  = (TH2F*)file1->Get(	"esposizionedgeo"	);
        TH2F * esposizionedgeoNaF  = (TH2F*)file1->Get(	"esposizionedgeoNaF"	);
        TH2F * esposizionedgeoAgl  = (TH2F*)file1->Get(	"esposizionedgeoAgl"	);

	LATcorr * EffpreSelDATA     = new LATcorr(file1,"EffpreSelDATA"  ,"Results");	
	LATcorr * LikelihoodLATcorr = new LATcorr(file1,"EffLikDATA" 	 ,"Results");
        LATcorr * DistanceLATcorr   = new LATcorr(file1,"EffDistDATA"	 ,"Results");

	cout<<"******* TOTAL LAT. CORRECTION *************"<<endl;

	
	TH1F * PreLATCorrR = (TH1F *)(( (TH2F *)EffpreSelDATA ->  LATcorrR_fit ) -> ProjectionX("",1,1)) -> Clone();
	
	PreLATCorrR -> Multiply( (TH1F *)( ( (TH2F *)EffpreSelDATA ->  LATcorrR_fit ) -> ProjectionX("",2,2))->Clone()	);
	PreLATCorrR -> Multiply( (TH1F *)( ( (TH2F *)EffpreSelDATA ->  LATcorrR_fit ) -> ProjectionX("",3,3))->Clone()	);		

	TH1F * PreLATCorrTOF = (TH1F *) PreLATCorrR -> Clone();
	TH1F * PreLATCorrNaF = (TH1F *) PreLATCorrR -> Clone();
	TH1F * PreLATCorrAgl = (TH1F *) PreLATCorrR -> Clone(); //preselections don't have difference in TOF,NaF,Agl

	TH1F *  TOTLATCorrR   = (TH1F *) PreLATCorrR -> Clone();
	TH1F *  TOTLATCorrTOF = (TH1F *) PreLATCorrR -> Clone();
	TH1F *  TOTLATCorrNaF = (TH1F *) PreLATCorrR -> Clone();
	TH1F *  TOTLATCorrAgl = (TH1F *) PreLATCorrR -> Clone();

	TOTLATCorrR    -> Multiply ( ( (TH1F *)LikelihoodLATcorr ->  LATcorrTOF_fit) );
	TOTLATCorrR    -> Multiply ( ( (TH1F *)DistanceLATcorr   ->  LATcorrTOF_fit) );	
	TOTLATCorrTOF  -> Multiply ( ( (TH1F *)LikelihoodLATcorr ->  LATcorrTOF_fit) );
	TOTLATCorrTOF  -> Multiply ( ( (TH1F *)DistanceLATcorr   ->  LATcorrTOF_fit) );
	TOTLATCorrNaF  -> Multiply ( ( (TH1F *)LikelihoodLATcorr ->  LATcorrNaF_fit) );
	TOTLATCorrNaF  -> Multiply ( ( (TH1F *)DistanceLATcorr   ->  LATcorrNaF_fit) );
	TOTLATCorrAgl  -> Multiply ( ( (TH1F *)LikelihoodLATcorr ->  LATcorrAgl_fit) );
	TOTLATCorrAgl  -> Multiply ( ( (TH1F *)DistanceLATcorr   ->  LATcorrAgl_fit) );

	
	//Only pres.
	TH1F * CorrezioneLATpre_pR = (TH1F *) Weighted_CorrLAT ( esposizionegeo_R , PreLATCorrR  	);
	
	TH1F * CorrezioneLATpre_dR = (TH1F *) Weighted_CorrLAT ( esposizionegeo_R , PreLATCorrR  	);

	//Full set
	TH1F * CorrezioneLAT_pR   = (TH1F *) Weighted_CorrLAT ( esposizionegeo_R   , TOTLATCorrTOF 	);
	TH1F * CorrezioneLAT_pTOF = (TH1F *) Weighted_CorrLAT ( esposizionepgeoTOF , TOTLATCorrTOF 	);
	TH1F * CorrezioneLAT_pNaF = (TH1F *) Weighted_CorrLAT ( esposizionepgeoNaF , TOTLATCorrNaF 	);
	TH1F * CorrezioneLAT_pAgl = (TH1F *) Weighted_CorrLAT ( esposizionepgeoAgl , TOTLATCorrAgl 	);

	TH1F * CorrezioneLAT_dR   = (TH1F *) Weighted_CorrLAT ( esposizionegeo_R   , TOTLATCorrTOF 	);
	TH1F * CorrezioneLAT_dTOF = (TH1F *) Weighted_CorrLAT ( esposizionedgeoTOF , TOTLATCorrTOF	);
	TH1F * CorrezioneLAT_dNaF = (TH1F *) Weighted_CorrLAT ( esposizionedgeoNaF , TOTLATCorrNaF	);
	TH1F * CorrezioneLAT_dAgl = (TH1F *) Weighted_CorrLAT ( esposizionedgeoAgl , TOTLATCorrAgl	);

	cout<<"*** Updating P1 file ****"<<endl;
        nomefile=percorso + "/Risultati/risultati/"+mese+"_"+frac+"_P1.root";
        file1 = TFile::Open(nomefile.c_str(),"UPDATE");
        if(!file1){
                nomefile=percorso + "/Risultati/"+mese+"/"+mese+"_"+frac+"_P1.root";
                file1 =TFile::Open(nomefile.c_str(),"UPDATE");
        }
	file1->cd("Results");
	
	PreLATCorrR    -> Write(     "PreLATCorr_LATcorrR_fit" 	);
 	PreLATCorrTOF  -> Write(     "PreLATCorr_LATcorrTOF_fit"); 	
	PreLATCorrNaF  -> Write(     "PreLATCorr_LATcorrNaF_fit"); 	
	PreLATCorrAgl  -> Write(     "PreLATCorr_LATcorrAgl_fit");       

                                                                
	TOTLATCorrR    -> Write(     "TOTLATCorr_LATcorrR_fit"  ); 	
	TOTLATCorrTOF  -> Write(     "TOTLATCorr_LATcorrTOF_fit"); 	
	TOTLATCorrNaF  -> Write(     "TOTLATCorr_LATcorrNaF_fit"); 	
	TOTLATCorrAgl  -> Write(     "TOTLATCorr_LATcorrAgl_fit");
	
	CorrezioneLATpre_pR -> Write(  "CorrezioneLATPrep_R"	);       
	CorrezioneLATpre_dR -> Write(  "CorrezioneLATPred_R"	);
	
	CorrezioneLAT_pR  -> Write(  "CorrezioneLATp_R"  	);      
        CorrezioneLAT_pTOF-> Write(  "CorrezioneLATp_TOF"	);      
	CorrezioneLAT_pNaF-> Write(  "CorrezioneLATp_NaF"       );
	CorrezioneLAT_pAgl-> Write(  "CorrezioneLATp_Agl"	);

	CorrezioneLAT_dR  -> Write(  "CorrezioneLATd_R"  	);
	CorrezioneLAT_dTOF-> Write(  "CorrezioneLATd_TOF"	);
	CorrezioneLAT_dNaF-> Write(  "CorrezioneLATd_NaF"	);
	CorrezioneLAT_dAgl-> Write(  "CorrezioneLATd_Agl"       );		

        file1->Write();
	file1->Close();



	TCanvas * c26 = new TCanvas("Latitude pile-up correction (R bins)");
	TCanvas * c26_bis = new TCanvas("Latitude pile-up correction (Beta bins)");


	
	c26->Divide(1,2);
	c26->cd(1);
	gPad->SetGridy();
        gPad->SetGridx();
	TGraphErrors *CorrLAT_tot_Spl=new TGraphErrors();
	for(int m=1;m<11;m++) {
		CorrLAT_tot_Spl->SetPoint(m-1,geomagC[m],TOTLATCorrTOF -> GetBinContent(m+1));
		CorrLAT_tot_Spl->SetPointError(m-1,0,TOTLATCorrTOF -> GetBinError(m+1));
	}
	TGraphErrors *CorrLAT_pre_Spl=new TGraphErrors();
	for(int m=1;m<11;m++) {
		CorrLAT_pre_Spl->SetPoint(m-1,geomagC[m],PreLATCorrR -> GetBinContent(m+1));
		CorrLAT_pre_Spl->SetPointError(m-1,0,PreLATCorrR -> GetBinError(m+1));
	}
	CorrLAT_tot_Spl->SetLineColor(2);
	CorrLAT_tot_Spl->SetMarkerColor(2);
	CorrLAT_tot_Spl->SetLineWidth(2);
	CorrLAT_tot_Spl->SetMarkerStyle(8);
	CorrLAT_pre_Spl->SetLineColor(4);
	CorrLAT_pre_Spl->SetMarkerColor(4);
	CorrLAT_pre_Spl->SetLineWidth(2);
	CorrLAT_pre_Spl->SetMarkerStyle(8);
	
	CorrLAT_tot_Spl->SetFillStyle(3002);
	CorrLAT_tot_Spl->GetXaxis()->SetTitle("Latitude");
        CorrLAT_tot_Spl->GetYaxis()->SetTitle("Eff. Corr. Factor");
	CorrLAT_tot_Spl->Draw("APC");
	CorrLAT_pre_Spl->Draw("PCsame");

	c26->cd(2);
        gPad->SetGridy();
        gPad->SetGridx();
	gPad->SetLogx();
	TGraphErrors *CorrLAT_totM1_Spl=new TGraphErrors();
	for(int i=0;i<43;i++){
			CorrLAT_totM1_Spl->SetPoint(i,R_cent[i],CorrezioneLAT_pR->GetBinContent(i+1));
		}
	CorrLAT_totM1_Spl->SetLineColor(2);
        CorrLAT_totM1_Spl->SetMarkerColor(2);
        CorrLAT_totM1_Spl->SetLineWidth(2);
        CorrLAT_totM1_Spl->SetMarkerStyle(8);
        CorrLAT_totM1_Spl->SetFillStyle(3002);
	CorrLAT_totM1_Spl->GetXaxis()->SetTitle("R [GV]");
        CorrLAT_totM1_Spl->GetYaxis()->SetTitle("Eff. Corr. Factor");
	CorrLAT_totM1_Spl->Draw("APC");
	TGraphErrors *CorrLAT_totM2_Spl=new TGraphErrors();
        for(int i=0;i<43;i++){
                CorrLAT_totM2_Spl->SetPoint(i,R_cent[i],CorrezioneLATpre_pR->GetBinContent(i+1));
        }
        CorrLAT_totM2_Spl->SetLineColor(4);
        CorrLAT_totM2_Spl->SetMarkerColor(4);
        CorrLAT_totM2_Spl->SetLineWidth(2);
        CorrLAT_totM2_Spl->SetMarkerStyle(8);
        CorrLAT_totM2_Spl->SetFillStyle(3002);
	CorrLAT_totM2_Spl->Draw("PCsame");

	c26_bis->Divide(1,3);
        c26_bis->cd(1);
        gPad->SetGridy();
        gPad->SetGridx();
	TGraphErrors * CorrLATp_TOF_Spl=new TGraphErrors();
	int point =0;
	for(int m=0;m<18;m++){
			if(CorrezioneLAT_pTOF->GetBinContent(m+1)>0)
			CorrLATp_TOF_Spl->SetPoint(point,Ekincent[m],CorrezioneLAT_pTOF->GetBinContent(m+1));		
			point++;
		}
	CorrLATp_TOF_Spl->SetLineColor(2);
        CorrLATp_TOF_Spl->SetMarkerColor(2);
        CorrLATp_TOF_Spl->SetLineWidth(2);
        CorrLATp_TOF_Spl->SetMarkerStyle(8);
        CorrLATp_TOF_Spl->SetFillStyle(3002);
        CorrLATp_TOF_Spl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        CorrLATp_TOF_Spl->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATp_TOF_Spl->Draw("APC");
	TGraphErrors * CorrLATd_TOF_Spl=new TGraphErrors();
        point=0;
	for(int m=0;m<18;m++){
                        if(CorrezioneLAT_dTOF->GetBinContent(m+1)>0)
			CorrLATd_TOF_Spl->SetPoint(point,Ekincent[m],CorrezioneLAT_dTOF->GetBinContent(m+1));
                	point++;
		}
	CorrLATd_TOF_Spl->SetLineColor(4);
        CorrLATd_TOF_Spl->SetMarkerColor(4);
        CorrLATd_TOF_Spl->SetLineWidth(2);
        CorrLATd_TOF_Spl->SetMarkerStyle(8);
        CorrLATd_TOF_Spl->SetFillStyle(3002);		
	CorrLATd_TOF_Spl->Draw("PCsame");
	 
	c26_bis->cd(2);
        gPad->SetGridy();
        gPad->SetGridx();
        TGraphErrors * CorrLATp_NaF_Spl=new TGraphErrors();
        for(int m=0;m<18;m++){
                        CorrLATp_NaF_Spl->SetPoint(m,EkincentNaF[m],CorrezioneLAT_pNaF->GetBinContent(m+1));
                }
        CorrLATp_NaF_Spl->SetLineColor(2);
        CorrLATp_NaF_Spl->SetMarkerColor(2);
        CorrLATp_NaF_Spl->SetLineWidth(2);
        CorrLATp_NaF_Spl->SetMarkerStyle(8);
        CorrLATp_NaF_Spl->SetFillStyle(3002);
        CorrLATp_NaF_Spl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        CorrLATp_NaF_Spl->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATp_NaF_Spl->Draw("APC");
        TGraphErrors * CorrLATd_NaF_Spl=new TGraphErrors();
        for(int m=0;m<18;m++){
                        CorrLATd_NaF_Spl->SetPoint(m,EkincentNaF[m],CorrezioneLAT_dNaF->GetBinContent(m+1));
                }
        CorrLATd_NaF_Spl->SetLineColor(4);
        CorrLATd_NaF_Spl->SetMarkerColor(4);
        CorrLATd_NaF_Spl->SetLineWidth(2);
        CorrLATd_NaF_Spl->SetMarkerStyle(8);
        CorrLATd_NaF_Spl->SetFillStyle(3002);
        CorrLATd_NaF_Spl->Draw("PCsame");

	c26_bis->cd(3);
        gPad->SetGridy();
        gPad->SetGridx();
        TGraphErrors * CorrLATp_Agl_Spl=new TGraphErrors();
        for(int m=0;m<18;m++){
                        CorrLATp_Agl_Spl->SetPoint(m,EkincentAgl[m],CorrezioneLAT_pAgl->GetBinContent(m+1));
                }
        CorrLATp_Agl_Spl->SetLineColor(2);
        CorrLATp_Agl_Spl->SetMarkerColor(2);
        CorrLATp_Agl_Spl->SetLineWidth(2);
        CorrLATp_Agl_Spl->SetMarkerStyle(8);
        CorrLATp_Agl_Spl->SetFillStyle(3002);
        CorrLATp_Agl_Spl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        CorrLATp_Agl_Spl->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATp_Agl_Spl->Draw("APC");
        TGraphErrors * CorrLATd_Agl_Spl=new TGraphErrors();
        for(int m=0;m<18;m++){
                        CorrLATd_Agl_Spl->SetPoint(m,EkincentAgl[m],CorrezioneLAT_dAgl->GetBinContent(m+1));
                }
        CorrLATd_Agl_Spl->SetLineColor(4);
        CorrLATd_Agl_Spl->SetMarkerColor(4);
        CorrLATd_Agl_Spl->SetLineWidth(2);
        CorrLATd_Agl_Spl->SetMarkerStyle(8);
        CorrLATd_Agl_Spl->SetFillStyle(3002);
        CorrLATd_Agl_Spl->Draw("PCsame");

	cout<<"*** Updating Results file ***"<<endl;
        nomefile=percorso + "/CodesforAnalysis/Final_plots/"+mese+".root";
        TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
        f_out->mkdir("DATA-driven Results/Latitude effect/Correction");
        f_out->cd("DATA-driven Results/Latitude effect/Correction");
        c26->Write();
        c26_bis->Write();
        f_out->Write();
        f_out->Close();
	
	return;


}


TH1 * ExposureTime(TH2 * esposizionegeo) {
	return (TH1 *) esposizionegeo -> ProjectionX("",0,10) -> Clone();
}

TH1 * Weighted_CorrLAT(TH2 * esposizionegeo, TH1 * LATcorr){
	TH2F * temp = (TH2F *)esposizionegeo -> Clone();
	for(int m=0;m<11;m++){
		for(int i=0; i< temp -> GetNbinsX(); i++){
			temp->SetBinContent(i+1,m,esposizionegeo->GetBinContent(i+1,m)*LATcorr -> GetBinContent(m+1));
		}
	}
	TH1F * temp2 =(TH1F *)temp -> ProjectionX("",0,10) -> Clone();
	TH1F * Exptime =(TH1F *) ExposureTime(esposizionegeo);

temp2 -> Divide ( Exptime );
	return (TH1 *)temp2 -> Clone();	
}
