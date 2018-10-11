void DrawQvsRBeta(TH2F * h6, TH2F *h9, TH2F * h6_, TH2F * h9_, std::string name, FileSaver finalResults){

   TCanvas * c1_ = new TCanvas("Protons MC");
   c1_->cd();
   PlotTH2F(gPad, h6, "Beta TOF",name.c_str(),"colz");
   TCanvas * c2_ = new TCanvas("Deutons MC");
   c2_->cd();
   PlotTH2F(gPad, h9, "Beta TOF",name.c_str(),"colz");
   TCanvas * c1_p = new TCanvas("Projections");
   c1_p->cd();
   {	  
	   TH1D * LowBetaP = h6->ProjectionY("0.5<#beta<0.6 Protons",h6->GetXaxis()->FindBin(0.5),h6->GetXaxis()->FindBin(0.6));	
	   TH1D * LowBetaD = h9->ProjectionY("0.5<#beta<0.6 Deuterons",h9->GetXaxis()->FindBin(0.5),h9->GetXaxis()->FindBin(0.6));	
	   LowBetaP->Scale(1/LowBetaP->Integral());
	   LowBetaD->Scale(1/LowBetaD->Integral());
	   PlotDistribution(gPad,(TH1F*)LowBetaP,name.c_str(), "Counts", 2, "same",1,1.2*LowBetaP->GetBinContent(LowBetaP->GetMaximumBin()));
	   PlotDistribution(gPad,(TH1F*)LowBetaD,name.c_str(), "Counts", 4, "same",1,1.2*LowBetaP->GetBinContent(LowBetaP->GetMaximumBin()));
	   TH1D * HighBetaP = h6->ProjectionY("0.7<#beta<0.8 Protons",h6->GetXaxis()->FindBin(0.7),h6->GetXaxis()->FindBin(0.8));	
	   TH1D * HighBetaD = h9->ProjectionY("0.7<#beta<0.8 Deuterons",h9->GetXaxis()->FindBin(0.7),h9->GetXaxis()->FindBin(0.8));	
	   HighBetaP->Scale(1/HighBetaP->Integral());
	   HighBetaD->Scale(1/HighBetaD->Integral());
	   HighBetaP->SetLineStyle(2);
	   HighBetaD->SetLineStyle(2);
	   PlotDistribution(gPad,(TH1F*)HighBetaP,name.c_str(), "Counts", 2, "same",1,1.2*LowBetaP->GetBinContent(LowBetaP->GetMaximumBin()));
	   PlotDistribution(gPad,(TH1F*)HighBetaD,name.c_str(), "Counts", 4, "same",1,1.2*LowBetaP->GetBinContent(LowBetaP->GetMaximumBin()));
   }
   finalResults.Add(c1_);
   finalResults.Add(c2_);
   finalResults.Add(c1_p);
   finalResults.writeObjsInFolder(("Cleaning/"+name+"vs Beta").c_str());

   TCanvas * c1__ = new TCanvas("Protons MC");
   c1__->cd();
   PlotTH2F(gPad, h6_, "Beta TOF",name.c_str(),"colz");
   TCanvas * c2__ = new TCanvas("Deutons MC");
   c2__->cd();
   PlotTH2F(gPad, h9_, "Beta TOF",name.c_str(),"colz");
    
   finalResults.Add(c1__);
   finalResults.Add(c2__);
   finalResults.writeObjsInFolder(("Cleaning/"+name+"vs R").c_str());


};

void DrawCleaning(FileSaver finalHistos,FileSaver finalResults){
   TH2F * h3d = (TH2F*)finalHistos.Get("MassTOFvsRupdown/MassTOFvsRupdown");
   TH2F * h4d = (TH2F*)finalHistos.Get("RvsChisquare_x/RvsChisquare_x");
   TH2F * h5d = (TH2F*)finalHistos.Get("RvsChisquare_y/RvsChisquare_y");

   TH2F * h6  = (TH2F*)finalHistos.Get("UtofQvsBeta_MCP/UtofQvsBeta_MCP"); 	
   TH2F * h7  = (TH2F*)finalHistos.Get("LtofQvsBeta_MCP/LtofQvsBeta_MCP"); 	
   TH2F * h8  = (TH2F*)finalHistos.Get("InnerQvsBeta_MCP/InnerQvsBeta_MCP"); 	

   TH2F * h9  = (TH2F*)finalHistos.Get("UtofQvsBeta_MCD/UtofQvsBeta_MCD"); 	
   TH2F * h10  = (TH2F*)finalHistos.Get("LtofQvsBeta_MCD/LtofQvsBeta_MCD"); 	
   TH2F * h11  = (TH2F*)finalHistos.Get("InnerQvsBeta_MCD/InnerQvsBeta_MCD"); 	

   TH2F * h6_  = (TH2F*)finalHistos.Get("UtofQvsR_MCP/UtofQvsR_MCP"); 	
   TH2F * h7_  = (TH2F*)finalHistos.Get("LtofQvsR_MCP/LtofQvsR_MCP"); 	
   TH2F * h8_  = (TH2F*)finalHistos.Get("InnerQvsR_MCP/InnerQvsR_MCP"); 	

   TH2F * h9_  = (TH2F*)finalHistos.Get("UtofQvsR_MCD/UtofQvsR_MCD"); 	
   TH2F * h10_  = (TH2F*)finalHistos.Get("LtofQvsR_MCD/LtofQvsR_MCD"); 	
   TH2F * h11_  = (TH2F*)finalHistos.Get("InnerQvsR_MCD/InnerQvsR_MCD"); 	


   
   cout<<h3d<<" "<<h4d<<endl;
   TH2F * chieffX = (TH2F *)h4d->Clone();
   chieffX->Reset();
   chieffX->SetName("chieffX");
   chieffX->SetTitle("chieffX");
  
   TH2F * chieffY = (TH2F *)h4d->Clone();
   chieffY->Reset();
   chieffY->SetName("chieffY");
   chieffY->SetTitle("chieffY");
 

   TCanvas *c2 = new TCanvas("ChiSquare Profiles", "ChiSquare Profiles",4586,195,884,613);
   SetCanvas(c2);
   h3d->GetXaxis()->SetTitle("|Rup-Rdown|/R");
   h3d->GetYaxis()->SetTitle("Mass [GeV/c^{2}]");
   h3d->Draw();

   TH1D * ProfilesX[100];

   for(int i=0;i<100;i++){
        TH1D * ProfileY = h4d->ProjectionY("",1,i+1);
        TH1F * Profile = (TH1F*)ProfileY->Clone();
        TH1D * ProfileX = h4d->ProjectionX("",1,i+1);
        ProfilesX[i] = (TH1D*) ProfileX->Clone();
        TH1D * Tot = h4d->ProjectionY();
        for(int j=0;j<100;j++){
                if(Tot->GetBinContent(j+1)>0){
                        chieffX->SetBinContent(i+1,j+1,Profile->GetBinContent(j+1)/Tot->GetBinContent(j+1));
                   }
                }
        }

   TH1D * ProfilesY[100];

   for(int i=0;i<100;i++){
        TH1D * ProfileY = h5d->ProjectionY("",1,i+1);
        TH1F * Profile = (TH1F*)ProfileY->Clone();
        TH1D * ProfileX = h5d->ProjectionX("",1,i+1);
        ProfilesY[i] = (TH1D*) ProfileX->Clone();
        TH1D * Tot = h5d->ProjectionY();
        for(int j=0;j<100;j++){
                if(Tot->GetBinContent(j+1)>0){
                        chieffY->SetBinContent(i+1,j+1,Profile->GetBinContent(j+1)/Tot->GetBinContent(j+1));
                   }
                }
        }
 
 
  
   TCanvas *c3 = new TCanvas("ChiSquare Profiles", "ChiSquare Profiles",4586,195,884,613);
//   SetCanvas(c3);
//   gStyle->SetPalette(55);
   gPad->SetLogy();
   TLegend * leg = new TLegend(0.67,0.2,0.94,0.4);
   leg->SetFillColor(0);
   leg->SetLineColor(0);

   int nColors = 256;//gStyle->GetNumberOfColors();
    int nHistos = 100;
   for(int i=0;i<100;i++) {
                int histocolor=(float)nColors / nHistos * i;
              //  ProfilesX[i]->SetLineColor(gStyle->GetColorPalette(histocolor));
                ProfilesX[i]->SetLineWidth(2);
                ProfilesX[i]->Scale(1/ProfilesX[i]->Integral());
                std::ostringstream buff;
                buff<<h4d->GetYaxis()->GetBinCenter(i+1);
                std::string bin = buff.str();
                leg->AddEntry(ProfilesX[i],("R: "+bin).c_str(), "L");
                ProfilesX[i]->GetXaxis()->SetTitle("#chi^{2} - Y");
                ProfilesX[i]->Draw("same");
        }
   leg->Draw("same");

   for(int i=0;i<100;i++) {
                int histocolor=(float)nColors / nHistos * i;
              //  ProfilesX[i]->SetLineColor(gStyle->GetColorPalette(histocolor));
                ProfilesY[i]->SetLineWidth(2);
                ProfilesY[i]->Scale(1/ProfilesX[i]->Integral());
                std::ostringstream buff;
                buff<<h4d->GetYaxis()->GetBinCenter(i+1);
                std::string bin = buff.str();
                leg->AddEntry(ProfilesX[i],("R: "+bin).c_str(), "L");
                ProfilesY[i]->GetXaxis()->SetTitle("#chi^{2} - Y");
                ProfilesY[i]->Draw("same");
        }
   leg->Draw("same");
  


   TCanvas *c4 = new TCanvas("ChiX ISOEff", "ChiX ISOEff",4586,195,884,613);
 //  SetCanvas(c4);
//  gStyle->SetPalette(55);
   chieffX->GetXaxis()->SetTitle("#chi^{2}(X)");
   chieffX->GetYaxis()->SetTitle("R [GV]");
   chieffX->Draw("cont1z");

   TCanvas *c5 = new TCanvas("ChiY ISOEff", "ChiY ISOEff",4586,195,884,613);
   //  SetCanvas(c4);
   //  gStyle->SetPalette(55);
   chieffY->GetXaxis()->SetTitle("#chi^{2}(Y)");
   chieffY->GetYaxis()->SetTitle("R [GV]");
   chieffY->Draw("cont1z");


   finalResults.Add(c2);
   finalResults.Add(c3);
   finalResults.Add(c4);
   finalResults.Add(c5);

   finalResults.writeObjsInFolder("Cleaning");

   DrawQvsRBeta(h6,h9,h6_,h9_,"Upper TOF Z",finalResults);
   DrawQvsRBeta(h7,h10,h7_,h10_,"Lower TOF Z",finalResults);
   DrawQvsRBeta(h8,h11,h8_,h11_,"Inner Z",finalResults);

}


void DrawBDT(FileSaver finalHistos,FileSaver finalResults){

   TH2F * h1 = (TH2F*)finalHistos.Get("RICHBDTvsMassNaFP/RICHBDTvsMassNaFP");
   TH2F * h2 = (TH2F*)finalHistos.Get("RICHBDTvsMassAglP/RICHBDTvsMassAglP");
   TH2F * h3 = (TH2F*)finalHistos.Get("RICHBDTvsMassNaFD/RICHBDTvsMassNaFD");
   TH2F * h4 = (TH2F*)finalHistos.Get("RICHBDTvsMassAglD/RICHBDTvsMassAglD"); 
   TH2F * h5 = (TH2F*)finalHistos.Get("RICHBDTvsBetaNaF/RICHBDTvsBetaNaF");
   TH2F * h6 = (TH2F*)finalHistos.Get("RICHBDTvsBetaAgl/RICHBDTvsBetaAgl"); 
   TH1F * Eff = (TH1F*) finalHistos.Get("BDTDiscr/BDTDiscr");

	TCanvas * c1 = new TCanvas("RICHBDTvsMassNaFP");
	c1->cd();
	PlotTH2F(gPad, h1, "Mass [GeV/c^{2}]","BDT discriminant","colz");
	TCanvas * c2 = new TCanvas("RICHBDTvsMassAglP");
	c2->cd();
	PlotTH2F(gPad, h2, "Mass [GeV/c^{2}]","BDT discriminant","colz");
	TCanvas * c3 = new TCanvas("RICHBDTvsMassNaFD");
	c3->cd();
	PlotTH2F(gPad, h3, "Mass [GeV/c^{2}]","BDT discriminant","colz");
	TCanvas * c4 = new TCanvas("RICHBDTvsMassAglD");
	c4->cd();
	PlotTH2F(gPad, h4, "Mass [GeV/c^{2}]","BDT discriminant","colz");
	TCanvas * c5 = new TCanvas("RICHBDTvsBetaNaF");
	c5->cd();
	PlotTH2F(gPad, h5, "#beta NaF","BDT discriminant","colz");
	TCanvas * c6 = new TCanvas("RICHBDTvsBetaAgl");
	c6->cd();
	PlotTH2F(gPad, h6, "#beta Agl","BDT discriminant","colz");


	TCanvas * c7 = new TCanvas("BDTcutEfficiency");
	c7->cd();                    
	TGraph * EffvsCut = new TGraph();
	for(int i=0;i<Eff->GetNbinsX();i++){
		EffvsCut->SetPoint(i+1,Eff->GetBinCenter(i+1),Eff->Integral(i+1,Eff->GetNbinsX())/Eff->Integral());
	}
	EffvsCut->SetMarkerStyle(8);
	EffvsCut->SetMarkerColor(1);
	EffvsCut->SetMarkerSize(1);                                  
	EffvsCut->Draw("APC");

	finalResults.Add(c1);
	finalResults.Add(c2);
	finalResults.Add(c3);
	finalResults.Add(c4);
	finalResults.Add(c5);
	finalResults.Add(c6);
	finalResults.Add(c7);

	finalResults.writeObjsInFolder("RichBDT");




}



void DrawMasses(FileSaver finalHistos,FileSaver finalResults){


   TH1F * h1 = (TH1F*)finalHistos.Get("MassTOF_noSel/MassTOF_noSel"	);
   TH1F * h2 = (TH1F*)finalHistos.Get("MassTOF_Sel/MassTOF_Sel"     );
   TH1F * h3 = (TH1F*)finalHistos.Get("MassTOF_Qual/MassTOF_Qual"    );
   TH1F * h4 = (TH1F*)finalHistos.Get("MassNaF_noSel/MassNaF_noSel"	); 
   TH1F * h5 = (TH1F*)finalHistos.Get("MassNaF_Sel/MassNaF_Sel"    );
   TH1F * h6 = (TH1F*)finalHistos.Get("MassNaF_Qual/MassNaF_Qual"   ); 
   TH1F * h7 = (TH1F*)finalHistos.Get("MassAgl_noSel/MassAgl_noSel");
   TH1F * h8 = (TH1F*)finalHistos.Get("MassAgl_Sel/MassAgl_Sel"  );
   TH1F * h9 = (TH1F*)finalHistos.Get("MassAgl_Qual/MassAgl_Qual");

   cout<<h1<<" "<<h2<<" "<<h3<<" "<<h4<<endl; 	
   cout<<h5<<" "<<h6<<" "<<h7<<" "<<h8<<" "<<h9<<endl; 	
 
   TCanvas * c1 = new TCanvas("Cleaning Mass TOF");
        c1->cd();	

	h1->SetFillColor(17);
        h1->SetFillStyle(3000);
	PlotDistribution(gPad,h1,"Mass [GeV/c^{2}]", "Counts", 1, "same",1,1.2*h1->GetBinContent(h1->GetMaximumBin()));
	h2->SetFillColor(15);
        h2->SetFillStyle(3000);
	PlotDistribution(gPad,h2,"Mass [GeV/c^{2}]", "Counts", 1, "same",1,1.2*h1->GetBinContent(h1->GetMaximumBin()));
	h3->SetFillColor(1);
        h3->SetFillStyle(3000);
	PlotDistribution(gPad,h3,"Mass [GeV/c^{2}]", "Counts", 1, "same",1,1.2*h1->GetBinContent(h1->GetMaximumBin()));

	finalResults.Add(c1);
	finalResults.writeObjsInFolder("Cleaning Mass");

   TCanvas * c2 = new TCanvas("Cleaning Mass NaF");
        c2->cd();	

	h4->SetFillColor(17);
        h4->SetFillStyle(3000);
	PlotDistribution(gPad,h4,"Mass [GeV/c^{2}]", "Counts", 1, "same",1,1.2*h4->GetBinContent(h1->GetMaximumBin()));
	h5->SetFillColor(15);
        h5->SetFillStyle(3000);
	PlotDistribution(gPad,h5,"Mass [GeV/c^{2}]", "Counts", 1, "same",1,1.2*h4->GetBinContent(h1->GetMaximumBin()));
	h6->SetFillColor(1);
        h6->SetFillStyle(3000);
	PlotDistribution(gPad,h6,"Mass [GeV/c^{2}]", "Counts", 1, "same",1,1.2*h4->GetBinContent(h1->GetMaximumBin()));

	finalResults.Add(c2);
	finalResults.writeObjsInFolder("Cleaning Mass");

   TCanvas * c3 = new TCanvas("Cleaning Mass Agl");
        c3->cd();	

	h7->SetFillColor(17);
        h7->SetFillStyle(3000);
	PlotDistribution(gPad,h7,"Mass [GeV/c^{2}]", "Counts", 1, "same",1,1.2*h7->GetBinContent(h1->GetMaximumBin()));
	h8->SetFillColor(15);
        h8->SetFillStyle(3000);
	PlotDistribution(gPad,h8,"Mass [GeV/c^{2}]", "Counts", 1, "same",1,1.2*h7->GetBinContent(h1->GetMaximumBin()));
	h9->SetFillColor(1);
        h9->SetFillStyle(3000);
	PlotDistribution(gPad,h9,"Mass [GeV/c^{2}]", "Counts", 1, "same",1,1.2*h7->GetBinContent(h1->GetMaximumBin()));

	finalResults.Add(c3);
	finalResults.writeObjsInFolder("Cleaning Mass");


}

TH1D * GetMeans(TH2F * h1){
	string nameh1 = h1->GetName();
	nameh1 = nameh1 + "_1";
	TH1D * h1_1 = (TH1D*)gDirectory->Get(nameh1.c_str());	
	return h1_1;

}


void DrawBetaRes(FileSaver finalHistos,FileSaver finalResults){

   TH2F * h1 = (TH2F*)finalHistos.Get("BetagenvsBetaMeasTOFP/BetagenvsBetaMeasTOFP");
   TH2F * h2 = (TH2F*)finalHistos.Get("BetagenvsBetaMeasNaFP/BetagenvsBetaMeasNaFP");
   TH2F * h3 = (TH2F*)finalHistos.Get("BetagenvsBetaMeasAglP/BetagenvsBetaMeasAglP");
   TH2F * h4 = (TH2F*)finalHistos.Get("BetagenvsBetaMeasTOFD/BetagenvsBetaMeasTOFD"); 
   TH2F * h5 = (TH2F*)finalHistos.Get("BetagenvsBetaMeasNaFD/BetagenvsBetaMeasNaFD");
   TH2F * h6 = (TH2F*)finalHistos.Get("BetagenvsBetaMeasAglD/BetagenvsBetaMeasAglD"); 
	
	TF1 * ideal = new TF1("ideal","x",0,1.5); 

	h1->FitSlicesY(0,50,h1->GetNbinsX(),0);
	h2->FitSlicesY(0,50,h2->GetNbinsX(),0);
	h3->FitSlicesY(0,50,h3->GetNbinsX(),0);
	h4->FitSlicesY(0,50,h4->GetNbinsX(),0);
	h5->FitSlicesY(0,50,h5->GetNbinsX(),0);
	h6->FitSlicesY(0,50,h6->GetNbinsX(),0);

	TH1D * h1_1 = GetMeans(h1);	
	TH1D * h2_1 = GetMeans(h2);	
	TH1D * h3_1 = GetMeans(h3);	
	TH1D * h4_1 = GetMeans(h4);	
	TH1D * h5_1 = GetMeans(h5);	
	TH1D * h6_1 = GetMeans(h6);	




	TCanvas * c1 = new TCanvas("Protons MC #beta shift");
	c1->cd();
//	PlotTH2F(gPad, h1, "#beta_{gen}","#beta_{meas}","colz");
	h1_1->Draw();
	ideal->Draw("same");
	TCanvas * c2 = new TCanvas("Deuterons");
	c2->cd();
//	PlotTH2F(gPad, h4, "#beta_{gen}","#beta_{meas}","colz");
	h4_1->Draw();
	ideal->Draw("same");
	finalResults.Add(c1);
	finalResults.Add(c2);
	finalResults.writeObjsInFolder("ResponseBeta/BetaTOF");


	TCanvas * c3 = new TCanvas("Protons MC #beta shift");
	c3->cd();
	//PlotTH2F(gPad, h2, "#beta_{gen}","#beta_{meas}","colz");
	h2_1->Draw();
	ideal->Draw("same");
	TCanvas * c4 = new TCanvas("Deuterons MC #beta shift");
	c4->cd();
	//PlotTH2F(gPad, h5, "#beta_{gen}","#beta_{meas}","colz");
	h5_1->Draw();
	ideal->Draw("same");
	finalResults.Add(c3);
	finalResults.Add(c4);
	finalResults.writeObjsInFolder("ResponseBeta/BetaNaF");

	TCanvas * c5 = new TCanvas("Protons MC #beta shift");
	c5->cd();
	//PlotTH2F(gPad, h3, "#beta_{gen}","#beta_{meas}","colz");
	h3_1->Draw();
	ideal->Draw("same");
	TCanvas * c6 = new TCanvas("Deuterons MC #beta shift");
	c6->cd();
	//PlotTH2F(gPad, h6, "#beta_{gen}","#beta_{meas}","colz");
	h6_1->Draw();
	ideal->Draw("same");
	finalResults.Add(c5);
	finalResults.Add(c6);
	finalResults.writeObjsInFolder("ResponseBeta/BetaAgl");


}
