#include "MacrosForPlots.h"

float GammaFromBeta (float beta){ return 1/(sqrt(1-beta*beta)); }
float EkFromBeta (float beta, float mass) {return GammaFromBeta(beta)*mass - mass;}


void DrawBins(TVirtualPad * c1, TH2F* histo, std::vector<float> orizontal_set,std::vector<float> vertical_set,int col_oriz, int col_vert) {
		c1->cd();
		std::vector< TLine *> oriz_lines;
		std::vector< TLine *> vert_lines;

		for(int i =0; i<orizontal_set.size(); i++) {
			TLine * oriz = new TLine(histo->GetXaxis()->GetBinLowEdge(1),orizontal_set[i],vertical_set[i],orizontal_set[i]);    
			TLine * vert = new TLine(vertical_set[i],histo->GetYaxis()->GetBinLowEdge(1),vertical_set[i],orizontal_set[i]);
			oriz_lines.push_back(oriz);
			vert_lines.push_back(vert);
		}
		for(int i =0; i<orizontal_set.size(); i++) {
			oriz_lines[i]->SetLineColor(col_oriz);
			vert_lines[i]->SetLineColor(col_vert);
	
			oriz_lines[i]->Draw("same");
			vert_lines[i]->Draw("same");
		}	
}

TH2F * RebinHisto(TH2F * histo,std::vector<float> orizontal_set,std::vector<float> vertical_set,std::string name){

	TH2F * new_histo = new TH2F(name.c_str(),name.c_str(),orizontal_set.size(),0,orizontal_set.size(),vertical_set.size(),0,vertical_set.size())	;
	
	for(int i=1;i<orizontal_set.size()-1;i++)
		for(int j=1;j<vertical_set.size()-1;j++){
			float bincontent=histo->Integral(histo->GetXaxis()->FindBin(vertical_set[j]),histo->GetXaxis()->FindBin(vertical_set[j+1]),histo->GetYaxis()->FindBin(orizontal_set[i]),histo->GetYaxis()->FindBin(orizontal_set[i+1]));		
			new_histo->SetBinContent(j+1,i+1,bincontent);
		}
	return new_histo;
}


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

TH1D * ReduceOnTimeHisto(TH1D *OnTime){
		std::string name = OnTime->GetTitle();
		TH1D * OnTime_red = new TH1D((name +"_red").c_str()      , (name+"_red").c_str()       , 2,0,2);
		OnTime_red   ->SetBinContent(1,OnTime->Integral(OnTime->FindBin(0),OnTime->FindBin(4.5)));
		OnTime_red   ->SetBinContent(2,OnTime->Integral(OnTime->FindBin(4.5),OnTime->FindBin(15)));
		return OnTime_red;

	}



void DrawBDT(FileSaver finalHistos,FileSaver finalResults){
	return;
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

TH1D * GetSTD(TH2F * h1){
	string nameh1 = h1->GetName();
	nameh1 = nameh1 + "_2";
	TH1D * h1_1 = (TH1D*)gDirectory->Get(nameh1.c_str());	
	return h1_1;

}

void DrawMassRes(FileSaver finalHistos,FileSaver finalResults){

	TH2F * h1 = (TH2F*)finalHistos.Get("MassvsBetaTOF/MassvsBetaTOF");
	TH2F * h2 = (TH2F*)finalHistos.Get("MassvsBetaNaF/MassvsBetaNaF");
	TH2F * h3 = (TH2F*)finalHistos.Get("MassvsBetaAgl/MassvsBetaAgl");

	cout<<h1<<" "<<h2<<" "<<h3<<endl;
	h1->FitSlicesY(0,20,h1->GetNbinsX(),5);
	h2->FitSlicesY(0,20,h2->GetNbinsX(),5);
	h3->FitSlicesY(0,20,h3->GetNbinsX(),5);

	TH1D * h1_1 = GetSTD(h1);	
	TH1D * h2_1 = GetSTD(h2);	
	TH1D * h3_1 = GetSTD(h3);	

	h1_1->Rebin(4);
	h2_1->Rebin(4);
	h3_1->Rebin(4);;

	h1_1->Scale(1/4.);
	h2_1->Scale(1/4.);
	h3_1->Scale(1/4.);;

	h1_1->Smooth();
	h2_1->Smooth();
	h3_1->Smooth();;


	TH2F * Frame = new TH2F("Frame","Frame", 1000,0,25,1000,0,1);
	TGraphErrors * ResoTOF = new TGraphErrors();
	TGraphErrors * ResoNaF = new TGraphErrors();
	TGraphErrors * ResoAgl = new TGraphErrors();
	int j=0;
	for(int i=0;i<h1_1->GetNbinsX();i++) if(h1_1->GetBinContent(i+1)>0 && h1_1->GetBinError(i+1)<0.05){ ResoTOF->SetPoint(j,EkFromBeta(h1_1->GetBinCenter(i+1),0.938) ,2*h1_1->GetBinContent(i+1)); j++;}
	j=0;                                                                                        
	for(int i=0;i<h1_1->GetNbinsX();i++) if(h2_1->GetBinContent(i+1)>0 && h2_1->GetBinError(i+1)<0.05){ ResoNaF->SetPoint(j,EkFromBeta(h2_1->GetBinCenter(i+1),0.938) ,2*h2_1->GetBinContent(i+1)); j++;}
	j=0;                                                                                        
	for(int i=0;i<h1_1->GetNbinsX();i++) if(h2_1->GetBinContent(i+1)>0 && h2_1->GetBinError(i+1)<0.05){ ResoAgl->SetPoint(j,EkFromBeta(h3_1->GetBinCenter(i+1),0.938) ,2*h3_1->GetBinContent(i+1)); j++;}

	
	TCanvas *c1 = new TCanvas("MassresoTOF");
	gPad->SetLogx();
	Frame->Draw();
	ResoTOF->SetLineColor(2);
        ResoNaF->SetLineColor(3);
        ResoAgl->SetLineColor(4);

	ResoTOF->SetLineWidth(2);
        ResoNaF->SetLineWidth(2);
        ResoAgl->SetLineWidth(2);

	ResoTOF->Draw("PLsame");
	ResoNaF->Draw("PLsame");
	ResoAgl->Draw("PLsame");	

	finalResults.Add(c1);
	finalResults.writeObjsInFolder("MassReso/TOF");
}

void DrawBetaSmear(FileSaver finalHistos,FileSaver finalResults){

   TH1F * htof 	  = (TH1F*)finalHistos.Get("BetaTOF_Z1/BetaTOF_Z1");
   TH1F * hnaf    = (TH1F*)finalHistos.Get("BetaNaF_Z1/BetaNaF_Z1");
   TH1F * hagl 	  = (TH1F*)finalHistos.Get("BetaAgl_Z1/BetaAgl_Z1");
   TH1F * htof_he = (TH1F*)finalHistos.Get("BetaTOF_Z2/BetaTOF_Z2");
   TH1F * hnaf_he = (TH1F*)finalHistos.Get("BetaNaF_Z2/BetaNaF_Z2");
   TH1F * hagl_he = (TH1F*)finalHistos.Get("BetaAgl_Z2/BetaAgl_Z2");
   TH1F * htof_mc = (TH1F*)finalHistos.Get("BetaTOF_P_MC/BetaTOF_P_MC"); 
   TH1F * hnaf_mc = (TH1F*)finalHistos.Get("BetaNaF_P_MC/BetaNaF_P_MC");
   TH1F * hagl_mc = (TH1F*)finalHistos.Get("BetaAgl_P_MC/BetaAgl_P_MC"); 
   TH1F * htof_mc_he = (TH1F*)finalHistos.Get("BetaTOF_He_MC/BetaTOF_He_MC"); 
   TH1F * hnaf_mc_he = (TH1F*)finalHistos.Get("BetaNaF_He_MC/BetaNaF_He_MC");
   TH1F * hagl_mc_he = (TH1F*)finalHistos.Get("BetaAgl_He_MC/BetaAgl_He_MC"); 
   TH1F * hagl_mc_s = (TH1F*)finalHistos.Get("BetaAgl_HE_MC_Smeared/BetaAgl_HE_MC_Smeared"); 
   TH1F * htof_mc_s = (TH1F*)finalHistos.Get("BetaTOF_HE_MC_Smeared/BetaTOF_HE_MC_Smeared"); 
 


	htof	->Scale(1/htof->Integral());
	hnaf	->Scale(1/hnaf->Integral());
	hagl	->Scale(1/hagl->Integral());
	htof_he ->Scale(1/htof_he->Integral());
	hnaf_he ->Scale(1/hnaf_he->Integral());
	hagl_he ->Scale(1/hagl_he->Integral());

	htof_mc ->Scale(1/htof_mc->Integral());
	hnaf_mc ->Scale(1/hnaf_mc->Integral());
	hagl_mc ->Scale(1/hagl_mc->Integral());
	htof_mc_he ->Scale(1/htof_mc_he->Integral());
	hnaf_mc_he ->Scale(1/hnaf_mc_he->Integral());
	hagl_mc_he ->Scale(1/hagl_mc_he->Integral());




	float norm1_P	=htof->GetBinContent(htof->GetMaximumBin())*0.6667;
	float norm2_P	=htof->GetBinContent(htof->GetMaximumBin())*0.332;
	float norm3_P	=htof->GetBinContent(htof->GetMaximumBin())*0.000095;
	float mean_P	=1;
	float mean2_P	=1;
	float mean3_P	=1.027;
	float sigma_P	=0.0311;
	float sigma2_P	=1.25*0.0311;
	float sigma3_P	=2.29*0.0311;

	cout<<norm1_P<<" "<<norm2_P<<endl;
	TF1 * modelTOF_P = new TF1("ModelTOF_P","gaus(0)+gaus(3)+gaus(6)",0.1,2);

	modelTOF_P->SetParLimits(0,0.7*	norm1_P 	,1.3*	norm1_P 	); 
	modelTOF_P->SetParLimits(1,0.95*	mean_P		,1.05*	mean_P		);
	modelTOF_P->SetParLimits(2,0.7*	sigma_P		,1.3*	sigma_P		);
	modelTOF_P->SetParLimits(3,0.7*	norm2_P 	,1.3*	norm2_P 	);
	modelTOF_P->SetParLimits(4,0.99*	mean2_P		,1.01*	mean2_P		);
	modelTOF_P->SetParLimits(5,0.5*	sigma2_P	,1.5*	sigma2_P	);
	modelTOF_P->SetParLimits(6,0.1*	norm3_P 	,4*	norm3_P 	);
	modelTOF_P->FixParameter(7,mean3_P);
	modelTOF_P->SetParLimits(8,0.5*	sigma3_P	,1.5*	sigma3_P	);

	htof->Fit("ModelTOF_P","");

	float norm1_He	=htof_he->GetBinContent(htof_he->GetMaximumBin())*0.6667;
	float norm2_He	=htof_he->GetBinContent(htof_he->GetMaximumBin())*0.332;
	float norm3_He	=htof_he->GetBinContent(htof_he->GetMaximumBin())*0.000095;
	float mean_He	=1;
	float mean2_He	=1;
	float mean3_He	=1.027;
	float sigma_He	=0.019;
	float sigma2_He	=1.25*0.019;
	float sigma3_He	=2.29*0.019;

	TF1 * modelTOF_He = new TF1("ModelTOF_He","gaus(0)+gaus(3)+gaus(6)",0.1,2);

	modelTOF_He->SetParLimits(0,0.7*	norm1_He 	,1.3*	norm1_He 	); 
	modelTOF_He->SetParLimits(1,0.95*	mean_He		,1.05*	mean_He		);
	modelTOF_He->SetParLimits(2,0.7*	sigma_He		,1.3*	sigma_He		);
	modelTOF_He->SetParLimits(3,0.7*	norm2_He 	,1.3*	norm2_He 	);
	modelTOF_He->SetParLimits(4,0.99*	mean2_He		,1.01*	mean2_He		);
	modelTOF_He->SetParLimits(5,0.5*	sigma2_He	,1.5*	sigma2_He	);
	modelTOF_He->SetParLimits(6,0.1*	norm3_He 	,4*	norm3_He 	);
	modelTOF_He->FixParameter(7,mean3_He);
	modelTOF_He->SetParLimits(8,0.5*	sigma3_He	,1.5*	sigma3_He	);

	htof_he->Fit("ModelTOF_He","");

	htof	->SetMarkerStyle(8	);
	hnaf	->SetMarkerStyle(8	);
	hagl	->SetMarkerStyle(8	);
	htof_he ->SetMarkerStyle(8	);
	hnaf_he ->SetMarkerStyle(8	);
	hagl_he ->SetMarkerStyle(8	);

	htof	->SetLineWidth(2	);
	hnaf	->SetLineWidth(2	);
	hagl	->SetLineWidth(2	);
	htof_he ->SetLineWidth(2	);
	hnaf_he ->SetLineWidth(2	);
	hagl_he ->SetLineWidth(2	);

	htof_mc	->SetLineWidth(2	);
	hnaf_mc	->SetLineWidth(2	);
	hagl_mc	->SetLineWidth(2	);
	htof_mc_he ->SetLineWidth(2	);
	hnaf_mc_he ->SetLineWidth(2	);
	hagl_mc_he ->SetLineWidth(2	);

	htof_mc	->SetLineColor(2	);
	hnaf_mc	->SetLineColor(2	);
	hagl_mc	->SetLineColor(2	);
	htof_mc_he ->SetLineColor(3	);
	hnaf_mc_he ->SetLineColor(3	);
	hagl_mc_he ->SetLineColor(3	);
	modelTOF_P->SetLineStyle(2);
	modelTOF_He->SetLineStyle(2);




	finalResults.Add(htof   	);
	finalResults.Add(hnaf   	);
	finalResults.Add(hagl   	);
	finalResults.Add(htof_he	);
	finalResults.Add(hnaf_he	);
	finalResults.Add(hagl_he	);
	finalResults.Add(htof_mc   	);
	finalResults.Add(hnaf_mc   	);
	finalResults.Add(hagl_mc   	);
	finalResults.Add(htof_mc_he	);
	finalResults.Add(hnaf_mc_he	);
	finalResults.Add(hagl_mc_he	);


	finalResults.Add(modelTOF_P);
	finalResults.Add(modelTOF_He);








	finalResults.writeObjsInFolder("BetaSmear");



}

void DrawBetaRes(FileSaver finalHistos,FileSaver finalResults){

   TH2F * h1 = (TH2F*)finalHistos.Get("BetagenvsBetaMeasTOFP/BetagenvsBetaMeasTOFP");
   TH2F * h2 = (TH2F*)finalHistos.Get("BetagenvsBetaMeasNaFP/BetagenvsBetaMeasNaFP");
   TH2F * h3 = (TH2F*)finalHistos.Get("BetagenvsBetaMeasAglP/BetagenvsBetaMeasAglP");
   TH2F * h4 = (TH2F*)finalHistos.Get("BetagenvsBetaMeasTOFD/BetagenvsBetaMeasTOFD"); 
   TH2F * h5 = (TH2F*)finalHistos.Get("BetagenvsBetaMeasNaFD/BetagenvsBetaMeasNaFD");
   TH2F * h6 = (TH2F*)finalHistos.Get("BetagenvsBetaMeasAglD/BetagenvsBetaMeasAglD"); 
   TH2F * h7 = (TH2F*)finalHistos.Get("BetagenvsBetaMeasTOFHe/BetagenvsBetaMeasTOFHe");
   TH2F * h8 = (TH2F*)finalHistos.Get("BetagenvsBetaMeasNaFHe/BetagenvsBetaMeasNaFHe");
   TH2F * h9 = (TH2F*)finalHistos.Get("BetagenvsBetaMeasAglHe/BetagenvsBetaMeasAglHe");
   
	TF1 * ideal = new TF1("ideal","x",0,50); 
	ideal->SetLineStyle(2);

	h1->FitSlicesY(0,50,h1->GetNbinsX(),0);
	h2->FitSlicesY(0,50,h2->GetNbinsX(),0);
	h3->FitSlicesY(0,50,h3->GetNbinsX(),0);
	h4->FitSlicesY(0,50,h4->GetNbinsX(),0);
	h5->FitSlicesY(0,50,h5->GetNbinsX(),0);
	h6->FitSlicesY(0,50,h6->GetNbinsX(),0);
	h7->FitSlicesY(0,50,h4->GetNbinsX(),0);
	h8->FitSlicesY(0,50,h5->GetNbinsX(),0);
	h9->FitSlicesY(0,50,h6->GetNbinsX(),0);


	TH1D * h1_1 = GetMeans(h1);	
	TH1D * h2_1 = GetMeans(h2);	
	TH1D * h3_1 = GetMeans(h3);	
	TH1D * h4_1 = GetMeans(h4);	
	TH1D * h5_1 = GetMeans(h5);	
	TH1D * h6_1 = GetMeans(h6);	
	TH1D * h7_1 = GetMeans(h7);	
	TH1D * h8_1 = GetMeans(h8);	
	TH1D * h9_1 = GetMeans(h9);	


        TF1 * fitfuncP  = new TF1("funcP","x*(1-[0]/x^[1])",0,1); 
	TF1 * fitfuncD = new TF1("funcD","x-0.4*(x-x*(1-[0]/x^[1])) ",0,1);
	TF1 * fitfuncHe = new TF1("funcHe","x-0.4*(x-x*(1-[0]/x^[1])) ",0,1);


	TCanvas * c1 = new TCanvas("Proton slowdown model");
	c1->cd();
	c1->SetLogz();
	((TH1F*)h1_1->Clone())->Fit("funcP","Q");

	PlotTH2F(gPad, h1, "#beta_{gen}","#beta_{meas}","colz");
	ideal->Draw("same");
	h1_1->SetMarkerStyle(8);
	h1_1->SetMarkerSize(1);
	h1_1->Draw("Psame");
	fitfuncP->Draw("same");

	TCanvas * c2 = new TCanvas("Deuterons slowdown model ");
	c2->cd();
	c2->SetLogz();
	
	fitfuncD->SetParameter(0,fitfuncP->GetParameter(0));
	fitfuncD->SetParameter(1,fitfuncP->GetParameter(1));
	cout<<"TOF slow down parameters; "<<endl;
	cout<<"par 0: "<<fitfuncP->GetParameter(0)<<endl;
	cout<<"par 1: "<<fitfuncP->GetParameter(1)<<endl;
	
	PlotTH2F(gPad, h4, "#beta_{gen}","#beta_{meas}","colz");
	h4_1->SetMarkerStyle(8);
	h4_1->SetMarkerSize(1);
	h4_1->Draw("Psame");
	ideal->Draw("same");
	fitfuncD->Draw("same");

	TCanvas * c1_ = new TCanvas("Helium slowdown model");
	c1_->cd();
	c1_->SetLogz();
	((TH1F*)h7_1->Clone())->Fit("funcHe","Q");

	PlotTH2F(gPad, h7, "#beta_{gen}","#beta_{meas}","colz");
	ideal->Draw("same");
	h7_1->SetMarkerStyle(8);
	h7_1->SetMarkerSize(1);
	h7_1->Draw("Psame");
	fitfuncHe->Draw("same");



	TCanvas * BinningTOF = new TCanvas("binningTOF");
	BinningTOF->Divide(2,0);

	std::vector<float> orizontal_set=Global.GetToFPBins().BetaBins ();
	std::vector<float> vertical_set_P=Global.GetToFPBins().BetaTOIBins () ;
	std::vector<float> vertical_set_D=Global.GetToFDBins().BetaTOIBins () ;

	BinningTOF->cd(1);
	gPad->SetLogz();
	PlotTH2F(gPad, h1, "#beta_{gen}","#beta_{meas}","colz");
	DrawBins(gPad,h1,orizontal_set,vertical_set_P,1,2);	
	ideal->Draw("same");
	fitfuncP->Draw("same");
	BinningTOF->cd(2);
	PlotTH2F(gPad, h4, "#beta_{gen}","#beta_{meas}","colz");
	gPad->SetLogz();
	DrawBins(gPad,h4,orizontal_set,vertical_set_D,1,4);	
	ideal->Draw("same");
	fitfuncD->Draw("same");


	TCanvas * MigrrectTOF = new TCanvas("Migr_rectTOF");
	MigrrectTOF->Divide(2,0);
	
	MigrrectTOF->cd(1);
	TH2F * Migrrect_P = RebinHisto(h1,orizontal_set,orizontal_set,"Migrrect_P");
	Migrrect_P->Draw("colz");
	ideal->Draw("same");
	MigrrectTOF->cd(2);
	TH2F * Migrrect_D = RebinHisto(h4,orizontal_set,orizontal_set,"Migrrect_D");
	Migrrect_D->Draw("colz");
	ideal->Draw("same");
	
	TCanvas * MigrmodTOF = new TCanvas("Migr_modTOF");
	MigrmodTOF->Divide(2,0);
	
	MigrmodTOF->cd(1);
	TH2F * Migrmod_P = RebinHisto(h1,orizontal_set,vertical_set_P,"Migrmod_P");
	Migrmod_P->Draw("colz");
	ideal->Draw("same");
	MigrmodTOF->cd(2);
	TH2F * Migrmod_D = RebinHisto(h4,orizontal_set,vertical_set_D,"Migrmod_D");
	Migrmod_D->Draw("colz");
	ideal->Draw("same");
	
	finalResults.Add(c1);
	finalResults.Add(c2);
	finalResults.Add(c1_);
	finalResults.Add(BinningTOF);
	finalResults.Add(MigrrectTOF);
	finalResults.Add(MigrmodTOF);
	finalResults.writeObjsInFolder("ResponseBeta/BetaTOF");





	TCanvas * c3 = new TCanvas("Proton slowdown model");
	c3->cd();
	c3->SetLogz();
	((TH1F*)h2_1->Clone())->Fit("funcP","Q");

	
	PlotTH2F(gPad, h2, "#beta_{gen}","#beta_{meas}","colz");
	ideal->Draw("same");
	h2_1->SetMarkerStyle(8);
	h2_1->SetMarkerSize(1);
	h2_1->Draw("Psame");
	
	TCanvas * c4 = new TCanvas("Deuteron slowdown modle");
	c4->cd();
        c4->SetLogz();
	fitfuncD->SetParameter(0,fitfuncP->GetParameter(0));
	fitfuncD->SetParameter(1,fitfuncP->GetParameter(1));

	cout<<"NaF slow down parameters; "<<endl;
	cout<<"par 0: "<<fitfuncP->GetParameter(0)<<endl;
	cout<<"par 1: "<<fitfuncP->GetParameter(1)<<endl;

	PlotTH2F(gPad, h5, "#beta_{gen}","#beta_{meas}","colz");
	h5_1->SetMarkerStyle(8);
	h5_1->SetMarkerSize(1);
	h5_1->Draw("Psame");
	ideal->Draw("same");
	fitfuncD->Draw("same");
	finalResults.Add(c3);
	finalResults.Add(c4);
	finalResults.writeObjsInFolder("ResponseBeta/BetaNaF");



	TCanvas * c5 = new TCanvas("Proton slowdown model");
	c5->cd();
	c5->SetLogz();
	
	((TH1F*)h3_1->Clone())->Fit("funcP","Q");

	PlotTH2F(gPad, h3, "#beta_{gen}","#beta_{meas}","colz");
	h3_1->SetMarkerStyle(8);
	h3_1->SetMarkerSize(1);
	h3_1->Draw("same");
	ideal->Draw("same");
	TCanvas * c6 = new TCanvas("Deuteron slowdown  model");
	c6->cd();
	c6->SetLogz();
	fitfuncD->SetParameter(0,fitfuncP->GetParameter(0));
	fitfuncD->SetParameter(1,fitfuncP->GetParameter(1));

	PlotTH2F(gPad, h6, "#beta_{gen}","#beta_{meas}","colz");
	h6_1->SetMarkerStyle(8);
	h6_1->SetMarkerSize(1);
	h6_1->Draw("Psame");
	ideal->Draw("same");
	fitfuncD->Draw("same");

	cout<<"Agl slow down parameters; "<<endl;
	cout<<"par 0: "<<fitfuncP->GetParameter(0)<<endl;
	cout<<"par 1: "<<fitfuncP->GetParameter(1)<<endl;

	finalResults.Add(c5);
	finalResults.Add(c6);
	finalResults.writeObjsInFolder("ResponseBeta/BetaAgl");
}


void DrawBetavsRig(FileSaver finalHistos,FileSaver finalResults){
	TH2F * RigvsBeta_TOFP = (TH2F*)finalHistos.Get("RigvsBeta_TOFP/RigvsBeta_TOFP");
	TH2F * RigvsBeta_NaFP = (TH2F*)finalHistos.Get("RigvsBeta_NaFP/RigvsBeta_NaFP");
	TH2F * RigvsBeta_AglP = (TH2F*)finalHistos.Get("RigvsBeta_AglP/RigvsBeta_AglP");
	TH2F * RigvsBeta_TOFD = (TH2F*)finalHistos.Get("RigvsBeta_TOFD/RigvsBeta_TOFD");
	TH2F * RigvsBeta_NaFD = (TH2F*)finalHistos.Get("RigvsBeta_NaFD/RigvsBeta_NaFD");
	TH2F * RigvsBeta_AglD = (TH2F*)finalHistos.Get("RigvsBeta_AglD/RigvsBeta_AglD");

	TH2F * RigvsBeta_TOF_data = (TH2F*)finalHistos.Get("RigvsBeta_TOF_data/RigvsBeta_TOF_data");
	TH2F * RigvsBeta_NaF_data = (TH2F*)finalHistos.Get("RigvsBeta_NaF_data/RigvsBeta_NaF_data");
	TH2F * RigvsBeta_Agl_data = (TH2F*)finalHistos.Get("RigvsBeta_Agl_data/RigvsBeta_Agl_data");
	TH2F * RigvsBeta_TOF_dataprim = (TH2F*)finalHistos.Get("RigvsBeta_TOF_dataprim/RigvsBeta_TOF_dataprim");
	TH2F * RigvsBeta_NaF_dataprim = (TH2F*)finalHistos.Get("RigvsBeta_NaF_dataprim/RigvsBeta_NaF_dataprim");
	TH2F * RigvsBeta_Agl_dataprim = (TH2F*)finalHistos.Get("RigvsBeta_Agl_dataprim/RigvsBeta_Agl_dataprim");

	TCanvas * c1 = new TCanvas("MC RvsBeta TOF");
	c1->cd();
	c1->SetLogz();
	RigvsBeta_TOFP->SetMarkerColor(7);
	RigvsBeta_TOFD->SetMarkerColor(13);
	RigvsBeta_TOFD->Draw();
	RigvsBeta_TOFP->Draw("same");
	
	std::vector<float> orizontal_set=Global.GetToFPBins().BetaBins ();
	std::vector<float> vertical_set_P=Global.GetToFPBins().RigBins () ;
	std::vector<float> vertical_set_D=Global.GetToFDBins().RigBins () ;
	
	DrawBins(c1,RigvsBeta_TOFP,orizontal_set,vertical_set_P,1,2);
	DrawBins(c1,RigvsBeta_TOFD,orizontal_set,vertical_set_D,1,4);

	TCanvas * c2 = new TCanvas("MC Binning TOF");
	c2->cd();
	c2->SetLogz();
	c2->Divide(2,1);
	
	c2->cd(1);
	TH2F * Rebinned_P = RebinHisto(RigvsBeta_TOFP,orizontal_set,vertical_set_P,"Rebinned_P");
	Rebinned_P->Draw("colz");

	c2->cd(2);
	TH2F * Rebinned_D = RebinHisto(RigvsBeta_TOFD,orizontal_set,vertical_set_D,"Rebinned_D");
	Rebinned_D->Draw("colz");


	finalResults.Add(c1);
	finalResults.Add(c2);
	finalResults.writeObjsInFolder("RigvsBeta/TOF_MC");

}



void DrawAcceptanceMatrix(FileSaver finalHistos,FileSaver finalResults){

	Variables * vars = new Variables();
	TH1F * Denominator = new TH1F("Denominator","MC spectrum",PRB.size(),0,PRB.size());
	TH1F * DenominatorW = new TH1F("DenominatorW","Reweighted",PRB.size(),0,PRB.size());
	
	float range = log(100)-log(0.5);
	
	Histogram Spectrum = vars->reweighter.getTo(); 


	for(int i=0;i<Denominator->GetNbinsX();i++){
		
		std::ostringstream ss;
		ss << fixed << setprecision(1) << PRB.RigBinCent(i);
		std::string s(ss.str());
		Denominator->SetBinContent(i+1,(log(PRB.RigBin(i+1))-log(PRB.RigBin(i)))/range );
		Denominator->GetXaxis()->SetBinLabel(i+1,s.c_str());
		DenominatorW->SetBinContent(i+1,Spectrum.integrate(PRB.RigBin(i),PRB.RigBin(i+1)));
		DenominatorW->GetXaxis()->SetBinLabel(i+1,s.c_str());
	}

	DenominatorW->Scale(1/DenominatorW->Integral());
	Denominator->Scale(1e6);
	DenominatorW->Scale(1e6);


	Denominator->SetLineWidth(2);
	DenominatorW->SetLineWidth(2);
	Denominator->SetLineColor(4);
	DenominatorW->SetLineColor(2);
	DenominatorW->GetXaxis()->LabelsOption("v");
	DenominatorW->GetXaxis()->SetTitle("R [GV]");
	DenominatorW->GetYaxis()->SetTitle("Gen. Events");


	TCanvas * c1 = new TCanvas("Denominator");
        c1->cd();
	gPad->SetGridx();
	DenominatorW->Draw();
	Denominator->Draw("same");
	finalResults.Add(c1);
        finalResults.writeObjsInFolder("Acceptance Matrix");




	std::string names[7] = {"N Tof Clusters","Q L tof","Q U tof","tof Coo chi","tof Time chi","N Tracks","On Time"};

	for(int i=0; i<7;i++) {
		TH2F * AcceptanceMatrix= (TH2F*)finalHistos.Get(("Acceptance Matrix v" + to_string(i) +"/Acceptance Matrix v" +to_string(i)).c_str());
		cout<<"Acc Matrix: "<<("Acceptance Matrix v" + to_string(i) +"/Acceptance Matrix v" +to_string(i)).c_str()<< AcceptanceMatrix<<endl; 
		AcceptanceMatrix->SetName(("tracker + " + names[i]).c_str());
		AcceptanceMatrix->SetTitle(("tracker + " + names[i]).c_str());
	
		for(int i=0;i<AcceptanceMatrix->GetNbinsX();i++){
			for(int j=0;j<AcceptanceMatrix->GetNbinsY();j++){
				if(DenominatorW->GetBinContent(i+1)>0 && (i+1)>6)
					AcceptanceMatrix->SetBinContent(i+1,j+1,AcceptanceMatrix->GetBinContent(i+1,j+1)/DenominatorW->GetBinContent(i+1));
				else AcceptanceMatrix->SetBinContent(i+1,j+1,AcceptanceMatrix->GetBinContent(i+1,j+1)/(2*Denominator->GetBinContent(i+1)));	
				}
		}

		for(int i=0;i<AcceptanceMatrix->GetNbinsX();i++){

			std::ostringstream ss;
			ss << fixed << setprecision(1) << PRB.RigBinCent(i);
			std::string s(ss.str());

			AcceptanceMatrix->GetXaxis()->SetBinLabel(i+1,s.c_str());
			AcceptanceMatrix->GetYaxis()->SetBinLabel(i+1,s.c_str());
		}
		AcceptanceMatrix->GetXaxis()->LabelsOption("v");

		AcceptanceMatrix->GetXaxis()->SetTitle("R_{gen} [GV]");
		AcceptanceMatrix->GetYaxis()->SetTitle("R_{meas} [GV]");

		AcceptanceMatrix->GetXaxis()->SetTitleOffset(1.35);
		AcceptanceMatrix->GetYaxis()->SetLabelSize(0.03);

		TCanvas * c2 = new TCanvas(("tracker + " + names[i]).c_str());
		c2->cd();
		gPad->SetGridx();
		gPad->SetGridy();
		c2->SetLogz();
		AcceptanceMatrix->Draw("colz");
		finalResults.Add(c2);

	}

        finalResults.writeObjsInFolder("Acceptance Matrix");

	TH2F * AcceptanceMatrix= (TH2F*)finalHistos.Get("Acceptance Matrix/Acceptance Matrix");

	for(int i=0;i<AcceptanceMatrix->GetNbinsX();i++){
			for(int j=0;j<AcceptanceMatrix->GetNbinsY();j++){
				if(DenominatorW->GetBinContent(i+1)>0 && (i+1)>6){
					AcceptanceMatrix->SetBinError(i+1,j+1, pow(AcceptanceMatrix->GetBinContent(i+1,j+1),0.5)/DenominatorW->GetBinContent(i+1));
					AcceptanceMatrix->SetBinContent(i+1,j+1,AcceptanceMatrix->GetBinContent(i+1,j+1)/DenominatorW->GetBinContent(i+1));
					}
				else{
					AcceptanceMatrix->SetBinError(i+1,j+1, pow(AcceptanceMatrix->GetBinContent(i+1,j+1),0.5)/(2*Denominator->GetBinContent(i+1)));
					AcceptanceMatrix->SetBinContent(i+1,j+1,AcceptanceMatrix->GetBinContent(i+1,j+1)/(2*Denominator->GetBinContent(i+1)));
					}
				}
		}

        TH1D*	Acceptance_gen  = AcceptanceMatrix->ProjectionX("projx",1,43);
	TH1D*   Acceptance_meas = AcceptanceMatrix->ProjectionY("projy",1,43);


	Acceptance_meas->GetXaxis()->LabelsOption("v");
	Acceptance_gen->GetXaxis()->LabelsOption("v");


	TCanvas * c4 = new TCanvas("Acceptance Matrix");
	for(int i=0;i<AcceptanceMatrix->GetNbinsX();i++){
		for(int j=0;j<AcceptanceMatrix->GetNbinsY();j++){
			if(DenominatorW->GetBinContent(i+1)>0 && (i+1)>6)
				AcceptanceMatrix->SetBinContent(i+1,j+1,AcceptanceMatrix->GetBinContent(i+1,j+1)/Denominator->GetBinContent(i+1));
			else AcceptanceMatrix->SetBinContent(i+1,j+1,AcceptanceMatrix->GetBinContent(i+1,j+1)/(2*Denominator->GetBinContent(i+1)));	
		}
	}

	for(int i=0;i<AcceptanceMatrix->GetNbinsX();i++){

		std::ostringstream ss;
		ss << fixed << setprecision(1) << PRB.RigBinCent(i);
		std::string s(ss.str());

		AcceptanceMatrix->GetXaxis()->SetBinLabel(i+1,s.c_str());
		AcceptanceMatrix->GetYaxis()->SetBinLabel(i+1,s.c_str());
	}
	AcceptanceMatrix->GetXaxis()->LabelsOption("v");

	AcceptanceMatrix->GetXaxis()->SetTitle("R_{gen} [GV]");
	AcceptanceMatrix->GetYaxis()->SetTitle("R_{meas} [GV]");

	AcceptanceMatrix->GetXaxis()->SetTitleOffset(1.35);
	AcceptanceMatrix->GetYaxis()->SetLabelSize(0.03);

	c4->cd();
	gPad->SetGridx();
	gPad->SetGridy();
	c4->SetLogz();
	AcceptanceMatrix->Draw("colz");
	finalResults.Add(c4);


	TCanvas * c3 = new TCanvas("Acceptance Projections");
        c3->cd();
	Acceptance_gen->SetLineColor(2);
	Acceptance_meas->SetLineColor(4);
	Acceptance_gen->Draw("P");
        Acceptance_meas->Draw("Psame");
	finalResults.Add(c3);
        finalResults.writeObjsInFolder("Acceptance Matrix");



}

TH1D * EvalRunningIntegral(TH1D * Distrib){
	TH1D * RunIntegral = (TH1D*) Distrib->Clone();
	RunIntegral->Clear();

	for(int i=0;i<RunIntegral->GetNbinsX();i++){
		RunIntegral->SetBinContent(i+1,Distrib->Integral(0,i+1));
//		RunIntegral->SetBinContent(i+1,Distrib->Integral(i+1,Distrib->GetNbinsX()));
		RunIntegral->SetBinError(i+1,0);
	}
	return RunIntegral;
}


void DrawVariables_app(FileSaver finalHistos,FileSaver finalResults, std::string app){

	TH2F * hd[7];
	TH2F * h[7];
	TCanvas * c[7];
	TCanvas * d[7];
	std::string names[7] = {"N Tof Clusters","Q L tof","Q U tof","tof Coo chi","tof Time chi","N Tracks","On Time"};

	hd[0]  =(TH2F*)finalHistos.Get( ("N Tof Clusters" +app+"/N Tof Clusters"+app).c_str());
	hd[1]  =(TH2F*)finalHistos.Get( ("Q L tof" +app+"/Q L tof"+app).c_str()	)	;
	hd[2]  =(TH2F*)finalHistos.Get( ("Q U tof"+app+"/Q U tof"+app).c_str()	)	 ;  
	hd[3]  =(TH2F*)finalHistos.Get( ("tof Coo chi"+app+"/tof Coo chi"+app).c_str());	   
	hd[4]  =(TH2F*)finalHistos.Get( ("tof Time chi"+app+"/tof Time chi"+app).c_str());		   
	hd[5]  =(TH2F*)finalHistos.Get( ("N Tracks"+app+"/N Tracks"+app).c_str()	)	 ;  
        hd[6]  =(TH2F*)finalHistos.Get( ("On Time"+app+"/On Time"+app).c_str()	)	;
                                         
	h[0]   =(TH2F*)finalHistos.Get( ("N Tof Clusters MC"+app+"/N Tof Clusters MC"+app).c_str()	);
	h[1]   =(TH2F*)finalHistos.Get( ("Q L tof MC"+app+"/Q L tof MC"+app).c_str()	)	;
	h[2]   =(TH2F*)finalHistos.Get( ("Q U tof MC"+app+"/Q U tof MC"+app).c_str()	)	 ;  
	h[3]   =(TH2F*)finalHistos.Get( ("tof Coo chi MC"+app+"/tof Coo chi MC"+app).c_str()	);	   
	h[4]   =(TH2F*)finalHistos.Get( ("tof Time chi MC"+app+"/tof Time chi MC"+app).c_str());		   
	h[5]   =(TH2F*)finalHistos.Get( ("N Tracks MC"+app+"/N Tracks MC"+app).c_str()	)	 ;  
        h[6]   =(TH2F*)finalHistos.Get( ("On Time MC"+app+"/On Time MC"+app).c_str()	)	;


	for(int i=0; i<7; i++) {
		
		cout<<h[i]<<" "<<hd[i]<<endl;	
		TH1D * HighDT = hd[i]->ProjectionX((names[i] + "DT_High").c_str(),15,40);
		TH1D * MedDT  = hd[i]->ProjectionX((names[i] + "DT_Med").c_str() ,8,15);
		TH1D * LowDT  = hd[i]->ProjectionX((names[i] + "DT_Low").c_str() ,1,7);

		TH1D * High = h[i]->ProjectionX((names[i] + "MC_High").c_str(),15,40);
		TH1D * Med  = h[i]->ProjectionX((names[i] + "MC_Med").c_str() ,8,15);
		TH1D * Low  = h[i]->ProjectionX((names[i] + "MC_Low").c_str() ,1,7);

		c[i] = new TCanvas(names[i].c_str());
		c[i]->Divide(3,1);
		d[i] = new TCanvas((names[i]+" R.I.").c_str());
		d[i]->Divide(3,1);
	
		c[i]->cd(1);
		gPad->SetLogy();


		LowDT->SetLineColor(1);
		LowDT->SetLineWidth(2);
		Low  ->SetLineColor(2);
		Low  ->SetLineWidth(2);
		LowDT->SetTitle("1-7 GV");	
		LowDT->Scale(1/LowDT->Integral());
		Low->Scale(1/Low->Integral());
	
		LowDT->Draw("hist");
		Low->Draw("same,hist");
		
		d[i]->cd(1);
		TH1D * LowDT_RI = EvalRunningIntegral(LowDT);
		TH1D * Low_RI = EvalRunningIntegral(Low);		
		LowDT_RI->Draw("hist");
		Low_RI->Draw("same,hist");

		c[i]->cd(2);
		gPad->SetLogy();

		MedDT->SetLineColor(1);
		MedDT->SetLineWidth(2);
		Med  ->SetLineColor(2);
		Med  ->SetLineWidth(2);
		MedDT->SetTitle("8-15 GV");	
		MedDT->Scale(1/MedDT->Integral());
		Med->Scale(1/Med->Integral());
	

		MedDT->Draw("hist");
		Med->Draw("same,hist");

		d[i]->cd(2);
		TH1D * MedDT_RI = EvalRunningIntegral(MedDT);
		TH1D * Med_RI = EvalRunningIntegral(Med);		
		MedDT_RI->Draw("hist");
		Med_RI->Draw("same,hist");

		c[i]->cd(3);
		gPad->SetLogy();


		HighDT->SetLineColor(1);
		HighDT->SetLineWidth(2);
		High  ->SetLineColor(2);
		High  ->SetLineWidth(2);
		HighDT->Scale(1/HighDT->Integral());
		High->Scale(1/High->Integral());
		HighDT->SetTitle("15-40 GV");	
		HighDT->Draw("hist");
		High->Draw("same,hist");

		d[i]->cd(3);
		TH1D * HighDT_RI = EvalRunningIntegral(HighDT);
		TH1D * High_RI = EvalRunningIntegral(High);		
		HighDT_RI->Draw("hist");
		High_RI->Draw("same,hist");


			finalResults.Add(c[i]);
			finalResults.Add(d[i]);

	}
	TCanvas * c2 = new TCanvas("NTracks vs QLtof");
	c2->cd();

       TH2F * NTrTracksvsQLtof = (TH2F*)finalHistos.Get(("NTrTracksvsQLtof"+app+"/NTrTracksvsQLtof"+app).c_str());
       TH2F * NTrTracksvsQLtof_MC = (TH2F*)finalHistos.Get(("NTrTracksvsQLtof MC"+app+"/NTrTracksvsQLtof MC"+app).c_str());
	cout<<NTrTracksvsQLtof<<" "<<NTrTracksvsQLtof_MC<<endl;

       TH1D * OneTrack      = NTrTracksvsQLtof->ProjectionY(   	"1 Track"	, NTrTracksvsQLtof->GetXaxis()->FindBin(0.5),NTrTracksvsQLtof->GetXaxis()->FindBin(1.5));	
       TH1D * OneTrack_MC   = NTrTracksvsQLtof_MC->ProjectionY(	"1 Track MC"	, NTrTracksvsQLtof->GetXaxis()->FindBin(0.5),NTrTracksvsQLtof->GetXaxis()->FindBin(1.5));
       TH1D * TwoTrack      = NTrTracksvsQLtof->ProjectionY(	"2 Track"		, NTrTracksvsQLtof->GetXaxis()->FindBin(1.5),NTrTracksvsQLtof->GetXaxis()->FindBin(2.5));	
       TH1D * TwoTrack_MC   = NTrTracksvsQLtof_MC->ProjectionY(	"2 Track MC"	, NTrTracksvsQLtof->GetXaxis()->FindBin(1.5),NTrTracksvsQLtof->GetXaxis()->FindBin(2.5));
       TH1D * ThreeTrack    = NTrTracksvsQLtof->ProjectionY(	"3 Track"		, NTrTracksvsQLtof->GetXaxis()->FindBin(2.5),NTrTracksvsQLtof->GetXaxis()->FindBin(3.5));	
       TH1D * ThreeTrack_MC = NTrTracksvsQLtof_MC->ProjectionY(	"3 Track MC"	, NTrTracksvsQLtof->GetXaxis()->FindBin(2.5),NTrTracksvsQLtof->GetXaxis()->FindBin(3.5));
       TH1D * FourTrack     = NTrTracksvsQLtof->ProjectionY(	"4 Track"		, NTrTracksvsQLtof->GetXaxis()->FindBin(3.5),NTrTracksvsQLtof->GetXaxis()->FindBin(4.5));	
       TH1D * FourTrack_MC  = NTrTracksvsQLtof_MC->ProjectionY(	"4 Track MC"	, NTrTracksvsQLtof->GetXaxis()->FindBin(3.5),NTrTracksvsQLtof->GetXaxis()->FindBin(4.5));
    
      OneTrack     ->SetTitle("1 Track"   ); 
      OneTrack_MC  ->SetTitle("1 Track MC"); 
      TwoTrack     ->SetTitle("2 Track"   ); 
      TwoTrack_MC  ->SetTitle("2 Track MC"); 
      ThreeTrack   ->SetTitle("3 Track"   ); 
      ThreeTrack_MC->SetTitle("3 Track MC"); 
      FourTrack    ->SetTitle("4 Track"   ); 
      FourTrack_MC ->SetTitle("4 Track MC"); 


       OneTrack      ->Scale(1/OneTrack     ->Integral() ); 
       OneTrack_MC   ->Scale(1/OneTrack_MC  ->Integral() ); 
       TwoTrack      ->Scale(1/TwoTrack    ->Integral()  ); 
       TwoTrack_MC   ->Scale(1/TwoTrack_MC ->Integral()  ); 
       ThreeTrack    ->Scale(1/ThreeTrack  ->Integral()  ); 
       ThreeTrack_MC ->Scale(1/ThreeTrack_MC->Integral() ); 
       FourTrack     ->Scale(1/FourTrack ->Integral()    ); 
       FourTrack_MC  ->Scale(1/FourTrack_MC->Integral() ); 

       OneTrack->SetLineColor(1);
       OneTrack_MC->SetLineColor(1);
       OneTrack->SetLineWidth(2);
       OneTrack_MC->SetLineWidth(2);
       OneTrack_MC->SetLineStyle(2);

      OneTrack->GetXaxis()->SetTitle("Q Ltof");
        OneTrack->Draw("hist");
       OneTrack_MC->Draw("hist,same");

       TwoTrack->SetLineColor(2);
       TwoTrack_MC->SetLineColor(2);
       TwoTrack->SetLineWidth(2);
       TwoTrack_MC->SetLineWidth(2);
       TwoTrack_MC->SetLineStyle(2);

       TwoTrack->Draw("hist,same");
       TwoTrack_MC->Draw("hist,same");

       ThreeTrack->SetLineColor(3);
       ThreeTrack_MC->SetLineColor(3);
       ThreeTrack->SetLineWidth(2);
       ThreeTrack_MC->SetLineWidth(2);
       ThreeTrack_MC->SetLineStyle(2);

       ThreeTrack->Draw("hist,same");
       ThreeTrack_MC->Draw("hist,same");

       FourTrack->SetLineColor(4);
       FourTrack_MC->SetLineColor(4);
       FourTrack->SetLineWidth(2);
       FourTrack_MC->SetLineWidth(2);
       FourTrack_MC->SetLineStyle(2);

       FourTrack->Draw("hist,same");
       FourTrack_MC->Draw("hist,same");


	finalResults.Add(c2);

	TCanvas * c3 = new TCanvas("NTracks vs QUtof");
	c3->cd();

       TH2F * NTrTracksvsQUtof = (TH2F*)finalHistos.Get(("NTrTracksvsQUtof"+app+"/NTrTracksvsQUtof"+app).c_str());
       TH2F * NTrTracksvsQUtof_MC = (TH2F*)finalHistos.Get(("NTrTracksvsQUtof MC"+app+"/NTrTracksvsQUtof MC"+app).c_str());
	cout<<NTrTracksvsQUtof<<" "<<NTrTracksvsQUtof_MC<<endl;

       TH1D * OneTrack2      = NTrTracksvsQUtof->ProjectionY("1 Track", NTrTracksvsQUtof->GetXaxis()->FindBin(0.5),NTrTracksvsQUtof->GetXaxis()->FindBin(1.5));	
       TH1D * OneTrack2_MC   = NTrTracksvsQUtof_MC->ProjectionY("1 Track MC", NTrTracksvsQUtof->GetXaxis()->FindBin(0.5),NTrTracksvsQUtof->GetXaxis()->FindBin(1.5));
       TH1D * TwoTrack2      = NTrTracksvsQUtof->ProjectionY("2 Track", NTrTracksvsQUtof->GetXaxis()->FindBin(1.5),NTrTracksvsQUtof->GetXaxis()->FindBin(2.5));	
       TH1D * TwoTrack2_MC   = NTrTracksvsQUtof_MC->ProjectionY("2 Track MC", NTrTracksvsQUtof->GetXaxis()->FindBin(1.5),NTrTracksvsQUtof->GetXaxis()->FindBin(2.5));
       TH1D * ThreeTrack2    = NTrTracksvsQUtof->ProjectionY("3 Track", NTrTracksvsQUtof->GetXaxis()->FindBin(2.5),NTrTracksvsQUtof->GetXaxis()->FindBin(3.5));	
       TH1D * ThreeTrack2_MC = NTrTracksvsQUtof_MC->ProjectionY("3 Track MC", NTrTracksvsQUtof->GetXaxis()->FindBin(2.5),NTrTracksvsQUtof->GetXaxis()->FindBin(3.5));
       TH1D * FourTrack2     = NTrTracksvsQUtof->ProjectionY("4 Track", NTrTracksvsQUtof->GetXaxis()->FindBin(3.5),NTrTracksvsQUtof->GetXaxis()->FindBin(4.5));	
       TH1D * FourTrack2_MC  = NTrTracksvsQUtof_MC->ProjectionY("4 Track MC", NTrTracksvsQUtof->GetXaxis()->FindBin(3.5),NTrTracksvsQUtof->GetXaxis()->FindBin(4.5));

      OneTrack2     ->SetTitle("1 Track"   ); 
      OneTrack2_MC  ->SetTitle("1 Track MC"); 
      TwoTrack2     ->SetTitle("2 Track"   ); 
      TwoTrack2_MC  ->SetTitle("2 Track MC"); 
      ThreeTrack2   ->SetTitle("3 Track"   ); 
      ThreeTrack2_MC->SetTitle("3 Track MC"); 
      FourTrack2    ->SetTitle("4 Track"   ); 
      FourTrack2_MC ->SetTitle("4 Track MC"); 


       OneTrack2      ->Scale(1/OneTrack2     ->Integral() ); 
       OneTrack2_MC   ->Scale(1/OneTrack2_MC  ->Integral() ); 
       TwoTrack2      ->Scale(1/TwoTrack2    ->Integral()  ); 
       TwoTrack2_MC   ->Scale(1/TwoTrack2_MC ->Integral()  ); 
       ThreeTrack2    ->Scale(1/ThreeTrack2  ->Integral()  ); 
       ThreeTrack2_MC ->Scale(1/ThreeTrack2_MC->Integral() ); 
       FourTrack2     ->Scale(1/FourTrack2 ->Integral()    ); 
       FourTrack2_MC  ->Scale(1/FourTrack2_MC->Integral() ); 

       OneTrack2->SetLineColor(1);
       OneTrack2_MC->SetLineColor(1);
       OneTrack2->SetLineWidth(2);
       OneTrack2_MC->SetLineWidth(2);
       OneTrack2_MC->SetLineStyle(2);

	OneTrack2->GetXaxis()->SetTitle("Q Utof");
       OneTrack2->Draw("hist");
       OneTrack2_MC->Draw("hist,same");

       TwoTrack2->SetLineColor(2);
       TwoTrack2_MC->SetLineColor(2);
       TwoTrack2->SetLineWidth(2);
       TwoTrack2_MC->SetLineWidth(2);
       TwoTrack2_MC->SetLineStyle(2);

       TwoTrack2->Draw("hist,same");
       TwoTrack2_MC->Draw("hist,same");

       ThreeTrack2->SetLineColor(3);
       ThreeTrack2_MC->SetLineColor(3);
       ThreeTrack2->SetLineWidth(2);
       ThreeTrack2_MC->SetLineWidth(2);
       ThreeTrack2_MC->SetLineStyle(2);

       ThreeTrack2->Draw("hist,same");
       ThreeTrack2_MC->Draw("hist,same");

       FourTrack2->SetLineColor(4);
       FourTrack2_MC->SetLineColor(4);
       FourTrack2->SetLineWidth(2);
       FourTrack2_MC->SetLineWidth(2);
       FourTrack2_MC->SetLineStyle(2);

       FourTrack2->Draw("hist,same");
       FourTrack2_MC->Draw("hist,same");


finalResults.Add(c3);
	TCanvas * c4 = new TCanvas("NTracks vs OnTime");
	c4->cd();

       TH2F * OnTimeVsNTrTracks = (TH2F*)finalHistos.Get(("OnTimeVsNTrTracks"+app+"/OnTimeVsNTrTracks"+app).c_str());
       TH2F * OnTimeVsNTrTracks_MC = (TH2F*)finalHistos.Get(("OnTimeVsNTrTracks MC"+app+"/OnTimeVsNTrTracks MC"+app).c_str());
	cout<<OnTimeVsNTrTracks<<" "<<OnTimeVsNTrTracks_MC<<endl;

       TH1D * OneTrack3      = OnTimeVsNTrTracks->ProjectionX("1 Track", OnTimeVsNTrTracks->GetYaxis()->FindBin(0.5),OnTimeVsNTrTracks->GetYaxis()->FindBin(1.5));	
       TH1D * OneTrack3_MC   = OnTimeVsNTrTracks_MC->ProjectionX("1 Track MC", OnTimeVsNTrTracks->GetYaxis()->FindBin(0.5),OnTimeVsNTrTracks->GetYaxis()->FindBin(1.5));
       TH1D * TwoTrack3      = OnTimeVsNTrTracks->ProjectionX("2 Track", OnTimeVsNTrTracks->GetYaxis()->FindBin(1.5),OnTimeVsNTrTracks->GetYaxis()->FindBin(2.5));	
       TH1D * TwoTrack3_MC   = OnTimeVsNTrTracks_MC->ProjectionX("2 Track MC", OnTimeVsNTrTracks->GetYaxis()->FindBin(1.5),OnTimeVsNTrTracks->GetYaxis()->FindBin(2.5));
       TH1D * ThreeTrack3    = OnTimeVsNTrTracks->ProjectionX("3 Track", OnTimeVsNTrTracks->GetYaxis()->FindBin(2.5),OnTimeVsNTrTracks->GetYaxis()->FindBin(3.5));	
       TH1D * ThreeTrack3_MC = OnTimeVsNTrTracks_MC->ProjectionX("3 Track MC", OnTimeVsNTrTracks->GetYaxis()->FindBin(2.5),OnTimeVsNTrTracks->GetYaxis()->FindBin(3.5));
       TH1D * FourTrack3     = OnTimeVsNTrTracks->ProjectionX("4 Track", OnTimeVsNTrTracks->GetYaxis()->FindBin(3.5),OnTimeVsNTrTracks->GetYaxis()->FindBin(4.5));	
       TH1D * FourTrack3_MC  = OnTimeVsNTrTracks_MC->ProjectionX("4 Track MC", OnTimeVsNTrTracks->GetYaxis()->FindBin(3.5),OnTimeVsNTrTracks->GetYaxis()->FindBin(4.5));

      OneTrack3     ->SetTitle("1 Track"   ); 
      OneTrack3_MC  ->SetTitle("1 Track MC"); 
      TwoTrack3     ->SetTitle("2 Track"   ); 
      TwoTrack3_MC  ->SetTitle("2 Track MC"); 
      ThreeTrack3   ->SetTitle("3 Track"   ); 
      ThreeTrack3_MC->SetTitle("3 Track MC"); 
      FourTrack3    ->SetTitle("4 Track"   ); 
      FourTrack3_MC ->SetTitle("4 Track MC"); 



	TH1D * OneTrack3_red     = ReduceOnTimeHisto(OneTrack3     	);  
        TH1D * OneTrack3_MC_red  = ReduceOnTimeHisto(OneTrack3_MC  	);
        TH1D * TwoTrack3_red     = ReduceOnTimeHisto(TwoTrack3     	);
        TH1D * TwoTrack3_MC_red  = ReduceOnTimeHisto(TwoTrack3_MC  	);
        TH1D * ThreeTrack3_red   = ReduceOnTimeHisto(ThreeTrack3   	);
        TH1D * ThreeTrack3_MC_red= ReduceOnTimeHisto(ThreeTrack3_MC	);
        TH1D * FourTrack3_red    = ReduceOnTimeHisto(FourTrack3    	);
        TH1D * FourTrack3_MC_red = ReduceOnTimeHisto(FourTrack3_MC 	);


       OneTrack3_red      ->Scale(1/OneTrack3_red     ->Integral() ); 
       OneTrack3_MC_red   ->Scale(1/OneTrack3_MC_red  ->Integral() ); 
       TwoTrack3_red      ->Scale(1/TwoTrack3_red     ->Integral()  ); 
       TwoTrack3_MC_red   ->Scale(1/TwoTrack3_MC_red  ->Integral()  ); 
       ThreeTrack3_red    ->Scale(1/ThreeTrack3_red   ->Integral()  ); 
       ThreeTrack3_MC_red ->Scale(1/ThreeTrack3_MC_red->Integral() ); 
       FourTrack3_red     ->Scale(1/FourTrack3_red    ->Integral()    ); 
       FourTrack3_MC_red  ->Scale(1/FourTrack3_MC_red ->Integral() ); 

       OneTrack3_red->SetLineColor(1);
       OneTrack3_MC_red->SetLineColor(1);
       OneTrack3_red->SetLineWidth(2);
       OneTrack3_MC_red->SetLineWidth(2);
       OneTrack3_MC_red->SetLineStyle(2);

	OneTrack3_red->GetXaxis()->SetTitle("On Time Cl.");
       OneTrack3_red->Draw("hist");
       OneTrack3_MC_red->Draw("hist,same");

       TwoTrack3_red->SetLineColor(2);
       TwoTrack3_MC_red->SetLineColor(2);
       TwoTrack3_red->SetLineWidth(2);
       TwoTrack3_MC_red->SetLineWidth(2);
       TwoTrack3_MC_red->SetLineStyle(2);

       TwoTrack3_red->Draw("hist,same");
       TwoTrack3_MC_red->Draw("hist,same");

       ThreeTrack3_red->SetLineColor(3);
       ThreeTrack3_MC_red->SetLineColor(3);
       ThreeTrack3_red->SetLineWidth(2);
       ThreeTrack3_MC_red->SetLineWidth(2);
       ThreeTrack3_MC_red->SetLineStyle(2);

       ThreeTrack3_red->Draw("hist,same");
       ThreeTrack3_MC_red->Draw("hist,same");

       FourTrack3_red->SetLineColor(4);
       FourTrack3_MC_red->SetLineColor(4);
       FourTrack3_red->SetLineWidth(2);
       FourTrack3_MC_red->SetLineWidth(2);
       FourTrack3_MC_red->SetLineStyle(2);

       FourTrack3_red->Draw("hist,same");
       FourTrack3_MC_red->Draw("hist,same");


	finalResults.Add(c4);


	finalResults.Add(c4);

	TCanvas * c5 = new TCanvas("NTofCl vs OnTime");
	c5->cd();

       TH2F * TofClustersVsOnTime = (TH2F*)finalHistos.Get(("TofClustersVsOnTime"+app+"/TofClustersVsOnTime"+app).c_str());
       TH2F * TofClustersVsOnTime_MC = (TH2F*)finalHistos.Get(("TofClustersVsOnTime MC"+app+"/TofClustersVsOnTime MC"+app).c_str());
	cout<<TofClustersVsOnTime<<" "<<TofClustersVsOnTime_MC<<endl;

	TH1D * OneClust      = TofClustersVsOnTime->ProjectionY(	"1 Cluster"	, TofClustersVsOnTime->GetXaxis()->FindBin(0.5),TofClustersVsOnTime->GetXaxis()->FindBin(1.5));	
	TH1D * OneClust_MC   = TofClustersVsOnTime_MC->ProjectionY(	"1 Cluster MC"	, TofClustersVsOnTime->GetXaxis()->FindBin(0.5),TofClustersVsOnTime->GetXaxis()->FindBin(1.5));
	TH1D * TwoClust      = TofClustersVsOnTime->ProjectionY(	"2 Cluster"	, TofClustersVsOnTime->GetXaxis()->FindBin(1.5),TofClustersVsOnTime->GetXaxis()->FindBin(2.5));	
	TH1D * TwoClust_MC   = TofClustersVsOnTime_MC->ProjectionY(	"2 Cluster MC"	, TofClustersVsOnTime->GetXaxis()->FindBin(1.5),TofClustersVsOnTime->GetXaxis()->FindBin(2.5));
	TH1D * ThreeClust    = TofClustersVsOnTime->ProjectionY(	"3 Cluster"	, TofClustersVsOnTime->GetXaxis()->FindBin(2.5),TofClustersVsOnTime->GetXaxis()->FindBin(3.5));	
	TH1D * ThreeClust_MC = TofClustersVsOnTime_MC->ProjectionY(	"3 Cluster MC"	, TofClustersVsOnTime->GetXaxis()->FindBin(2.5),TofClustersVsOnTime->GetXaxis()->FindBin(3.5));
	TH1D * FourClust     = TofClustersVsOnTime->ProjectionY(	"4 Cluster"	, TofClustersVsOnTime->GetXaxis()->FindBin(3.5),TofClustersVsOnTime->GetXaxis()->FindBin(4.5));	
	TH1D * FourClust_MC  = TofClustersVsOnTime_MC->ProjectionY(	"4 Cluster MC"	, TofClustersVsOnTime->GetXaxis()->FindBin(3.5),TofClustersVsOnTime->GetXaxis()->FindBin(4.5));

	OneClust     ->SetTitle("1 Cluster"   ); 
	OneClust_MC  ->SetTitle("1 Cluster MC"); 
	TwoClust     ->SetTitle("2 Cluster"   ); 
	TwoClust_MC  ->SetTitle("2 Cluster MC"); 
	ThreeClust   ->SetTitle("3 Cluster"   ); 
	ThreeClust_MC->SetTitle("3 Cluster MC"); 
	FourClust    ->SetTitle("4 Cluster"   ); 
	FourClust_MC ->SetTitle("4 Cluster MC"); 



	TH1D * OneClust_red     = ReduceOnTimeHisto(OneClust     	);  
	TH1D * OneClust_MC_red  = ReduceOnTimeHisto(OneClust_MC  	);
	TH1D * TwoClust_red     = ReduceOnTimeHisto(TwoClust     	);
	TH1D * TwoClust_MC_red  = ReduceOnTimeHisto(TwoClust_MC  	);
	TH1D * ThreeClust_red   = ReduceOnTimeHisto(ThreeClust   	);
	TH1D * ThreeClust_MC_red= ReduceOnTimeHisto(ThreeClust_MC	);
	TH1D * FourClust_red    = ReduceOnTimeHisto(FourClust    	);
	TH1D * FourClust_MC_red = ReduceOnTimeHisto(FourClust_MC 	);

     	if(OneClust_red          ->Integral() > 0)  OneClust_red      ->Scale(1/OneClust_red     ->Integral() ); 
	if(OneClust_MC_red       ->Integral() > 0)  OneClust_MC_red   ->Scale(1/OneClust_MC_red  ->Integral() ); 
	if(TwoClust_red          ->Integral() > 0)  TwoClust_red      ->Scale(1/TwoClust_red     ->Integral()  ); 
	if(TwoClust_MC_red       ->Integral() > 0)  TwoClust_MC_red   ->Scale(1/TwoClust_MC_red  ->Integral()  ); 
	if(ThreeClust_red        ->Integral() > 0)  ThreeClust_red    ->Scale(1/ThreeClust_red   ->Integral()  ); 
	if(ThreeClust_MC_red     ->Integral() > 0)  ThreeClust_MC_red ->Scale(1/ThreeClust_MC_red->Integral() ); 
	if(FourClust_red         ->Integral() > 0)  FourClust_red     ->Scale(1/FourClust_red    ->Integral()    ); 
        if(FourClust_MC_red      ->Integral() > 0)  FourClust_MC_red  ->Scale(1/FourClust_MC_red ->Integral() ); 

       OneClust_red->SetLineColor(1);
       OneClust_MC_red->SetLineColor(1);
       OneClust_red->SetLineWidth(2);
       OneClust_MC_red->SetLineWidth(2);
       OneClust_MC_red->SetLineStyle(2);

       OneClust_red->GetXaxis()->SetTitle("On Time Cl.");
       OneClust_red->Draw("hist");
       OneClust_MC_red->Draw("hist,same");

       TwoClust_red->SetLineColor(2);
       TwoClust_MC_red->SetLineColor(2);
       TwoClust_red->SetLineWidth(2);
       TwoClust_MC_red->SetLineWidth(2);
       TwoClust_MC_red->SetLineStyle(2);

       TwoClust_red->Draw("hist,same");
       TwoClust_MC_red->Draw("hist,same");

       ThreeClust_red->SetLineColor(3);
       ThreeClust_MC_red->SetLineColor(3);
       ThreeClust_red->SetLineWidth(2);
       ThreeClust_MC_red->SetLineWidth(2);
       ThreeClust_MC_red->SetLineStyle(2);

       ThreeClust_red->Draw("hist,same");
       ThreeClust_MC_red->Draw("hist,same");

       FourClust_red->SetLineColor(4);
       FourClust_MC_red->SetLineColor(4);
       FourClust_red->SetLineWidth(2);
       FourClust_MC_red->SetLineWidth(2);
       FourClust_MC_red->SetLineStyle(2);

       FourClust_red->Draw("hist,same");
       FourClust_MC_red->Draw("hist,same");


	finalResults.Add(c5);


	TCanvas * c6 = new TCanvas("NTofCl vs QUtof");
	c6->cd();

       TH2F * TofClustersVsQUtof = (TH2F*)finalHistos.Get(("TofClustersVsQUtof"+app+"/TofClustersVsQUtof"+app).c_str());
       TH2F * TofClustersVsQUtof_MC = (TH2F*)finalHistos.Get(("TofClustersVsQUtof MC"+app+"/TofClustersVsQUtof MC"+app).c_str());
	cout<<TofClustersVsQUtof<<" "<<TofClustersVsQUtof_MC<<endl;

       TH1D * OneClust2      = TofClustersVsQUtof->ProjectionY("1 Cluster", TofClustersVsQUtof->GetXaxis()->FindBin(0.5),TofClustersVsQUtof->GetXaxis()->FindBin(1.5));	
       TH1D * OneClust2_MC   = TofClustersVsQUtof_MC->ProjectionY("1 Cluster MC", TofClustersVsQUtof->GetXaxis()->FindBin(0.5),TofClustersVsQUtof->GetXaxis()->FindBin(1.5));
       TH1D * TwoClust2      = TofClustersVsQUtof->ProjectionY("2 Cluster", TofClustersVsQUtof->GetXaxis()->FindBin(1.5),TofClustersVsQUtof->GetXaxis()->FindBin(2.5));	
       TH1D * TwoClust2_MC   = TofClustersVsQUtof_MC->ProjectionY("2 Cluster MC", TofClustersVsQUtof->GetXaxis()->FindBin(1.5),TofClustersVsQUtof->GetXaxis()->FindBin(2.5));
       TH1D * ThreeClust2    = TofClustersVsQUtof->ProjectionY("3 Cluster", TofClustersVsQUtof->GetXaxis()->FindBin(2.5),TofClustersVsQUtof->GetXaxis()->FindBin(3.5));	
       TH1D * ThreeClust2_MC = TofClustersVsQUtof_MC->ProjectionY("3 Cluster MC", TofClustersVsQUtof->GetXaxis()->FindBin(2.5),TofClustersVsQUtof->GetXaxis()->FindBin(3.5));
       TH1D * FourClust2     = TofClustersVsQUtof->ProjectionY("4 Cluster", TofClustersVsQUtof->GetXaxis()->FindBin(3.5),TofClustersVsQUtof->GetXaxis()->FindBin(4.5));	
       TH1D * FourClust2_MC  = TofClustersVsQUtof_MC->ProjectionY("4 Cluster MC", TofClustersVsQUtof->GetXaxis()->FindBin(3.5),TofClustersVsQUtof->GetXaxis()->FindBin(4.5));

	OneClust2     ->SetTitle("1 Cluster"   ); 
	OneClust2_MC  ->SetTitle("1 Cluster MC"); 
	TwoClust2     ->SetTitle("2 Cluster"   ); 
	TwoClust2_MC  ->SetTitle("2 Cluster MC"); 
	ThreeClust2   ->SetTitle("3 Cluster"   ); 
	ThreeClust2_MC->SetTitle("3 Cluster MC"); 
	FourClust2    ->SetTitle("4 Cluster"   ); 
	FourClust2_MC ->SetTitle("4 Cluster MC"); 




       if(OneClust2          ->Integral() > 0)OneClust2      ->Scale(1/OneClust2     ->Integral() ); 
       if(OneClust2_MC       ->Integral() > 0)OneClust2_MC   ->Scale(1/OneClust2_MC  ->Integral() ); 
       if(TwoClust2          ->Integral() > 0)TwoClust2      ->Scale(1/TwoClust2    ->Integral()  ); 
       if(TwoClust2_MC       ->Integral() > 0)TwoClust2_MC   ->Scale(1/TwoClust2_MC ->Integral()  ); 
       if(ThreeClust2        ->Integral() > 0)ThreeClust2    ->Scale(1/ThreeClust2  ->Integral()  ); 
       if(ThreeClust2_MC     ->Integral() > 0)ThreeClust2_MC ->Scale(1/ThreeClust2_MC->Integral() ); 
       if(FourClust2         ->Integral() > 0)FourClust2     ->Scale(1/FourClust2 ->Integral()    ); 
       if(FourClust2_MC      ->Integral() > 0)FourClust2_MC  ->Scale(1/FourClust2_MC->Integral() ); 

       FourClust2->SetLineColor(4);
       FourClust2_MC->SetLineColor(4);
       FourClust2->SetLineWidth(2);
       FourClust2_MC->SetLineWidth(2);
       FourClust2_MC->SetLineStyle(2);

       FourClust2->GetXaxis()->SetTitle("Q UTof");
       FourClust2->Draw("hist");
       FourClust2_MC->Draw("hist,same");

       OneClust2->SetLineColor(1);
       OneClust2_MC->SetLineColor(1);
       OneClust2->SetLineWidth(2);
       OneClust2_MC->SetLineWidth(2);
       OneClust2_MC->SetLineStyle(2);

       OneClust2->Draw("hist,same");
       OneClust2_MC->Draw("hist,same");

       TwoClust2->SetLineColor(2);
       TwoClust2_MC->SetLineColor(2);
       TwoClust2->SetLineWidth(2);
       TwoClust2_MC->SetLineWidth(2);
       TwoClust2_MC->SetLineStyle(2);

       TwoClust2->Draw("hist,same");
       TwoClust2_MC->Draw("hist,same");

       ThreeClust2->SetLineColor(3);
       ThreeClust2_MC->SetLineColor(3);
       ThreeClust2->SetLineWidth(2);
       ThreeClust2_MC->SetLineWidth(2);
       ThreeClust2_MC->SetLineStyle(2);

       ThreeClust2->Draw("hist,same");
       ThreeClust2_MC->Draw("hist,same");

     	finalResults.Add(c6);

	TCanvas * c7 = new TCanvas("NTofCl vs QLtof");
	c7->cd();

       TH2F * TofClustersVsQLtof = (TH2F*)finalHistos.Get(("TofClustersVsQLtof"+app+"/TofClustersVsQLtof"+app).c_str());
       TH2F * TofClustersVsQLtof_MC = (TH2F*)finalHistos.Get(("TofClustersVsQLtof MC"+app+"/TofClustersVsQLtof MC"+app).c_str());
	cout<<TofClustersVsQLtof<<" "<<TofClustersVsQLtof_MC<<endl;

       TH1D * OneClust3      = TofClustersVsQLtof->ProjectionY("1 Cluster", TofClustersVsQLtof->GetXaxis()->FindBin(0.5),TofClustersVsQLtof->GetXaxis()->FindBin(1.5));	
       TH1D * OneClust3_MC   = TofClustersVsQLtof_MC->ProjectionY("1 Cluster MC", TofClustersVsQLtof->GetXaxis()->FindBin(0.5),TofClustersVsQLtof->GetXaxis()->FindBin(1.5));
       TH1D * TwoClust3      = TofClustersVsQLtof->ProjectionY("2 Cluster", TofClustersVsQLtof->GetXaxis()->FindBin(1.5),TofClustersVsQLtof->GetXaxis()->FindBin(2.5));	
       TH1D * TwoClust3_MC   = TofClustersVsQLtof_MC->ProjectionY("2 Cluster MC", TofClustersVsQLtof->GetXaxis()->FindBin(1.5),TofClustersVsQLtof->GetXaxis()->FindBin(2.5));
       TH1D * ThreeClust3    = TofClustersVsQLtof->ProjectionY("3 Cluster", TofClustersVsQLtof->GetXaxis()->FindBin(2.5),TofClustersVsQLtof->GetXaxis()->FindBin(3.5));	
       TH1D * ThreeClust3_MC = TofClustersVsQLtof_MC->ProjectionY("3 Cluster MC", TofClustersVsQLtof->GetXaxis()->FindBin(2.5),TofClustersVsQLtof->GetXaxis()->FindBin(3.5));
       TH1D * FourClust3     = TofClustersVsQLtof->ProjectionY("4 Cluster", TofClustersVsQLtof->GetXaxis()->FindBin(3.5),TofClustersVsQLtof->GetXaxis()->FindBin(4.5));	
       TH1D * FourClust3_MC  = TofClustersVsQLtof_MC->ProjectionY("4 Cluster MC", TofClustersVsQLtof->GetXaxis()->FindBin(3.5),TofClustersVsQLtof->GetXaxis()->FindBin(4.5));

	OneClust3     ->SetTitle("1 Cluster"   ); 
	OneClust3_MC  ->SetTitle("1 Cluster MC"); 
	TwoClust3     ->SetTitle("2 Cluster"   ); 
	TwoClust3_MC  ->SetTitle("2 Cluster MC"); 
	ThreeClust3   ->SetTitle("3 Cluster"   ); 
	ThreeClust3_MC->SetTitle("3 Cluster MC"); 
	FourClust3    ->SetTitle("4 Cluster"   ); 
	FourClust3_MC ->SetTitle("4 Cluster MC"); 


       if(OneClust3          ->Integral() > 0)OneClust3      ->Scale(1/OneClust3     ->Integral() ); 
       if(OneClust3_MC       ->Integral() > 0)OneClust3_MC   ->Scale(1/OneClust3_MC  ->Integral() ); 
       if(TwoClust3          ->Integral() > 0)TwoClust3      ->Scale(1/TwoClust3    ->Integral()  ); 
       if(TwoClust3_MC       ->Integral() > 0)TwoClust3_MC   ->Scale(1/TwoClust3_MC ->Integral()  ); 
       if(ThreeClust3        ->Integral() > 0)ThreeClust3    ->Scale(1/ThreeClust3  ->Integral()  ); 
       if(ThreeClust3_MC     ->Integral() > 0)ThreeClust3_MC ->Scale(1/ThreeClust3_MC->Integral() ); 
       if(FourClust3         ->Integral() > 0)FourClust3     ->Scale(1/FourClust3 ->Integral()    ); 
       if(FourClust3_MC      ->Integral() > 0)FourClust3_MC  ->Scale(1/FourClust3_MC->Integral() ); 


       FourClust3->SetLineColor(4);
       FourClust3_MC->SetLineColor(4);
       FourClust3->SetLineWidth(2);
       FourClust3_MC->SetLineWidth(2);
       FourClust3_MC->SetLineStyle(2);

       FourClust3->GetXaxis()->SetTitle("Q LTof");
       FourClust3->Draw("hist");
       FourClust3_MC->Draw("hist,same");


       OneClust3->SetLineColor(1);
       OneClust3_MC->SetLineColor(1);
       OneClust3->SetLineWidth(2);
       OneClust3_MC->SetLineWidth(2);
       OneClust3_MC->SetLineStyle(2);

       OneClust3->Draw("hist,same");
       OneClust3_MC->Draw("hist,same");

       TwoClust3->SetLineColor(2);
       TwoClust3_MC->SetLineColor(2);
       TwoClust3->SetLineWidth(2);
       TwoClust3_MC->SetLineWidth(2);
       TwoClust3_MC->SetLineStyle(2);

       TwoClust3->Draw("hist,same");
       TwoClust3_MC->Draw("hist,same");

       ThreeClust3->SetLineColor(3);
       ThreeClust3_MC->SetLineColor(3);
       ThreeClust3->SetLineWidth(2);
       ThreeClust3_MC->SetLineWidth(2);
       ThreeClust3_MC->SetLineStyle(2);

       ThreeClust3->Draw("hist,same");
       ThreeClust3_MC->Draw("hist,same");

	finalResults.Add(c7);

	TCanvas * c8 = new TCanvas("NTofCl vs NTrTracks");
	c8->cd();
	
	TH2F * TofClustersVsNTrtracks = (TH2F*)finalHistos.Get(("TofClustersVsNTrtracks"+app+"/TofClustersVsNTrtracks"+app).c_str());
        TH2F * TofClustersVsNTrtracks_MC = (TH2F*)finalHistos.Get(("TofClustersVsNTrtracks MC"+app+"/TofClustersVsNTrtracks MC"+app).c_str());
	cout<<TofClustersVsNTrtracks<<" "<<TofClustersVsNTrtracks_MC<<endl;

       TH1D * OneClust4      = TofClustersVsNTrtracks->ProjectionY("1 Cluster", TofClustersVsNTrtracks->GetXaxis()->FindBin(0.5),TofClustersVsNTrtracks->GetXaxis()->FindBin(1.5));	
       TH1D * OneClust4_MC   = TofClustersVsNTrtracks_MC->ProjectionY("1 Cluster MC", TofClustersVsNTrtracks->GetXaxis()->FindBin(0.5),TofClustersVsNTrtracks->GetXaxis()->FindBin(1.5));
       TH1D * TwoClust4      = TofClustersVsNTrtracks->ProjectionY("2 Cluster", TofClustersVsNTrtracks->GetXaxis()->FindBin(1.5),TofClustersVsNTrtracks->GetXaxis()->FindBin(2.5));	
       TH1D * TwoClust4_MC   = TofClustersVsNTrtracks_MC->ProjectionY("2 Cluster MC", TofClustersVsNTrtracks->GetXaxis()->FindBin(1.5),TofClustersVsNTrtracks->GetXaxis()->FindBin(2.5));
       TH1D * ThreeClust4    = TofClustersVsNTrtracks->ProjectionY("3 Cluster", TofClustersVsNTrtracks->GetXaxis()->FindBin(2.5),TofClustersVsNTrtracks->GetXaxis()->FindBin(3.5));	
       TH1D * ThreeClust4_MC = TofClustersVsNTrtracks_MC->ProjectionY("3 Cluster MC", TofClustersVsNTrtracks->GetXaxis()->FindBin(2.5),TofClustersVsNTrtracks->GetXaxis()->FindBin(3.5));
       TH1D * FourClust4     = TofClustersVsNTrtracks->ProjectionY("4 Cluster", TofClustersVsNTrtracks->GetXaxis()->FindBin(3.5),TofClustersVsNTrtracks->GetXaxis()->FindBin(4.5));	
       TH1D * FourClust4_MC  = TofClustersVsNTrtracks_MC->ProjectionY("4 Cluster MC", TofClustersVsNTrtracks->GetXaxis()->FindBin(3.5),TofClustersVsNTrtracks->GetXaxis()->FindBin(4.5));

       if(OneClust4          ->Integral() > 0)OneClust4      ->Scale(1/OneClust4     ->Integral() ); 
       if(OneClust4_MC       ->Integral() > 0)OneClust4_MC   ->Scale(1/OneClust4_MC  ->Integral() ); 
       if(TwoClust4          ->Integral() > 0)TwoClust4      ->Scale(1/TwoClust4    ->Integral()  ); 
       if(TwoClust4_MC       ->Integral() > 0)TwoClust4_MC   ->Scale(1/TwoClust4_MC ->Integral()  ); 
       if(ThreeClust4        ->Integral() > 0)ThreeClust4    ->Scale(1/ThreeClust4  ->Integral()  ); 
       if(ThreeClust4_MC     ->Integral() > 0)ThreeClust4_MC ->Scale(1/ThreeClust4_MC->Integral() ); 
       if(FourClust4         ->Integral() > 0)FourClust4     ->Scale(1/FourClust4 ->Integral()    ); 
       if(FourClust4_MC      ->Integral() > 0)FourClust4_MC  ->Scale(1/FourClust4_MC->Integral() ); 
	
	OneClust4     ->SetTitle("1 Cluster"   ); 
	OneClust4_MC  ->SetTitle("1 Cluster MC"); 
	TwoClust4     ->SetTitle("2 Cluster"   ); 
	TwoClust4_MC  ->SetTitle("2 Cluster MC"); 
	ThreeClust4   ->SetTitle("3 Cluster"   ); 
	ThreeClust4_MC->SetTitle("3 Cluster MC"); 
	FourClust4    ->SetTitle("4 Cluster"   ); 
	FourClust4_MC ->SetTitle("4 Cluster MC"); 


       OneClust4    ->Rebin(2); 
       OneClust4_MC  ->Rebin(2); 
       TwoClust4     ->Rebin(2); 
       TwoClust4_MC  ->Rebin(2); 
       ThreeClust4   ->Rebin(2); 
       ThreeClust4_MC->Rebin(2); 
       FourClust4    ->Rebin(2); 
       FourClust4_MC ->Rebin(2); 


       FourClust4->SetLineColor(4);
       FourClust4_MC->SetLineColor(4);
       FourClust4->SetLineWidth(2);
       FourClust4_MC->SetLineWidth(2);
       FourClust4_MC->SetLineStyle(2);
       
       FourClust4->GetXaxis()->SetTitle("N Tr. Tracks");
       FourClust4->Draw("hist");
       FourClust4_MC->Draw("hist,same");


       OneClust4->SetLineColor(1);
       OneClust4_MC->SetLineColor(1);
       OneClust4->SetLineWidth(2);
       OneClust4_MC->SetLineWidth(2);
       OneClust4_MC->SetLineStyle(2);

       OneClust4->Draw("hist,same");
       OneClust4_MC->Draw("hist,same");

       TwoClust4->SetLineColor(2);
       TwoClust4_MC->SetLineColor(2);
       TwoClust4->SetLineWidth(2);
       TwoClust4_MC->SetLineWidth(2);
       TwoClust4_MC->SetLineStyle(2);

       TwoClust4->Draw("hist,same");
       TwoClust4_MC->Draw("hist,same");

       ThreeClust4->SetLineColor(3);
       ThreeClust4_MC->SetLineColor(3);
       ThreeClust4->SetLineWidth(2);
       ThreeClust4_MC->SetLineWidth(2);
       ThreeClust4_MC->SetLineStyle(2);

       ThreeClust4->Draw("hist,same");
       ThreeClust4_MC->Draw("hist,same");

       finalResults.Add(c8);



       finalResults.writeObjsInFolder(("Variables"+app).c_str());





}
 
void DrawVariables(FileSaver finalHistos,FileSaver finalResults){
	DrawVariables_app(finalHistos,finalResults,"");
	DrawVariables_app(finalHistos,finalResults," L1");

}

