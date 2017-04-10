using namespace std;

TH2F* CreateFrame (float xmin,float xmax,float ymin, float ymax,std::string Xaxis,std::string Yaxis){
	
		TH2F * Frame = new TH2F("Frame","Frame",1000,xmin,xmax,1000,ymin,ymax);
		Frame->SetStats(false);
		Frame->SetTitle("");

		Frame->SetStats(false);
		Frame->SetTitle("");

		Frame->GetYaxis()->SetMoreLogLabels();
		Frame->GetXaxis()->SetMoreLogLabels();

		Frame->GetYaxis()->SetNoExponent();
		Frame->GetXaxis()->SetNoExponent();

		Frame->GetYaxis()->SetRangeUser(ymin,ymax);
		Frame->GetYaxis()->SetRangeUser(ymin,ymax);

		Frame->GetXaxis()->SetTitle(Xaxis.c_str());
		Frame->GetYaxis()->SetTitle(Yaxis.c_str());	

		Frame->GetXaxis()->SetTitleFont(33);
                Frame->GetYaxis()->SetTitleFont(33); 

                Frame->GetXaxis()->SetLabelFont(33);
                Frame->GetYaxis()->SetLabelFont(33); 

		Frame->GetXaxis()->SetLabelSize(40);
                Frame->GetYaxis()->SetLabelSize(40); 

		Frame->GetXaxis()->SetTitleSize(40);
		Frame->GetYaxis()->SetTitleSize(40);	

		Frame->GetXaxis()->CenterTitle();
		Frame->GetYaxis()->CenterTitle();

		return Frame;	
	
}


void PlotDistribution(TVirtualPad * c, TH1F * Distribution, std::string Xaxis, std::string Yaxis, int color, std::string options, float ymin=-1,float ymax=-1,int rebin=1){
	c -> cd();
	gPad-> SetLogy();
	gPad-> SetTickx();
	gPad-> SetTicky();

	Distribution->SetLineColor(color);
	Distribution->SetLineWidth(3);

	Distribution->SetStats(false);
	Distribution->SetTitle("");


	TF1 * fit = Distribution->GetFunction("fitfunc");

	if(ymin>0&&ymax>0) Distribution->GetYaxis()->SetRangeUser(ymin,ymax);

	Distribution->Rebin(rebin);

	Distribution->GetXaxis()->SetTitle(Xaxis.c_str());
	Distribution->GetYaxis()->SetTitle(Yaxis.c_str());	

	Distribution->GetXaxis()->SetTitleSize(0.045);
	Distribution->GetYaxis()->SetTitleSize(0.045);	

	Distribution->GetXaxis()->CenterTitle();
	Distribution->GetYaxis()->CenterTitle();

	Distribution->GetXaxis()->SetTitleFont(32);
	Distribution->GetYaxis()->SetTitleFont(32); 

	Distribution->GetXaxis()->SetLabelFont(32);
	Distribution->GetYaxis()->SetLabelFont(32); 


	Distribution->Draw((options+",hist").c_str());

	if(fit){
		fit->SetLineColor(color);
		fit->SetLineWidth(5);
		fit->Draw("same");

	}

	return;
}

void PlotFunction(TVirtualPad * c, TF1 * Function, std::string Xaxis, std::string Yaxis, int color, std::string options, float xmin=-1,float xmax=-1,float ymin=-1,float ymax=-1,std::string legendname=""){
	c -> cd();
	gPad-> SetTickx();
	gPad-> SetTicky();

	Function->SetLineColor(color);
	Function->SetLineWidth(5);

	Function->SetRange(xmin,xmax);

	if(options==""){	
		TH2F * Frame = CreateFrame(xmin,xmax,ymin,ymax,Xaxis,Yaxis);

		TLegend * leg = new TLegend(0.6,0.95,0.95,0.7);
		leg->SetName("leg");

		Frame->Draw("same");
		Function->Draw("same");
		if(legendname=="")leg->AddEntry(Function,Function->GetName(),"l");
		else leg->AddEntry(Function,legendname.c_str(),"l");
		leg->SetLineWidth(3);
		leg->SetFillColor(0);
		leg->Draw("same");

	}

	else {
		Function->Draw(options.c_str());
		TLegend * leg = (TLegend*) gPad->FindObject("leg");
		if(legendname=="") leg->AddEntry(Function,Function->GetName(),"l");
		else leg->AddEntry(Function,legendname.c_str(),"l");
		leg->SetLineWidth(3);
		leg->SetFillColor(0);
		leg->Draw("same");
	}

	return;
}

void PlotFunction(TVirtualPad * c, TSpline3 * Function, std::string Xaxis, std::string Yaxis, int color, std::string options, float xmin=-1,float xmax=-1,float ymin=-1,float ymax=-1,std::string legendname=""){
	c -> cd();
	gPad-> SetTickx();
	gPad-> SetTicky();

	Function->SetLineColor(color);
	Function->SetLineWidth(5);


	if(options==""){	
		TH2F * Frame = CreateFrame(xmin,xmax,ymin,ymax,Xaxis,Yaxis);

		TLegend * leg = new TLegend(0.6,0.95,0.95,0.7);
		leg->SetName("leg");

		Frame->Draw("same");
		Function->Draw("same");
		if(legendname=="") leg->AddEntry(Function,Function->GetName(),"l");
		else leg->AddEntry(Function,legendname.c_str(),"l");
		leg->SetLineWidth(3);
		leg->SetFillColor(0);
		leg->Draw("same");

	}

	else {
		Function->Draw(options.c_str());
		TLegend * leg = (TLegend*) gPad->FindObject("leg");
		if(legendname=="") leg->AddEntry(Function,Function->GetName(),"l");
		else leg->AddEntry(Function,legendname.c_str(),"l");
		leg->SetLineWidth(3);
		leg->SetFillColor(0);
		leg->Draw("same");
	}


	return;
}

void PlotTH1FintoGraph(TVirtualPad * c, Binning bins, TH1F * Values, std::string Xaxis, std::string Yaxis, int color,bool Ekin=false, std::string options="same", float xmin=-1,float xmax=-1,float ymin=-1,float ymax=-1,std::string legendname=""){

	c->cd();
	c->SetTopMargin(0.1);
	c->SetBottomMargin(0.1);
	TGraphErrors * Graph = new TGraphErrors();
	int a=0;
	for(int i=0;i<Values->GetNbinsX();i++){
		if(Values->GetBinContent(i+1)!=0){
			if(Ekin) Graph->SetPoint(a,bins.EkPerMassBinCent(i), Values->GetBinContent(i+1));
			else Graph->SetPoint(a,bins.GetBinCenter(i), Values->GetBinContent(i+1));
			Graph->SetPointError(a,0,Values->GetBinError(i+1));
			a++;
		}
	}	

	Graph->SetLineColor(color);
	Graph->SetLineWidth(5);
	Graph->SetLineColor(color);
	Graph->SetMarkerSize(3);
	Graph->SetMarkerStyle(8);
	Graph->SetMarkerColor(color);

	if(options==""){	
		TH2F * Frame = CreateFrame(xmin,xmax,ymin,ymax,Xaxis,Yaxis);

		TLegend * leg = new TLegend(0.6,0.95,0.95,0.7);
		leg->SetName("leg");


		Frame->Draw("same");
		Graph->Draw("ePsame");
		if(legendname=="") leg->AddEntry(Graph,Values->GetName(),"ep");             
		else leg->AddEntry(Graph,legendname.c_str(),"ep");
		leg->SetLineWidth(3);
		leg->SetFillColor(0);
		leg->Draw("same");		

	}

	else {
		Graph->Draw(options.c_str());
		TLegend * leg = (TLegend*) gPad->FindObject("leg");
		if(legendname=="")leg->AddEntry(Graph,Values->GetName(),"ep");             
		else leg->AddEntry(Graph,legendname.c_str(),"ep");
		leg->SetLineWidth(3);
		leg->SetFillColor(0);
		leg->Draw("same");
	}
	return;
}



void PlotTH1FRatiointoGraph(TVirtualPad * c, Binning bins, TH1F * Values1, TH1F * Values2, std::string Xaxis, std::string Yaxis, int color,bool Ekin=false, std::string options="same", float xmin=-1,float xmax=-1,float ymin=-1,float ymax=-1,std::string legendname1="",std::string legendname2=""){


	TH1F * Ratio1=(TH1F*) Values1->Clone();
	Ratio1->Sumw2();
	Ratio1->Divide(Values1);

	TH1F * Ratio2=(TH1F*) Values2->Clone();
	Ratio2->Sumw2();
	Ratio2->Divide(Values1);

	std::string yaxis = Yaxis + " (ratio)";

	PlotTH1FintoGraph(c,bins, Ratio1,Xaxis,yaxis,2,Ekin,options,xmin,xmax,ymin,ymax,legendname1);
	PlotTH1FintoGraph(c,bins, Ratio2,Xaxis,yaxis,4,Ekin,"ePsame",xmin,xmax,ymin,ymax,legendname2);	

	return;


}
