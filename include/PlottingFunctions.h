using namespace std;

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
		TH2F * Frame = new TH2F("Frame","Frame",100,xmin,xmax,100,ymin,ymax);
		Frame->SetStats(false);
		Frame->SetTitle("");

		TLegend * leg = new TLegend(0.6,0.95,0.95,0.7);
		leg->SetName("leg");

		Frame->GetYaxis()->SetMoreLogLabels();
		Frame->GetXaxis()->SetMoreLogLabels();

		Frame->GetYaxis()->SetNoExponent();
		Frame->GetXaxis()->SetNoExponent();

		Frame->GetYaxis()->SetRangeUser(ymin,ymax);
		Frame->GetYaxis()->SetRangeUser(ymin,ymax);

		Frame->GetXaxis()->SetTitle(Xaxis.c_str());
		Frame->GetYaxis()->SetTitle(Yaxis.c_str());	

		Frame->GetXaxis()->SetTitleSize(0.045);
		Frame->GetYaxis()->SetTitleSize(0.045);	

		Frame->GetXaxis()->CenterTitle();
		Frame->GetYaxis()->CenterTitle();

		Frame->GetXaxis()->SetTitleFont(32);
		Frame->GetYaxis()->SetTitleFont(32); 

		Frame->GetXaxis()->SetLabelFont(32);
		Frame->GetYaxis()->SetLabelFont(32); 

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
		TH2F * Frame = new TH2F("Frame","Frame",100,xmin,xmax,100,ymin,ymax);
		Frame->SetStats(false);
		Frame->SetTitle("");

		TLegend * leg = new TLegend(0.6,0.95,0.95,0.7);
		leg->SetName("leg");


		Frame->GetYaxis()->SetMoreLogLabels();
		Frame->GetXaxis()->SetMoreLogLabels();

		Frame->GetYaxis()->SetNoExponent();
		Frame->GetXaxis()->SetNoExponent();

		Frame->GetYaxis()->SetRangeUser(ymin,ymax);
		Frame->GetYaxis()->SetRangeUser(ymin,ymax);

		Frame->GetXaxis()->SetTitle(Xaxis.c_str());
		Frame->GetYaxis()->SetTitle(Yaxis.c_str());	

		Frame->GetXaxis()->SetTitleSize(0.045);
		Frame->GetYaxis()->SetTitleSize(0.045);	

		Frame->GetXaxis()->CenterTitle();
		Frame->GetYaxis()->CenterTitle();

		Frame->GetXaxis()->SetTitleFont(32);
		Frame->GetYaxis()->SetTitleFont(32); 

		Frame->GetXaxis()->SetLabelFont(32);
		Frame->GetYaxis()->SetLabelFont(32); 

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

	c -> cd();
	gPad-> SetTickx();
	gPad-> SetTicky();
	gPad-> SetGridx();
	gPad-> SetGridy();
	gPad-> SetLogx();

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
		TH2F * Frame = new TH2F("Frame","Frame",100,xmin,xmax,100,ymin,ymax);

		TLegend * leg = new TLegend(0.6,0.95,0.95,0.7);
		leg->SetName("leg");

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

		Frame->GetXaxis()->SetTitleSize(0.045);
		Frame->GetYaxis()->SetTitleSize(0.045);	

		Frame->GetXaxis()->CenterTitle();
		Frame->GetYaxis()->CenterTitle();

		Frame->GetXaxis()->SetTitleFont(32);
		Frame->GetYaxis()->SetTitleFont(32); 

		Frame->GetXaxis()->SetLabelFont(32);
		Frame->GetYaxis()->SetLabelFont(32); 

		Frame->Draw();
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

