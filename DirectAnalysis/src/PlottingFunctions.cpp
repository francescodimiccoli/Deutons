#include "PlottingFunctions.h"

using namespace std;

void SetCanvas(TCanvas *c1){
 /*  gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
  c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetTopMargin(0.07666099);
   c1->SetBottomMargin(0.120954);
   c1->SetFrameBorderMode(0);
   c1->SetFrameLineWidth(4);
   c1->SetFrameBorderMode(0);
*/
}

TH1F * CreateHisto(std::string name,  Binning Bins,bool IsEkin){
	float x[Bins.size()+1];
	if(IsEkin) for(int i=0;i<Bins.size()+1;i++) x[i]=Bins.EkPerMasTOIBins()[i];
	else for(int i=0;i<Bins.size()+1;i++) { x[i]=Bins.RigTOIBins()[i]; cout<<x[i]<<endl; }
	
	TH1F * hist = new TH1F(name.c_str(),name.c_str(),Bins.size(),x);
	if(IsEkin) hist->GetXaxis()->SetTitle("Ekin [GeV/n]");
	else hist->GetXaxis()->SetTitle("R [GV]");
	
	return hist;
}


TH1F * ConvertBinnedHisto(TH1F * histo,std::string name,  Binning Bins,bool IsEkin){
	TH1F * plot = CreateHisto(name,Bins,IsEkin);
	for(int bin=0;bin<histo->GetNbinsX();bin++){
		plot->SetBinContent(bin+1,histo->GetBinContent(bin+1));
		plot->SetBinError(bin+1,histo->GetBinError(bin+1));
	}
	return plot;
}


TH2F* CreateFrame (TVirtualPad * c,float xmin,float xmax,float ymin, float ymax,std::string Xaxis,std::string Yaxis){
	
		TH2F * Frame = new TH2F("Frame","Frame",1e3,xmin,xmax,1e3,ymin,ymax);
		Frame->SetStats(false);
		Frame->SetTitle("");

		Frame->SetStats(false);
		Frame->SetTitle("");

		Frame->SetLineWidth(1);

		//Frame->GetYaxis()->SetMoreLogLabels();
		//Frame->GetXaxis()->SetMoreLogLabels();

		Frame->GetYaxis()->SetNoExponent();
		Frame->GetXaxis()->SetNoExponent();

		Frame->GetYaxis()->SetRangeUser(ymin,ymax);
		Frame->GetYaxis()->SetRangeUser(ymin,ymax);

		Frame->GetXaxis()->SetTitle(Xaxis.c_str());
		Frame->GetYaxis()->SetTitle(Yaxis.c_str());	

		Frame->GetXaxis()->SetTitleFont(21);
                Frame->GetYaxis()->SetTitleFont(21); 

                Frame->GetXaxis()->SetLabelFont(21);
                Frame->GetYaxis()->SetLabelFont(21); 

		Int_t wtopx,wtopy;
		UInt_t ww,wh;
		c->cd();
		gPad->GetCanvas()->GetCanvasPar(wtopx,wtopy,ww,wh);

		cout<<"CANVAS: "<<wtopx<<" "<<wtopy<<" "<<(int)ww<<" "<<(int)wh<<" "<<gPad->GetCanvas()->GetWh()<<endl;
		cout<<endl;
		int charwidthX = 25;
		int charwidthY = 25;
		
		Frame->GetXaxis()->SetLabelSize(0.035*gPad->GetCanvas()->GetWh()/(float)gPad->YtoPixel(gPad->GetY1()));
                Frame->GetYaxis()->SetLabelSize(0.035*gPad->GetCanvas()->GetWh()/(float)gPad->YtoPixel(gPad->GetY1())); 
                                                                                           
                Frame->GetXaxis()->SetTitleSize(0.035*gPad->GetCanvas()->GetWh()/(float)gPad->YtoPixel(gPad->GetY1()));
                Frame->GetYaxis()->SetTitleSize(0.035*gPad->GetCanvas()->GetWh()/(float)gPad->YtoPixel(gPad->GetY1()));	

		Frame->GetXaxis()->SetTitleOffset(1/(gPad->GetCanvas()->GetWh()/(float)gPad->YtoPixel(gPad->GetY1())));
                Frame->GetYaxis()->SetTitleOffset(1/(gPad->GetCanvas()->GetWh()/(float)gPad->YtoPixel(gPad->GetY1())));	


		Frame->GetXaxis()->CenterTitle();
		Frame->GetYaxis()->CenterTitle();

		return Frame;	
	
}


void PlotDistribution(TVirtualPad * c, TH1F * Distribution, std::string Xaxis, std::string Yaxis, int color, std::string options, float ymin,float ymax,float thickline,std::string legendname,bool filled ,bool dots,bool skipleg,int rebin){
	c -> cd();
//	gPad-> SetLogy();
	gPad-> SetTickx();
	gPad-> SetTicky();

	gPad->SetFrameLineWidth(5);
	Distribution->SetLineColor(color);
	Distribution->SetLineWidth(thickline);

//	Distribution->SetFillColor(color);
//	Distribution->SetFillStyle(3002);

	Distribution->SetStats(false);
	Distribution->SetTitle("");


	TF1 * fit = Distribution->GetFunction("fitfunc");

	if(ymin!=0&&ymax!=-1) Distribution->GetYaxis()->SetRangeUser(ymin,ymax);

	Distribution->Rebin(rebin);

	Distribution->GetXaxis()->SetTitle(Xaxis.c_str());
	Distribution->GetYaxis()->SetTitle(Yaxis.c_str());	

	Distribution->GetXaxis()->SetTitleSize(0.045);
	Distribution->GetYaxis()->SetTitleSize(0.045);	

	Distribution->GetXaxis()->SetLabelSize(0.045);
	Distribution->GetYaxis()->SetLabelSize(0.035);	

	Distribution->GetXaxis()->CenterTitle();
	Distribution->GetYaxis()->CenterTitle();

	Distribution->GetXaxis()->SetTitleFont(22);
	Distribution->GetYaxis()->SetTitleFont(22); 

	Distribution->GetXaxis()->SetLabelFont(22);
	Distribution->GetYaxis()->SetLabelFont(22); 

	Distribution->GetYaxis()->SetTitleOffset(0.8);


	Distribution->SetMarkerColor(color);
	Distribution->SetMarkerStyle(8);
	Distribution->SetMarkerSize(1.5);
	if(dots) Distribution->SetMarkerSize(thickline);



	if(filled){
		Distribution->SetFillColor(color);
		Distribution->SetFillStyle(3002);
		Distribution->SetLineStyle(2);
	}
	if(dots) Distribution->Draw((options).c_str());
	Distribution->Draw((options+",hist").c_str());
 	
	if(!skipleg){
		TLegend * leg = (TLegend*) gPad->FindObject("leg");
		if(leg){
			if(legendname=="") leg->AddEntry(Distribution,Distribution->GetName());
			else {
				if(dots) leg->AddEntry(Distribution,legendname.c_str(),"p");
				 leg->AddEntry(Distribution,legendname.c_str(),"l");
			}

			leg->SetTextFont(32);
			leg->Draw("same");
		}

		else{
			TLegend *leg = new TLegend(0.5510683,0.4323081,0.8697238,0.8649037,NULL,"brNDC");
			leg->SetBorderSize(0);
			leg->SetTextFont(32);
			leg->SetLineColor(1);
			leg->SetLineStyle(1);
			leg->SetLineWidth(0);
			leg->SetLineColor(0);
			leg->SetFillColor(0);
			leg->SetFillStyle(0);
			leg->SetName("leg");
			if(legendname=="") leg->AddEntry(Distribution,Distribution->GetName());
			else {
				if(dots) leg->AddEntry(Distribution,legendname.c_str(),"p"); 
				leg->AddEntry(Distribution,legendname.c_str(),"l"); 
			}
			leg->Draw("same"); 
		}
	}
	if(fit){
		fit->SetLineColor(color);
		fit->SetLineWidth(5);
		fit->Draw("same");

	}

	return;
}

void PlotTH1F(TVirtualPad * c, TH1F * Distribution, std::string Xaxis, std::string Yaxis, int color, std::string options,std::string legendname, float thickness, bool skipleg){

	c -> cd();
	gPad-> SetLogz();
	gPad-> SetTickx();
	gPad-> SetTicky();


	Distribution->SetStats(false);
	Distribution->SetTitle("");

	Distribution->SetLineColor(color);
	Distribution->SetMarkerColor(color);

	Distribution->SetLineWidth(thickness);
	Distribution->SetMarkerSize(thickness);


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

	Distribution->Draw((options+"same").c_str());

	if(!skipleg){
		TLegend * leg = (TLegend*) gPad->FindObject("leg");
		if(leg){
			if(legendname=="") leg->AddEntry(Distribution,Distribution->GetName());
				else leg->AddEntry(Distribution,legendname.c_str(),(options).c_str());

			leg->SetTextFont(32);
			leg->Draw("same");
		}

		else{
			leg = new TLegend(0.7, 0.7,0.95,0.95);
			leg->SetName("leg");
			if(legendname=="") leg->AddEntry(Distribution,Distribution->GetName());
				else leg->AddEntry(Distribution,legendname.c_str(),(options).c_str()); 
			leg->SetTextFont(32);
			leg->SetLineWidth(5);
			leg->SetFillColor(0);
			leg->Draw("same"); 
		}
	}


	return;
}

void PlotTH2F(TVirtualPad * c, TH2F * Distribution, std::string Xaxis, std::string Yaxis, std::string options){

	c -> cd();
	gPad-> SetTickx();
	gPad-> SetTicky();


	Distribution->SetStats(false);
	Distribution->SetTitle("");
	Distribution->SetLineWidth(4);
	
	Distribution->GetXaxis()->SetTitle(Xaxis.c_str());
	Distribution->GetYaxis()->SetTitle(Yaxis.c_str());	

	Distribution->GetXaxis()->CenterTitle();
	Distribution->GetYaxis()->CenterTitle();

	Distribution->GetXaxis()->SetTitleFont(23);
	Distribution->GetYaxis()->SetTitleFont(23); 

	Distribution->GetXaxis()->SetLabelFont(23);
	Distribution->GetYaxis()->SetLabelFont(23); 

	Distribution->GetXaxis()->SetTitleSize(40);
	Distribution->GetYaxis()->SetTitleSize(40);	
	Distribution->GetXaxis()->SetLabelSize(25);
	Distribution->GetYaxis()->SetLabelSize(25);	


	Distribution->Draw((options+"same").c_str());

	return;
}

void PlotGraph(TVirtualPad * c,TGraphErrors * graph,std::string Xaxis, std::string Yaxis, int color,std::string options, float xmin,float xmax,float ymin,float ymax,std::string legendname,int dotstyle,bool skipleg) {

	c -> cd();
	gPad-> SetTickx();
	gPad-> SetTicky();
	gPad-> SetGridx();
	gPad-> SetGridy();


	graph->SetLineColor(color);
	graph->SetMarkerColor(color);
	graph->SetMarkerStyle(dotstyle);
	graph->SetMarkerSize(3);
	graph->SetLineWidth(5);

	if(xmin!=-1) graph->GetXaxis()->SetRangeUser(xmin,xmax);
	if(ymin!=-1) graph->GetYaxis()->SetRangeUser(ymin,ymax);
	TH2F * Frame = (TH2F*) gPad->FindObject("Frame");

	if(Frame){
		TLegend * leg = (TLegend*) gPad->FindObject("leg");
		leg->SetName("leg");

		graph->Draw((options+"same").c_str());	
		if(legendname=="")leg->AddEntry(graph,graph->GetName(),"l");
		else leg->AddEntry(graph,legendname.c_str(),"l");
		leg->SetLineWidth(3);
		leg->SetFillColor(0);
		leg->Draw("same");
	}
	else{

		TH2F * Frame = CreateFrame(gPad,xmin,xmax,ymin,ymax,Xaxis,Yaxis);
		Frame->Draw();
		graph->Draw((options+"same").c_str());

		TLegend * leg = new TLegend(0.6,0.95,0.95,0.7);
		leg->SetName("leg");

		if(legendname=="") leg->AddEntry(graph,graph->GetName(),"P");             
		else if(!skipleg) leg->AddEntry(graph,legendname.c_str(),"P");
		leg->SetLineWidth(3);
		leg->SetFillColor(0);
		leg->SetTextFont(32);
		leg->Draw("same");		

	}

	return;

}




void PlotFunction(TVirtualPad * c, TF1 * Function, std::string Xaxis, std::string Yaxis, int color, std::string options, float xmin,float xmax,float ymin,float ymax,std::string legendname){
	c -> cd();
	gPad-> SetTickx();
	gPad-> SetTicky();

	Function->SetLineColor(color);
	Function->SetLineWidth(5);

	Function->SetRange(xmin,xmax);

	if(options==""){	
		TH2F * Frame = CreateFrame(gPad,xmin,xmax,ymin,ymax,Xaxis,Yaxis);

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

void PlotFunction(TVirtualPad * c, TSpline3 * Function, std::string Xaxis, std::string Yaxis, int color, std::string options, float xmin,float xmax,float ymin,float ymax,std::string legendname){
	c -> cd();
	gPad-> SetTickx();
	gPad-> SetTicky();

	Function->SetLineColor(color);
	Function->SetLineWidth(5);


	if(options==""){	
		TH2F * Frame = CreateFrame(gPad,xmin,xmax,ymin,ymax,Xaxis,Yaxis);

		TLegend * leg = new TLegend(0.6,0.95,0.95,0.7);
		leg->SetName("leg");

		Frame->Draw("same");
		Function->Draw("same");
		if(legendname=="") leg->AddEntry(Function,Function->GetName(),"l");
		else leg->AddEntry(Function,legendname.c_str(),"l");
		leg->SetLineWidth(3);
		leg->SetFillColor(0);
		leg->SetTextFont(32);
		leg->Draw("same");

	}

	else {
		Function->Draw(options.c_str());
		TLegend * leg = (TLegend*) gPad->FindObject("leg");
		if(legendname=="") leg->AddEntry(Function,Function->GetName(),"l");
		else leg->AddEntry(Function,legendname.c_str(),"l");
		leg->SetLineWidth(3);
		leg->SetFillColor(0);
		 leg->SetTextFont(32);
		leg->Draw("same");
	}


	return;
}

void PlotTH1FintoGraph(TVirtualPad * c, Binning bins, TH1F * Values, std::string Xaxis, std::string Yaxis, int color,bool Ekin, std::string options, float xmin,float xmax,float ymin,float ymax,std::string legendname,int dotstyle,bool skipleg, bool cleanhigherrors, TH1F* ErrUp, TH1F* ErrDw){

	c->cd();
	c->SetTopMargin(0.1);
	c->SetBottomMargin(0.1);
	gPad->SetFrameLineWidth(5);
	gPad->SetTickx();
	gPad->SetTicky();
	
	TGraphAsymmErrors * Graph = new TGraphAsymmErrors();
	int a=0;
	for(int i=0;i<Values->GetNbinsX();i++){
		if(Values->GetBinContent(i+1)!=0){
			if(Ekin) Graph->SetPoint(a,bins.EkPerMassBinCent(i), Values->GetBinContent(i+1));
			else Graph->SetPoint(a,bins.GetBinCenter(i), Values->GetBinContent(i+1));
			if(!ErrUp&&!ErrDw){
				Graph->SetPointError(a,0,0,Values->GetBinError(i+1),Values->GetBinError(i+1));
				if(cleanhigherrors) if(Values->GetBinError(i+1)/Values->GetBinContent(i+1)>0.05) Graph->RemovePoint(a);
			}
			else{
				Graph->SetPointError(a,0,0,ErrUp->GetBinContent(i+1),ErrDw->GetBinContent(i+1));	
			}
			a++;
		}
	}	

	Graph->SetLineColor(color);
	Graph->SetFillColor(color);
	Graph->SetFillStyle(3001);
	Graph->SetLineWidth(1.5);
	Graph->SetLineColor(color);
	Graph->SetMarkerSize(1.5);
	Graph->SetMarkerStyle(dotstyle);
	Graph->SetMarkerColor(color);

	TH2F * Frame = (TH2F*) gPad->FindObject("Frame");	
		
	if(Frame){
		Graph->Draw(options.c_str());
		TLegend * leg = (TLegend*) gPad->FindObject("leg");
		if(legendname=="")leg->AddEntry(Graph,Values->GetName(),"P");
		else if(!skipleg) leg->AddEntry(Graph,legendname.c_str(),"P");
		leg->SetLineWidth(3);
		leg->SetFillColor(0);
		leg->SetTextFont(32);
		leg->Draw("same");

	}

	else{

		TH2F * Frame = CreateFrame(gPad,xmin,xmax,ymin,ymax,Xaxis,Yaxis);
		Frame->Draw();
		Graph->Draw((options+"same").c_str());

		TLegend * leg = new TLegend(0.6,0.9,0.95,0.6);
		leg->SetName("leg");
		leg->SetBorderSize(0);
		leg->SetTextFont(32);
		leg->SetLineColor(1);
		leg->SetLineStyle(1);
		leg->SetLineWidth(0);
		leg->SetLineColor(0);
		leg->SetFillColor(0);
		leg->SetFillStyle(0);

		if(legendname=="") leg->AddEntry(Graph,Values->GetName(),"P");             
		else if(!skipleg) leg->AddEntry(Graph,legendname.c_str(),"P");
		leg->SetLineWidth(3);
		leg->SetFillColor(0);
		leg->SetTextFont(32);
		leg->Draw("same");		

	}
	
	return;
}







void PlotMergedRanges(TVirtualPad * c, TH1F * ValuesTOF, TH1F* ValuesNaF, TH1F* ValuesAgl, std::string Xaxis, std::string Yaxis, int color,bool Ekin, std::string options, float xmin,float xmax,float ymin,float ymax,std::string legendname,int dotstyle,bool skipleg, bool cleanhigherrors){

	c->cd();
	c->SetTopMargin(0.1);
	c->SetBottomMargin(0.1);
	gPad->SetFrameLineWidth(5);
	gPad->SetTickx();
	gPad->SetTicky();
	
	TGraphErrors * Graph = new TGraphErrors();
	int a=0;
	for(int i=0;i<18;i++){
		   if(ValuesTOF->GetBinContent(i+1)!=0){
			if(Ekin) Graph->SetPoint(a,Global.GetToFDBins().EkPerMassBinCent(i), ValuesTOF->GetBinContent(i+1));
			else Graph->SetPoint(a,Global.GetToFDBins().GetBinCenter(i), ValuesTOF->GetBinContent(i+1));
			Graph->SetPointError(a,0,ValuesTOF->GetBinError(i+1));
			if(cleanhigherrors) if(ValuesTOF->GetBinError(i+1)/ValuesTOF->GetBinContent(i+1)>0.05) Graph->RemovePoint(a);
			a++;
		}
	}	

	for(int i=1;i<17;i++){
		if(ValuesNaF->GetBinContent(i+1)!=0){
			if(Ekin) Graph->SetPoint(a,Global.GetNaFDBins().EkPerMassBinCent(i), ValuesNaF->GetBinContent(i+1));
			else Graph->SetPoint(a,Global.GetNaFDBins().GetBinCenter(i), ValuesNaF->GetBinContent(i+1));
			Graph->SetPointError(a,0,ValuesNaF->GetBinError(i+1));
			if(cleanhigherrors) if(ValuesNaF->GetBinError(i+1)/ValuesNaF->GetBinContent(i+1)>0.05) Graph->RemovePoint(a);
			a++;
		}
	}	
	for(int i=5;i<17;i++){
		if(ValuesAgl->GetBinContent(i+1)!=0){
			if(Ekin) Graph->SetPoint(a,Global.GetAglDBins().EkPerMassBinCent(i), ValuesAgl->GetBinContent(i+1));
			else Graph->SetPoint(a,Global.GetAglDBins().GetBinCenter(i), ValuesAgl->GetBinContent(i+1));
			Graph->SetPointError(a,0,ValuesAgl->GetBinError(i+1));
			if(cleanhigherrors) if(ValuesAgl->GetBinError(i+1)/ValuesAgl->GetBinContent(i+1)>0.05) Graph->RemovePoint(a);
			a++;
		}
	}	
	Graph->SetLineColor(color);
	Graph->SetFillColor(color);
	Graph->SetFillStyle(3001);
	Graph->SetLineWidth(3);
	Graph->SetLineColor(color);
	Graph->SetMarkerSize(1.5);
	Graph->SetMarkerStyle(dotstyle);
	Graph->SetMarkerColor(color);

	TH2F * Frame = (TH2F*) gPad->FindObject("Frame");	
		
	if(Frame){
		Graph->Draw(options.c_str());
		TLegend * leg = (TLegend*) gPad->FindObject("leg");
		if(legendname=="")leg->AddEntry(Graph,"","P");
		else if(!skipleg) leg->AddEntry(Graph,legendname.c_str(),"P");
		leg->SetLineWidth(3);
		leg->SetFillColor(0);
		leg->SetTextFont(32);
		leg->Draw("same");

	}

	else{

		TH2F * Frame = CreateFrame(gPad,xmin,xmax,ymin,ymax,Xaxis,Yaxis);
		Frame->Draw();
		Graph->Draw((options+"same").c_str());

		TLegend * leg = new TLegend(0.6,0.9,0.95,0.6);
		leg->SetName("leg");
		leg->SetBorderSize(0);
		leg->SetTextFont(32);
		leg->SetLineColor(1);
		leg->SetLineStyle(1);
		leg->SetLineWidth(0);
		leg->SetLineColor(0);
		leg->SetFillColor(0);
		leg->SetFillStyle(0);

		if(legendname=="") leg->AddEntry(Graph,"","P");             
		else if(!skipleg) leg->AddEntry(Graph,legendname.c_str(),"P");
		leg->SetLineWidth(3);
		leg->SetFillColor(0);
		leg->SetTextFont(32);
		leg->Draw("same");		

	}
	
	return;

}


void PlotTH1FRatiointoGraph(TVirtualPad * c, Binning bins, TH1F * Values1, TH1F * Values2, std::string Xaxis, std::string Yaxis, int color,bool Ekin, std::string options, float xmin,float xmax,float ymin,float ymax,std::string legendname1,std::string legendname2){


	TH1F * Ratio1=(TH1F*) Values1->Clone();
	Ratio1->Sumw2();
	Ratio1->Divide(Values1);

	TH1F * Ratio2=(TH1F*) Values2->Clone();
	Ratio2->Sumw2();
	Ratio2->Divide(Values1);

	std::string yaxis = Yaxis + " (ratio)";

	PlotTH1FintoGraph(c,bins, Ratio1,Xaxis,yaxis,2,Ekin,options,xmin,xmax,ymin,ymax,legendname1,8,true);
	PlotTH1FintoGraph(c,bins, Ratio2,Xaxis,yaxis,4,Ekin,"ePsame",xmin,xmax,ymin,ymax,legendname2,8,true);	

	return;


}

void PlotRatioWithSplineintoGraph(TVirtualPad * c, Binning bins, TH1F * Values1,TSpline3 * Spline, std::string Xaxis, std::string Yaxis, int color,bool Ekin, std::string options, float xmin,float xmax,float ymin,float ymax,std::string legendname1,int dotstyle){
	TH1F * Ratio1=(TH1F*) Values1->Clone();
        Ratio1->Sumw2();
	TH1F * Errors = (TH1F*) Ratio1->Clone();
	
	
	if(Ekin){
		for(int i=0;i<Values1->GetNbinsX();i++){
			Ratio1->SetBinContent(i+1,Values1->GetBinContent(i+1)/Spline->Eval(bins.EkPerMassBinCent(i)));
			Ratio1->SetBinError(i+1,1.0*Values1->GetBinError(i+1)/Spline->Eval(bins.EkPerMassBinCent(i)));			
		}

		for(int i=0;i<Errors->GetNbinsX();i++){
			Errors->SetBinContent(i+1,1);
			Errors->SetBinError(i+1,1.0*Values1->GetBinError(i+1)/Spline->Eval(bins.EkPerMassBinCent(i)));			
		}
	}
	else{
		for(int i=0;i<Values1->GetNbinsX();i++){
			Ratio1->SetBinContent(i+1,Values1->GetBinContent(i+1)/Spline->Eval(bins.RigTOIBinsCent()[i]));
			Ratio1->SetBinError(i+1,1.0*Values1->GetBinError(i+1)/Spline->Eval(bins.RigTOIBinsCent()[i]));			
		}

		for(int i=0;i<Errors->GetNbinsX();i++){
			Errors->SetBinContent(i+1,1);
			Errors->SetBinError(i+1,1.0*Values1->GetBinError(i+1)/Spline->Eval(bins.RigTOIBinsCent()[i]));			
		}
	}
	PlotTH1FintoGraph(c,bins, Ratio1,Xaxis,Yaxis,color,Ekin,"Psame",xmin,xmax,ymin,ymax,legendname1,dotstyle);
//	PlotTH1FintoGraph(c,bins, Errors,Xaxis,Yaxis,color,Ekin,"3same",xmin,xmax,ymin,ymax,legendname1,dotstyle,true);

}
void PlotRatioWithSplineiAvg(TVirtualPad * c, Binning bins, TH1F * Values1,TSpline3 * Spline, std::string Xaxis, std::string Yaxis, int color,bool drawError , bool Ekin, std::string options, float xmin,float xmax,float ymin,float ymax,std::string legendname1,int dotstyle){
	TH1F * Ratio1=(TH1F*) Values1->Clone();
        Ratio1->Sumw2();
	TH1F * Errors_dw = CreateHisto("Errors_dw", bins,Ekin);
	TH1F * Errors_up = CreateHisto("Errors_up", bins,Ekin);
	TH1F * Compare=(TH1F*) Values1->Clone();
	
	if(Ekin){
		for(int i=0;i<Values1->GetNbinsX();i++){
			Ratio1->SetBinContent(i+1,Values1->GetBinContent(i+1)/((Spline->Eval(bins.EkPerMassBinCent(i))+Values1->GetBinContent(i+1))/2)	);
			Compare->SetBinContent(i+1,Spline->Eval(bins.EkPerMassBinCent(i))/((Spline->Eval(bins.EkPerMassBinCent(i))+Values1->GetBinContent(i+1))/2)  );
			Ratio1->SetBinError(i+1,1.0*Values1->GetBinError(i+1)/((Spline->Eval(bins.EkPerMassBinCent(i))+Values1->GetBinContent(i+1))/2)	);			
			Compare->SetBinError(i+1,0);			
		
		}
		for(int i=0;i<Errors_up->GetNbinsX();i++){
			if(Values1->GetBinError(i+1)>0 && Spline->Eval(bins.RigTOIBinsCent()[i])>0){
				float err = pow((pow(Values1->GetBinError(i+1)/Values1->GetBinContent(i+1),2)+0.01*0.01),0.5)*Values1->GetBinContent(i+1);
				Errors_dw->SetBinContent(i+1,1-err/((Spline->Eval(bins.EkPerMassBinCent(i))+Values1->GetBinContent(i+1))/2));
				Errors_dw->SetBinError(i+1,0);			
				Errors_up->SetBinContent(i+1,1+err/((Spline->Eval(bins.EkPerMassBinCent(i))+Values1->GetBinContent(i+1))/2));
				Errors_up->SetBinError(i+1,0);			
			}
			else  {
				Errors_dw->SetBinContent(i+1,0);
				Errors_up->SetBinContent(i+1,2);
				}		
	
		}
	}
	else{
		for(int i=0;i<Values1->GetNbinsX();i++){
			Ratio1->SetBinContent(i+1,Values1->GetBinContent(i+1)/((Spline->Eval(bins.RigTOIBinsCent()[i])+Values1->GetBinContent(i+1))/2)	);
			Ratio1->SetBinError(i+1,1.0*Values1->GetBinError(i+1)/((Spline->Eval(bins.RigTOIBinsCent()[i])+Values1->GetBinContent(i+1))/2)	);			
			Compare->SetBinContent(i+1,Spline->Eval(bins.RigTOIBinsCent()[i])/((Spline->Eval(bins.RigTOIBinsCent()[i])+Values1->GetBinContent(i+1))/2)	);
			Compare->SetBinError(i+1,0	);			
		}

		for(int i=0;i<Errors_up->GetNbinsX();i++){
			if(Values1->GetBinError(i+1)>0 && Spline->Eval(bins.RigTOIBinsCent()[i])>0){
				float err = pow((pow(Values1->GetBinError(i+1)/Values1->GetBinContent(i+1),2)+0.01*0.01),0.5)*Values1->GetBinContent(i+1);
				Errors_dw->SetBinContent(i+1,1-err/((Spline->Eval(bins.RigTOIBinsCent()[i])+Values1->GetBinContent(i+1))/2));
				Errors_dw->SetBinError(i+1,0	);			
				Errors_up->SetBinContent(i+1,1+err/((Spline->Eval(bins.RigTOIBinsCent()[i])+Values1->GetBinContent(i+1))/2));
				Errors_up->SetBinError(i+1,0	);			
			}
			else  {
				Errors_dw->SetBinContent(i+1,0);
				Errors_up->SetBinContent(i+1,2);
				}		
		}
	}
	PlotTH1FintoGraph(c,bins, Ratio1,Xaxis,Yaxis,color,Ekin,"Psame",xmin,xmax,ymin,ymax,legendname1,dotstyle);
	PlotTH1FintoGraph(c,bins, Compare,Xaxis,Yaxis,1,Ekin,"Psame",xmin,xmax,ymin,ymax,legendname1,dotstyle);
	if(drawError){
		Errors_dw->Draw("hist,same");
		Errors_up->Draw("hist,same");
	}
}
