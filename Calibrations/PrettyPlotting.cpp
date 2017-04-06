#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "../include/binning.h"
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include <TFile.h>
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "../include/GlobalBinning.h"

#include "../Ntuple-making/Commonglobals.cpp"
#include "../Ntuple-making/Variables.hpp"
#include "../include/filesaver.h"

#include "../include/FitError.h"
#include "../include/Resolution.h"
#include "../include/PlottingFunctions.h"


int colorbase=55;


void PlotDistrib(FileSaver Plots,Resolution * Reso,std::string Xaxis,std::string Yaxis, float xmin,float xmax,float ymin,float ymax){

	TLegend * leg = new TLegend(0.4, 0.7,0.95,0.95);
        for(int i=Reso->GetBinning().size()-1; i>0;i--){
                leg->AddEntry(Reso->Get_Histo(i),("R:" + to_string(Reso->GetBinning().GetBinLowEdge(i))).c_str());
        }


	TCanvas * c1 = new TCanvas("Distributions");
	c1->cd();
	for(int i=Reso->GetBinning().size()-1; i>0;i--)
		PlotDistribution(gPad, Reso->Get_Histo(i),Xaxis.c_str(),Yaxis.c_str(),colorbase + i,"same",ymin,ymax);
	
	leg->Draw("same");

	TCanvas * c2 = new TCanvas("Fits");
	c2->cd();
	for(int i=Reso->GetBinning().size()-1; i>0;i--){
		TF1 * fit = Reso->Get_Histo(i)->GetFunction("fitfunc");
		fit->SetName(to_string(i).c_str());
		if(i==Reso->GetBinning().size()-1) 
			PlotFunction(gPad,fit,Xaxis.c_str(),Yaxis.c_str(),colorbase + i,"",xmin,xmax,ymin,ymax);
		else PlotFunction(gPad,fit,Xaxis.c_str(),Yaxis.c_str(),colorbase + i,"same",xmin,xmax,ymin,ymax);
	}

	leg->Draw("same");

	Plots.Add(c1);
	Plots.Add(c2);
	Plots.writeObjsInFolder(Reso->GetName().c_str());

	delete c1,c2;

	return;
}

TH1F * GetResolution(Resolution * Reso){

	TH1F * resolution = (TH1F*)Reso->Get_Sigmas()->Clone();
	TH1F * means = new TH1F("means","means",Reso->GetBinning().size(),0,Reso->GetBinning().size());

	for(int i=0;i<Reso->GetBinning().size();i++){
		means->SetBinContent(i+1,Reso->GetBinning().GetBinCenter(i));
		means->SetBinError(i+1,0);
	}

	resolution->Sumw2();

	resolution->Multiply(means); //sigma(1/x)/(1/x)=sigma(x)/x

	return resolution;
}


void PlotFitResults(FileSaver Plots,Resolution * Reso,std::string Xaxis, float xmin=-1,float xmax=-1,float ymin=-1,float ymax=-1){

	TCanvas * c1 = new TCanvas("Sigma");
	c1->cd();
	PlotTH1FintoGraph(gPad,Reso->GetBinning(), Reso->Get_Sigmas(), Xaxis.c_str(), "#sigma",2,false,"",xmin,xmax,ymin,ymax);
	
	TCanvas * c2 = new TCanvas("Mean Shift");
	c2->cd();
	PlotTH1FintoGraph(gPad,Reso->GetBinning(), Reso->Get_Means(), Xaxis.c_str(), "Mean Shift (rec/gen)",2,false,"",xmin,xmax,ymin,ymax);

	TCanvas * c3 = new TCanvas("Resolution");
        c3->cd();
        PlotTH1FintoGraph(gPad,Reso->GetBinning(), Reso->Get_Resolutions(), Xaxis.c_str(), "Resolution (#sigma/mean)",2,false,"",xmin,xmax,ymin,ymax);


	Plots.Add(c1);
	Plots.Add(c2);
	Plots.Add(c3);
	Plots.writeObjsInFolder(Reso->GetName().c_str());

	return;
}


void PlotFitResults(FileSaver Plots,Resolution * Reso1,Resolution * Reso2,std::string Xaxis, float xmin=-1,float xmax=-1,float ymin=-1,float ymax=-1){

	TCanvas * c1 = new TCanvas("Sigma");
	
	c1->cd();
	PlotTH1FintoGraph(gPad,Reso1->GetBinning(), Reso1->Get_Sigmas(), Xaxis.c_str(),  "#sigma",2,false,"",xmin,xmax,ymin,ymax,"protons MC");
	PlotTH1FintoGraph(gPad,Reso2->GetBinning(), Reso2->Get_Sigmas(), Xaxis.c_str(),  "#sigma",4,false,"ePsame",xmin,xmax,ymin,ymax,"deuterons MC");
	TObject * Model1 = Reso1->Get_Model();
	TObject * Model2 = Reso2->Get_Model();
	string tf1="TF1";
	
	if(Model1->ClassName()==tf1) {
		TF1 * model1 = (TF1*) Reso1->Get_Model();
		PlotFunction(gPad,model1, Xaxis.c_str(),  "#sigma",2,"same",xmin,xmax,ymin,ymax,"protons param.");
	}
	else {
		TSpline3 * model1 = (TSpline3*) Reso1->Get_Model();
		PlotFunction(gPad,model1, Xaxis.c_str(),  "#sigma",2,"same",xmin,xmax,ymin,ymax,"protons param.");
	}
	

	if(Model2->ClassName()==tf1) {
		TF1 * model2 = (TF1*) Reso2->Get_Model();
		PlotFunction(gPad,model2, Xaxis.c_str(),  "#sigma",4,"same",xmin,xmax,ymin,ymax,"deuterons param.");
	}
	else {
		TSpline3 * model2 = (TSpline3*) Reso2->Get_Model();
		PlotFunction(gPad,model2, Xaxis.c_str(),  "#sigma",4,"same",xmin,xmax,ymin,ymax,"deuterons param.");
	}


	TCanvas * c2 = new TCanvas("Mean Shift");
	c2->cd();
	PlotTH1FintoGraph(gPad,Reso1->GetBinning(), Reso1->Get_Means(), Xaxis.c_str(), "Mean Shift (rec/gen)",2,false,"",xmin,xmax,ymin,ymax,"protons MC");
	PlotTH1FintoGraph(gPad,Reso2->GetBinning(), Reso2->Get_Means(), Xaxis.c_str(), "Mean Shift (rec/gen)",4,false,"ePsame",xmin,xmax,ymin,ymax,"deuterons MC");


        TCanvas * c3 = new TCanvas("Resolution");
        c3->cd();
        PlotTH1FintoGraph(gPad,Reso1->GetBinning(), Reso1->Get_Resolutions(), Xaxis.c_str(), "Resolution (#sigma/mean)",2,false,"",xmin,xmax,ymin,ymax,"protons MC");
	PlotTH1FintoGraph(gPad,Reso2->GetBinning(), Reso2->Get_Resolutions(), Xaxis.c_str(), "Resolution (#sigma/mean)",4,false,"ePsame",xmin,xmax,ymin,ymax,"deuterons MC");



	Plots.Add(c1);
	Plots.Add(c2);
	Plots.Add(c3);
	Plots.writeObjsInFolder(Reso1->GetName().c_str());
	return;
}


int main(int argc, char * argv[]){

	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	string INPUT(argv[1]);
	string OUTPUT(argv[2]);	

	FileSaver finalHistos;
	FileSaver Plots;
	
        finalHistos.setName(INPUT.c_str());
	Plots.setName(OUTPUT.c_str());
		
	
	bool checkfile = finalHistos.CheckFile();	


	cout<<"****************************** BINS ***************************************"<<endl;

	SetBins();	

	PRB.Print();
	DRB.Print();

	cout<<"**TOF**"<<endl;
	ToFDB.Print();
	ToFPB.Print();

	cout<<"**NaF**"<<endl;
	NaFDB.Print();
	NaFPB.Print();

	cout<<"**Agl**"<<endl;
	AglDB.Print();
	AglPB.Print();

	ToFDB.UseBetaEdges();
	ToFPB.UseBetaEdges();
	NaFDB.UseBetaEdges();
	NaFPB.UseBetaEdges();
	AglDB.UseBetaEdges();
	AglPB.UseBetaEdges();

	DRB.UseREdges();
	PRB.UseREdges();
	
	cout<<endl;
	
	cout<<"****************************** PLOTTING ***************************************"<<endl;


	Resolution * RigidityResolution_P = new Resolution(finalHistos,"RvsR Resolution (P)",PRB);
	Resolution * RigidityResolution_D = new Resolution(finalHistos,"RvsR Resolution (D)",PRB);

	PlotDistrib(Plots,RigidityResolution_P,"1/R_{gen}-1/R_{meas}","Normalized Distr.",-0.5,1.5,1e-4,1);
	PlotFitResults(Plots,RigidityResolution_P,RigidityResolution_D,"R_{gen} [GV]", 0.4,100,-0.1,1);


	Resolution * RigidityTOFResolution_P = new Resolution(finalHistos,"RvsBetaTOF Resolution (P)",ToFPB);
	Resolution * RigidityTOFResolution_D = new Resolution(finalHistos,"RvsBetaTOF Resolution (D)",ToFDB);

	PlotFitResults(Plots,RigidityTOFResolution_P,RigidityTOFResolution_D,"#beta_{gen}",0.45,1,-0.1,1);

        Resolution * RigidityNaFResolution_P = new Resolution(finalHistos,"RvsBetaNaF Resolution (P)",NaFPB);
	Resolution * RigidityNaFResolution_D = new Resolution(finalHistos,"RvsBetaNaF Resolution (D)",NaFDB);

	PlotFitResults(Plots,RigidityNaFResolution_P,RigidityNaFResolution_D,"#beta_{gen}",0.7,1,-0.1,1);

        Resolution * RigidityAglResolution_P = new Resolution(finalHistos,"RvsBetaAgl Resolution (P)",AglPB);
	Resolution * RigidityAglResolution_D = new Resolution(finalHistos,"RvsBetaAgl Resolution (D)",AglDB);

	PlotFitResults(Plots,RigidityAglResolution_P,RigidityAglResolution_D,"#beta_{gen}",0.94,1,-0.1,1);


	Resolution * BetaTOF_RResolution_P = new Resolution(finalHistos,"BetaTOFvsR Resolution (P)",PRB);
	Resolution * BetaTOF_RResolution_D = new Resolution(finalHistos,"BetaTOFvsR Resolution (D)",PRB);

	PlotFitResults(Plots,BetaTOF_RResolution_P,BetaTOF_RResolution_D,"R_{gen}",0.4,100,-0.05,0.1);

        Resolution * BetaNaF_RResolution_P = new Resolution(finalHistos,"BetaNaFvsR Resolution (P)",PRB);
	Resolution * BetaNaF_RResolution_D = new Resolution(finalHistos,"BetaNaFvsR Resolution (D)",PRB);

	PlotFitResults(Plots,BetaNaF_RResolution_P,BetaNaF_RResolution_D,"R_{gen}",0.4,100,-0.005,0.018);

        Resolution * BetaAgl_RResolution_P = new Resolution(finalHistos,"BetaAglvsR Resolution (P)",PRB);
	Resolution * BetaAgl_RResolution_D = new Resolution(finalHistos,"BetaAglvsR Resolution (D)",PRB);

	PlotFitResults(Plots,BetaAgl_RResolution_P,BetaAgl_RResolution_D,"R_{gen}",0.4,100,-0.0005,0.0018);


	Resolution * BetaTOFResolution_P = new Resolution(finalHistos,"BetaTOFvsBeta Resolution (P)",ToFPB);
	Resolution * BetaTOFResolution_D = new Resolution(finalHistos,"BetaTOFvsBeta Resolution (D)",ToFDB);

	PlotFitResults(Plots,BetaTOFResolution_P,BetaTOFResolution_D,"#beta_{gen}",0.45,1,-0.05,0.1);

        Resolution * BetaNaFResolution_P = new Resolution(finalHistos,"BetaNaFvsBeta Resolution (P)",NaFPB);
	Resolution * BetaNaFResolution_D = new Resolution(finalHistos,"BetaNaFvsBeta Resolution (D)",NaFDB);

	PlotFitResults(Plots,BetaNaFResolution_P,BetaNaFResolution_D,"#beta_{gen}",0.7,1,-0.005,0.018);

        Resolution * BetaAglResolution_P = new Resolution(finalHistos,"BetaAglvsBeta Resolution (P)",AglPB);
	Resolution * BetaAglResolution_D = new Resolution(finalHistos,"BetaAglvsBeta Resolution (D)",AglDB);

	PlotFitResults(Plots,BetaAglResolution_P,BetaAglResolution_D,"#beta_{gen}",0.94,1,-0.0005,0.0018);

	Resolution * EdepUTOFResolution_P= new Resolution(finalHistos,"EdepUTOFvsBeta Resolution (P)",ToFPB);
	Resolution * EdepUTOFResolution_D= new Resolution(finalHistos,"EdepUTOFvsBeta Resolution (D)",ToFDB);

	PlotFitResults(Plots,EdepUTOFResolution_P,EdepUTOFResolution_D,"#beta_{gen}",0.45,1,0,0.7);

	Resolution * EdepLTOFResolution_P= new Resolution(finalHistos,"EdepLTOFvsBeta Resolution (P)",ToFPB);
	Resolution * EdepLTOFResolution_D= new Resolution(finalHistos,"EdepLTOFvsBeta Resolution (D)",ToFDB);

	PlotFitResults(Plots,EdepLTOFResolution_P,EdepLTOFResolution_D,"#beta_{gen}",0.45,1,0,0.7);
	
	Resolution * EdepTrackResolution_P= new Resolution(finalHistos,"EdepTrackvsBeta Resolution (P)",ToFPB);
	Resolution * EdepTrackResolution_D= new Resolution(finalHistos,"EdepTrackvsBeta Resolution (D)",ToFDB);

	PlotFitResults(Plots,EdepTrackResolution_P,EdepTrackResolution_D,"#beta_{gen}",0.45,1,0,10);
	


	
	return 0;
	
}






