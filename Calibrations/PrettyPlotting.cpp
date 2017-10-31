#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "../DirectAnalysis/include/binning.h"
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
#include "../DirectAnalysis/include/GlobalBinning.h"

#include "../Ntuple-making/Commonglobals.cpp"
#include "../DirectAnalysis/include/Variables.hpp"
#include "../DirectAnalysis/include/Cuts.h"
#include "../DirectAnalysis/include/filesaver.h"

#include "../DirectAnalysis/include/FitError.h"
#include "../DirectAnalysis/include/Resolution.h"
#include "../DirectAnalysis/include/PlottingFunctions.h"


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



void PlotFitResults(FileSaver Plots,Resolution * Reso1,Resolution * Reso2,std::string Xaxis, float xmin=-1,float xmax=-1,float ymin=-1,float ymax=-1,std::string entry1="protons MC",std::string entry2="deuterons MC"){

	TCanvas * c1 = new TCanvas(" Sigma");
	c1->SetCanvasSize(2000,1500);
	
	TPad * c1_up = new TPad("upperPad", "upperPad",0.0,0.3,1.0,1.0);
	c1_up->Draw();

	TPad * c1_do = new TPad("lowerPad", "lowerPad",0.0,0.0,1.0,0.3);
        c1_do->Draw();

	PlotTH1FintoGraph(c1_up,Reso1->GetBinning(), Reso1->Get_Sigmas(), Xaxis.c_str(),  "#sigma",2,false,"",xmin,xmax,ymin,ymax,entry1.c_str());
	PlotTH1FintoGraph(c1_up,Reso2->GetBinning(), Reso2->Get_Sigmas(), Xaxis.c_str(),  "#sigma",4,false,"ePsame",xmin,xmax,ymin,ymax,entry2.c_str());
	TObject * Model1 = Reso1->Get_SigmasModel();
	TObject * Model2 = Reso2->Get_SigmasModel();
	string tf1="TF1";
	
	if(Model1->ClassName()==tf1) {
		TF1 * model1 = (TF1*) Reso1->Get_SigmasModel();
		PlotFunction(c1_up,model1, Xaxis.c_str(),  "#sigma",2,"same",xmin,xmax,ymin,ymax,"Modelization");
	}
	else {
		TSpline3 * model1 = (TSpline3*) Reso1->Get_SigmasModel();
		PlotFunction(c1_up,model1, Xaxis.c_str(),  "#sigma",2,"same",xmin,xmax,ymin,ymax,"Modelization");
	}
	

	if(Model2->ClassName()==tf1) {
		TF1 * model2 = (TF1*) Reso2->Get_SigmasModel();
		PlotFunction(c1_up,model2, Xaxis.c_str(),  "#sigma",4,"same",xmin,xmax,ymin,ymax,"Modelization");
	}
	else {
		TSpline3 * model2 = (TSpline3*) Reso2->Get_SigmasModel();
		PlotFunction(c1_up,model2, Xaxis.c_str(),  "#sigma",4,"same",xmin,xmax,ymin,ymax,"Modelization");
	}

	PlotTH1FRatiointoGraph(c1_do,Reso1->GetBinning(), Reso1->Get_Sigmas(),Reso2->Get_Sigmas(), Xaxis.c_str(),  "#sigma",2,false,"",xmin,xmax,0.1,2.5,entry1.c_str(),entry2.c_str());


	TCanvas * c2 = new TCanvas(" Mean");
	c2->SetCanvasSize(2000,1500);
	
	TPad * c2_up = new TPad("upperPad", "upperPad",0.0,0.3,1.0,1.0);
	c2_up->Draw();

	TPad * c2_do = new TPad("lowerPad", "lowerPad",0.0,0.0,1.0,0.3);
        c2_do->Draw();


	PlotTH1FintoGraph(c2_up,Reso1->GetBinning(), Reso1->Get_Means(), Xaxis.c_str(), "Mean Shift (rec/gen)",2,false,"",xmin,xmax,ymin,ymax,entry1.c_str());
	PlotTH1FintoGraph(c2_up,Reso2->GetBinning(), Reso2->Get_Means(), Xaxis.c_str(), "Mean Shift (rec/gen)",4,false,"ePsame",xmin,xmax,ymin,ymax,entry2.c_str());
	
	Model1 = Reso1->Get_MeansModel();
	Model2 = Reso2->Get_MeansModel();
	
	if(Model1->ClassName()==tf1) {
		TF1 * model1 = (TF1*) Reso1->Get_MeansModel();
		PlotFunction(c2_up,model1, Xaxis.c_str(),  "#sigma",2,"same",xmin,xmax,ymin,ymax,"Modelization");
	}
	else {
		TSpline3 * model1 = (TSpline3*) Reso1->Get_MeansModel();
		PlotFunction(c2_up,model1, Xaxis.c_str(),  "#sigma",2,"same",xmin,xmax,ymin,ymax,"Modelization");
	}
	

	if(Model2->ClassName()==tf1) {
		TF1 * model2 = (TF1*) Reso2->Get_MeansModel();
		PlotFunction(c2_up,model2, Xaxis.c_str(),  "#sigma",4,"same",xmin,xmax,ymin,ymax,"Modelization");
	}
	else {
		TSpline3 * model2 = (TSpline3*) Reso2->Get_MeansModel();
		PlotFunction(c2_up,model2, Xaxis.c_str(),  "#sigma",4,"same",xmin,xmax,ymin,ymax,"Modelization");
	}

	PlotTH1FRatiointoGraph(c2_do,Reso1->GetBinning(), Reso1->Get_Means(),Reso2->Get_Means(), Xaxis.c_str(),  "Means",2,false,"",xmin,xmax,0.2,2.5,entry1.c_str(),entry2.c_str());
        PlotFunction(c2_do,Reso1->ModelMeansRatio(Reso2), Xaxis.c_str(),  "#sigma (ratio)",4,"same",xmin,xmax,ymin,ymax,"Calibration");


        TCanvas * c3 = new TCanvas(" Resol.");
        c3->cd();
        PlotTH1FintoGraph(gPad,Reso1->GetBinning(), Reso1->Get_Resolutions(), Xaxis.c_str(), "Resolution (#sigma/mean)",2,false,"",xmin,xmax,ymin,ymax,entry1.c_str());
	PlotTH1FintoGraph(gPad,Reso2->GetBinning(), Reso2->Get_Resolutions(), Xaxis.c_str(), "Resolution (#sigma/mean)",4,false,"ePsame",xmin,xmax,ymin,ymax,entry2.c_str());



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

	PResB.Print();

	cout<<"**TOF**"<<endl;
	ToFResB.Print();

	cout<<"**NaF**"<<endl;
	NaFResB.Print();

	cout<<"**Agl**"<<endl;
	AglResB.Print();

	ToFResB.UseBetaEdges();
	NaFResB.UseBetaEdges();
	AglResB.UseBetaEdges();

	PResB.UseREdges();
	
	cout<<endl;
	
	cout<<"****************************** PLOTTING ***************************************"<<endl;

	/*
	Resolution * RigidityResolution_P = new Resolution(finalHistos,"RvsR Resolution (P)",PResB);
        Resolution * RigidityResolution_D = new Resolution(finalHistos,"RvsR Resolution (D)",PResB);

        PlotDistrib(Plots,RigidityResolution_P,"1/R_{gen}-1/R_{meas}","Normalized Distr.",-0.5,1.5,1e-4,1);
        PlotFitResults(Plots,RigidityResolution_P,RigidityResolution_D,"R_{gen} [GV]", 0.4,100,-0.1,1);

	

	Resolution * BetaTOFResolution_P = new Resolution(finalHistos,"BetaTOFvsBeta Resolution (P)",ToFResB);
        Resolution * BetaTOFResolution_D = new Resolution(finalHistos,"BetaTOFvsBeta Resolution (D)",ToFResB);

        PlotFitResults(Plots,BetaTOFResolution_P,BetaTOFResolution_D,"#beta_{gen}",0.45,1,-0.05,0.5);

	

	Resolution * BetaNaFResolution_P = new Resolution(finalHistos,"BetaNaFvsBeta Resolution (P)",NaFResB);
        Resolution * BetaNaFResolution_D = new Resolution(finalHistos,"BetaNaFvsBeta Resolution (D)",NaFResB);

        PlotFitResults(Plots,BetaNaFResolution_P,BetaNaFResolution_D,"#beta_{gen}",0.7,1,-0.005,0.018);

	

	Resolution * BetaAglResolution_P = new Resolution(finalHistos,"BetaAglvsBeta Resolution (P)",AglResB);
        Resolution * BetaAglResolution_D = new Resolution(finalHistos,"BetaAglvsBeta Resolution (D)",AglResB);

        PlotFitResults(Plots,BetaAglResolution_P,BetaAglResolution_D,"#beta_{gen}",0.94,1,-0.0005,0.0018);

	

	Resolution * EdepUTOFResolution_P= new Resolution(finalHistos,"EdepUTOFvsBeta Resolution (P)",ToFResB);
        Resolution * EdepUTOFResolution_D= new Resolution(finalHistos,"EdepUTOFvsBeta Resolution (D)",ToFResB);

        PlotFitResults(Plots,EdepUTOFResolution_P,EdepUTOFResolution_D,"#beta_{gen}",0.45,1,0,0.7);

        

	Resolution * EdepLTOFResolution_P= new Resolution(finalHistos,"EdepLTOFvsBeta Resolution (P)",ToFResB);
        Resolution * EdepLTOFResolution_D= new Resolution(finalHistos,"EdepLTOFvsBeta Resolution (D)",ToFResB);

        PlotFitResults(Plots,EdepLTOFResolution_P,EdepLTOFResolution_D,"#beta_{gen}",0.45,1,0,0.7);
        
        

	Resolution * EdepTrackResolution_P= new Resolution(finalHistos,"EdepTrackvsBeta Resolution (P)",ToFResB);
        Resolution * EdepTrackResolution_D= new Resolution(finalHistos,"EdepTrackvsBeta Resolution (D)",ToFResB);

        PlotFitResults(Plots,EdepTrackResolution_P,EdepTrackResolution_D,"#beta_{gen}",0.45,1,0,13);

	

	Resolution * EdepUTOFMC_P = new Resolution(finalHistos,"EdepUTOFvsBeta Measured MC",ToFResB);
	Resolution * EdepUTOFDT_P = new Resolution(finalHistos,"EdepUTOFvsBeta Measured DT",ToFResB);

	PlotFitResults(Plots,EdepUTOFMC_P,EdepUTOFDT_P,"#beta_{TOF}",0.45,1,0,0.7,"protons MC","ISS data");

	

	Resolution * EdepLTOFMC_P = new Resolution(finalHistos,"EdepLTOFvsBeta Measured MC",ToFResB);
	Resolution * EdepLTOFDT_P = new Resolution(finalHistos,"EdepLTOFvsBeta Measured DT",ToFResB);

	PlotFitResults(Plots,EdepLTOFMC_P,EdepLTOFDT_P,"#beta_{TOF}",0.45,1,0,0.7,"protons MC","ISS data");

	

	Resolution * EdepTrackMC_P = new Resolution(finalHistos,"EdepTrackvsBeta Measured MC",ToFResB);
	Resolution * EdepTrackDT_P = new Resolution(finalHistos,"EdepTrackvsBeta Measured DT",ToFResB);

	PlotFitResults(Plots,EdepTrackMC_P,EdepTrackDT_P,"#beta_{TOF}",0.45,1,0,13,"protons MC","ISS data");

	*/

	Resolution * EdepTRDMC_P = new Resolution(finalHistos,"EdepTRDvsBeta Measured MC",ToFResB);
	Resolution * EdepTRDDT_P = new Resolution(finalHistos,"EdepTRDvsBeta Measured DT",ToFResB);

	PlotFitResults(Plots,EdepTRDMC_P,EdepTRDDT_P,"#beta_{TOF}",0.45,1,0,0.1,"protons MC","ISS data");



/*	
	Resolution * QUTOFMC_P = new Resolution(finalHistos,"QUTOFvsBeta Measured MC",ToFResB);
        Resolution * QUTOFDT_P = new Resolution(finalHistos,"QUTOFvsBeta Measured DT",ToFResB);

        PlotFitResults(Plots,QUTOFMC_P,QUTOFDT_P,"#beta_{TOF}",0.45,1,0,2,"protons MC","ISS data");

        

        Resolution * QLTOFMC_P = new Resolution(finalHistos,"QLTOFvsBeta Measured MC",ToFResB);
        Resolution * QLTOFDT_P = new Resolution(finalHistos,"QLTOFvsBeta Measured DT",ToFResB);

        PlotFitResults(Plots,QLTOFMC_P,QLTOFDT_P,"#beta_{TOF}",0.45,1,0,2,"protons MC","ISS data");

        

        Resolution * QInnerMC_P = new Resolution(finalHistos,"QInnervsBeta Measured MC",ToFResB);
        Resolution * QInnerDT_P = new Resolution(finalHistos,"QInnervsBeta Measured DT",ToFResB);

        PlotFitResults(Plots,QInnerMC_P,QInnerDT_P,"#beta_{TOF}",0.45,1,0,2,"protons MC","ISS data");

*/	
	return 0;
	
}






