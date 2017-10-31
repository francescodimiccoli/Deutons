#ifndef RESOLUTION_H
#define RESOLUTION_H

#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "../include/binning.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include <TVector3.h>
#include "TMath.h"
#include <TFile.h>
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include "TGraphErrors.h"
#include "TFractionFitter.h"
#include "TRandom3.h"
#include "../include/GlobalBinning.h"
#include "string.h"

#include "filesaver.h"
#include "FitError.h"
#include "Variables.hpp"
#include "Cuts.h"

using namespace std;

class Resolution{

	private:
		std::vector<TH1F*> Histos;
		Binning bins;
		std::string var;
		std::string cut;
		std::string discr_var;
		std::string basename;
		int nbins;
		float xmin,xmax;
		TH1F * Sigmas_Histo;
		TH1F * Means_Histo;
		TH1F * Resolutions_Histo;
		TObject * ModelSigmas;
		TObject * ModelMeans;
		double FitRangeEdges[2];
		bool fixedwindow = false;

	public:
		//standard constructor
		Resolution(std::string Basename,Binning Bins, std::string Cut, int Nbins, float Xmin, float Xmax,std::string Var="", std::string Discr_var=""){
			for(int i=0;i<Bins.size();i++){
				TH1F * Histo = new TH1F((Basename + to_string(i)).c_str(),(Basename + to_string(i)).c_str(),Nbins,Xmin,Xmax);
				Histos.push_back(Histo);	
			}
			basename=Basename;	
			var=Var;
			bins=Bins;
			xmin=Xmin;
			xmax=Xmax;
			nbins=Nbins;
			if(Cut.size()>0) cut=Cut; else cut="1>0";
			discr_var=Discr_var;	

			Sigmas_Histo      = new TH1F((Basename + "_sigmas").c_str(),(Basename + "_sigmas").c_str(),Bins.size(),0,Bins.size()); 		
			Means_Histo 	  = new TH1F((Basename + "_means" ).c_str(),(Basename + "_means" ).c_str(),Bins.size(),0,Bins.size());
			Resolutions_Histo = new TH1F((Basename + "_reso"  ).c_str(),(Basename + "_reso"  ).c_str(),Bins.size(),0,Bins.size());
		
		}

		//reading constructors
		Resolution(FileSaver File,std::string Basename, Binning Bins){
			std::string path = Basename + "/";
			for(int i=0;i<Bins.size();i++){
				TH1F * Histo = (TH1F*)  File.Get((path+Basename + to_string(i)).c_str());
				Histos.push_back(Histo);
				basename=Basename;
				bins=Bins;
			}
			nbins=Bins.size();

			Sigmas_Histo      = (TH1F*)  File.Get((path+ "Fit Results/" + Basename+ "_sigmas" ).c_str());
			Means_Histo       = (TH1F*)  File.Get((path+ "Fit Results/" + Basename+ "_means"  ).c_str());	
			Resolutions_Histo = (TH1F*)  File.Get((path+ "Fit Results/" + Basename+ "_reso"  ).c_str());	
			ModelSigmas 	  = (TObject *) File.Get((path+ "Fit Results/" + Basename+ "_ModelSigmas").c_str());			
			ModelMeans        = (TObject *) File.Get((path+ "Fit Results/" + Basename+ "_ModelMeans").c_str());	


			if(!Sigmas_Histo||!Means_Histo||!Resolutions_Histo){					
				Sigmas_Histo      = new TH1F((Basename + "_sigmas").c_str(),(Basename + "_sigmas").c_str(),Bins.size(),0,Bins.size());	
				Means_Histo 	  = new TH1F((Basename + "_means" ).c_str(),(Basename + "_means" ).c_str(),Bins.size(),0,Bins.size());
				Resolutions_Histo = new TH1F((Basename + "_reso"  ).c_str(),(Basename + "_reso"  ).c_str(),Bins.size(),0,Bins.size());
			}
		}

		void Fill(TTree * tree);	
		void Fill(TTree * tree, Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars));
		void FillEventByEvent(float var, float discr_var, bool cut);
		void Save(FileSaver finalhisto,bool recreate=false);
		void Normalize();
		void Eval_Resolution(std::vector<float> ExpectationValues={-1},float low_limit=0.05,float high_limit=0.85);
		void SetFixedFitWindow(float min, float max){FitRangeEdges[0]=min; FitRangeEdges[1]=max; fixedwindow=true; return;}

	bool CheckHistos();
	std::string GetName() {return basename;}
	Binning GetBinning(){return bins;}


	TH1F * Get_Means() {return Means_Histo;}
	TH1F * Get_Sigmas() {return Sigmas_Histo;}
	TH1F * Get_Resolutions() {return Resolutions_Histo;}

	TH1F * Get_Histo(int i) { return Histos[i];}
	TObject * Get_SigmasModel(){ return ModelSigmas;}
	TObject * Get_MeansModel(){ return ModelMeans;}



	TSpline3* ModelSigmasWithSpline()   {return ModelWithSpline(Sigmas_Histo,(basename+"_ModelSigmas").c_str(),bins);};	
	TF1* 	  ModelSigmasWithPoly()     {return ModelWithPoly(Sigmas_Histo,(basename+"_ModelSigmas").c_str(),bins);};
	
	TSpline3* ModelMeansWithSpline()    {return ModelWithSpline(Means_Histo,(basename+"_ModelMeans").c_str(),bins);};	
	TF1* 	  ModelMeansWithPoly()	    {return ModelWithPoly(Means_Histo,(basename+"_ModelMeans").c_str(),bins);};

	TF1 * ModelMeansRatio(Resolution * Reso1);

};


void Resolution::Fill(TTree * tree){

	cout<<basename.c_str()<<" Filling ..."<< endl;
	
	std::string histo = "htemp("+to_string(nbins)+","+to_string(xmin)+","+to_string(xmax)+")";
	cout<<histo.c_str()<<endl;

	for(int i=0;i<bins.size();i++){
		cout<<(" Filling Bin " + to_string(i)).c_str()<< endl;
		std::string discr_cut = discr_var + ">" + to_string(bins.GetBinLowEdge(i)) + "&&" + discr_var + "<=" + to_string(bins.GetBinLowEdge(i+1));   
		cout<<(var + ">>" + histo).c_str()<<" " <<(cut + "&&" + discr_cut).c_str()<<endl;
		tree->Draw((var + ">>" + histo).c_str(), (cut + "&&" + discr_cut).c_str() );
	
		Histos[i] = (TH1F*) gDirectory->Get("htemp");
		Histos[i]->SetName((basename + to_string(i)).c_str());
	}
	
	return;
}


void Resolution::Fill(TTree * tree, Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars) ){

	cout<<basename.c_str()<<" Filling ..."<< endl;
	vars->ReadBranches(tree);
	for(int i=0;i<tree->GetEntries();i++){

		UpdateProgressBar(i, tree->GetEntries());
		tree->GetEvent(i);
                vars->Update();
		FillEventByEvent( var(vars), discr_var(vars),ApplyCuts(cut,vars));
	}
	return;
}


void Resolution::FillEventByEvent(float var, float discr_var, bool CUT){
	
	int kbin;
	kbin = 	bins.GetBin(discr_var);
	if(CUT&&kbin>0){ 
		Histos[kbin]->Fill(var);	
	}
	return;	

}




void Resolution::Save(FileSaver finalhisto,bool recreate){

	for(int i=0;i<bins.size();i++) 
		finalhisto.Add(Histos[i]);
	finalhisto.writeObjsInFolder(basename.c_str(),recreate);
	return;
}



void Resolution::Normalize(){

	for(int i=0;i<Histos.size();i++){
			float integral=Histos[i]->Integral();
			if(integral>30){
			Histos[i]->Sumw2();
			Histos[i]->Scale(1/integral);
			}
			else
				for(int j=0;j<Histos[i]->GetNbinsX();j++){
					Histos[i]->SetBinContent(j+1,0);
					Histos[i]->SetBinError(j+1,0);

				}
		}
	return;

}


void Resolution::Eval_Resolution(std::vector<float> ExpectationValues,float low_limit,float high_limit){

	for(int i=0;i<Histos.size();i++){
		double quantiles[2] = {low_limit,high_limit};
		if(!fixedwindow) Histos[i]->GetQuantiles(2,FitRangeEdges,quantiles);

		TF1 * fitfunc = new TF1("fitfunc","gaus",-1,1);
		fitfunc->SetParameter(0,Histos[i]->GetBinContent(Histos[i]->GetMaximumBin()));
		fitfunc->SetParameter(1,Histos[i]->GetMean());
		fitfunc->SetParameter(2,Histos[i]->GetRMS());
		
		//scanning for best fit window
                float window=0.05;
		float step=window/2;
		float chi[5]={999999};
		cout<<"RANGE::: "<<quantiles[0]<<": "<<FitRangeEdges[0]<<", "<<quantiles[1]<<": "<<FitRangeEdges[1]<<endl;
		if(Histos[i]->Integral()>0){ 
			for(int j=0;j<5;j++){
				cout<<FitRangeEdges[0]-(window-j*step)*FitRangeEdges[0]<<" "<<FitRangeEdges[1]+(window-j*step)*FitRangeEdges[1]<<endl;
				Histos[i]->Fit("fitfunc","","",FitRangeEdges[0]-(window-j*step)*FitRangeEdges[0],FitRangeEdges[1]+(window-j*step)*FitRangeEdges[1]);
				chi[j]=fitfunc->GetChisquare();
				}
		}
		float j_minchisquare=0;
		float temp=chi[0];
		for(int j=0;j<5;j++){
			if(chi[j]<temp) {temp=chi[j]; j_minchisquare=(float)j;}
			cout<<j_minchisquare<<" "<<FitRangeEdges[0]-(window-j*step)*FitRangeEdges[0]<<" "<<FitRangeEdges[1]+(window-j*step)*FitRangeEdges[1]<<endl;
		}
		//final fit
		Histos[i]->Fit("fitfunc","","",FitRangeEdges[0]-(window-j_minchisquare*step)*FitRangeEdges[0],FitRangeEdges[1]+(window-j_minchisquare*step)*FitRangeEdges[1]);


		float mean;
		if(ExpectationValues[0]!=-1) mean = ExpectationValues[i]; 
		else mean  = bins.GetBinCenter(i);
		
		float meanshift   = fitfunc->GetParameter(1)/(1/mean);
		float meanshifterr= fitfunc->GetParError(1)/(1/mean);
		
		float sigma 	  = fitfunc->GetParameter(2);
		float sigmaerr    = fitfunc->GetParError(2);	

		if(Histos[i]->GetEntries()>0){	
			Sigmas_Histo->SetBinContent(i+1,sigma);
			Sigmas_Histo->SetBinError(i+1,sigmaerr);		

			Means_Histo->SetBinContent(i+1,meanshift);
			Means_Histo->SetBinError(i+1,meanshifterr);

			Resolutions_Histo->SetBinContent(i+1,sigma*mean);
			Resolutions_Histo->SetBinError(i+1,sigmaerr*mean);

		}
				
	}	
	return;
}



TF1 * Resolution::ModelMeansRatio(Resolution * Reso1){
	
	TH1F * Means1 = (TH1F *)Means_Histo->Clone();
	TH1F * Means2 = (TH1F *)Reso1->Means_Histo->Clone();

	Means1->Sumw2();
	Means2->Sumw2();
		
	Means2->Divide(Means1);
	

	TF1 * model = ModelWithPoly(Means2,(basename+"_ModelMeansRatio").c_str(),bins);
	return model;

}



bool Resolution::CheckHistos(){
	bool check=true;
	for(int i=0;i<bins.size();i++){
		if(!Histos[i])
			cout<<"**** BIN "<<i<<" : Histogram not found on File"<<endl;
		check=(check && Histos[i]);
	}		
	return check;
}

bool ReadCalibration(){ 
 
        cout<<"****************** CALIB. READING **************************"<<endl; 
        string nomecal=("/storage/gpfs_ams/ams/users/fdimicco/Deutons/DirectAnalysis/include/CalibTRD.root"); 
        FileSaver Calibration;   
        Calibration.setName(nomecal.c_str()); 
        bool checkfile = Calibration.CheckFile(); 
 
        if(checkfile) cout<<"calibration file found"<<endl; 
        else { cout<<"calibration file not found"<<endl; return false;} 
 
        Resolution * EdepTRDMC_P = new Resolution(Calibration,"EdepTRDvsBeta Measured MC",ToFResB); 
        Resolution * EdepTRDDT_P = new Resolution(Calibration,"EdepTRDvsBeta Measured DT",ToFResB); 
 
 
        EdepTRDbeta =(TSpline3 *) EdepTRDDT_P ->Get_MeansModel();  
        cout<<"TRD spline: "<<EdepTRDbeta<<endl; 
 
        return checkfile; 
 
} 
 


#endif
