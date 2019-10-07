#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include <TVector3.h>
#include "TMath.h"
#include <TFile.h>
#include <string>
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TChain.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TFractionFitter.h"
#include "TRandom3.h"
#include "fstream"
#include "TCanvas.h"

#include "InputFileReader.h"
#include "DBarReader.h"

#include "filesaver.h"
#include "binning.h"
#include "Globals.h"
#include "Variables.hpp"
#include "ParallelFiller.h"
#include "Cuts.h"

struct TFit {
   TH1F * Templ_P ;
   TH1F * Templ_D ;
   TH1F * Data;
   float wheightP,wheightD;
   TFractionFitter *Tfit;
   int Tfit_outcome;
};

void Do_TemplateFIT(TFit * Fit);

std::string GetOutFileName(std::string listname){
	std::ifstream infile(listname);	
	std::string line;
	std::getline(infile, line);
	line.erase(line.begin(),line.end()-15);
	
	std::cout<<"OUTPUT FILE: "<<line<<std::endl;
	
	return line;
};


float EvalChi2(TH1F * Data, TH1F * ModelP, TH1F* ModelD){
	int dof=0;
	float chi=0;
	TH1F * Model = (TH1F*) ModelP->Clone();
	Model->Add(ModelD);
	for(int i=0;i<Data->GetNbinsX();i++) {
		{
			if((pow(Data->GetBinError(i+1),2)+pow(Model->GetBinError(i+1),2))>0){
				chi+= pow(Data->GetBinContent(i+1) - Model->GetBinContent(i+1),2)/(pow(Data->GetBinError(i+1),2)+pow(Model->GetBinError(i+1),2));
				dof++;
			}
		}
	}
	return chi/dof;
}



float SmearBeta(float Beta, float Ri,float sigma,float tailfactor,float mean){

	float time = 1.2/(Beta*3e-4);

	float tailcontrolfactor=1;	//migration tail fixing
//	if(Beta>0.93) tailcontrolfactor = tailfactor;
	tailcontrolfactor = 1 + tailfactor*(atan(300*(Beta-mean))/3.1415926 + 0.5);
	float smeartime = (20 + Rand->Gaus(0,(float) tailcontrolfactor*sigma));
	time = time + smeartime;
	return 1.2/(time*3e-4);
}



int main(int argc, char * argv[])
{
	int comb=10;
	FRAC = 5;
	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	TH1::SetDefaultSumw2();     	
	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	string INPUT1(argv[1]);
	string INPUT2(argv[2]);
	string INPUT3(argv[3]);
	
	TChain * chain_RTI  	= InputFileReader(INPUT1.c_str(),"RTI");
	TChain * chainDT    	= InputFileReader(INPUT1.c_str(),"Event");
	TChain * chainMC    	= InputFileReader(INPUT2.c_str(),"Event");
	TChain * chainDT_Cpct    = InputFileReader(INPUT1.c_str(),"Compact");
	TChain * chainMC_Cpct    = InputFileReader(INPUT2.c_str(),"Compact");

	std::string outname = GetOutFileName(INPUT1.c_str());

	cout<<"****************************** BINS ***************************************"<<endl;
	SetUpUsualBinning();

	
	Variables * vars = new Variables();

	cout<<"****************************** ANALYIS ******************************************"<<endl;
	
	 FileSaver finalHistos;
    	 finalHistos.setName((INPUT3).c_str());

	 FileSaver finalResults;
    	 finalResults.setName((INPUT3 + "_Results").c_str());

	DBarReader readerMC(chainMC, true ,chain_RTI,chainMC_Cpct);
	DBarReader readerDT(chainDT, false,chain_RTI,chainDT_Cpct);

	TH1F * MassDistribution_data;
        TH1F * MassDistribution_orig;
        std::vector<std::vector<TH1F *>> MassDistribution_smea;
	TH1F * MassDistribution_origD;
        std::vector<std::vector<TH1F *>> MassDistribution_smeaD;


	
	
	if(finalHistos.CheckFile()){
		
		TFile * f = finalHistos.GetFile();
		MassDistribution_data  	= (TH1F*)f->Get("Example of Mass Distribution data");
		MassDistribution_orig	= (TH1F*)f->Get("Example of Mass Distribution Or");
		for(int i=0;i<comb;i++) {
			MassDistribution_smea.push_back(std::vector<TH1F *>());
			for(int j=0;j<comb;j++) MassDistribution_smea[i].push_back((TH1F*)f->Get(("P Templ./Example of Mass Distribution Smear_" + to_string(i) + to_string(j)).c_str()));
		}
		MassDistribution_origD	= (TH1F*)f->Get("Example of Mass Distribution Or D");
		for(int i=0;i<comb;i++) {
			MassDistribution_smeaD.push_back(std::vector<TH1F *>());
			for(int j=0;j<comb;j++) MassDistribution_smeaD[i].push_back((TH1F*)f->Get(("P Templ./Example of Mass Distribution SmearD_" + to_string(i) + to_string(j)).c_str()));
		}

		cout<<"*************file opening**************"<<endl;
		for(int j=0;j<comb;j++)for(int i=0;i<comb;i++) cout<<MassDistribution_smeaD[i][j]<<endl;
		cout<<"*************file opening**************"<<endl;
		for(int j=0;j<comb;j++)for(int i=0;i<comb;i++) cout<<MassDistribution_smea[i][j]<<endl;

		cout<<"*********** pre - scaling ************"<<endl;
		for(int j=0;j<comb;j++)for(int i=0;i<comb;i++) MassDistribution_smea[i][j]->Scale(MassDistribution_data->GetBinContent(MassDistribution_data->GetMaximumBin())/MassDistribution_smea[i][j]->GetBinContent(MassDistribution_smea[i][j]->GetMaximumBin()));
		for(int j=0;j<comb;j++)for(int i=0;i<comb;i++) 
		{
			TH1F * tmp = (TH1F *) MassDistribution_data->Clone();
			tmp->Add(MassDistribution_smea[i][j],-1);
			MassDistribution_smeaD[i][j]->Scale(tmp->GetBinContent(tmp->GetMaximumBin())/MassDistribution_smeaD[i][j]->GetBinContent(MassDistribution_smeaD[i][j]->GetMaximumBin()));
		}

		cout<<"************ fits *******************"<<endl;
		TFit * Fits[comb][comb];
		for(int i=0; i<comb; i++)
			for(int j=0; j<comb; j++)
				Fits[i][j]= new TFit;

		for(int i=0; i<comb; i++)
			for(int j=0; j<comb; j++){
				Fits[i][j]->Templ_P = (TH1F*) MassDistribution_smea[i][j]->Clone();
				Fits[i][j]->Templ_D = (TH1F*) MassDistribution_smeaD[i][j]->Clone();
				Fits[i][j]->Data    = (TH1F*) MassDistribution_data->Clone();
			//	Do_TemplateFIT(Fits[i][j]);
			}

		TH2F * ChiSquare = new TH2F("ChiSquare","ChiSquare",comb,0,comb,comb,0,comb);
		Int_t minx,miny,minz;
		for(int j=0;j<comb;j++) for(int i=0;i<comb;i++) ChiSquare->SetBinContent(i+1,j+1,EvalChi2(MassDistribution_data,Fits[i][j]->Templ_P,Fits[i][j]->Templ_D));
	
		ChiSquare->GetMinimumBin(minx,miny,minz);
		TH1F * BestP = (TH1F*) Fits[minx-1][miny-1]->Templ_P->Clone();
		TH1F * BestD = (TH1F*) Fits[minx-1][miny-1]->Templ_D->Clone();
		TH1F * Sum = (TH1F*) BestP->Clone();
		Sum->Add(BestD);
		Sum->SetName("best model");
		Sum->SetTitle("best model");


		BestP->SetLineColor(2);
		BestD->SetLineColor(4);
		BestP->SetLineWidth(2);
		BestD->SetLineColor(4);
		MassDistribution_data->SetLineColor(1);
		MassDistribution_data->SetLineWidth(3);

		finalResults.Add(ChiSquare);
		finalResults.Add(MassDistribution_data);
		finalResults.Add(MassDistribution_orig);
		finalResults.Add(MassDistribution_origD);
		finalResults.Add(BestP);
		finalResults.Add(BestD);
		finalResults.Add(Sum);

		finalResults.writeObjsInFolder("");	
		for(int j=0;j<comb;j++)for(int i=0;i<comb;i++) finalResults.Add(Fits[i][j]->Templ_P);
		finalResults.writeObjsInFolder("P Templ.");	
		for(int j=0;j<comb;j++)for(int i=0;i<comb;i++) finalResults.Add(Fits[i][j]->Templ_D);
		finalResults.writeObjsInFolder("D Templ.");	


	}
	
	/*
	else{
	MassDistribution_data = new TH1F("Example of Mass Distribution data","Example of Mass Distribution data",80,0,6);
	MassDistribution_orig = new TH1F("Example of Mass Distribution Or","Example of Mass Distribution Or",80,0,6);
	for(int j=0;j<comb;j++) for(int i=0;i<comb;i++) MassDistribution_smea[i][j] = new TH1F(("Example of Mass Distribution Smear_"+to_string(i)+ to_string(j)).c_str(),("Example of Mass Distribution Smear_"+to_string(i)+ to_string(j)).c_str(),80,0,6);
	MassDistribution_origD = new TH1F("Example of Mass Distribution Or D","Example of Mass Distribution Or D",80,0,6);
	for(int j=0;j<comb;j++) for(int i=0;i<comb;i++) MassDistribution_smeaD[i][j] = new TH1F(("Example of Mass Distribution SmearD_"+to_string(i)+ to_string(j)).c_str(),("Example of Mass Distribution SmearD_"+to_string(i)+ to_string(j)).c_str(),80,0,6);



	for(int i=0;i<readerMC.GetTreeEntries()/FRAC;i++){
		UpdateProgressBar(i, readerMC.GetTreeEntries()/FRAC);
		//if(i%(int)FRAC!=0) continue;
		readerMC.FillVariables(i,vars);
		vars->Update();
		if(ApplyCuts("IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&IsOnlyFromToF",vars)){
			float betasmear= SmearBeta(vars->Beta,vars->R,90,0,0.92);
			if(betasmear>0.83 && betasmear<0.85){
				MassDistribution_orig->Fill(vars->R/betasmear * pow((1-pow(betasmear,2)),0.5),vars->mcweight);	
				MassDistribution_origD->Fill((1.857/0.938)*vars->R/betasmear * pow((1-pow(betasmear,2)),0.5),vars->mcweight);	
				}
			for(int j=0;j<comb;j++) for(int i=0;i<comb;i++){
				betasmear = SmearBeta(vars->Beta,vars->R,90,0.9+0.1*i,0.87+0.01*j);
				if(betasmear>0.83 && betasmear<0.85){
					MassDistribution_smea[i][j]->Fill(vars->R/betasmear * pow((1-pow(betasmear,2)),0.5),vars->mcweight);	
					MassDistribution_smeaD[i][j]->Fill((1.857/0.938)*vars->R/betasmear * pow((1-pow(betasmear,2)),0.5),vars->mcweight);	
					}
				}
		}
	}

	for(int i=0;i<readerDT.GetCompactEntries()/FRAC;i++){
		UpdateProgressBar(i, readerDT.GetCompactEntries()/FRAC);
		//if(i%(int)FRAC!=0) continue;
		readerDT.FillCompact(i,vars);
		vars->Update();
		if(ApplyCuts("IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&IsOnlyFromToF",vars)){
			float betasmear=vars->Beta;
			if(betasmear>0.83 && betasmear<0.85)
				MassDistribution_data->Fill(vars->R/betasmear * pow((1-pow(betasmear,2)),0.5));	

		}
	}
	finalHistos.Add(MassDistribution_data);
	finalHistos.Add(MassDistribution_orig);
	finalHistos.Add(MassDistribution_origD);
	finalHistos.writeObjsInFolder("");	
	for(int j=0;j<comb;j++)for(int i=0;i<comb;i++) finalHistos.Add(MassDistribution_smea[i][j]);
	finalHistos.writeObjsInFolder("P Templ.");	
	for(int j=0;j<comb;j++)for(int i=0;i<comb;i++) finalHistos.Add(MassDistribution_smeaD[i][j]);
	finalHistos.writeObjsInFolder("P Templ.");	
	
	}

	*/
	return 0;
}


void Do_TemplateFIT(TFit * Fit){
         TObjArray *Tpl;
         Tpl = new TObjArray(2);
         Tpl -> Add( Fit ->  Templ_P );
         Tpl -> Add( Fit ->  Templ_D );
         Fit -> Tfit = new TFractionFitter(Fit -> Data, Tpl ,"q");
         float highP=1;
         float lowP=0.6;
         float highD=0.25;
         float lowD=0.00001;
         float min=0.5;
         float max=2.5;
         Fit -> Tfit -> SetRangeX(Fit -> Data -> FindBin(min), Fit -> Data -> FindBin(max));

         Fit -> Tfit_outcome = Fit -> Tfit -> Fit();

         if(Fit -> Tfit_outcome==0){
                TH1F * Result = (TH1F *) Fit-> Tfit -> GetPlot();
                float itot= Result->Integral();
                double w1,e1 = 0;
                double w2,e2 = 0;
                Fit -> Tfit ->GetResult(0,w1,e1);
                Fit -> Tfit ->GetResult(1,w2,e2);
                float i1 = Fit-> Templ_P  ->Integral(Fit->Templ_P -> FindBin(min), Fit->Templ_P -> FindBin(max));
                float i2 = Fit-> Templ_D  ->Integral(Fit->Templ_D -> FindBin(min), Fit->Templ_D -> FindBin(max));
                Fit ->wheightP= w1*itot/i1;
                Fit ->wheightD= w2*itot/i2;
                cout<<w1<<" "<<w2<<endl;

                Fit ->  Templ_P  -> Scale(Fit ->wheightP);
                Fit ->  Templ_D  -> Scale(Fit ->wheightD);
         }
        else{
                Fit ->wheightP= 0;
                Fit ->wheightD= 0;
        }


        return;
}

