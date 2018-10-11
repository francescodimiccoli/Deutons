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
#include "TRandom3.h"
#include "../include/Globals.h"
#include "TKey.h"
#include "TFractionFitter.h"
#include "TText.h"
#include "TPaveLabel.h"

#include "../include/Variables.hpp"
#include "../include/Cuts.h"
#include "../include/filesaver.h"
#include "../include/TemplateFITbetasmear.h"

#include "../include/FitError.h"
#include "../include/PlottingFunctions.h"

std::string Convert (float number){
    std::ostringstream buff;
    buff<<number;
    std::string output= buff.str();
    output.erase(4,output.end()-output.begin()-4);	
    return output;   
}


int colorbase = 55;
float minscan =1.65;
std::vector<TH1F*> GetListOfTemplates(TFile * file,std::string path){

        std::vector<TH1F*> Templates;

        if(file->GetDirectory(path.c_str())){
                TList * list = file->GetDirectory(path.c_str())->GetListOfKeys();


                TIter next(list);
                TKey * key;
                while((key = (TKey*)next())){
                        Templates.push_back((TH1F *)file->Get((path + "/" + key->GetName()).c_str()));
                }
        }
        return Templates;


}

void DrawFits(TFile * file, std::string basename, Binning Bins, FileSaver Plots);
void DrawBranching(TFile * file, std::string basename, Binning Bins, FileSaver Plots);
void DrawFragments(TFile * file, std::string basename, Binning Bins, FileSaver Plots);
void DrawFragmentsCheck(TFile * file, std::string basename, Binning Bins, FileSaver Plots);


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

	SetUpUsualBinning();

	cout<<"****************************** PLOTTING ***************************************"<<endl;	 

	DrawFits(finalHistos.GetFile(),"HeContTOF",ToFDB,Plots);
	DrawFits(finalHistos.GetFile(),"HeContNaF",NaFDB,Plots);
	DrawFits(finalHistos.GetFile(),"HeContAgl",AglDB,Plots);

	DrawBranching(finalHistos.GetFile(),"HeContCheck9",ToFDB,Plots);
	DrawBranching(finalHistos.GetFile(),"HeContTOF",ToFDB,Plots);
	DrawBranching(finalHistos.GetFile(),"HeContNaF",NaFDB,Plots);
	DrawBranching(finalHistos.GetFile(),"HeContAgl",AglDB,Plots);

/*        for(int i=0;i<10;i++){
                 std::string basename = "HeContCheck" + to_string(i);
		 DrawBranching(finalHistos.GetFile(),basename.c_str(),ToFDB,Plots);
	} 
   
*/
	DrawFragments(finalHistos.GetFile(),"HeContTOF",ToFDB,Plots);
	DrawFragments(finalHistos.GetFile(),"HeContNaF",NaFDB,Plots);
	DrawFragments(finalHistos.GetFile(),"HeContAgl",AglDB,Plots);

	DrawFragmentsCheck(finalHistos.GetFile(),"HeContCheck",ToFDB,Plots);

}


void DrawFragmentsCheck(TFile * file, std::string basename, Binning Bins, FileSaver Plots){
		std::string pathdatacheck[10];
		std::vector<std::vector<TH1F*>> Data;
		TH1F * SumData[10];
		TH1F * SumRatio[10];

		TCanvas * c3 = new TCanvas("Fragments Test");
                c3->SetCanvasSize(2000,1500);
		TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.35,1.0,1.0);	
		pad1->Draw();
		TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.35);
		pad2->Draw();
		for(int i=0;i<10;i++) { 
			pathdatacheck[i]  = (basename + to_string(i)+"/Fragments/");
			Data.push_back(GetListOfTemplates(file, pathdatacheck[i]));	
			cout<<"check: "<<Data[i].size()<<endl;
			for(int bin =0; bin < Bins.size();bin++){
				Data[i][bin]->Rebin(3);
				if(bin==0) SumData[i] = (TH1F *)Data[i][0]->Clone();
				else SumData[i]->Add(Data[i][bin]);
			}
		}
		float P[10];
		float D[10];
		float T[10];

		for(int i=0;i<10;i++){
			P[i]=SumData[i]->Integral(SumData[i]->FindBin(0.5),SumData[i]->FindBin(1.2));
			D[i]=SumData[i]->Integral(SumData[i]->FindBin(1.6),SumData[i]->FindBin(2.2));
			T[i]=SumData[i]->Integral(SumData[i]->FindBin(2.7),SumData[i]->FindBin(3.7));
		}


		for(int i=0;i<10;i++) {
			SumData[i]->Scale(1/SumData[i]->Integral());
		}
		for(int i=0;i<10;i++) {
			SumRatio[i]=(TH1F*)SumData[i]->Clone();
			SumRatio[i]->Divide(SumData[9]);
		}
		for(int i=0;i<10;i++){ 
			PlotDistribution(pad1, SumData[i] ,"Reconstructed Mass [GeV/c^{2}]","Distribution Width",i,"ePsame",-SumData[i]->GetBinContent(SumData[i]->GetMaximumBin())*2.33,SumData[i]->GetBinContent(SumData[i]->GetMaximumBin())*2.33,3,("Cut: qL1>"+to_string(minscan+i*0.05)));
			PlotDistribution(pad2, SumRatio[i] ,"Reconstructed Mass [GeV/c^{2}]","Ratio",i,"ePsame",0.,2.5,2,("Cut: qL1>"+to_string(minscan+i*0.05)));
		}
	
		TCanvas * c4 = new TCanvas("Fragments Test Ratios");
                c3->SetCanvasSize(2000,1500);
		TGraphErrors * TOverP = new TGraphErrors();
		TGraphErrors * DOverP = new TGraphErrors();
		TGraphErrors * DOverT = new TGraphErrors();
		TGraphErrors * Palone = new TGraphErrors();
		TGraphErrors * Dalone = new TGraphErrors();
		TGraphErrors * Talone = new TGraphErrors();



		for(int i=0;i<10;i++) { TOverP->SetPoint(i,(minscan+i*0.05),(T[i]/P[i])/(T[0]/P[0]));
				       TOverP->SetPointError(i,0,pow(1/P[i]+1/T[i],0.5)*(T[i]/P[i])/(T[0]/P[0]));
				     } 
		for(int i=0;i<10;i++) { DOverP->SetPoint(i,(minscan+i*0.05),(D[i]/P[i])/(D[0]/P[0]));
				       DOverP->SetPointError(i,0,pow(1/P[i]+1/D[i],0.5)*(D[i]/P[i])/(D[0]/P[0]));
				     } 
		for(int i=0;i<10;i++) { DOverT->SetPoint(i,(minscan+i*0.05),(D[i]/T[i])/(D[0]/T[0]));
				       DOverT->SetPointError(i,0,pow(1/T[i]+1/D[i],0.5)*(D[i]/T[i])/(D[0]/T[0]));
				     } 
			for(int i=0;i<10;i++) { Palone->SetPoint(i,(minscan+i*0.05),(P[i]/T[i]));
				       Palone->SetPointError(i,0,pow(1/P[i],0.5)*P[i]/T[i]);
				     } 
		for(int i=0;i<10;i++) { Dalone->SetPoint(i,(minscan+i*0.05),(D[i]/T[i]));
				       Dalone->SetPointError(i,0,pow(1/D[i],0.5)*D[i]/T[i]);
				     } 
		for(int i=0;i<10;i++) { Talone->SetPoint(i,(minscan+i*0.05),(T[i]/T[0]));
				       Talone->SetPointError(i,0,pow(1/T[i],0.5)*T[i]/T[0]);
				     } 

		PlotGraph(gPad,TOverP,"Cut Value","Fragments Ratio",2,"PL",1.5,2.3,0.5,2.5,"T over P");;
		PlotGraph(gPad,DOverP,"Cut Value","Fragments Ratio",1,"PLsame",1.5,2.3,0.5,2.5,"D over P");;
		PlotGraph(gPad,DOverT,"Cut Value","Fragments Ratio",4,"PLsame",1.5,2.3,0.5,125,"D over T");;
	
		TCanvas * c5 = new TCanvas("Tritium over L1 Cut");
                c5->SetCanvasSize(2000,1500);
		
		PlotGraph(gPad,Talone,"Cut Value","Counts (over first)",3,"PLsame",1.5,2.3,1e-6,1.5,"He->T");;	

		TCanvas * c6 = new TCanvas("Fragments Ratio Evolution");
                c6->SetCanvasSize(2000,1500);
		PlotGraph(gPad,Palone,"Cut Value","Counts (over tritium)",2,"PLsame",1.5,2.3,1e-6,1e2,"He->P/He->T");;	
		PlotGraph(gPad,Dalone,"Cut Value","Counts (over tritium)",4,"PLsame",1.5,2.3,1e-6,1e2,"He->D/He->T");;	

		Plots.Add(c3);
		Plots.Add(c4);
		Plots.Add(c5);
               	Plots.Add(c6);
		Plots.writeObjsInFolder("FragmentationCheck");
}

void DrawBranching(TFile * file, std::string basename, Binning Bins, FileSaver Plots){
	std::string pathdata    = (basename + "/Fragments/");
        std::string pathmc      = (basename + "/FragmentsMC/");
	std::string pathq2data  = (basename + "/Q2Templates/");
        std::string pathmassPmc    = (basename + "/PMassMC/");
 	std::string pathmassDmc    = (basename + "/DMassMC/");
	std::string pathmassTmc    = (basename + "/TMassMC/");
	std::string pathq2mc    = (basename + "/Q2TemplatesMC/");


	std::vector<TH1F*> Data=GetListOfTemplates(file, pathdata);
	std::vector<TH1F*> MC=GetListOfTemplates(file, pathmc);			
	std::vector<TH1F*> Q2Data=GetListOfTemplates(file, pathq2data);
	std::vector<TH1F*> Q2MC=GetListOfTemplates(file, pathq2mc);			
	std::vector<TH1F*> PMassMC=GetListOfTemplates(file, pathmassPmc);			
	std::vector<TH1F*> DMassMC=GetListOfTemplates(file, pathmassDmc);			
	std::vector<TH1F*> TMassMC=GetListOfTemplates(file, pathmassTmc);			


	TH1F * SumData ;
	TH1F * SumMC; 
	TH1F * SumQ2Data ;
	TH1F * SumQ2MC; 
	TH1F * SumPMassMC; 
	TH1F * SumDMassMC; 
	TH1F * SumTMassMC; 


	for(int bin =2; bin < Bins.size();bin++){

		Data[bin]->Rebin(3);
		MC[bin]->Rebin(3);
		Q2Data[bin]->Rebin(3);
		Q2MC[bin]->Rebin(3);
		PMassMC[bin]->Rebin(3);
		DMassMC[bin]->Rebin(3);
		TMassMC[bin]->Rebin(3);
	
		if(bin==2){
		 SumData = (TH1F *)Data[2]->Clone();
		 SumMC   = (TH1F *)MC[2]   ->Clone();
		 SumQ2Data = (TH1F *)Q2Data[2]->Clone();
		 SumQ2MC   = (TH1F *)Q2MC[2]   ->Clone();
		 SumPMassMC   = (TH1F *)PMassMC[2]   ->Clone();
		 SumDMassMC   = (TH1F *)DMassMC[2]   ->Clone();
		 SumTMassMC   = (TH1F *)TMassMC[2]   ->Clone();
		
		}
		else {
			SumData->Add(Data[bin]);
			SumMC->Add(MC[bin]);
			SumQ2Data->Add(Q2Data[bin]);
			SumQ2MC->Add(Q2MC[bin]);
			SumPMassMC ->Add(PMassMC[bin]);
                        SumDMassMC ->Add(DMassMC[bin]);
			SumTMassMC ->Add(TMassMC[bin]);
		}

	}

	TCanvas * c4 = new TCanvas("Branching Ratio DT");
        c4->SetCanvasSize(2000,1500);
	c4->cd();
	SumData->Scale(1/SumData->Integral());
	
	float normP = (SumData->GetBinContent(SumData->FindBin(1.1))+SumData->GetBinContent(SumData->FindBin(1.1)+1)+SumData->GetBinContent(SumData->FindBin(1.1)-1) +SumData->GetBinContent(SumData->FindBin(1.1)+2)+SumData->GetBinContent(SumData->FindBin(1.1)-2))/5;
	float normD = (SumData->GetBinContent(SumData->FindBin(1.9))+SumData->GetBinContent(SumData->FindBin(1.9)+1)+SumData->GetBinContent(SumData->FindBin(1.9)-1))/3;
	float normT = (SumData->GetBinContent(SumData->FindBin(3))+SumData->GetBinContent(SumData->FindBin(3)+1)+SumData->GetBinContent(SumData->FindBin(3)-1))/3;

	SumPMassMC ->Scale(normP/ SumPMassMC -> GetBinContent(SumPMassMC ->FindBin(1.1)));
	SumDMassMC ->Scale(normD/ SumDMassMC -> GetBinContent(SumDMassMC ->FindBin(1.9)));
	SumTMassMC ->Scale(normT/ SumTMassMC -> GetBinContent(SumTMassMC ->FindBin(3)));

	PlotDistribution(gPad, SumData ,"Reconstructed Mass [GeV/c^{2}]","Distribution Width",1,"ePsame",-SumData->GetBinContent(SumData->GetMaximumBin())*2.33,SumData->GetBinContent(SumData->GetMaximumBin())*2.33,3,"Fragments");   
	PlotDistribution(gPad, SumPMassMC ,"Reconstructed Mass [GeV/c^{2}]","Distribution Width",2,"same",-SumData->GetBinContent(SumData->GetMaximumBin())*2.33,SumData->GetBinContent(SumData->GetMaximumBin())*2.33,3,"P mass Template");   
	PlotDistribution(gPad, SumDMassMC ,"Reconstructed Mass [GeV/c^{2}]","Distribution Width",4,"same",-SumData->GetBinContent(SumData->GetMaximumBin())*2.33,SumData->GetBinContent(SumData->GetMaximumBin())*2.33,3,"D mass Template");   
	PlotDistribution(gPad, SumTMassMC ,"Reconstructed Mass [GeV/c^{2}]","Distribution Width",3,"same",-SumData->GetBinContent(SumData->GetMaximumBin())*2.33,SumData->GetBinContent(SumData->GetMaximumBin())*2.33,3,"T mass Template");   
		
	std::string Pfract = "He->P: " + Convert(SumPMassMC->Integral()/SumData->Integral() *100) + "%"; 
	std::string Dfract = "He->D: " + Convert(SumDMassMC->Integral()/SumData->Integral() *100) + "%"; 
	std::string Tfract = "He->T: " + Convert(SumTMassMC->Integral()/SumData->Integral() *100) + "%"; 

	{
	TPad *newpad=new TPad("newpad","a transparent pad",0.5,0.2,0.8,0.6);
	newpad->SetFillStyle(4000);
	newpad->Draw();
	newpad->cd();
	
	TText *t = new TText(.1,.7,Pfract.c_str());
	t->SetTextAlign(11);
	t->SetTextColor(kRed+1);
	t->SetTextFont(43);
	t->SetTextSize(40);
	t->Draw();
	TText *f = new TText(.1,.5,Dfract.c_str());
	f->SetTextAlign(11);
	f->SetTextColor(kBlue+1);
	f->SetTextFont(43);
	f->SetTextSize(40);
	f->Draw("same");
	TText *g = new TText(.1,.3,Tfract.c_str());
	g->SetTextAlign(11);
	g->SetTextColor(kGreen+1);
	g->SetTextFont(43);
	g->SetTextSize(40);
	g->Draw("same");
	}

	Plots.Add(c4);	
	Plots.writeObjsInFolder((basename + "/BranchingRatio").c_str());

	TCanvas * c5 = new TCanvas("Branching Ratio MC");
        c5->SetCanvasSize(2000,1500);
	c5->cd();
	SumMC->Scale(1/SumMC->Integral());
	
	normP = (SumMC->GetBinContent(SumMC->FindBin(0.93))+SumMC->GetBinContent(SumMC->FindBin(0.93)+1)+SumMC->GetBinContent(SumMC->FindBin(0.93)-1) +SumMC->GetBinContent(SumMC->FindBin(0.93)+2)+SumMC->GetBinContent(SumMC->FindBin(0.93)-2))/5;
	normD = (SumMC->GetBinContent(SumMC->FindBin(1.9))+SumMC->GetBinContent(SumMC->FindBin(1.9)+1)+SumMC->GetBinContent(SumMC->FindBin(1.9)-1))/3;
	normT = (SumMC->GetBinContent(SumMC->FindBin(3))+SumMC->GetBinContent(SumMC->FindBin(3)+1)+SumMC->GetBinContent(SumMC->FindBin(3)-1))/3;

	SumPMassMC ->Scale(normP/ SumPMassMC -> GetBinContent(SumPMassMC ->FindBin(0.93)));
	SumDMassMC ->Scale(normD/ SumDMassMC -> GetBinContent(SumDMassMC ->FindBin(1.9)));
	SumTMassMC ->Scale(normT/ SumTMassMC -> GetBinContent(SumTMassMC ->FindBin(3)));

	PlotDistribution(gPad, SumMC ,"Reconstructed Mass [GeV/c^{2}]","Distribution Width",3,"same",-SumMC->GetBinContent(SumMC->GetMaximumBin())*2.33,SumMC->GetBinContent(SumMC->GetMaximumBin())*2.33,3,"He MC Fragments");   
	PlotDistribution(gPad, SumPMassMC ,"Reconstructed Mass [GeV/c^{2}]","Distribution Width",2,"same",-SumMC->GetBinContent(SumMC->GetMaximumBin())*2.33,SumMC->GetBinContent(SumMC->GetMaximumBin())*2.33,3,"P mass Template");   
	PlotDistribution(gPad, SumDMassMC ,"Reconstructed Mass [GeV/c^{2}]","Distribution Width",4,"same",-SumMC->GetBinContent(SumMC->GetMaximumBin())*2.33,SumMC->GetBinContent(SumMC->GetMaximumBin())*2.33,3,"D mass Template");   
	PlotDistribution(gPad, SumTMassMC ,"Reconstructed Mass [GeV/c^{2}]","Distribution Width",3,"same",-SumMC->GetBinContent(SumMC->GetMaximumBin())*2.33,SumMC->GetBinContent(SumMC->GetMaximumBin())*2.33,3,"T mass Template");   
		
	Pfract = "He->P: " + Convert(SumPMassMC->Integral()/SumMC->Integral() *100) + "%"; 
	Dfract = "He->D: " + Convert(SumDMassMC->Integral()/SumMC->Integral() *100) + "%"; 
	Tfract = "He->T: " + Convert(SumTMassMC->Integral()/SumMC->Integral() *100) + "%"; 

	{
	TPad *newpad=new TPad("newpad","a transparent pad",0.5,0.2,0.8,0.6);
	newpad->SetFillStyle(4000);
	newpad->Draw();
	newpad->cd();
	
	
	TText *t = new TText(.1,.7,Pfract.c_str());
	t->SetTextAlign(11);
	t->SetTextColor(kRed+1);
	t->SetTextFont(43);
	t->SetTextSize(40);
	t->Draw();
	TText *f = new TText(.1,.5,Dfract.c_str());
	f->SetTextAlign(11);
	f->SetTextColor(kBlue+1);
	f->SetTextFont(43);
	f->SetTextSize(40);
	f->Draw("same");
	TText *g = new TText(.1,.3,Tfract.c_str());
	g->SetTextAlign(11);
	g->SetTextColor(kGreen+1);
	g->SetTextFont(43);
	g->SetTextSize(40);
	g->Draw("same");
	}

	Plots.Add(c5);	
	Plots.writeObjsInFolder((basename + "/BranchingRatio").c_str());


}


void DrawFragments(TFile * file, std::string basename, Binning Bins, FileSaver Plots){
	std::string pathdata    = (basename + "/Fragments/");
        std::string pathmc      = (basename + "/FragmentsMC/");
	std::string pathq2data  = (basename + "/Q2Templates/");
        std::string pathq2mc    = (basename + "/Q2TemplatesMC/");


	std::vector<TH1F*> Data=GetListOfTemplates(file, pathdata);
	std::vector<TH1F*> MC=GetListOfTemplates(file, pathmc);			
	std::vector<TH1F*> Q2Data=GetListOfTemplates(file, pathq2data);
	std::vector<TH1F*> Q2MC=GetListOfTemplates(file, pathq2mc);			


	TH1F * SumData ;
	TH1F * SumMC; 
	TH1F * SumQ2Data ;
	TH1F * SumQ2MC; 


	for(int bin =2; bin < Bins.size();bin++){
		TCanvas * c3 = new TCanvas(("Fragments Bin "+to_string(bin)).c_str());
                c3->SetCanvasSize(2000,1500);
		gPad->SetLogy();

		Data[bin]->Rebin(3);
		MC[bin]->Rebin(3);
		Q2Data[bin]->Rebin(3);
		Q2MC[bin]->Rebin(3);
	

		//if(Data[bin]->Integral()>0) Data[bin]->Scale(1/Data[bin]->Integral());
		//if(MC[bin]->Integral()>0)   MC[bin]->Scale(1/MC[bin]->Integral());
	
		if(bin==2){
		 SumData = (TH1F *)Data[2]->Clone();
		 SumMC   = (TH1F *)MC[2]   ->Clone();
		 SumQ2Data = (TH1F *)Q2Data[2]->Clone();
		 SumQ2MC   = (TH1F *)Q2MC[2]   ->Clone();
		}
		else {
			SumData->Add(Data[bin]);
			SumMC->Add(MC[bin]);
			SumQ2Data->Add(Q2Data[bin]);
			SumQ2MC->Add(Q2MC[bin]);
		}	

		Data[bin]->Scale(1/Q2Data[bin]->Integral());
		MC[bin]->Scale(1/Q2MC[bin]->Integral());

		PlotDistribution(gPad, Data[bin] ,"Reconstructed Mass [GeV/c^{2}]","Distribution Width",1,"ePsame",1e-7,Data[bin]->GetBinContent(Data[bin]->GetMaximumBin())*2.33,3,"Data-Driven");	
		PlotDistribution(gPad, MC[bin] ,"Reconstructed Mass [GeV/c^{2}]","Distribution Width",3,"same",1e-7,Data[bin]->GetBinContent(Data[bin]->GetMaximumBin())*2.33,6,"Helium MC");	
		
		Plots.Add(c3);	
	}

	TCanvas * c4 = new TCanvas("Fragments Tot");
        c4->SetCanvasSize(2000,1500);
	c4->cd();
	gPad->SetLogy();

	SumData->Scale(1/SumQ2Data->Integral());
	SumMC->Scale(1/(1.135*SumQ2MC->Integral())); //1.135 He3/He4 correction

	PlotDistribution(gPad, SumData ,"Reconstructed Mass [GeV/c^{2}]","Distribution Width",1,"ePsame",1e-7,SumData->GetBinContent(SumData->GetMaximumBin())*2.33,3,"Data-Driven");   
        PlotDistribution(gPad, SumMC ,"Reconstructed Mass [GeV/c^{2}]","Distribution Width",3,"same",1e-7,SumData->GetBinContent(SumData->GetMaximumBin())*2.33,6,"Helium MC");	

	float D_DT=SumData->Integral(SumData->FindBin(1.6),SumData->FindBin(2.2));
        float T_DT=SumData->Integral(SumData->FindBin(2.7),SumData->FindBin(3.7));
	float D_MC=SumMC->Integral(SumMC->FindBin(1.6),SumMC->FindBin(2.2));
        float T_MC=SumMC->Integral(SumMC->FindBin(2.7),SumMC->FindBin(3.7));

	std::string MC_Cross_Section = "He fragm. Probability: " + Convert(SumMC->Integral()*100) + "%"; 
	std::string DT_Cross_Section = "He fragm. Probability: " + Convert(SumData->Integral()*100) + "%"; 
	std::string MC_ToverD = "He->T / He->D: " + Convert(T_MC/D_MC*100) + "%"; 
	std::string DT_ToverD = "He->T / He->D: " + Convert(T_DT/D_DT*100) + "%"; 



	TPad *newpad=new TPad("newpad","a transparent pad",0.5,0.2,0.8,0.6);
	newpad->SetFillStyle(4000);
	newpad->Draw();
	newpad->cd();
	
	TText *t = new TText(.1,.7,MC_Cross_Section.c_str());
	t->SetTextAlign(11);
	t->SetTextColor(kGreen+1);
	t->SetTextFont(43);
	t->SetTextSize(40);
	t->Draw();
	TText *f = new TText(.1,.5,MC_ToverD.c_str());
	f->SetTextAlign(11);
	f->SetTextColor(kGreen+1);
	f->SetTextFont(43);
	f->SetTextSize(40);
	f->Draw("same");
	
	c4->cd();
	TPad *newpad2=new TPad("newpad2","another pad",0.2,0.2,0.5,0.6);
	newpad2->SetFillStyle(4000);
	newpad2->Draw();
	newpad2->cd();
	
	TText *t2 = new TText(.1,.7,DT_Cross_Section.c_str());
	t2->SetTextAlign(11);
	t2->SetTextColor(kBlack);
	t2->SetTextFont(43);
	t2->SetTextSize(40);
	t2->Draw();
	TText *f2 = new TText(.1,.5,DT_ToverD.c_str());
	f2->SetTextAlign(11);
	f2->SetTextColor(kBlack);
	f2->SetTextFont(43);
	f2->SetTextSize(40);
	f2->Draw("same");




	Plots.Add(c4);	
	Plots.writeObjsInFolder((basename + "/Fragments").c_str());
}

void DrawFits(TFile * file, std::string basename, Binning Bins, FileSaver Plots){

	std::string pathdata  = (basename + "/L1Distrib./");
        std::string pathtempl1= (basename + "/Q1Templates/");
        std::string pathtempl2= (basename + "/Q2Templates/");
	std::string pathfits  = (basename + "/Fits/");

	std::vector<TH1F*> Data=GetListOfTemplates(file, pathdata);
	std::vector<TH1F*> Templates1=GetListOfTemplates(file, pathtempl1);		
	std::vector<TH1F*> Templates2=GetListOfTemplates(file, pathtempl2);		
	std::vector<TH1F*> Fits=GetListOfTemplates(file, pathfits);		


	cout<<Data.size()<<" "<<Templates1.size()<<" "<<Templates2.size()<<endl;
	for(int bin =0; bin < Bins.size();bin++){
		TCanvas * c3 = new TCanvas(("Bin "+to_string(bin)).c_str());
		c3->SetCanvasSize(2000,1500);
		TH1F * Sum = (TH1F*) Templates1[bin]->Clone();
		Sum->Add(Templates2[bin]);

		PlotDistribution(gPad, Templates1[bin] ,"Layer 1 Charge","Counts",2,"same",10,Data[bin]->GetBinContent(Data[bin]->GetMaximumBin())*1.33,6,"Layer 2 Q=1 Template");
		PlotDistribution(gPad, Templates2[bin] ,"Layer 1 Charge","Counts",3,"same",10,Data[bin]->GetBinContent(Data[bin]->GetMaximumBin())*1.33,6,"Layer 2 Q=2 Template");
		PlotDistribution(gPad, Sum 	       ,"Layer 1 Charge","Counts",1,"same",10,Data[bin]->GetBinContent(Data[bin]->GetMaximumBin())*1.33,2,"Sum of contributions");
		PlotDistribution(gPad, Fits[bin]       ,"Layer 1 Charge","Counts",6,"same",10,Data[bin]->GetBinContent(Data[bin]->GetMaximumBin())*1.33,2,"Fraction Fit");
	
		PlotDistribution(gPad, Data[bin],"Layer 1 Charge","Counts",1,"ePsame",1,1e5,2,"ISS data",false,true);


		Plots.Add(c3);
	}

	Plots.writeObjsInFolder((basename + "/Charge Fits").c_str());




}



















