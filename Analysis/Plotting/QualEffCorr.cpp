#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "../../include/binning.h"
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
#include "../../include/GlobalBinning.h"
#include "TKey.h"
#include "TFractionFitter.h"

#include "../../Ntuple-making/Commonglobals.cpp"
#include "../../include/Variables.hpp"
#include "../../include/Cuts.h"
#include "../../include/filesaver.h"

#include "../../include/FitError.h"
#include "../../include/Resolution.h"
#include "../../include/Efficiency.h"

#include "../../include/PlottingFunctions.h"

#include "../../include/EffCorrTemplate.h"

void PlotCorrection(FileSaver Plots, EffCorrTemplate* Correction, Binning bins);
std::vector<TH1F*> GetListOfTemplates(TFile * file,std::string path);
void PlotFitLatitudes(FileSaver finalHistos, FileSaver Plots, EffCorrTemplate* Correction,Binning bins,std::string lat);
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

        cout<<"**TOF**"<<endl;
        ToFDB.Print();

        cout<<"**NaF**"<<endl;
        NaFDB.Print();

        cout<<"**Agl**"<<endl;
        AglDB.Print();

        ToFDB.UseBetaEdges();
        NaFDB.UseBetaEdges();
        AglDB.UseBetaEdges();

        PRB.UseREdges();


        cout<<endl;

	cout<<"**************************** PLOTTING ***************************************"<<endl;

	EffCorrTemplate* DistCorr_TOF = new EffCorrTemplate(finalHistos,"DistanceCorrTOF","Quality Eff. Corr",ToFDB,"ControlSample","ControlSample&DistanceCut","","");	
	EffCorrTemplate* DistCorr_NaF = new EffCorrTemplate(finalHistos,"DistanceCorrNaF","Quality Eff. Corr",NaFDB,"ControlSample&IsFromNaF","ControlSample&IsFromNaF&DistanceCut","","",true);	
	EffCorrTemplate* DistCorr_Agl = new EffCorrTemplate(finalHistos,"DistanceCorrAgl","Quality Eff. Corr",AglDB,"ControlSample&IsFromAgl","ControlSample&IsFromAgl&DistanceCut","","",true);	
	EffCorrTemplate* LikCorr_TOF = new EffCorrTemplate(finalHistos,"LikelihoodCorrTOF","Quality Eff. Corr",ToFDB,"ControlSample&DistanceCut","ControlSample&DistanceCut&LikelihoodCut","","");	
	EffCorrTemplate* LikCorr_NaF = new EffCorrTemplate(finalHistos,"LikelihoodCorrNaF","Quality Eff. Corr",NaFDB,"ControlSample&DistanceCut&IsFromNaF","ControlSample&IsFromNaF&DistanceCut&LikelihoodCut","","",true);	
	EffCorrTemplate* LikCorr_Agl = new EffCorrTemplate(finalHistos,"LikelihoodCorrAgl","Quality Eff. Corr",AglDB,"ControlSample&DistanceCut&IsFromAgl","ControlSample&IsFromAgl&DistanceCut&LikelihoodCut","","",true);	

	PlotCorrection(Plots,DistCorr_TOF,ToFDB);
	PlotCorrection(Plots,DistCorr_NaF,NaFDB);
	PlotCorrection(Plots,DistCorr_Agl,AglDB);

	PlotCorrection(Plots,LikCorr_TOF,ToFDB);
	PlotCorrection(Plots,LikCorr_NaF,NaFDB);
	PlotCorrection(Plots,LikCorr_Agl,AglDB);

	PlotFitLatitudes(finalHistos,Plots,DistCorr_TOF,ToFDB,"Glob");
	PlotFitLatitudes(finalHistos,Plots,DistCorr_NaF,NaFDB,"Glob");
	PlotFitLatitudes(finalHistos,Plots,DistCorr_Agl,AglDB,"Glob");

	PlotFitLatitudes(finalHistos,Plots,LikCorr_TOF,ToFDB,"Glob");
	PlotFitLatitudes(finalHistos,Plots,LikCorr_NaF,NaFDB,"Glob");
	PlotFitLatitudes(finalHistos,Plots,LikCorr_Agl,AglDB,"Glob");

	return 0;
}

void PlotCorrection(FileSaver Plots, EffCorrTemplate* Correction, Binning bins){
	
	TCanvas * c3 = new TCanvas("Latitude Corrections"); 
        c3->SetCanvasSize(2000,1500);

	for(int lat=0;lat<10;lat++)
		PlotTH1FintoGraph(gPad,bins, (TH1F*)Correction->GetCorrectionLat_P(lat),"Kinetic Energy [GeV/nucl.]", "Efficiency Correction Factor",55+5*lat,true,"Psame",0.5*bins.EkPerMassBin(0),1.3*bins.EkPerMassBin(bins.size()-1),0.1,1.45,("Lat zone " + to_string(lat)).c_str());

	TCanvas * c4 = new TCanvas("Global Corrections"); 
        c4->SetCanvasSize(2000,1500);

	PlotTH1FintoGraph(gPad,bins, (TH1F*)Correction->GetGlobCorrection_P(),"Kinetic Energy [GeV/nucl.]", "Efficiency Correction Factor",2,true,"Psame",0.5*bins.EkPerMassBin(0),1.3*bins.EkPerMassBin(bins.size()-1),0.1,1.45,"Protons Global Correction");
	PlotTH1FintoGraph(gPad,bins, (TH1F*)Correction->GetGlobCorrection_D(),"Kinetic Energy [GeV/nucl.]", "Efficiency Correction Factor",4,true,"Psame",0.5*bins.EkPerMassBin(0),1.3*bins.EkPerMassBin(bins.size()-1),0.1,1.45,"Deuterons Global Correction");
	
	Plots.Add(c3);
	Plots.Add(c4);
	Plots.writeObjsInFolder(Correction->GetName().c_str());
	
}

void PlotFitLatitudes(FileSaver finalHistos, FileSaver Plots, EffCorrTemplate* Correction,Binning bins,std::string lat){
	
	std::vector<TH1F *> HistosBefore =GetListOfTemplates(finalHistos.GetFile(),(Correction->GetDirPath()+"/"+Correction->GetName()+"/TFit/Before/"+lat).c_str());
	std::vector<TH1F *> HistosAfter  =GetListOfTemplates(finalHistos.GetFile(),(Correction->GetDirPath()+"/"+Correction->GetName()+"/TFit/After/"+lat).c_str());
	cout<<(Correction->GetDirPath()+"/"+Correction->GetName()+"/TFit/Before/Glob").c_str()<<endl;

	TCanvas * c1 = new TCanvas("Global Fits - Before");
	c1->SetCanvasSize(2000,1500);
	c1->Divide(6,3);
	TCanvas * c2 = new TCanvas("Global Fits - After");
	c2->SetCanvasSize(2000,1500);
	c2->Divide(6,3);
	
	for(int bin =0; bin<bins.size();bin++){
		c1->cd(bin+1);
		gPad->SetLogy();
		float maxscale=HistosBefore[(bin+1)*3-1]->GetBinContent(HistosBefore[(bin+1)*3-1]->GetMaximumBin())*1.3;
		PlotDistribution(gPad, HistosBefore[(bin+1)*3-3],"[GeV/c^2]","",2,"same",1,maxscale,10,"Protons MC Template");
		PlotDistribution(gPad, HistosBefore[(bin+1)*3-2],"[GeV/c^2]","",4,"same",1,maxscale,10,"Deuterons MC Template");
		PlotDistribution(gPad, HistosBefore[(bin+1)*3-1],"[GeV/c^2]","",1,"ePsame",1,maxscale,2.5,"ISS data",false,true);				
		c2->cd(bin+1);
		gPad->SetLogy();
		PlotDistribution(gPad, HistosAfter[(bin+1)*3-3],"[GeV/c^2]","",2,"same",1,maxscale,10,"Protons MC Template");
		PlotDistribution(gPad, HistosAfter[(bin+1)*3-2],"[GeV/c^2]","",4,"same",1,maxscale,10,"Deuterons MC Template");
		PlotDistribution(gPad, HistosAfter[(bin+1)*3-1],"[GeV/c^2]","",1,"ePsame",1,maxscale,2.5,"ISS data",false,true);				
	}

	Plots.Add(c1);
        Plots.Add(c2);
        Plots.writeObjsInFolder((Correction->GetName()+"/Fits/"+lat).c_str());	
	return;
}

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
	
