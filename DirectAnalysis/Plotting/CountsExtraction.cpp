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
#include "TPaveLabel.h"
#include "TGraph.h"
#include "../include/Variables.hpp"
#include "../include/Cuts.h"
#include "../include/filesaver.h"
#include "../include/TemplateFITbetasmear.h"
#include "RangeMerger.h"
#include "../include/PlottingFunctions.h"
#include <sstream>

void CheckSlowDown(FileSaver Plots);


std::string Convert (float number){
    std::ostringstream buff;
    buff<<number;
    std::string output= buff.str();
    output.erase(4,output.end()-output.begin()-4);	
    return output;   
}




int colorbase = 55;
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

void NormalizeResiduals(TH1F *ratio){
	for (int i=0;i<ratio->GetNbinsX();i++){
		if(ratio->GetBinError(i+1)>0) { ratio->SetBinContent(i+1,ratio->GetBinContent(i+1)/ratio->GetBinError(i+1));
					        ratio->SetBinError(i+1,1);
						}	
		else  ratio->SetBinContent(i+1,0);	
	}
}
void CleanResiduals(TH1F *Ratio,TH1F *MC,TH1F * Datas){
	for (int i=0;i<Ratio->GetNbinsX();i++){
		if(Datas->GetBinContent(i+1)==0||MC->GetBinContent(i+1)==0)
			Ratio->SetBinContent(i+1,-9999);
	}
}

void DrawParameters(FileSaver finalHistos,FileSaver Plots,std::string path, Binning bins, std::string title,std::string x, std::string x_udm,float xmin, float xmax,float y1min,float y1max,float y2min,float y2max);
void DrawFits(TemplateFIT * FIT,FileSaver finalHistos, FileSaver Plots, bool IsFitNoise = false,bool IsSmearingCheck = false);

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
	cout<<"****************************** CHECK SLOWDOWN ***************************************"<<endl;
	CheckSlowDown(Plots); 
        cout<<"****************************** PLOTTING FITS ***************************************"<<endl;

	//	TemplateFIT * SmearingCheck = new TemplateFIT(finalHistos,"SmearingCheck",PRB);

	TemplateFIT *ToFfits= new TemplateFIT(finalHistos,"TOFDfits",Global.GetToFDBins(),false,11,50,100,1);
	TemplateFIT *NaFfits= new TemplateFIT(finalHistos,"NaFDfits",Global.GetNaFDBins(),true,11,250,150,0);
	TemplateFIT *Aglfits= new TemplateFIT(finalHistos,"AglDfits",Global.GetAglDBins(),true,11,10,25,1);

	TemplateFIT *ToFfits_P= new TemplateFIT(finalHistos,"TOFPfits",Global.GetToFPBins(),false,11,50,100,1);
	TemplateFIT *NaFfits_P= new TemplateFIT(finalHistos,"NaFPfits",Global.GetNaFPBins(),true,11,250,150,0);
	TemplateFIT *Aglfits_P= new TemplateFIT(finalHistos,"AglPfits",Global.GetAglPBins(),true,11,10,25,1);


	//	DrawFits(SmearingCheck,finalHistos,Plots,false,true);
	DrawFits(ToFfits,finalHistos,Plots);
	DrawFits(NaFfits,finalHistos,Plots);
	DrawFits(Aglfits,finalHistos,Plots);
//
	DrawFits(ToFfits_P,finalHistos,Plots);
	DrawFits(NaFfits_P,finalHistos,Plots);
	DrawFits(Aglfits_P,finalHistos,Plots);


	cout<<"****************************** RESULTS ***************************************"<<endl;

	TCanvas * c4 = new TCanvas("Deuteron Counts");

	//std::string pathresCheck = (SmearingCheck->GetName() + "/Fit Results/");
	std::string pathresTOF   = (ToFfits->GetName() + "/Fit Results/");
	std::string pathresNaF   = (NaFfits->GetName() + "/Fit Results/");
	std::string pathresAgl   = (Aglfits->GetName() + "/Fit Results/");
		
	TH1F * DCountsTOF = (TH1F*) finalHistos.Get((pathresTOF+"Deuteron Counts").c_str());
	TH1F * DCountsNaF = (TH1F*) finalHistos.Get((pathresNaF+"Deuteron Counts").c_str());
	TH1F * DCountsAgl = (TH1F*) finalHistos.Get((pathresAgl+"Deuteron Counts").c_str());
	TH1F * PCountsTOF = (TH1F*) finalHistos.Get((pathresTOF+"Proton Counts").c_str());
	TH1F * PCountsNaF = (TH1F*) finalHistos.Get((pathresNaF+"Proton Counts").c_str());
	TH1F * PCountsAgl = (TH1F*) finalHistos.Get((pathresAgl+"Proton Counts").c_str());
	TH1F * TCountsTOF = (TH1F*) finalHistos.Get((pathresTOF+"Tritium Counts").c_str());
	TH1F * TCountsNaF = (TH1F*) finalHistos.Get((pathresNaF+"Tritium Counts").c_str());
	TH1F * TCountsAgl = (TH1F*) finalHistos.Get((pathresAgl+"Tritium Counts").c_str());
	
	TH1F * DCountsPrimTOF = (TH1F*) finalHistos.Get((pathresTOF+"Primary Deuteron Counts").c_str());
	TH1F * DCountsPrimNaF = (TH1F*) finalHistos.Get((pathresNaF+"Primary Deuteron Counts").c_str());
	TH1F * DCountsPrimAgl = (TH1F*) finalHistos.Get((pathresAgl+"Primary Deuteron Counts").c_str());
	TH1F * PCountsPrimTOF = (TH1F*) finalHistos.Get((pathresTOF+"Primary Proton Counts").c_str());
	TH1F * PCountsPrimNaF = (TH1F*) finalHistos.Get((pathresNaF+"Primary Proton Counts").c_str());
	TH1F * PCountsPrimAgl = (TH1F*) finalHistos.Get((pathresAgl+"Primary Proton Counts").c_str());
	
	TH1F * PCountsPrimTOF_rigbins = (TH1F*) finalHistos.Get("TOFPCounts/TOFPCounts/TOFPCounts_before");
	TH1F * PCountsPrimNaF_rigbins = (TH1F*) finalHistos.Get("NaFPCounts/NaFPCounts/NaFPCounts_before");
	TH1F * PCountsPrimAgl_rigbins = (TH1F*) finalHistos.Get("AglPCounts/AglPCounts/AglPCounts_before");
	

	
	PlotTH1FintoGraph(gPad,Global.GetToFDBins(), DCountsTOF,"Kinetic Energy [GeV/nucl.]", "Counts",4,true,"Psame",0.1,10,10,2*DCountsTOF->GetBinContent(DCountsTOF->GetMaximumBin()),"Deuteron Counts (TOF)",8);
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(), DCountsNaF,"Kinetic Energy [GeV/nucl.]", "Counts",4,true,"Psame",0.1,10,10,7e4,"Deuteron Counts (NaF)",22);
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(), DCountsAgl,"Kinetic Energy [GeV/nucl.]", "Counts",4,true,"Psame",0.1,10,10,7e4,"Deuteron Counts (Agl)",29);

	PlotTH1FintoGraph(gPad,Global.GetToFDBins(), DCountsTOF,"Kinetic Energy [GeV/nucl.]", "Counts",4,true,"Psame",0.1,10,10,7e4,"Primary Counts (TOF)",4);
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(), DCountsNaF,"Kinetic Energy [GeV/nucl.]", "Counts",4,true,"Psame",0.1,10,10,7e4,"Primary Counts (NaF)",26);
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(), DCountsAgl,"Kinetic Energy [GeV/nucl.]", "Counts",4,true,"Psame",0.1,10,10,7e4,"Primary Counts (Agl)",30);


	TCanvas * c4_ = new TCanvas("Proton Counts");

	PlotTH1FintoGraph(gPad,Global.GetToFDBins(), PCountsTOF,"Kinetic Energy [GeV/nucl.]", "Counts",2,true,"Psame",0.1,10,10,2*PCountsTOF->GetBinContent(DCountsTOF->GetMaximumBin()),"Proton Counts (TOF)",8);
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(), PCountsNaF,"Kinetic Energy [GeV/nucl.]", "Counts",2,true,"Psame",0.1,10,10,7e4,"Proton Counts (NaF)",22);
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(), PCountsAgl,"Kinetic Energy [GeV/nucl.]", "Counts",2,true,"Psame",0.1,10,10,7e4,"Proton Counts (Agl)",29);

	PlotTH1FintoGraph(gPad,Global.GetToFDBins(), PCountsPrimTOF,"Kinetic Energy [GeV/nucl.]", "Counts",2,true,"Psame",0.1,10,10,7e4,"Primary Counts (TOF)",4);
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(), PCountsPrimNaF,"Kinetic Energy [GeV/nucl.]", "Counts",2,true,"Psame",0.1,10,10,7e4,"Primary Counts (NaF)",26);
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(), PCountsPrimAgl,"Kinetic Energy [GeV/nucl.]", "Counts",2,true,"Psame",0.1,10,10,7e4,"Primary Counts (Agl)",30);

	PlotTH1FintoGraph(gPad,Global.GetToFDBins(), PCountsPrimTOF_rigbins,"Kinetic Energy [GeV/nucl.]", "Counts",1,true,"Psame",0.1,10,10,7e4,"RigBin Counts (TOF)",4);
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(), PCountsPrimNaF_rigbins,"Kinetic Energy [GeV/nucl.]", "Counts",1,true,"Psame",0.1,10,10,7e4,"RigBin Counts (NaF)",26);
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(), PCountsPrimAgl_rigbins,"Kinetic Energy [GeV/nucl.]", "Counts",1,true,"Psame",0.1,10,10,7e4,"RigBin Counts (Agl)",30);


	Plots.Add(c4);
	Plots.Add(c4_);
	Plots.writeObjsInFolder("Results");




	TCanvas * c5 = new TCanvas("D/P Raw Counts ratio");
	
	TH1F * RatioTOF = (TH1F*)DCountsTOF->Clone();
	RatioTOF->Divide(PCountsTOF);
	TH1F * RatioNaF = (TH1F*)DCountsNaF->Clone();
	RatioNaF->Divide(PCountsNaF);
	TH1F * RatioAgl = (TH1F*)DCountsAgl->Clone();
	RatioAgl->Divide(PCountsAgl);

	PlotTH1FintoGraph(gPad,Global.GetToFDBins(), RatioTOF,"Kinetic Energy [GeV/nucl.]", "Counts ratio",2,true,"Psame",0.1,10,1e-3,1,"D/P Counts ratio (TOF)",8);
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(), RatioNaF,"Kinetic Energy [GeV/nucl.]", "Counts ratio",2,true,"Psame",0.1,10,1e-3,1,"D/P Counts ratio (NaF)",22);
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(), RatioAgl,"Kinetic Energy [GeV/nucl.]", "Counts ratio",2,true,"Psame",0.1,10,1e-3,1,"D/P Counts ratio (Agl)",29);


	Plots.Add(c5);
	Plots.writeObjsInFolder("Results");

	TCanvas * c5_ = new TCanvas("T/D Raw Counts ratio");
	
	TH1F * RatioTOF_T = (TH1F*)TCountsTOF->Clone();
	RatioTOF_T->Divide(DCountsTOF);
	TH1F * RatioNaF_T = (TH1F*)TCountsNaF->Clone();
	RatioNaF_T->Divide(DCountsNaF);
	TH1F * RatioAgl_T = (TH1F*)TCountsAgl->Clone();
	RatioAgl_T->Divide(DCountsAgl);

	for(int i=0;i<RatioTOF_T->GetNbinsX();i++) RatioTOF_T->SetBinError(i+1,RatioTOF_T->GetBinError(i+1)/2);
	for(int i=0;i<RatioNaF_T->GetNbinsX();i++) RatioNaF_T->SetBinError(i+1,RatioNaF_T->GetBinError(i+1)/2);
	for(int i=0;i<RatioAgl_T->GetNbinsX();i++) RatioAgl_T->SetBinError(i+1,RatioAgl_T->GetBinError(i+1)/2);


	PlotTH1FintoGraph(gPad,Global.GetToFDBins(), RatioTOF_T,"Kinetic Energy [GeV/nucl.]", "Counts ratio",3,true,"Psame",0.1,10,1e-6,0.1,"T/D Counts ratio (TOF)",8);
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(), RatioNaF_T,"Kinetic Energy [GeV/nucl.]", "Counts ratio",3,true,"Psame",0.1,10,1e-6,0.1,"T/D Counts ratio (NaF)",22);
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(), RatioAgl_T,"Kinetic Energy [GeV/nucl.]", "Counts ratio",3,true,"Psame",0.1,10,1e-6,0.1,"T/D Counts ratio (Agl)",29);


	Plots.Add(c5_);
	Plots.writeObjsInFolder("Results");


	TCanvas * c6 = new TCanvas("Uncertainty Break-down");

	TH1F * StatErrTOFP = (TH1F*) finalHistos.Get((pathresTOF+"StatErrorP").c_str());
       	TH1F * StatErrTOFD = (TH1F*) finalHistos.Get((pathresTOF+"StatErrorD").c_str());
       	TH1F * StatErrTOFT = (TH1F*) finalHistos.Get((pathresTOF+"StatErrorT").c_str());
	TH1F * SystErrTOF = (TH1F*) finalHistos.Get((pathresTOF+"SystError").c_str());
	TH1F * HeCoErrTOF = (TH1F*) finalHistos.Get((pathresTOF+"HeliumContamination/HeContError").c_str());
	
	TH1F * StatErrNaFP = (TH1F*) finalHistos.Get((pathresNaF+"StatErrorP").c_str());
       	TH1F * StatErrNaFD = (TH1F*) finalHistos.Get((pathresNaF+"StatErrorD").c_str());
       	TH1F * StatErrNaFT = (TH1F*) finalHistos.Get((pathresNaF+"StatErrorT").c_str());
        TH1F * SystErrNaF = (TH1F*) finalHistos.Get((pathresNaF+"SystError").c_str());
	TH1F * HeCoErrNaF = (TH1F*) finalHistos.Get((pathresNaF+"HeliumContamination/HeContError").c_str());

	TH1F * StatErrAglP = (TH1F*) finalHistos.Get((pathresAgl+"StatErrorP").c_str());
       	TH1F * StatErrAglD = (TH1F*) finalHistos.Get((pathresAgl+"StatErrorD").c_str());
       	TH1F * StatErrAglT = (TH1F*) finalHistos.Get((pathresAgl+"StatErrorT").c_str());
        TH1F * SystErrAgl = (TH1F*) finalHistos.Get((pathresAgl+"SystError").c_str());
	TH1F * HeCoErrAgl = (TH1F*) finalHistos.Get((pathresAgl+"HeliumContamination/HeContError").c_str());
	
	TH1F * TotErrTOF  = (TH1F *)StatErrTOFP->Clone();
	TH1F * TotErrNaF  = (TH1F *)StatErrNaFP->Clone();
	TH1F * TotErrAgl  = (TH1F *)StatErrAglP->Clone();

	HeCoErrTOF->Divide(DCountsTOF);
	HeCoErrNaF->Divide(DCountsNaF);
	HeCoErrAgl->Divide(DCountsAgl);


	
	for(int bin=0;bin<DCountsTOF->GetNbinsX();bin++) {
		if(DCountsTOF->GetBinContent(bin+1)>0){
			StatErrTOFD->SetBinContent(bin+1,StatErrTOFD->GetBinContent(bin+1)*PCountsTOF->GetBinContent(bin+1)/DCountsTOF->GetBinContent(bin+1));
			StatErrTOFD->SetBinError(bin+1,0);
			TotErrTOF->SetBinContent(bin+1,DCountsTOF->GetBinError(bin+1)/DCountsTOF->GetBinContent(bin+1));
			TotErrTOF->SetBinError(bin+1,0);
		}
		if(DCountsNaF->GetBinContent(bin+1)>0){	
			StatErrNaFD->SetBinContent(bin+1,StatErrNaFD->GetBinContent(bin+1)*PCountsNaF->GetBinContent(bin+1)/DCountsNaF->GetBinContent(bin+1));
			StatErrNaFD->SetBinError(bin+1,0);
			TotErrNaF->SetBinContent(bin+1,DCountsNaF->GetBinError(bin+1)/DCountsNaF->GetBinContent(bin+1));
			TotErrNaF->SetBinError(bin+1,0);
		}
		if(DCountsAgl->GetBinContent(bin+1)>0){	
			StatErrAglD->SetBinContent(bin+1,StatErrAglD->GetBinContent(bin+1)*PCountsAgl->GetBinContent(bin+1)/DCountsAgl->GetBinContent(bin+1));
                        StatErrAglD->SetBinError(bin+1,0);
			TotErrAgl->SetBinContent(bin+1,DCountsAgl->GetBinError(bin+1)/DCountsAgl->GetBinContent(bin+1));
			TotErrAgl->SetBinError(bin+1,0);
		}	

	}	



	


/*	c6->cd(1);
	PlotDistribution(gPad,TotErrTOF ,"TOF Range Bin","Relative error",2,"same",1e-4,1.1,10,"T. Fit Total Error");
	PlotDistribution(gPad,SystErrTOF,"TOF Range Bin","Relative error",4,"same",1e-4,1.1,4,"T. Fit Systematic Error");
	PlotDistribution(gPad,StatErrTOFD,"TOF Range Bin","Relative error",1,"same",1e-4,1.1,4,"T. Fit Statistical Error");
	PlotDistribution(gPad,HeCoErrTOF,"TOF Range Bin","Relative error",3,"same",1e-4,1.1,4,"Helium Fragm. Error");
	c6->cd(2);
	PlotDistribution(gPad,TotErrNaF ,"NaF Range Bin","Relative error",2,"same",1e-4,1.1,10,"T. Fit Total Error");
	PlotDistribution(gPad,SystErrNaF,"NaF Range Bin","Relative error",4,"same",1e-4,1.1,4,"T. Fit Systematic Error");
	PlotDistribution(gPad,StatErrNaFD,"NaF Range Bin","Relative error",1,"same",1e-4,1.1,4,"T. Fit Statistical Error");
	PlotDistribution(gPad,HeCoErrNaF,"NaF Range Bin","Relative error",3,"same",1e-4,1.1,4,"Helium Fragm. Error");
	c6->cd(3);
	PlotDistribution(gPad,TotErrAgl ,"Agl Range Bin","Relative error",2,"same",1e-4,1.1,10,"T. Fit Total Error");
	PlotDistribution(gPad,SystErrAgl,"Agl Range Bin","Relative error",4,"same",1e-4,1.1,4,"T. Fit Systematic Error");
	PlotDistribution(gPad,StatErrAglD,"Agl Range Bin","Relative error",1,"same",1e-4,1.1,4,"T. Fit Statistical Error");
	PlotDistribution(gPad,HeCoErrAgl,"Agl Range Bin","Relative error",3,"same",1e-4,1.1,4,"Helium Fragm. Error");
*/


	TH1F * MergedTOT  =  Global.MergeSubDResult_D(TotErrTOF ,TotErrNaF,TotErrAgl);
	TH1F * MergedSyst =  Global.MergeSubDResult_D(SystErrTOF,SystErrNaF,SystErrAgl);
	TH1F * MergedStat =  Global.MergeSubDResult_D(StatErrTOFD,StatErrNaFD,StatErrAglD);
	TH1F * MergedHeCo =  Global.MergeSubDResult_D(HeCoErrTOF,HeCoErrNaF,HeCoErrAgl);

	TH1F * MergedTOT_plot =ConvertBinnedHisto(MergedTOT ,"Total",Global.GetGlobalDBins());
	TH1F * MergedSyst_plot=ConvertBinnedHisto(MergedSyst,"Fit Syst. Error",Global.GetGlobalDBins());
	TH1F * MergedStat_plot=ConvertBinnedHisto(MergedStat,"Fit Stat. Error",Global.GetGlobalDBins());
	TH1F * MergedHeCo_plot=ConvertBinnedHisto(MergedHeCo,"He Fragm. Error",Global.GetGlobalDBins());

	MergedTOT_plot ->Draw("hist,same"); 
	MergedSyst_plot ->Draw("hist,same"); 
	MergedStat_plot ->Draw("hist,same"); 
	MergedHeCo_plot ->Draw("hist,same"); 


	Plots.Add(c6);
	Plots.writeObjsInFolder("Results");




	TCanvas * c7 = new TCanvas("Chi Square");
	c7->SetCanvasSize(5000,1000);
	c7->Divide(3,1);

	TH1F * BestChiSquareTOF         = (TH1F*) finalHistos.Get((pathresTOF+"Best ChiSquare").c_str());
        TH1F * OriginalChiSquareTOF     = (TH1F*) finalHistos.Get((pathresTOF+"Original ChiSquare").c_str());
        TH1F * BestChiSquareNaF         = (TH1F*) finalHistos.Get((pathresNaF+"Best ChiSquare").c_str());
        TH1F * OriginalChiSquareNaF     = (TH1F*) finalHistos.Get((pathresNaF+"Original ChiSquare").c_str());
        TH1F * BestChiSquareAgl         = (TH1F*) finalHistos.Get((pathresAgl+"Best ChiSquare").c_str());
        TH1F * OriginalChiSquareAgl     = (TH1F*) finalHistos.Get((pathresAgl+"Original ChiSquare").c_str());

	c7->cd(1);
	PlotDistribution(gPad,BestChiSquareTOF ,"TOF Range Bin","#chi^2 of T. Fit",2,"Psame",1e-1,30,3,"Best #chi^{2} mod. Template");
	PlotDistribution(gPad,OriginalChiSquareTOF,"TOF Range Bin","#chi^2 of T. Fit",4,"Psame",1e-1,30,3,"Original MC Templates");
	
	c7->cd(2);
	PlotDistribution(gPad,BestChiSquareNaF ,"NaF Range Bin","#chi^2 of T. Fit",2,"Psame",1e-1,35,3,"Best #chi^{2} mod. Template");
	PlotDistribution(gPad,OriginalChiSquareNaF,"NaF Range Bin","#chi^2 of T. Fit",4,"Psame",1e-1,35,3,"Original MC Templates");

	c7->cd(3);
	PlotDistribution(gPad,BestChiSquareAgl ,"Agl Range Bin","#chi^2 of T. Fit",2,"Psame",1e-1,35,3,"Best #chi^{2} mod. Template");
	PlotDistribution(gPad,OriginalChiSquareAgl,"Agl Range Bin","#chi^2 of T. Fit",4,"Psame",1e-1,35,3,"Original MC Templates");
	
	Plots.Add(c7);
	Plots.writeObjsInFolder("Results");



	cout<<"****************************** RESULTS RIGIDITY ***************************************"<<endl;
	std::string pathresTOF_P   = (ToFfits_P->GetName() + "/Fit Results/");
	std::string pathresNaF_P   = (NaFfits_P->GetName() + "/Fit Results/");
	std::string pathresAgl_P   = (Aglfits_P->GetName() + "/Fit Results/");
	

	SetUpRigTOIBinning();


	TH1F * PCountsPrim_rigTOF = (TH1F*) finalHistos.Get((pathresTOF_P+"Proton Counts").c_str());
	TH1F * PCountsPrim_rigNaF = (TH1F*) finalHistos.Get((pathresNaF_P+"Proton Counts").c_str());
	TH1F * PCountsPrim_rigAgl = (TH1F*) finalHistos.Get((pathresAgl_P+"Proton Counts").c_str());

	TCanvas * d4 = new TCanvas("Rigidity Counts");
	d4->cd();
	gPad->SetLogx();
	gPad->SetLogy();

	for(int i=0; i<Global.GetToFPBins().size();i++) PCountsPrim_rigTOF->SetBinContent(i+1,PCountsPrim_rigTOF->GetBinContent(i+1)/(Global.GetToFPBins().RigTOIBins()[i+1]-Global.GetToFPBins().RigTOIBins()[i])); 
	for(int i=0; i<Global.GetNaFPBins().size();i++) PCountsPrim_rigNaF->SetBinContent(i+1,PCountsPrim_rigNaF->GetBinContent(i+1)/(Global.GetNaFPBins().RigTOIBins()[i+1]-Global.GetNaFPBins().RigTOIBins()[i])); 
	for(int i=0; i<Global.GetAglPBins().size();i++) PCountsPrim_rigAgl->SetBinContent(i+1,PCountsPrim_rigAgl->GetBinContent(i+1)/(Global.GetAglPBins().RigTOIBins()[i+1]-Global.GetAglPBins().RigTOIBins()[i])); 

	PlotTH1FintoGraph(gPad,Global.GetToFPBins(), PCountsPrim_rigTOF,"R [GV]", "Counts density [GV^{-1}]",2,false,"Psame",0.5,20,10,10*PCountsPrim_rigTOF->GetBinContent(PCountsPrim_rigTOF->GetMaximumBin()),"P Counts (TOF)",4);
	PlotTH1FintoGraph(gPad,Global.GetNaFPBins(), PCountsPrim_rigNaF,"R [GV]", "Counts density [GV^{-1}]",2,false,"Psame",0.5,20,10,10*PCountsPrim_rigTOF->GetBinContent(PCountsPrim_rigTOF->GetMaximumBin()),"P Counts (NaF)",26);
	PlotTH1FintoGraph(gPad,Global.GetAglPBins(), PCountsPrim_rigAgl,"R [GV]", "Counts density [GV^{-1}]",2,false,"Psame",0.5,20,10,10*PCountsPrim_rigTOF->GetBinContent(PCountsPrim_rigTOF->GetMaximumBin()),"P Counts (Agl)",30);

	for(int i=0; i<Global.GetToFDBins().size();i++) DCountsTOF->SetBinContent(i+1,DCountsTOF->GetBinContent(i+1)/(Global.GetToFDBins().RigTOIBins()[i+1]-Global.GetToFDBins().RigTOIBins()[i])); 
	for(int i=0; i<Global.GetNaFDBins().size();i++) DCountsNaF->SetBinContent(i+1,DCountsNaF->GetBinContent(i+1)/(Global.GetNaFDBins().RigTOIBins()[i+1]-Global.GetNaFDBins().RigTOIBins()[i])); 
	for(int i=0; i<Global.GetAglDBins().size();i++) DCountsAgl->SetBinContent(i+1,DCountsAgl->GetBinContent(i+1)/(Global.GetAglDBins().RigTOIBins()[i+1]-Global.GetAglDBins().RigTOIBins()[i])); 

	PlotTH1FintoGraph(gPad,Global.GetToFDBins(), DCountsTOF,"R [GV]", "Counts density [GV^{-1}]",4,false,"Psame",0.5,20,10,7e4,"D Counts (TOF)",4);
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(), DCountsNaF,"R [GV]", "Counts density [GV^{-1}]",4,false,"Psame",0.5,20,10,7e4,"D Counts (NaF)",26);
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(), DCountsAgl,"R [GV]", "Counts density [GV^{-1}]",4,false,"Psame",0.5,20,10,7e4,"D Counts (Agl)",30);

	Plots.Add(d4);
	Plots.writeObjsInFolder("Results Rig");


	
	cout<<"****************************** RESULTS PARAMETERS ***************************************"<<endl;


	
//	DrawParameters(finalHistos,Plots,pathresTOF,Global.GetToFDBins(),"Parameters TOF","Measured #beta","[ps]",0.45,0.9,-200,200,-40,180);
//	DrawParameters(finalHistos,Plots,pathresNaF,Global.GetNaFDBins(),"Parameters NaF","Measured #beta","[rad/10^{4}]",0.7,0.98,-1000,1000,-1000,2000);
//	DrawParameters(finalHistos,Plots,pathresAgl,Global.GetAglDBins(),"Parameters Agl","Measured #beta","[rad/10^{4}]",0.95,1.005,-230,230,-150,400);


/*
	TCanvas * c8 = new TCanvas("Chi Square - Smearing Check");
        c8->SetCanvasSize(5000,1000);

	TH1F * BestChiSquareCheck         = (TH1F*) finalHistos.Get((pathresCheck+"Best ChiSquare").c_str());
        TH1F * OriginalChiSquareCheck      = (TH1F*) finalHistos.Get((pathresCheck+"Original ChiSquare").c_str());

	c8->cd();
	PlotDistribution(gPad,BestChiSquareCheck ,"Rig. Range Bin","#chi^2 of T. Fit",2,"Psame",1e-2,130,3,"Best #chi^{2} mod. Template");
        PlotDistribution(gPad,OriginalChiSquareCheck,"Rig. Range Bin","#chi^2 of T. Fit",4,"Psame",1e-2,130,3,"Original MC Templates");

	Plots.Add(c8);
        Plots.writeObjsInFolder("Results");

	DrawParameters(finalHistos,Plots,pathresCheck,PRB,"Parameters Check","Measured Rig [GV]","[ps]",0.5,120,-150,150,-40,180);
*/






	return 0;
}



void DrawParameters(FileSaver finalHistos,FileSaver Plots,std::string path, Binning bins, std::string title, std::string x,std::string x_udm, float xmin, float xmax,float y1min,float y1max,float y2min,float y2max){

	TCanvas * c7 = new TCanvas(title.c_str());

	TH1F * ShiftBest = (TH1F*) finalHistos.Get((path+"Best Fit Shift").c_str());
        TH1F * SigmaBest = (TH1F*) finalHistos.Get((path+"Best Fit Sigma").c_str());
	

	TPad * c7_up = new TPad("upperPad", "upperPad",0.0,0.5,1.0,1.0);
	c7_up->Draw();

	TPad * c7_do = new TPad("lowerPad", "lowerPad",0.0,0.0,1.0,0.5);
	c7_do->Draw();

	
	PlotTH1FintoGraph(c7,bins, ShiftBest, x.c_str(),  ("Mean shift "+x_udm).c_str(),2,false,"ep",xmin,xmax,y1min,y1max,"Best #chi^{2} Shift",8,false,false);
	PlotTH1FintoGraph(c7_do,bins, SigmaBest, x.c_str(),  ("Additive #sigma "+x_udm).c_str(),4,false,"ep",xmin,xmax,y2min,y2max,"Best #chi^{2} #sigma");
		

	Plots.Add(c7);
	Plots.writeObjsInFolder("Results");
	return;
}

TH1F * CreateNoiseD( TH1F * NoiseP, TH1F * SignP, TH1F * SignD ){
	int binpeakP = SignP->FindBin(1.0);
	int binpeakD = SignD->FindBin(2.0);
	TH1F * NoiseD = (TH1F *) NoiseP->Clone();
	for(int i=0;i<NoiseD->GetNbinsX();i++){
		if(i+1-(binpeakD-binpeakP)<1) NoiseD->SetBinContent(i+1,0);
		NoiseD->SetBinContent(i+1,NoiseP->GetBinContent(i+1-(binpeakD-binpeakP)));
	}
	NoiseD->Scale(SignD->Integral()/SignP->Integral());
	return NoiseD;
}

void DrawFits(TemplateFIT * FIT,FileSaver finalHistos,FileSaver Plots,bool IsFitNoise, bool IsSmearingCheck){

	std::string pathdata  = (FIT->GetName() + "/Fit Results/Data");
	std::string pathtemplP= (FIT->GetName() + "/Fit Results/ScaledTemplatesP");
	std::string pathtemplD= (FIT->GetName() + "/Fit Results/ScaledTemplatesD");	
	std::string pathtemplHe=(FIT->GetName() + "/Fit Results/ScaledTemplatesHe");	
	std::string pathtemplNoise=(FIT->GetName() + "/Fit Results/ScaledTemplatesNoise");	

	std::string pathfit   = (FIT->GetName() + "/Fit Results/FractionFits");
	std::string pathres   = (FIT->GetName() + "/Fit Results/");
	
	TFile * infile = finalHistos.GetFile();	

	TCanvas * CollectiveView = new TCanvas("Summary 1 ");
	CollectiveView->Divide(FIT->GetBinning().size()/3,3);
	CollectiveView->SetCanvasSize(3000,2000);

	TCanvas * CollectiveView2 = new TCanvas("Summary 2 ");
	CollectiveView2->Divide(FIT->GetBinning().size()/3,3);
	CollectiveView2->SetCanvasSize(3000,2000);


	for(int i=0; i<FIT->GetBinning().size();i++){
			
		std::string pathbinP    = pathtemplP + "/Bin"+to_string(i);
		std::string pathbinD    = pathtemplD + "/Bin"+to_string(i);
		std::string pathbinHe   = pathtemplHe+ "/Bin"+to_string(i);
		std::string pathbinNoise   = pathtemplNoise+ "/Bin"+to_string(i);
		std::string pathbindata = pathdata   + "/Bin"+to_string(i);
		std::string pathbinfit  = pathfit    + "/Bin"+to_string(i);

		std::vector<TH1F*> TemplatesP=GetListOfTemplates(infile, pathbinP);
		std::vector<TH1F*> TemplatesD=GetListOfTemplates(infile, pathbinD);		
		std::vector<TH1F*> TemplatesHe=GetListOfTemplates(infile, pathbinHe);		
		std::vector<TH1F*> TemplatesNoise=GetListOfTemplates(infile, pathbinNoise);		

		std::vector<TH1F*> Datas     =GetListOfTemplates(infile, pathbindata);
		std::vector<TH1F*> Fits      =GetListOfTemplates(infile, pathbinfit);

		//cout<<pathbinHe<<" "<<TemplatesHe.size()<<" "<<TemplatesHe[0]<<endl;
	

		TCanvas * c1 = new TCanvas("Modified Templates");
		c1->SetCanvasSize(2000,1500);
		std::string bintitle = ("R: ");// + Convert(FIT->GetBinning().RigTOIBins()[i]) + " GV");
		TPaveLabel* title = new TPaveLabel(0.35,0.94,0.6,0.97,bintitle.c_str(),"brndc");
		title->SetFillColor(0);
		title->SetLineColor(0);
		title->SetShadowColor(0);
		title->SetTextColor(1);
		title->SetTextFont(22);
		title->SetTextSize(2.5);

		for(int j=TemplatesP.size()-1;j>=0;j--){
			if(j==0) PlotDistribution(gPad, TemplatesP[j],"Reconstructed Mass [GeV/c^2]","Counts",1,"same",1,TemplatesP[j]->GetBinContent(TemplatesP[j]->GetMaximumBin())*1.13,3);
			else     PlotDistribution(gPad, TemplatesP[j],"Reconstructed Mass [GeV/c^2]","Counts",colorbase + j,"same",1,TemplatesP[j]->GetBinContent(TemplatesP[j]->GetMaximumBin())*1.13,3,"",false,false,true);		
		}
		c1->cd();
		
		Plots.Add(c1);
		Plots.writeObjsInFolder((FIT->GetName()+"/Fits/Bin"+to_string(i)).c_str());	


		TCanvas * c2 = new TCanvas("Modified T. Fits");
                c2->SetCanvasSize(2000,1500);

		PlotDistribution(gPad, TemplatesP[1], "Reconstructed Mass [GeV/c^2]","Counts",1,"same",3,Datas[0]->GetBinContent(Datas[0]->GetMaximumBin())*1.13,3,"Original Protons MC Template");
		if(!IsSmearingCheck) PlotDistribution(gPad, TemplatesD[0], "Reconstructed Mass [GeV/c^2]","Counts",4,"same",1,1e5,3,"Original Deuterons MC Template");
		if(!IsSmearingCheck) PlotDistribution(gPad, TemplatesHe[0],"Reconstructed Mass [GeV/c^2]","Counts",3,"same",1,1e5,3,"Original He Fragm. MC Template");

		for(int j=TemplatesP.size()-1;j>=1;j--){
                        PlotDistribution(gPad, TemplatesP[j],"Reconstructed Mass [GeV/c^2]","Counts",2,"same",1,1e5,1,"",false,false,true);
			if(!IsSmearingCheck) PlotDistribution(gPad, TemplatesD[j],"Reconstructed Mass [GeV/c^2]","Counts",4,"same",1,1e5,1,"",false,false,true);
                	//PlotDistribution(gPad, TemplatesHe[j],"Reconstructed Mass [GeV/c^2]","Counts",3,"same",1,1e5,1,"",false,false,true);
		}
		PlotDistribution(gPad, Datas[0],"Reconstructed Mass [GeV/c^2]","Counts",1,"ePsame",1,1e5,3,"ISS data",false,true);
		c2->cd();
	
		Plots.Add(c2);
		 Plots.writeObjsInFolder((FIT->GetName()+"/Fits/Bin"+to_string(i)).c_str());
		
		
		TCanvas * c3 = new TCanvas("Template Fits");
                c3->SetCanvasSize(2000,1500);
		title->Draw();		
		TPad * c3_up = new TPad("upperPad", "upperPad",0.0,0.3,1.0,1.0);
	        c3_up->Draw();

        	TPad * c3_do = new TPad("lowerPad", "lowerPad",0.0,0.0,1.0,0.3);
        	c3_do->Draw();

		c3_up->cd();
		gPad->SetLogy();
		
		//TemplatesHe[1]->Scale(100);
		TH1F * Sum = (TH1F*) TemplatesP[1]->Clone();
		Sum->Add(TemplatesD[1]);
		Sum->Add(TemplatesHe[1]);
		if(IsFitNoise) Sum->Add(TemplatesNoise[1]);

		//PlotDistribution(gPad, TemplatesP[0] ,"Reconstructed Mass [GeV/c^2]","Counts",2,"same",1e-1,Datas[0]->GetBinContent(Datas[0]->GetMaximumBin())*1.33,2,"Original Protons MC Template");
		//if(!IsSmearingCheck) PlotDistribution(gPad, TemplatesD[0] ,"Reconstructed Mass [GeV/c^2]","Counts",4,"same",1e-1,1e5,2,"Original Deuterons MC Template");
	//	if(!IsFitNoise) if(!IsSmearingCheck) PlotDistribution(gPad, TemplatesHe[0],"Reconstructed Mass [GeV/c^2]","Counts",3,"same",1e-1,1e5,2,"Original He Fragm. MC Template");
		PlotDistribution(gPad, TemplatesP[1] ,"Reconstructed Mass [GeV/c^{2}]","Weighted Counts",2,"hist,esame",1,Datas[0]->GetBinContent(Datas[0]->GetMaximumBin())*1.33,2,"Best #chi^{2} Protons MC Template");
		if(!IsSmearingCheck) PlotDistribution(gPad, TemplatesD[1] ,"Reconstructed Mass [GeV/c^{2}]","Weighted Counts",4,"hist,esame",1,1e5,2,"Best #chi^{2} Deuterons MC Template");
		if(!IsSmearingCheck) PlotDistribution(gPad, TemplatesHe[1],"Reconstructed Mass [GeV/c^{2}]","Weighted Counts",3,"hist,esame",1,1e5,2,"Best #chi^{2} Tritium MC Template");
		if(IsFitNoise) {
			TH1F * NoiseD = CreateNoiseD(TemplatesNoise[1],TemplatesP[1],TemplatesD[1]);
			PlotDistribution(gPad, TemplatesNoise[2],"Reconstructed Mass [GeV/c^{2}]","Weighted Counts",kRed-9,"esame",1e-1,1e5,1,"Noise P Template",true);
			PlotDistribution(gPad, NoiseD,"Reconstructed Mass [GeV/c^{2}]","Counts",kBlue-7,"esame",1e-1,1e5,1,"Noise D Template",true);
		}		
		
		PlotDistribution(gPad, Datas[0],"Reconstructed Mass [GeV/c^{2}]","Weighted Counts",1,"ePsame",1,1e5,1,"ISS data",false,true);
		if(Fits.size()>0) { //PlotDistribution(gPad, Fits[0],"Reconstructed Mass [Gev/c^{2}]","counts",6,"same",1,1e5,4,"Fraction Fit");
				     PlotDistribution(gPad, Sum,"Reconstructed Mass [Gev/c^{2}]","counts",1,"hist,esame",1,1e5,2,"Sum of Contributions");	
		}
		c3_up->cd();
		title->Draw();
	
		c3_do->cd();
		TH1F * Ratio = (TH1F*) Datas[0]->Clone();
		TH1F * Line = (TH1F*) Datas[0]->Clone();
		TH1F * LineUp = (TH1F*) Datas[0]->Clone();
		TH1F * LineDo = (TH1F*) Datas[0]->Clone();

		Ratio->Add(Sum,-1);
		NormalizeResiduals(Ratio);	
		CleanResiduals(Ratio,Sum,Datas[0]);
		for(int i=0;i<Line->GetNbinsX();i++){Line->SetBinContent(i+1,0); Line->SetBinError(i+1,0);}	
		for(int i=0;i<Line->GetNbinsX();i++){LineUp->SetBinContent(i+1,2); Line->SetBinError(i+1,0);}	
		for(int i=0;i<Line->GetNbinsX();i++){LineDo->SetBinContent(i+1,-2); Line->SetBinError(i+1,0);}	

		PlotDistribution(gPad, Ratio,"","",1,"ePsame",-10,10,2,"Weighted Residuals",false,true);	
		PlotDistribution(gPad, Line,"","",1,"Lsame",-10,10,4,"",false,false,true);
		PlotDistribution(gPad, LineUp,"","",1,"Lsame",-10,10,2,"",false,false,true);
		PlotDistribution(gPad, LineDo,"","",1,"Lsame",-10,10,2,"",false,false,true);
		

			Plots.Add(c3);
                Plots.writeObjsInFolder((FIT->GetName()+"/Fits/Bin"+to_string(i)).c_str());

		CollectiveView->cd(i);
		c3_up->DrawClone();

		CollectiveView2->cd(i);
	

		TCanvas * c4 = new TCanvas("ChiSquare");
		gPad->SetLogz();	
		TH2F * Chi = (TH2F*) infile->Get((FIT->GetName()+"/Fit Results/Spreads/ChiSquare/ChiSquare Bin "+to_string(i)).c_str());
                Chi->GetZaxis()->SetRangeUser(0.2,100);
                PlotTH2F(gPad, Chi, "#sigma_{smear}/#sigma_{MC} (%)","shift (%)", "colz");	


		Plots.Add(c4);
                Plots.writeObjsInFolder((FIT->GetName()+"/Fits/Bin"+to_string(i)).c_str());


	}
	 Plots.Add(CollectiveView);
         Plots.writeObjsInFolder((FIT->GetName()+"/Fits/").c_str());	
	 Plots.Add(CollectiveView2);
         Plots.writeObjsInFolder((FIT->GetName()+"/Fits/").c_str());	
	
	return;

}


void CheckSlowDown(FileSaver Plots){
	

	Binning BinningRIGD(deuton);
	Binning BinningRIGP(proton);

        BinningRIGD.setBinsFromRigidity(25, 0.8, 5,ResponseTOF,0.00347548,5.8474);
        BinningRIGP.setBinsFromRigidity(25, 0.8, 5,ResponseTOF,0.00347548,5.8474);
	
	TGraph * SlowdownRigP = new TGraph();
	TGraph * SlowdownRigD = new TGraph();

	for(int i=0;i<BinningRIGP.size();i++){
		SlowdownRigP->SetPoint(i,BinningRIGP.RigTOIBins()[i],(BinningRIGP.RigBins()[i]-BinningRIGP.RigTOIBins()[i]));
		SlowdownRigD->SetPoint(i,BinningRIGP.RigTOIBins()[i],(BinningRIGD.RigBins()[i]-BinningRIGD.RigTOIBins()[i]));
	}
	
	TCanvas * c1 = new TCanvas("Slowdown Rig");
	c1->cd();
	SlowdownRigP->SetLineColor(2);
	SlowdownRigP->SetLineWidth(4);
	SlowdownRigP->GetXaxis()->SetTitle("Rig[GV]");
	SlowdownRigD->SetLineColor(4);
	SlowdownRigD->SetLineWidth(4);
	SlowdownRigD->Draw("AL");
	SlowdownRigP->Draw("Lsame");
	

	Binning BinningBetaD(deuton);
	Binning BinningBetaP(proton);


	BinningBetaD.setBinsFromBeta(20, 0.55, 0.9,ResponseTOF,0.00347548,5.8474);
        BinningBetaP.setBinsFromBeta(20, 0.55, 0.9,ResponseTOF,0.00347548,5.8474);
	

	TGraph * SlowdownBetaP = new TGraph();
	TGraph * SlowdownBetaD = new TGraph();

	for(int i=0;i<BinningBetaP.size();i++){
		SlowdownBetaP->SetPoint(i,BinningBetaP.BetaTOIBins()[i],BinningBetaP.BetaBins()[i]-BinningBetaP.BetaTOIBins()[i]);
		SlowdownBetaD->SetPoint(i,BinningBetaP.BetaTOIBins()[i],BinningBetaD.BetaBins()[i]-BinningBetaD.BetaTOIBins()[i]);
	}
	
	TCanvas * c2 = new TCanvas("Slowdown Beta");
	c2->cd();
	SlowdownBetaP->SetLineColor(2);
	SlowdownBetaP->SetLineWidth(4);
	SlowdownBetaD->SetLineColor(4);
	SlowdownBetaD->SetLineWidth(4);
	SlowdownBetaD->GetXaxis()->SetTitle("Beta");
	
	SlowdownBetaD->Draw("AL");
	SlowdownBetaP->Draw("Lsame");
	


	Plots.Add(c1);
	Plots.Add(c2);
	Plots.writeObjsInFolder("EnergyLoss");	
}


