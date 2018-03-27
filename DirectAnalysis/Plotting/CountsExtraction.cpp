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

#include "../include/Variables.hpp"
#include "../include/Cuts.h"
#include "../include/filesaver.h"
#include "../include/TemplateFITbetasmear.h"

#include "../include/FitError.h"
#include "../include/PlottingFunctions.h"
#include <sstream>

std::string Convert (float number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();   
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

        cout<<"****************************** PLOTTING FITS ***************************************"<<endl;

//	TemplateFIT * SmearingCheck = new TemplateFIT(finalHistos,"SmearingCheck",PRB);
	TemplateFIT * ToFfits= new TemplateFIT(finalHistos,"TOFfits",ToFDB);
	TemplateFIT * NaFfits= new TemplateFIT(finalHistos,"NaFfits",NaFDB);
	TemplateFIT * Aglfits= new TemplateFIT(finalHistos,"Aglfits",AglDB);

//	DrawFits(SmearingCheck,finalHistos,Plots,false,true);
	DrawFits(ToFfits,finalHistos,Plots);
	DrawFits(NaFfits,finalHistos,Plots);
	DrawFits(Aglfits,finalHistos,Plots);
//


	cout<<"****************************** RESULTS ***************************************"<<endl;

	TCanvas * c4 = new TCanvas("Deuteron Counts");
	c4->SetCanvasSize(2000,1500);

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
	
	TH1F * DCountsPrimTOF = (TH1F*) finalHistos.Get((pathresTOF+"Primary Deuteron Counts").c_str());
	TH1F * DCountsPrimNaF = (TH1F*) finalHistos.Get((pathresNaF+"Primary Deuteron Counts").c_str());
	TH1F * DCountsPrimAgl = (TH1F*) finalHistos.Get((pathresAgl+"Primary Deuteron Counts").c_str());
	
	PlotTH1FintoGraph(gPad,ToFDB, DCountsTOF,"Kinetic Energy [GeV/nucl.]", "Counts",4,true,"Psame",0.1,10,10,7e4,"Deuteron Counts (TOF)",8);
	PlotTH1FintoGraph(gPad,NaFDB, DCountsNaF,"Kinetic Energy [GeV/nucl.]", "Counts",4,true,"Psame",0.1,10,10,7e4,"Deuteron Counts (NaF)",22);
	PlotTH1FintoGraph(gPad,AglDB, DCountsAgl,"Kinetic Energy [GeV/nucl.]", "Counts",4,true,"Psame",0.1,10,10,7e4,"Deuteron Counts (Agl)",29);

	PlotTH1FintoGraph(gPad,ToFDB, DCountsPrimTOF,"Kinetic Energy [GeV/nucl.]", "Counts",4,true,"Psame",0.1,10,10,7e4,"Primary Counts (TOF)",4);
	PlotTH1FintoGraph(gPad,NaFDB, DCountsPrimNaF,"Kinetic Energy [GeV/nucl.]", "Counts",4,true,"Psame",0.1,10,10,7e4,"Primary Counts (NaF)",26);
	PlotTH1FintoGraph(gPad,AglDB, DCountsPrimAgl,"Kinetic Energy [GeV/nucl.]", "Counts",4,true,"Psame",0.1,10,10,7e4,"Primary Counts (Agl)",30);


	Plots.Add(c4);
	Plots.writeObjsInFolder("Results");


	TCanvas * c5 = new TCanvas("D/P Raw Counts ratio");
        c5->SetCanvasSize(2000,1500);
	
	TH1F * RatioTOF = (TH1F*)DCountsTOF->Clone();
	RatioTOF->Divide(PCountsTOF);
	TH1F * RatioNaF = (TH1F*)DCountsNaF->Clone();
	RatioNaF->Divide(PCountsNaF);
	TH1F * RatioAgl = (TH1F*)DCountsAgl->Clone();
	RatioAgl->Divide(PCountsAgl);

	PlotTH1FintoGraph(gPad,ToFDB, RatioTOF,"Kinetic Energy [GeV/nucl.]", "Counts ratio",2,true,"Psame",0.1,10,1e-3,1e-1,"d/P Counts ratio (TOF)",8);
	PlotTH1FintoGraph(gPad,NaFDB, RatioNaF,"Kinetic Energy [GeV/nucl.]", "Counts ratio",2,true,"Psame",0.1,10,1e-3,1e-1,"d/P Counts ratio (NaF)",22);
	PlotTH1FintoGraph(gPad,AglDB, RatioAgl,"Kinetic Energy [GeV/nucl.]", "Counts ratio",2,true,"Psame",0.1,10,1e-3,1e-1,"d/P Counts ratio (Agl)",29);


	Plots.Add(c5);
	Plots.writeObjsInFolder("Results");

	
	TCanvas * c6 = new TCanvas("Uncertainty Break-down");
        c6->SetCanvasSize(5000,1000);
	c6->Divide(3,1);

	TH1F * StatErrTOF = (TH1F*) finalHistos.Get((pathresTOF+"StatError").c_str());
        TH1F * SystErrTOF = (TH1F*) finalHistos.Get((pathresTOF+"SystError").c_str());
	TH1F * HeCoErrTOF = (TH1F*) finalHistos.Get((pathresTOF+"HeliumContamination/HeContError").c_str());
	TH1F * StatErrNaF = (TH1F*) finalHistos.Get((pathresNaF+"StatError").c_str());
        TH1F * SystErrNaF = (TH1F*) finalHistos.Get((pathresNaF+"SystError").c_str());
	TH1F * HeCoErrNaF = (TH1F*) finalHistos.Get((pathresNaF+"HeliumContamination/HeContError").c_str());
	TH1F * StatErrAgl = (TH1F*) finalHistos.Get((pathresAgl+"StatError").c_str());
        TH1F * SystErrAgl = (TH1F*) finalHistos.Get((pathresAgl+"SystError").c_str());
	TH1F * HeCoErrAgl = (TH1F*) finalHistos.Get((pathresAgl+"HeliumContamination/HeContError").c_str());
	
	TH1F * TotErrTOF  = (TH1F *)StatErrTOF->Clone();
	TH1F * TotErrNaF  = (TH1F *)StatErrNaF->Clone();
	TH1F * TotErrAgl  = (TH1F *)StatErrAgl->Clone();
	
	for(int bin=0;bin<DCountsTOF->GetNbinsX();bin++) {
		if(DCountsTOF->GetBinContent(bin+1)>0){
			TotErrTOF->SetBinContent(bin+1,DCountsTOF->GetBinError(bin+1)/DCountsTOF->GetBinContent(bin+1));
			TotErrTOF->SetBinError(bin+1,0);
		}
		if(DCountsNaF->GetBinContent(bin+1)>0){	
			TotErrNaF->SetBinContent(bin+1,DCountsNaF->GetBinError(bin+1)/DCountsNaF->GetBinContent(bin+1));
			TotErrNaF->SetBinError(bin+1,0);
		}
		if(DCountsAgl->GetBinContent(bin+1)>0){	
			TotErrAgl->SetBinContent(bin+1,DCountsAgl->GetBinError(bin+1)/DCountsAgl->GetBinContent(bin+1));
			TotErrAgl->SetBinError(bin+1,0);
		}	

	}	


	c6->cd(1);
	PlotDistribution(gPad,TotErrTOF ,"TOF Range Bin","Relative error",2,"same",1e-4,1.1,10,"T. Fit Total Error");
	PlotDistribution(gPad,SystErrTOF,"TOF Range Bin","Relative error",4,"same",1e-4,1.1,4,"T. Fit Systematic Error");
	PlotDistribution(gPad,StatErrTOF,"TOF Range Bin","Relative error",1,"same",1e-4,1.1,4,"T. Fit Statistical Error");
	PlotDistribution(gPad,HeCoErrTOF,"TOF Range Bin","Relative error",3,"same",1e-4,1.1,4,"Helium Fragm. Error");
	c6->cd(2);
	PlotDistribution(gPad,TotErrNaF ,"NaF Range Bin","Relative error",2,"same",1e-4,1.1,10,"T. Fit Total Error");
	PlotDistribution(gPad,SystErrNaF,"NaF Range Bin","Relative error",4,"same",1e-4,1.1,4,"T. Fit Systematic Error");
	PlotDistribution(gPad,StatErrNaF,"NaF Range Bin","Relative error",1,"same",1e-4,1.1,4,"T. Fit Statistical Error");
	PlotDistribution(gPad,HeCoErrNaF,"NaF Range Bin","Relative error",3,"same",1e-4,1.1,4,"Helium Fragm. Error");
	c6->cd(3);
	PlotDistribution(gPad,TotErrAgl ,"Agl Range Bin","Relative error",2,"same",1e-4,1.1,10,"T. Fit Total Error");
	PlotDistribution(gPad,SystErrAgl,"Agl Range Bin","Relative error",4,"same",1e-4,1.1,4,"T. Fit Systematic Error");
	PlotDistribution(gPad,StatErrAgl,"Agl Range Bin","Relative error",1,"same",1e-4,1.1,4,"T. Fit Statistical Error");
	PlotDistribution(gPad,HeCoErrAgl,"Agl Range Bin","Relative error",3,"same",1e-4,1.1,4,"Helium Fragm. Error");

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



	
	DrawParameters(finalHistos,Plots,pathresTOF,ToFDB,"Parameters TOF","Measured #beta","[ps]",0.45,0.9,-100,100,-40,180);
	DrawParameters(finalHistos,Plots,pathresNaF,NaFDB,"Parameters NaF","Measured #beta","[rad/10^{4}]",0.7,0.98,-1000,1000,-1000,2000);
	DrawParameters(finalHistos,Plots,pathresAgl,AglDB,"Parameters Agl","Measured #beta","[rad/10^{4}]",0.95,1.005,-230,230,-150,400);


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
	c7->SetCanvasSize(2000,1500);

	TH1F * ShiftBest = (TH1F*) finalHistos.Get((path+"Best Fit Shift").c_str());
        TH1F * SigmaBest = (TH1F*) finalHistos.Get((path+"Best Fit Sigma").c_str());
	

	TPad * c7_up = new TPad("upperPad", "upperPad",0.0,0.5,1.0,1.0);
	c7_up->Draw();

	TPad * c7_do = new TPad("lowerPad", "lowerPad",0.0,0.0,1.0,0.5);
	c7_do->Draw();

	PlotTH1FintoGraph(c7_up,bins, ShiftBest, x.c_str(),  ("Mean shift "+x_udm).c_str(),2,false,"ep",xmin,xmax,y1min,y1max,"Best #chi^{2} Shift");
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
	std::string pathtrans = (FIT->GetName() + "/Fit Results/TrasnferFunctions");
	std::string pathres   = (FIT->GetName() + "/Fit Results/");
	
	TFile * infile = finalHistos.GetFile();	

	TCanvas * CollectiveView = new TCanvas("Summary 1 ");
	CollectiveView->Divide(FIT->GetBinning().size()/3,3);
	CollectiveView->SetCanvasSize(3000,2000);

	TCanvas * CollectiveView2 = new TCanvas("Summary 2 ");
	CollectiveView2->Divide(FIT->GetBinning().size()/3,3);
	CollectiveView2->SetCanvasSize(3000,2000);



	for(int i=1; i<FIT->GetBinning().size();i++){
		
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
		std::vector<TH1F*> Transfer  =GetListOfTemplates(infile, pathtrans);

		//cout<<pathbinHe<<" "<<TemplatesHe.size()<<" "<<TemplatesHe[0]<<endl;
	

		TCanvas * c1 = new TCanvas("Modified Templates");
		c1->SetCanvasSize(2000,1500);
		std::string bintitle = ("Kin. Energy: " + Convert(FIT->GetBinning().EkPerMassBinCent(i)) + " GeV/n");
		TPaveLabel* title = new TPaveLabel(0.35,0.94,0.6,0.97,bintitle.c_str(),"brndc");
		title->SetFillColor(0);
		title->SetLineColor(0);
		title->SetShadowColor(0);
		title->SetTextColor(1);
		title->SetTextFont(22);
		title->SetTextSize(2.5);

		for(int j=TemplatesP.size()-1;j>=0;j--){
			if(j==0) PlotDistribution(gPad, TemplatesP[j],"Reconstructed Mass [GeV/c^2]","Counts",1,"same",1,TemplatesP[j]->GetBinContent(TemplatesP[j]->GetMaximumBin())*1.13,10);
			else     PlotDistribution(gPad, TemplatesP[j],"Reconstructed Mass [GeV/c^2]","Counts",colorbase + j,"same",1,TemplatesP[j]->GetBinContent(TemplatesP[j]->GetMaximumBin())*1.13,7,"",false,false,true);		
		}
		c1->cd();
		
		Plots.Add(c1);
		Plots.writeObjsInFolder((FIT->GetName()+"/Fits/Bin"+to_string(i)).c_str());	


		TCanvas * c2 = new TCanvas("Modified T. Fits");
                c2->SetCanvasSize(2000,1500);

		PlotDistribution(gPad, TemplatesP[0], "Reconstructed Mass [GeV/c^2]","Counts",2,"same",1,Datas[0]->GetBinContent(Datas[0]->GetMaximumBin())*1.13,10,"Original Protons MC Template");
		if(!IsSmearingCheck) PlotDistribution(gPad, TemplatesD[0], "Reconstructed Mass [GeV/c^2]","Counts",4,"same",1,1e5,10,"Original Deuterons MC Template");
		if(!IsSmearingCheck) PlotDistribution(gPad, TemplatesHe[0],"Reconstructed Mass [GeV/c^2]","Counts",3,"same",1,1e5,10,"Original He Fragm. MC Template");

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
		PlotDistribution(gPad, TemplatesP[1] ,"Reconstructed Mass [GeV/c^{2}]","Weighted Counts",2,"esame",8,Datas[0]->GetBinContent(Datas[0]->GetMaximumBin())*1.33,10,"Best #chi^{2} Protons MC Template");
		if(!IsSmearingCheck) PlotDistribution(gPad, TemplatesD[1] ,"Reconstructed Mass [GeV/c^{2}]","Weighted Counts",4,"esame",8,1e5,10,"Best #chi^{2} Deuterons MC Template");
		if(!IsSmearingCheck) PlotDistribution(gPad, TemplatesHe[1],"Reconstructed Mass [GeV/c^{2}]","Weighted Counts",3,"esame",8,1e5,10,"Best #chi^{2} Tritium MC Template");
		if(IsFitNoise) {
			TH1F * NoiseD = CreateNoiseD(TemplatesNoise[1],TemplatesP[1],TemplatesD[1]);
			PlotDistribution(gPad, TemplatesNoise[1],"Reconstructed Mass [GeV/c^{2}]","Weighted Counts",kRed-9,"esame",1e-1,1e5,10,"Noise P Template",true);
			PlotDistribution(gPad, NoiseD,"Reconstructed Mass [GeV/c^{2}]","Counts",kBlue-7,"esame",1e-1,1e5,10,"Noise D Template",true);
		}		
		
		PlotDistribution(gPad, Datas[0],"Reconstructed Mass [GeV/c^{2}]","Weighted Counts",1,"ePsame",8,1e5,3,"ISS data",false,true);
		if(Fits.size()>0) { PlotDistribution(gPad, Fits[0],"Reconstructed Mass [Gev/c^{2}]","counts",6,"same",1,1e5,4,"Fraction Fit");
				     PlotDistribution(gPad, Sum,"Reconstructed Mass [Gev/c^{2}]","counts",1,"esame",8,1e5,6,"Sum of Contributions");	
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
		PlotDistribution(gPad, Line,"","",1,"Lsame",-10,10,7,"",false,false,true);
		PlotDistribution(gPad, LineUp,"","",1,"Lsame",-10,10,2,"",false,false,true);
		PlotDistribution(gPad, LineDo,"","",1,"Lsame",-10,10,2,"",false,false,true);
		

			Plots.Add(c3);
                Plots.writeObjsInFolder((FIT->GetName()+"/Fits/Bin"+to_string(i)).c_str());

		CollectiveView->cd(i);
		c3_up->DrawClone();

		CollectiveView2->cd(i);
		//PlotDistribution(gPad, Ratio,"","",1,"ePsame",-10,10,2,"",false,true,true);	
		//PlotDistribution(gPad, Line,"","",1,"Lsame",-10,10,7,"",false,false,true);
	
		TCanvas * c5 = new TCanvas("Transfer Functions");
                c5->SetCanvasSize(2000,1500);

		for(int j=0;j<Transfer.size();j++){
			Transfer[j]->Smooth(3);
			PlotDistribution(gPad, Transfer[j],"Reconstructed Mass [GeV/c^{2}]","Prim. / (Prim. + Sec.)",j,"hist,same",0,1,7,("Bin. "+to_string(j)).c_str(),false,false);
		}

		Plots.Add(c5);
                Plots.writeObjsInFolder((FIT->GetName()+"/Fits").c_str());

	

		TCanvas * c4 = new TCanvas("ChiSquare");
                c4->SetCanvasSize(2000,1500);
		gPad->SetLogz();	
		TH2F * Chi = (TH2F*) infile->Get((FIT->GetName()+"/Fit Results/Spreads/ChiSquare/ChiSquare Bin "+to_string(i)).c_str());
                Chi->GetZaxis()->SetRangeUser(0.2,100);
                PlotTH2F(gPad, Chi, "Additive #sigma","Mean shift", "colz");	


		Plots.Add(c4);
                Plots.writeObjsInFolder((FIT->GetName()+"/Fits/Bin"+to_string(i)).c_str());

		TCanvas * c6 = new TCanvas("OverCutoff Events");
                c6->SetCanvasSize(2000,1500);
		gPad->SetLogy();
	
		TH1F * OverCutoffP = (TH1F *) TemplatesP[1]->Clone();
		OverCutoffP->Multiply(Transfer[i]);
		TH1F * OverCutoffD = (TH1F *) TemplatesD[1]->Clone();
		OverCutoffD->Multiply(Transfer[i]);
		TH1F * OverCutoffHe = (TH1F *) TemplatesHe[1]->Clone();
		OverCutoffHe->Multiply(Transfer[i]);
		
		TH1F * NoCutoffP = (TH1F *) TemplatesP[1]->Clone();

		NoCutoffP->Scale(
			OverCutoffP->GetBinContent(OverCutoffP->GetMaximumBin())/
			NoCutoffP->GetBinContent(NoCutoffP->GetMaximumBin()) );

		PlotDistribution(gPad, NoCutoffP,"Reconstructed Mass [GeV/c^2]","Primary Counts",2,"same",1e-1,Datas[1]->GetBinContent(Datas[1]->GetMaximumBin())*1.13,3,"Best #chi^{2} Protons MC Template");
		PlotDistribution(gPad, OverCutoffP,"Reconstructed Mass [GeV/c^2]","Counts",2,"same",1e-1,Datas[0]->GetBinContent(Datas[0]->GetMaximumBin())*1.13,10,"Best #chi^{2} Protons MC (Cutoff filtered)");
		if(!IsSmearingCheck) PlotDistribution(gPad, OverCutoffD,"Reconstructed Mass [GeV/c^2]","Counts",4,"same",1e-1,Datas[0]->GetBinContent(Datas[0]->GetMaximumBin())*1.13,10,"Best #chi^{2} Deutons MC (Cutoff filtered)");
		if(!IsSmearingCheck) 
			if(!IsFitNoise) PlotDistribution(gPad, OverCutoffHe,"Reconstructed Mass [GeV/c^2]","Counts",3,"same",1e-1,Datas[0]->GetBinContent(Datas[0]->GetMaximumBin())*1.13,10,"Best #chi^{2} Deutons MC (Cutoff filtered)");
			else {
                        TH1F * NoiseD = CreateNoiseD(OverCutoffHe,OverCutoffP,OverCutoffD);
                        PlotDistribution(gPad, OverCutoffHe,"Reconstructed Mass [GeV/c^2]","Counts",kRed-9,"same",1e-1,1e5,10,"Noise P Template",true);
                        PlotDistribution(gPad, NoiseD,"Reconstructed Mass [GeV/c^2]","Counts",kBlue-7,"same",1e-1,1e5,10,"Noise D Template",true);
                }
		
		PlotDistribution(gPad, Datas[1],"Reconstructed Mass [GeV/c^2]","Primary Counts",1,"ePsame",1e-1,Datas[1]->GetBinContent(Datas[1]->GetMaximumBin())*1.13,3,"ISS data",false,true);

		Plots.Add(c6);
                Plots.writeObjsInFolder((FIT->GetName()+"/Fits/Bin"+to_string(i)).c_str());


	}
	 Plots.Add(CollectiveView);
         Plots.writeObjsInFolder((FIT->GetName()+"/Fits/").c_str());	
	 Plots.Add(CollectiveView2);
         Plots.writeObjsInFolder((FIT->GetName()+"/Fits/").c_str());	
	
	return;

}
