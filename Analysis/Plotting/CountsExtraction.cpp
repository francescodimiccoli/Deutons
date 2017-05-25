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
#include "../../include/TemplateFIT.h"

#include "../../include/FitError.h"
#include "../../include/Resolution.h"
#include "../../include/PlottingFunctions.h"

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

void DrawParameters(FileSaver finalHistos,FileSaver Plots,std::string path, Binning bins, std::string title, std::string x_udm,float xmin, float xmax,float y1min,float y1max,float y2min,float y2max);
void DrawFits(TemplateFIT * FIT,FileSaver finalHistos, FileSaver Plots);

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

	TemplateFIT * ToFfits= new TemplateFIT(finalHistos,"TOFfits",ToFDB);
	TemplateFIT * NaFfits= new TemplateFIT(finalHistos,"NaFfits",NaFDB);
	TemplateFIT * Aglfits= new TemplateFIT(finalHistos,"Aglfits",AglDB);

	DrawFits(ToFfits,finalHistos,Plots);
	DrawFits(NaFfits,finalHistos,Plots);
	DrawFits(Aglfits,finalHistos,Plots);



	cout<<"****************************** RESULTS ***************************************"<<endl;

	TCanvas * c4 = new TCanvas("Deuteron Counts");
	c4->SetCanvasSize(2000,1500);

	std::string pathresTOF   = (ToFfits->GetName() + "/Fit Results/");
	std::string pathresNaF   = (NaFfits->GetName() + "/Fit Results/");
	std::string pathresAgl   = (Aglfits->GetName() + "/Fit Results/");
		
	TH1F * DCountsTOF = (TH1F*) finalHistos.Get((pathresTOF+"Deuteron Counts").c_str());
	TH1F * PCountsTOF = (TH1F*) finalHistos.Get((pathresTOF+"Proton Counts").c_str());
	TH1F * DCountsNaF = (TH1F*) finalHistos.Get((pathresNaF+"Deuteron Counts").c_str());
	TH1F * PCountsNaF = (TH1F*) finalHistos.Get((pathresNaF+"Proton Counts").c_str());
	TH1F * DCountsAgl = (TH1F*) finalHistos.Get((pathresAgl+"Deuteron Counts").c_str());
	TH1F * PCountsAgl = (TH1F*) finalHistos.Get((pathresAgl+"Proton Counts").c_str());
	
	
	DCountsTOF->SetName("Deuteron Counts (TOF)");
	PCountsTOF->SetName("Proton Counts (TOF)");
	DCountsNaF->SetName("Deuteron Counts (NaF)");
	PCountsNaF->SetName("Proton Counts (NaF)");
	DCountsAgl->SetName("Deuteron Counts (Agl)");
	PCountsAgl->SetName("Proton Counts (Agl)");


	PlotTH1FintoGraph(gPad,ToFDB, DCountsTOF,"Kinetic Energy [GeV/nucl.]", "Counts",4,true,"ePsame",0.1,10,10,4e4);
	PlotTH1FintoGraph(gPad,NaFDB, DCountsNaF,"Kinetic Energy [GeV/nucl.]", "Counts",4,true,"ePsame",0.1,10,10,4e4);
	PlotTH1FintoGraph(gPad,AglDB, DCountsAgl,"Kinetic Energy [GeV/nucl.]", "Counts",4,true,"ePsame",0.1,10,10,4e4);

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

	PlotTH1FintoGraph(gPad,ToFDB, RatioTOF,"Kinetic Energy [GeV/nucl.]", "Counts ratio",2,true,"ePsame",0.1,10,1e-3,1e-1,"d/P Counts ratio (TOF)");
	PlotTH1FintoGraph(gPad,NaFDB, RatioNaF,"Kinetic Energy [GeV/nucl.]", "Counts ratio",2,true,"ePsame",0.1,10,1e-3,1e-1,"d/P Counts ratio (NaF)");
	PlotTH1FintoGraph(gPad,AglDB, RatioAgl,"Kinetic Energy [GeV/nucl.]", "Counts ratio",2,true,"ePsame",0.1,10,1e-3,1e-1,"d/P Counts ratio (Agl)");


	Plots.Add(c5);
	Plots.writeObjsInFolder("Results");

	
	TCanvas * c6 = new TCanvas("Uncertainty Break-down");
        c6->SetCanvasSize(5000,1000);
	c6->Divide(3,1);

	TH1F * StatErrTOF = (TH1F*) finalHistos.Get((pathresTOF+"StatError").c_str());
        TH1F * SystErrTOF = (TH1F*) finalHistos.Get((pathresTOF+"SystError").c_str());
	TH1F * HeCoErrTOF = (TH1F*) finalHistos.Get((pathresTOF+"HeContTOF_Eff").c_str());
	TH1F * StatErrNaF = (TH1F*) finalHistos.Get((pathresNaF+"StatError").c_str());
        TH1F * SystErrNaF = (TH1F*) finalHistos.Get((pathresNaF+"SystError").c_str());
	TH1F * HeCoErrNaF = (TH1F*) finalHistos.Get((pathresNaF+"HeContNaF_Eff").c_str());
	TH1F * StatErrAgl = (TH1F*) finalHistos.Get((pathresAgl+"StatError").c_str());
        TH1F * SystErrAgl = (TH1F*) finalHistos.Get((pathresAgl+"SystError").c_str());
	TH1F * HeCoErrAgl = (TH1F*) finalHistos.Get((pathresAgl+"HeContAgl_Eff").c_str());


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
	
	DrawParameters(finalHistos,Plots,pathresTOF,ToFDB,"Parameters TOF","[ps]",0.45,0.9,-100,100,-40,180);
	DrawParameters(finalHistos,Plots,pathresNaF,NaFDB,"Parameters NaF","[rad/10^{4}]",0.7,0.98,-1000,1000,-1000,2000);
	DrawParameters(finalHistos,Plots,pathresAgl,AglDB,"Parameters Agl","[rad/10^{4}]",0.95,1.005,-230,230,-150,400);




	return 0;
}

void DrawParameters(FileSaver finalHistos,FileSaver Plots,std::string path, Binning bins, std::string title, std::string x_udm, float xmin, float xmax,float y1min,float y1max,float y2min,float y2max){

	TCanvas * c7 = new TCanvas(title.c_str());
	c7->SetCanvasSize(2000,1500);

	TH1F * ShiftBest = (TH1F*) finalHistos.Get((path+"Best Fit Shift").c_str());
        TH1F * SigmaBest = (TH1F*) finalHistos.Get((path+"Best Fit Sigma").c_str());
	

	TPad * c7_up = new TPad("upperPad", "upperPad",0.0,0.5,1.0,1.0);
	c7_up->Draw();

	TPad * c7_do = new TPad("lowerPad", "lowerPad",0.0,0.0,1.0,0.5);
	c7_do->Draw();

	PlotTH1FintoGraph(c7_up,bins, ShiftBest, "Measured #beta",  ("Mean shift "+x_udm).c_str(),2,false,"ep",xmin,xmax,y1min,y1max,"Best #chi^{2} Shift");
	PlotTH1FintoGraph(c7_do,bins, SigmaBest, "Measured #beta",  ("Additive #sigma "+x_udm).c_str(),4,false,"ep",xmin,xmax,y2min,y2max,"Best #chi^{2} #sigma");
		

	Plots.Add(c7);
	Plots.writeObjsInFolder("Results");
	return;
}


void DrawFits(TemplateFIT * FIT,FileSaver finalHistos,FileSaver Plots){

	std::string pathdata  = (FIT->GetName() + "/Fit Results/Data");
	std::string pathtemplP= (FIT->GetName() + "/Fit Results/ScaledTemplatesP");
	std::string pathtemplD= (FIT->GetName() + "/Fit Results/ScaledTemplatesD");	
	std::string pathfit   = (FIT->GetName() + "/Fit Results/FractionFits");
	std::string pathres   = (FIT->GetName() + "/Fit Results/");
	
	TFile * infile = finalHistos.GetFile();	


	for(int i=1; i<FIT->GetBinning().size();i++){
		
		std::string pathbinP    = pathtemplP + "/Bin"+to_string(i);
		std::string pathbinD    = pathtemplD + "/Bin"+to_string(i);
		std::string pathbindata = pathdata   + "/Bin"+to_string(i);
		std::string pathbinfit  = pathfit    + "/Bin"+to_string(i);

		std::vector<TH1F*> TemplatesP=GetListOfTemplates(infile, pathbinP);
		std::vector<TH1F*> TemplatesD=GetListOfTemplates(infile, pathbinD);		
		std::vector<TH1F*> Datas     =GetListOfTemplates(infile, pathbindata);
		std::vector<TH1F*> Fits      =GetListOfTemplates(infile, pathbinfit);

	

		TCanvas * c1 = new TCanvas("Modified Templates");
		c1->SetCanvasSize(2000,1500);

		for(int j=TemplatesP.size()-1;j>=0;j--){
			if(j==0) PlotDistribution(gPad, TemplatesP[j],"Reconstructed Mass [GeV/c^2]","Counts",1,"same",1,TemplatesP[j]->GetBinContent(TemplatesP[j]->GetMaximumBin())*1.13,10);
			else     PlotDistribution(gPad, TemplatesP[j],"Reconstructed Mass [GeV/c^2]","Counts",colorbase + j,"same",1,TemplatesP[j]->GetBinContent(TemplatesP[j]->GetMaximumBin())*1.13,7,"",false,false,true);		
		}
	
		Plots.Add(c1);
		Plots.writeObjsInFolder((FIT->GetName()+"/Fits/Bin"+to_string(i)).c_str());	


		TCanvas * c2 = new TCanvas("Modified T. Fits");
                c2->SetCanvasSize(2000,1500);

		PlotDistribution(gPad, TemplatesP[0],"Reconstructed Mass [GeV/c^2]","Counts",2,"same",1,Datas[0]->GetBinContent(Datas[0]->GetMaximumBin())*1.13,10,"Original Protons MC Template");
		PlotDistribution(gPad, TemplatesD[0],"Reconstructed Mass [GeV/c^2]","Counts",4,"same",1,1e5,10,"Original Deuterons MC Template");
		for(int j=TemplatesP.size()-1;j>=1;j--){
                        PlotDistribution(gPad, TemplatesP[j],"Reconstructed Mass [GeV/c^2]","Counts",2,"same",1,1e5,1,"",false,false,true);
			PlotDistribution(gPad, TemplatesD[j],"Reconstructed Mass [GeV/c^2]","Counts",4,"same",1,1e5,1,"",false,false,true);
                }
		PlotDistribution(gPad, Datas[0],"Reconstructed Mass [GeV/c^2]","Counts",1,"ePsame",1,1e5,3,"ISS data",false,true);
	


		Plots.Add(c2);
                Plots.writeObjsInFolder((FIT->GetName()+"/Fits/Bin"+to_string(i)).c_str());
		
		
		TCanvas * c3 = new TCanvas("Template Fits");
                c3->SetCanvasSize(2000,1500);

		PlotDistribution(gPad, TemplatesP[0],"Reconstructed Mass [GeV/c^2]","Counts",2,"same",1,Datas[0]->GetBinContent(Datas[0]->GetMaximumBin())*1.13,2,"Original Protons MC Template");
		PlotDistribution(gPad, TemplatesD[0],"Reconstructed Mass [GeV/c^2]","Counts",4,"same",1,1e5,2,"Original Deuterons MC Template");
		PlotDistribution(gPad, TemplatesP[1],"Reconstructed Mass [GeV/c^2]","Counts",2,"same",1,Datas[0]->GetBinContent(Datas[0]->GetMaximumBin())*1.13,10,"Best #chi^{2} Protons MC Template");
		PlotDistribution(gPad, TemplatesD[1],"Reconstructed Mass [GeV/c^2]","Counts",4,"same",1,1e5,10,"Best #chi^{2} Deuterons MC Template");
		
		if(Fits.size()>0) PlotDistribution(gPad, Fits[0],"Reconstructed Mass [GeV/c^2]","Counts",6,"same",1,1e5,4,"Fraction Fit");
		PlotDistribution(gPad, Datas[0],"Reconstructed Mass [GeV/c^2]","Counts",1,"ePsame",1,1e5,3,"ISS data",false,true);


		Plots.Add(c3);
                Plots.writeObjsInFolder((FIT->GetName()+"/Fits/Bin"+to_string(i)).c_str());
	

		TCanvas * c4 = new TCanvas("ChiSquare");
                c4->SetCanvasSize(2000,1500);
		gPad->SetLogz();	
		TH2F * Chi = (TH2F*) infile->Get((FIT->GetName()+"/Fit Results/Spreads/ChiSquare/ChiSquare Bin "+to_string(i)).c_str());
                Chi->GetZaxis()->SetRangeUser(0.5,100);
                PlotTH2F(gPad, Chi, "#sigma deformation (%)","Mean shift (%)", "colz");	


		Plots.Add(c4);
                Plots.writeObjsInFolder((FIT->GetName()+"/Fits/Bin"+to_string(i)).c_str());
	

		

	

	}

	return;

}
