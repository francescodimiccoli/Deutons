#ifndef PLOTTER_H
#define PLOTTER_H

#include "InputFileReader.h"
#include "DBarReader.h"

#include "filesaver.h"
#include "binning.h"
#include "Globals.h"
#include "Variables.hpp"
#include "ParallelFiller.h"

#include "HistoBooker.h"
#include "PlottingFunctions.h"
#include "MacrosForPlots.h"

class Plotter{

	typedef void (*PlottingFunction) (FileSaver finalhistos, FileSaver finalresults);

	
	private:
	HistoBooker Booker;
	HistoBooker BookerMC;
	FileSaver finalHistos;
	FileSaver finalResults;
        std::vector<PlottingFunction>  plottingfunctions;
 
	public:
	Plotter(FileSaver finalhistos,FileSaver finalresults) {finalHistos = finalhistos; finalResults = finalresults;}
	void FillAllAnalyses(TChain * chainDT,TChain *chainMC);
	void BookMassAnalysis();
	void BookCleaningCutsAnalysis();
	void BookRichBDTAnalysis();
	void BookBetaResMatrixAnalysis();
	void BookGenAcceptanceAnalysis();
	void BookTrackingEfficiencyAnalysis();
	void DoAllAnalyses();
};

void Plotter::FillAllAnalyses(TChain * chainDT,TChain *chainMC){
	Booker.FillEverything(DBarReader(chainDT, false));
	Booker.SaveEverything(finalHistos);
	
	BookerMC.FillEverything(DBarReader(chainMC, true));
        BookerMC.SaveEverything(finalHistos);
}

void Plotter::BookMassAnalysis(){
	Booker.BookSingleHisto("MassTOF_noSel" ,100,0,5,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&TofBetaSafetyCut",GetRecMassTOF);
        Booker.BookSingleHisto("MassTOF_Sel"   ,100,0,5,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&TofBetaSafetyCut",GetRecMassTOF);
        Booker.BookSingleHisto("MassTOF_Qual"  ,100,0,5,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime&TofBetaSafetyCut",GetRecMassTOF);
        Booker.BookSingleHisto("MassNaF_noSel" ,100,0,5,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&NafBetaSafetyCut&IsFromNaF_nosel",GetRecMassRICH);
        Booker.BookSingleHisto("MassNaF_Sel"   ,100,0,5,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&NafBetaSafetyCut&IsFromNaF",GetRecMassRICH);
        Booker.BookSingleHisto("MassNaF_Qual"  ,100,0,5,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&NafBetaSafetyCut&IsFromNaF&RICHBDTCut",GetRecMassRICH);
        Booker.BookSingleHisto("MassAgl_noSel" ,100,0,5,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&AglBetaSafetyCut&IsFromAgl_nosel",GetRecMassRICH);
        Booker.BookSingleHisto("MassAgl_Sel"   ,100,0,5,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&AglBetaSafetyCut&IsFromAgl",GetRecMassRICH);
        Booker.BookSingleHisto("MassAgl_Qual"  ,100,0,5,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&AglBetaSafetyCut&IsFromAgl&RICHBDTCut",GetRecMassRICH);

	BookerMC.BookSingleHisto("MassNaF_noSel_MCP",100,0,5,"IsPurePMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&NafBetaSafetyCut&IsFromNaF_nosel",GetRecMassRICH);
        BookerMC.BookSingleHisto("MassNaF_Sel_MCP",100,0,5,"IsPurePMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&NafBetaSafetyCut&IsFromNaF",GetRecMassRICH);
        BookerMC.BookSingleHisto("MassNaF_Qual_MCP",100,0,5,"IsPurePMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&NafBetaSafetyCut&IsFromNaF&RICHBDTCut",GetRecMassRICH);
        BookerMC.BookSingleHisto("MassAgl_noSel_MCP",100,0,5,"IsPurePMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&AglBetaSafetyCut&IsFromAgl_nosel",GetRecMassRICH);
        BookerMC.BookSingleHisto("MassAgl_Sel_MCP",100,0,5,"IsPurePMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&AglBetaSafetyCut&IsFromAgl",GetRecMassRICH);
        BookerMC.BookSingleHisto("MassAgl_Qual_MCP",100,0,5,"IsPurePMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&AglBetaSafetyCut&IsFromAgl&RICHBDTCut",GetRecMassRICH);
        BookerMC.BookSingleHisto("MassNaF_noSel_MCD",100,0,5,"IsPureDMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&NafBetaSafetyCut&IsFromNaF_nosel",GetRecMassRICH);
        BookerMC.BookSingleHisto("MassNaF_Sel_MCD",100,0,5,"IsPureDMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&NafBetaSafetyCut&IsFromNaF",GetRecMassRICH);
        BookerMC.BookSingleHisto("MassNaF_Qual_MCD",100,0,5,"IsPureDMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&NafBetaSafetyCut&IsFromNaF&RICHBDTCut",GetRecMassRICH);
        BookerMC.BookSingleHisto("MassAgl_noSel_MCD",100,0,5,"IsPureDMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&AglBetaSafetyCut&IsFromAgl_nosel",GetRecMassRICH);
        BookerMC.BookSingleHisto("MassAgl_Sel_MCD",100,0,5,"IsPureDMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&AglBetaSafetyCut&IsFromAgl",GetRecMassRICH);
        BookerMC.BookSingleHisto("MassAgl_Qual_MCD",100,0,5,"IsPureDMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&AglBetaSafetyCut&IsFromAgl&RICHBDTCut",GetRecMassRICH);
        BookerMC.BookSingleHisto("MassAgl_SelNoPID_MCD",100,0,5,"IsDeutonMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&AglBetaSafetyCut&IsFromAgl",GetRecMassRICH);
        BookerMC.BookSingleHisto("MassAgl_SelNoPID_MCP",100,0,5,"IsProtonMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&&AglBetaSafetyCut&IsFromAgl",GetRecMassRICH);

	plottingfunctions.push_back(DrawMasses);
}

void Plotter::BookCleaningCutsAnalysis(){
	Booker.BookSingleScatter("MassTOFvsRupdown",100,0,1,50,1,5,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning",GetRupdown,GetRecMassTOF);
        Booker.BookSingleScatter("RvsChisquare_x",100,0,20,100,0,35,"IsPhysTrig&IsMinimumBias&IsLooseCharge1",GetChisquareX,GetRigidity);
        Booker.BookSingleScatter("RvsChisquare_y",100,0,20,100,0,35,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX",GetChisquareY,GetRigidity);
	Booker.BookSingleScatter("UtofQvsR",100,0,35,100,0,3,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsGoodChiSquareY",GetRigidity,GetUtofQ);
	Booker.BookSingleScatter("LtofQvsR",100,0,35,100,0,3,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsGoodChiSquareY",GetRigidity,GetLtofQ);
	Booker.BookSingleScatter("InnerQvsR",100,0,35,100,0,3,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsGoodChiSquareY",GetRigidity,GetInnerQ);

	BookerMC.BookSingleScatter("RvsChisquare_x_MCP",100,0,20,100,0,35,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsPurePMC",GetChisquareX,GetRigidity);
        BookerMC.BookSingleScatter("RvsChisquare_y_MCP",100,0,20,100,0,35,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsPurePMC",GetChisquareY,GetRigidity);
	BookerMC.BookSingleScatter("RvsChisquare_x_MCD",100,0,20,100,0,35,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsPureDMC",GetChisquareX,GetRigidity);
        BookerMC.BookSingleScatter("RvsChisquare_y_MCD",100,0,20,100,0,35,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsPureDMC",GetChisquareY,GetRigidity);

	BookerMC.BookSingleScatter("UtofQvsR_MCP",100,0,35,100,0,3,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsGoodChiSquareY&IsPurePMC",GetRigidity,GetUtofQ);
	BookerMC.BookSingleScatter("LtofQvsR_MCP",100,0,35,100,0,3,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsGoodChiSquareY&IsPurePMC",GetRigidity,GetLtofQ);
	BookerMC.BookSingleScatter("InnerQvsR_MCP",100,0,35,100,0,3,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsGoodChiSquareY&IsPurePMC",GetRigidity,GetInnerQ);

	BookerMC.BookSingleScatter("UtofQvsR_MCD",100,0,35,100,0,3,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsGoodChiSquareY&IsPureDMC",GetRigidity,GetUtofQ);
	BookerMC.BookSingleScatter("LtofQvsR_MCD",100,0,35,100,0,3,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsGoodChiSquareY&IsPureDMC",GetRigidity,GetLtofQ);
	BookerMC.BookSingleScatter("InnerQvsR_MCD",100,0,35,100,0,3,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsGoodChiSquareY&IsPureDMC",GetRigidity,GetInnerQ);

	BookerMC.BookSingleScatter("UtofQvsBeta_MCP",100,0,1,100,0,3,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsGoodChiSquareY&IsPurePMC",GetBetaTOF,GetUtofQ);
	BookerMC.BookSingleScatter("LtofQvsBeta_MCP",100,0,1,100,0,3,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsGoodChiSquareY&IsPurePMC",GetBetaTOF,GetLtofQ);
	BookerMC.BookSingleScatter("InnerQvsBeta_MCP",100,0,1,100,0,3,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsGoodChiSquareY&IsPurePMC",GetBetaTOF,GetInnerQ);

	BookerMC.BookSingleScatter("UtofQvsBeta_MCD",100,0,1,100,0,3,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsGoodChiSquareY&IsPureDMC",GetBetaTOF,GetUtofQ);
	BookerMC.BookSingleScatter("LtofQvsBeta_MCD",100,0,1,100,0,3,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsGoodChiSquareY&IsPureDMC",GetBetaTOF,GetLtofQ);
	BookerMC.BookSingleScatter("InnerQvsBeta_MCD",100,0,1,100,0,3,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsGoodChiSquareX&IsGoodChiSquareY&IsPureDMC",GetBetaTOF,GetInnerQ);

	plottingfunctions.push_back(DrawCleaning);	
}

void Plotter::BookRichBDTAnalysis(){
	BookerMC.BookSingleScatter(	"RICHBDTvsMassNaFP",100,0,4.5,100,-0.5,0.6,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&NafBetaSafetyCut&IsurePMC",GetRecMassRICH,GetRICHBDT,GetRICHBDT);
        BookerMC.BookSingleScatter(	"RICHBDTvsMassAglP",100,0,4.5,100,-0.5,0.6,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&AglBetaSafetyCut&IsPurePMC",GetRecMassRICH,GetRICHBDT,GetRICHBDT);
        BookerMC.BookSingleScatter(	"RICHBDTvsMassNaFD",100,0,4.5,100,-0.5,0.6,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&NafBetaSafetyCut&IsPureDMC",GetRecMassRICH,GetRICHBDT,GetRICHBDT);
        BookerMC.BookSingleScatter(	"RICHBDTvsMassAglD",100,0,4.5,100,-0.5,0.6,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&AglBetaSafetyCut&IsPureDMC",GetRecMassRICH,GetRICHBDT,GetRICHBDT);
	BookerMC.BookSingleHisto(	"BDTDiscrNaF_MCP",100,-0.5,1.,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&NafBetaSafetyCut&IsPurePMC",GetRICHBDT);
        BookerMC.BookSingleHisto(	"BDTDiscrNaF_MCD",100,-0.5,1.,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&NafBetaSafetyCut&IsPureDMC",GetRICHBDT);
	BookerMC.BookSingleHisto(	"BDTDiscrAgl_MCP",100,-0.5,1.,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&AglBetaSafetyCut&IsPurePMC",GetRICHBDT);
        BookerMC.BookSingleHisto(	"BDTDiscrAgl_MCD",100,-0.5,1.,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&AglBetaSafetyCut&IsPureDMC",GetRICHBDT);

	Booker.BookSingleScatter(	"RICHBDTvsBetaNaF",100,0.8,1.1,100,-0.5,0.6,"IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsFromNaF",GetBetaRICH,GetRICHBDT,GetRICHBDT);
        Booker.BookSingleScatter(	"RICHBDTvsBetaAgl",100,0.95,1.05,100,-0.5,0.6,"IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsFromAgl",GetBetaRICH,GetRICHBDT,GetRICHBDT);
	Booker.BookSingleHisto(		"BDTDiscr",100,-0.5,1.,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&NafBetaSafetyCut",GetRICHBDT);
	Booker.BookSingleHisto(		"BDTDiscr",100,-0.5,1.,"IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&AglBetaSafetyCut",GetRICHBDT);
	
	plottingfunctions.push_back(DrawBDT);

}

void Plotter::BookBetaResMatrixAnalysis(){
	
	BookerMC.BookSingleScatter("BetagenvsBetaMeasTOFP",300,0.3,1,300,0.3,1,"IsPurePMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime"          ,GetBetaGen,GetBetaTOF,GetBetaTOF);
        BookerMC.BookSingleScatter("BetagenvsBetaMeasNaFP",300,0.7,1,300,0.7,1,"IsPurePMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",GetBetaGen,GetBetaRICH,GetBetaRICH);
        BookerMC.BookSingleScatter("BetagenvsBetaMeasAglP",300,0.9,1,300,0.9,1,"IsPurePMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",GetBetaGen,GetBetaRICH,GetBetaRICH);
        BookerMC.BookSingleScatter("BetagenvsBetaMeasTOFD",300,0.3,1,300,0.3,1,"IsPureDMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime"                    ,GetBetaGen,GetBetaTOF,GetBetaTOF);
        BookerMC.BookSingleScatter("BetagenvsBetaMeasNaFD",300,0.7,1,300,0.7,1,"IsPureDMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"          ,GetBetaGen,GetBetaRICH,GetBetaRICH);
        BookerMC.BookSingleScatter("BetagenvsBetaMeasAglD",300,0.9,1,300,0.9,1,"IsPureDMC&IsPhysTrig&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"          ,GetBetaGen,GetBetaRICH,GetBetaRICH);


	plottingfunctions.push_back(DrawBetaRes);
}


void Plotter::BookGenAcceptanceAnalysis(){
	BookerMC.BookSingleHisto("ProtonMC_GenSpectrum",1000,0,100,"IsProtonMC",GetGenMomentum);
	BookerMC.BookSingleHisto("DeutonMC_GenSpectrum",1000,0,100,"IsDeutonMC",GetGenMomentum);
}


void Plotter::BookTrackingEfficiencyAnalysis(){
	BookerMC.BookSingleScatter(         "EdepEcalvsR_P",5000,0,100,5000,0,100,"IsProtonMC&IsPositive&IsMinimumBias&IsGoodTOFStandaloneQ1&IsExtrapolInsideL8",GetRigidity,GetEdepECAL,GetEdepECAL);
	Booker.BookSingleScatter(           "EdepEcalvsR",5000,0,100,5000,0,100,"IsPositive&IsMinimumBias&IsGoodTOFStandaloneQ1&IsExtrapolInsideL8",GetRigidity,GetEdepECAL,GetEdepECAL);
	BookerMC.BookSingleScatter(           "ProxyvsR_MC",5000,0,100,5000,0,100,"IsProtonMC&IsPositive&IsMinimumBias&IsGoodTOFStandaloneQ1&IsExtrapolInsideL8",GetRigidity,GetMomentumProxy,GetMomentumProxy);
	Booker.BookSingleScatter(           "ProxyvsR",5000,0,100,5000,0,100,"IsPositive&IsMinimumBias&IsGoodTOFStandaloneQ1&IsExtrapolInsideL8",GetRigidity,GetMomentumProxy,GetMomentumProxy);

}


void Plotter::DoAllAnalyses(){

	for(int i =0; i<plottingfunctions.size();i++) plottingfunctions[i](finalHistos,finalResults);
}




#endif 
