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
	void FillAllAnalyses(TChain * chain_RTI, TChain * chainDT,TChain *chainMC,TChain * chainDT_Cpct,TChain *chainMC1_Cpct,TChain *chainMC2_Cpct,TChain *chainMC3_Cpct,Variables * vars);
	void BookMassAnalysis();
	void BookCleaningCutsAnalysis();
	void BookRichBDTAnalysis();
	void BookBetaResMatrixAnalysis();
	void BookGenAcceptanceAnalysis();
	void BookTrackingEfficiencyAnalysis();
	void BookRigvsBetaAnalysis();
	void BookAcceptanceMatrixAnalysis();
	void BookCutVariablesAnalysis();	
	void BookMassResoAnalysis();
	void BookSmearingAnalysis();
	void BookHeliumRICHVariablesAnalysis();
	
	void DoAllAnalyses();
};

void Plotter::FillAllAnalyses(TChain * chain_RTI, TChain * chainDT,TChain *chainMC,TChain * chainDT_Cpct,TChain *chainMC1_Cpct,TChain *chainMC2_Cpct,TChain *chainMC3_Cpct,Variables * vars){
	Booker.FillEverything(DBarReader(chainDT, false,chain_RTI,chainDT_Cpct),vars);
	Booker.SaveEverything(finalHistos);
	
	BookerMC.FillEverything(DBarReader(chainMC, true,chain_RTI,chainMC1_Cpct),DBarReader(chainMC, true,chain_RTI,chainMC2_Cpct),DBarReader(chainMC, true,chain_RTI,chainMC3_Cpct),vars);
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

}

void Plotter::BookRichBDTAnalysis(){

	BookerMC.BookSingleHisto(	"BDTDiscrNaF_MCP",100,-0.5,1.,"IsPhysTrig&IsBaseline&IsLooseCharge1&IsCleaning&IsFromNaF&NafBetaSafetyCut&IsProtonMC",GetRICHBDT);
        BookerMC.BookSingleHisto(	"BDTDiscrNaF_MCD",100,-0.5,1.,"IsPhysTrig&IsBaseline&IsLooseCharge1&IsCleaning&IsFromNaF&NafBetaSafetyCut&IsDeutonMC",GetRICHBDT);
	BookerMC.BookSingleHisto(	"BDTDiscrAgl_MCP",100,-0.5,1.,"IsPhysTrig&IsBaseline&IsLooseCharge1&IsCleaning&IsFromNaF&AglBetaSafetyCut&IsProtonMC",GetRICHBDT);
        BookerMC.BookSingleHisto(	"BDTDiscrAgl_MCD",100,-0.5,1.,"IsPhysTrig&IsBaseline&IsLooseCharge1&IsCleaning&IsFromNaF&AglBetaSafetyCut&IsDeutonMC",GetRICHBDT);

	Booker.BookSingleHisto(		"BDTDiscr",100,-0.5,1.,"IsPhysTrig&IsBaseline&IsLooseCharge1&IsCleaning&IsFromNaF&NafBetaSafetyCut",GetRICHBDT);
	Booker.BookSingleHisto(		"BDTDiscr",100,-0.5,1.,"IsPhysTrig&IsBaseline&IsLooseCharge1&IsCleaning&IsFromAgl&AglBetaSafetyCut",GetRICHBDT);
	Booker.BookSingleHisto(		"BDTDiscr_He",100,-0.5,1.,"IsPhysTrig&IsBaselineHe&IsLooseCharge2&IsCleaningHe&IsFromNaF&NafBetaSafetyCut",GetRICHBDT);
	Booker.BookSingleHisto(		"BDTDiscr_He",100,-0.5,1.,"IsPhysTrig&IsBaselineHe&IsLooseCharge2&IsCleaningHe&IsFromAgl&AglBetaSafetyCut",GetRICHBDT);

	Booker.BookSingleScatter(		"RICHBDTvsBetaNaF_Z1",100,0.95,1.05,100,-0.5,0.6,"IsPhysTrig&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&IsHighEn",GetInverseBetaRICH,GetRICHBDT,GetRICHBDT);
	Booker.BookSingleScatter(		"RICHBDTvsBetaAgl_Z1",100,0.95,1.05,100,-0.5,0.6,"IsPhysTrig&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&IsHighEn",GetInverseBetaRICH,GetRICHBDT,GetRICHBDT);
	Booker.BookSingleScatter(		"RICHBDTvsBetaNaF_Z2",100,0.95,1.05,100,-0.5,0.6,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsHighEn",GetInverseBetaRICH,GetRICHBDT,GetRICHBDT);
	Booker.BookSingleScatter(		"RICHBDTvsBetaAgl_Z2",100,0.95,1.05,100,-0.5,0.6,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsHighEn",GetInverseBetaRICH,GetRICHBDT,GetRICHBDT);
	
	Booker.BookSingleScatter(	"RICHBDTvsMass_Z1",100,0,4.5,100,-0.5,0.6,"IsPhysTrig&IsBaseline&IsLooseCharge1&IsCleaning&IsFromNaF&NafBetaSafetyCut",GetRecMassRICH,GetRICHBDT,GetRICHBDT);
        Booker.BookSingleScatter(	"RICHBDTvsMass_Z1",100,0,4.5,100,-0.5,0.6,"IsPhysTrig&IsBaseline&IsLooseCharge1&IsCleaning&IsFromAgl&AglBetaSafetyCut",GetRecMassRICH,GetRICHBDT,GetRICHBDT);
        Booker.BookSingleScatter(	"RICHBDTvsMass_Z2",100,0,4.5,100,-0.5,0.6,"IsPhysTrig&IsBaselineHe&IsLooseCharge2&IsCleaningHe&IsFromNaF&NafBetaSafetyCut",GetRecMassRICH,GetRICHBDT,GetRICHBDT);
        Booker.BookSingleScatter(	"RICHBDTvsMass_Z2",100,0,4.5,100,-0.5,0.6,"IsPhysTrig&IsBaselineHe&IsLooseCharge2&IsCleaningHe&IsFromAgl&AglBetaSafetyCut",GetRecMassRICH,GetRICHBDT,GetRICHBDT);


	plottingfunctions.push_back(DrawBDT);

}

void Plotter::BookBetaResMatrixAnalysis(){
	
	

	BookerMC.BookSingleScatter("BetagenvsBetaMeasTOFP",300,0.3,1,300,0.3,1,"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime",GetBetaGen ,GetBetaTOF ,GetBetaTOF);
        BookerMC.BookSingleScatter("BetagenvsBetaMeasNaFP",300,0.7,1,300,0.7,1,"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF" ,GetBetaGen ,GetBetaRICH,GetBetaRICH);
        BookerMC.BookSingleScatter("BetagenvsBetaMeasAglP",300,0.9,1,300,0.9,1,"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl" ,GetBetaGen ,GetBetaRICH,GetBetaRICH);
        BookerMC.BookSingleScatter("BetagenvsBetaMeasTOFD",300,0.3,1,300,0.3,1,"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime",GetBetaGen ,GetBetaTOF ,GetBetaTOF);
        BookerMC.BookSingleScatter("BetagenvsBetaMeasNaFD",300,0.7,1,300,0.7,1,"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF" ,GetBetaGen ,GetBetaRICH,GetBetaRICH);
        BookerMC.BookSingleScatter("BetagenvsBetaMeasAglD",300,0.9,1,300,0.9,1,"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl" ,GetBetaGen ,GetBetaRICH,GetBetaRICH);

	BookerMC.BookSingleScatter("BetagenvsBetaMeasTOFHe",300,0.3,1,300,0.3,1,"IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsGoodTime",GetBetaGen ,GetBetaTOF ,GetBetaTOF);
        BookerMC.BookSingleScatter("BetagenvsBetaMeasNaFHe",300,0.7,1,300,0.7,1,"IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF" ,GetBetaGen ,GetBetaRICH,GetBetaRICH);
        BookerMC.BookSingleScatter("BetagenvsBetaMeasAglHe",300,0.9,1,300,0.9,1,"IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl" ,GetBetaGen ,GetBetaRICH,GetBetaRICH);

	BookerMC.BookSingleScatter("BetaSlowvsBetaMeasTOFP",300,0.3,1,300,-0.5,0.5,"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime",GetBetaTOF,GetBetaSlowTOF,GetBetaSlowTOF);
        BookerMC.BookSingleScatter("BetaSlowvsBetaMeasNaFP",300,0.7,1,300,-0.05,0.05,"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF",GetBetaRICH,GetBetaSlowRICH,GetBetaSlowRICH);
        BookerMC.BookSingleScatter("BetaSlowvsBetaMeasAglP",300,0.9,1,300,-0.05,0.05,"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl",GetBetaRICH,GetBetaSlowRICH,GetBetaSlowRICH);
        BookerMC.BookSingleScatter("BetaSlowvsBetaMeasTOFD",300,0.3,1,300,-0.5,0.5,"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime",GetBetaTOF,GetBetaSlowTOF,GetBetaSlowTOF);
        BookerMC.BookSingleScatter("BetaSlowvsBetaMeasNaFD",300,0.7,1,300,-0.05,0.05,"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF",GetBetaRICH,GetBetaSlowRICH,GetBetaSlowRICH);
        BookerMC.BookSingleScatter("BetaSlowvsBetaMeasAglD",300,0.9,1,300,-0.05,0.05,"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl",GetBetaRICH,GetBetaSlowRICH,GetBetaSlowRICH);
	BookerMC.BookSingleScatter("RSlowvsRMeas",1000,0.3,100,300,-3,3,  "IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime"  ,GetRigidity,GetRigSlow,GetRigSlow);




	plottingfunctions.push_back(DrawBetaRes);
}

void Plotter::BookSmearingAnalysis(){
	Booker.BookSingleHisto(		"BetaTOF_Z1",100,0.5 ,1.5 ,"IsPhysTrig&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF&IsHighEn",GetInverseBetaTOF);
	Booker.BookSingleHisto(		"BetaNaF_Z1",100,0.95,1.05,"IsPhysTrig&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut&IsHighEn",GetInverseBetaRICH);
	Booker.BookSingleHisto(		"BetaAgl_Z1",100,0.95,1.05,"IsPhysTrig&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut&IsHighEn",GetInverseBetaRICH);
	
	Booker.BookSingleHisto(		"BetaTOF_Z2",100,0.5 ,1.5 ,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsGoodTime&QualityTOF&IsHighEn",GetInverseBetaTOF);
	Booker.BookSingleHisto(		"BetaNaF_Z2",100,0.95,1.05,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsHighEn",GetInverseBetaRICH);
	Booker.BookSingleHisto(		"BetaAgl_Z2",100,0.95,1.05,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsHighEn",GetInverseBetaRICH);
	
       	BookerMC.BookSingleHisto(	"BetaTOF_P_MC",100,0.5 ,1.5 ,"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF&IsHighEn",GetInverseBetaTOF);
	BookerMC.BookSingleHisto(	"BetaNaF_P_MC",100,0.95,1.05,"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut&IsHighEn",GetInverseBetaRICH);
	BookerMC.BookSingleHisto(	"BetaAgl_P_MC",100,0.95,1.05,"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut&IsHighEn",GetInverseBetaRICH);

       	BookerMC.BookSingleHisto(	"BetaTOF_D_MC",100,0.5 ,1.5 ,"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF&IsHighEn",GetInverseBetaTOF);
	BookerMC.BookSingleHisto(	"BetaNaF_D_MC",100,0.95,1.05,"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut&IsHighEn",GetInverseBetaRICH);
	BookerMC.BookSingleHisto(	"BetaAgl_D_MC",100,0.95,1.05,"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut&IsHighEn",GetInverseBetaRICH);

       	BookerMC.BookSingleHisto(	"BetaTOF_He_MC",100,0.5 ,1.5 ,"IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsGoodTime&QualityTOF&IsHighEn",GetInverseBetaTOF);
	BookerMC.BookSingleHisto(	"BetaNaF_He_MC",100,0.95,1.05,"IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsHighEn",GetInverseBetaRICH);
	BookerMC.BookSingleHisto(	"BetaAgl_He_MC",100,0.95,1.05,"IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsHighEn",GetInverseBetaRICH);


	BookerMC.BookSingleHisto(	"BetaTOF_HE_MC_Smeared",100,0.5,1.5,"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&QualityTOF&IsHighEn",GetInverseBetaTOF_Smear);
	BookerMC.BookSingleHisto(	"BetaAgl_HE_MC_Smeared",100,0.95,1.05,"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut&IsHighEn",GetInverseBetaRICH_Smear);

	plottingfunctions.push_back(DrawBetaSmear);

}	



void Plotter::BookHeliumRICHVariablesAnalysis(){

	Booker.BookSingleScatter("RICHTOThitsvsBetaAgl_Z2",100,0.95,1.05,100,0,25,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsHighEn",GetInverseBetaRICH,GetNRichTOThits,GetNRichTOThits);
	Booker.BookSingleScatter("RICHUSEDhitsvsBetaAgl_Z2",100,0.95,1.05,100,0,25,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsHighEn",GetInverseBetaRICH,GetNRichUSEDhits,GetNRichUSEDhits);
	Booker.BookSingleScatter("RICHPMTsvsBetaAgl_Z2",100,0.95,1.05,100,0,25,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsHighEn",GetInverseBetaRICH,GetNRichPMTs,GetNRichPMTs);
	Booker.BookSingleScatter("RICHProbvsBetaAgl_Z2",100,0.95,1.05,100,0,1,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsHighEn",GetInverseBetaRICH,GetRICHProb,GetRICHProb);
	Booker.BookSingleScatter("RICHCollovsBetaAgl_Z2",100,0.95,1.05,100,0,1,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsHighEn",GetInverseBetaRICH,GetRICHCollovertotal,GetRICHCollovertotal);
	Booker.BookSingleScatter("RICHTOFConsvsBetaAgl_Z2",100,0.95,1.05,100,-1,1,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsHighEn",GetInverseBetaRICH,GetRICHTOFBetaConsistency,GetRICHTOFBetaConsistency);

	Booker.BookSingleScatter("RICHTOThitsvsBetaNaF_Z2",100,0.95,1.05,100,0,25,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsHighEn",GetInverseBetaRICH,GetNRichTOThits,GetNRichTOThits);
	Booker.BookSingleScatter("RICHUSEDhitsvsBetaNaF_Z2",100,0.95,1.05,100,0,25,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsHighEn",GetInverseBetaRICH,GetNRichUSEDhits,GetNRichUSEDhits);
	Booker.BookSingleScatter("RICHPMTsvsBetaNaF_Z2",100,0.95,1.05,100,0,25,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsHighEn",GetInverseBetaRICH,GetNRichPMTs,GetNRichPMTs);
	Booker.BookSingleScatter("RICHProbvsBetaNaF_Z2",100,0.95,1.05,100,0,1,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsHighEn",GetInverseBetaRICH,GetRICHProb,GetRICHProb);
	Booker.BookSingleScatter("RICHCollovsBetaNaF_Z2",100,0.95,1.05,100,0,1,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsHighEn",GetInverseBetaRICH,GetRICHCollovertotal,GetRICHCollovertotal);
	Booker.BookSingleScatter("RICHTOFConsvsBetaNaF_Z2",100,0.95,1.05,100,-1,1,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsHighEn",GetInverseBetaRICH,GetRICHTOFBetaConsistency,GetRICHTOFBetaConsistency);

	Booker.BookSingleScatter("RICHTOThitsvsMassNaF_Z2",100,0.5,7.05,100,0,25,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsBinNaF",GetRecMassRICH,GetNRichTOThits,GetNRichTOThits);
	Booker.BookSingleScatter("RICHUSEDhitsvsMassNaF_Z2",100,0.5,7.05,100,0,25,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsBinNaF",GetRecMassRICH,GetNRichUSEDhits,GetNRichUSEDhits);
	Booker.BookSingleScatter("RICHPMTsvsMassNaF_Z2",100,0.5,7.05,100,0,25,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsBinNaF",GetRecMassRICH,GetNRichPMTs,GetNRichPMTs);
	Booker.BookSingleScatter("RICHProbvsMassNaF_Z2",100,0.5,7.05,100,0,1,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsBinNaF",GetRecMassRICH,GetRICHProb,GetRICHProb);
	Booker.BookSingleScatter("RICHCollovsMassNaF_Z2",100,0.5,7.05,100,0,1,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsBinNaF",GetRecMassRICH,GetRICHCollovertotal,GetRICHCollovertotal);
	Booker.BookSingleScatter("RICHTOFConsvsMassNaF_Z2",100,0.5,7.05,100,-1,1,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsBinNaF",GetRecMassRICH,GetRICHTOFBetaConsistency,GetRICHTOFBetaConsistency);
	Booker.BookSingleScatter("RICHTotPhelNaF_Z2",100,0.5,7.05,100,0,20,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&IsBinNaF",GetRecMassRICH,GetRICHTotPhel,GetRICHTotPhel);

	Booker.BookSingleScatter("RICHTOThitsvsMassAgl_Z2",100,0.5,7.05,100,0,25,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsBinAgl",GetRecMassRICH,GetNRichTOThits,GetNRichTOThits);
	Booker.BookSingleScatter("RICHUSEDhitsvsMassAgl_Z2",100,0.5,7.05,100,0,25,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsBinAgl",GetRecMassRICH,GetNRichUSEDhits,GetNRichUSEDhits);
	Booker.BookSingleScatter("RICHPMTsvsMassAgl_Z2",100,0.5,7.05,100,0,25,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsBinAgl",GetRecMassRICH,GetNRichPMTs,GetNRichPMTs);
	Booker.BookSingleScatter("RICHProbvsMassAgl_Z2",100,0.5,7.05,100,0,1,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsBinAgl",GetRecMassRICH,GetRICHProb,GetRICHProb);
	Booker.BookSingleScatter("RICHCollovsMassAgl_Z2",100,0.5,7.05,100,0,1,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsBinAgl",GetRecMassRICH,GetRICHCollovertotal,GetRICHCollovertotal);
	Booker.BookSingleScatter("RICHTOFConsvsMassAgl_Z2",100,0.5,7.05,100,-1,1,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsBinAgl",GetRecMassRICH,GetRICHTOFBetaConsistency,GetRICHTOFBetaConsistency);
	Booker.BookSingleScatter("RICHTotPhelAgl_Z2",100,0.5,7.05,100,0,20,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&IsBinAgl",GetRecMassRICH,GetRICHTotPhel,GetRICHTotPhel);


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

void Plotter::BookRigvsBetaAnalysis(){
	BookerMC.BookSingleScatter("RigvsBeta_TOFP",600,0,50,100,0.4,1.05,"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF",GetRigidity,GetBetaTOF);
	BookerMC.BookSingleScatter("RigvsBeta_TOFD",600,0,50,100,0.4,1.05,"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF",GetRigidity,GetBetaTOF);
	BookerMC.BookSingleScatter("RigvsBeta_TOFHe",600,0,50,100,0.4,1.05,"IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsGoodTime",GetRigidity,GetBetaTOF);

	BookerMC.BookSingleScatter("RigvsBeta_NaFP",600,1,50,100,0.7,1.05,"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",GetRigidity,GetBetaRICH);
	BookerMC.BookSingleScatter("RigvsBeta_NaFD",600,1,50,100,0.7,1.05,"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",GetRigidity,GetBetaRICH);
	BookerMC.BookSingleScatter("RigvsBeta_NaFHe",600,1,50,100,0.7,1.05,"IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF",GetRigidity,GetBetaRICH);

	BookerMC.BookSingleScatter("RigvsBeta_AglP",600,0,50,100,0.9,1.05,"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",GetRigidity,GetBetaRICH);
	BookerMC.BookSingleScatter("RigvsBeta_AglD",600,0,50,100,0.9,1.05,"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",GetRigidity,GetBetaRICH);
	BookerMC.BookSingleScatter("RigvsBeta_AglHe",600,0,50,100,0.9,1.05,"IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl",GetRigidity,GetBetaRICH);

	
	Booker.BookSingleScatter("RigvsBeta_TOF_data",600,0,50,100,0.4,1.05,"IsPhysTrig&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF",GetRigidity,GetBetaTOF);
	Booker.BookSingleScatter("RigvsBeta_TOF_dataprim",600,50,6,100,0.4,1.05,"IsPhysTrig&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF&IsPrimary",GetRigidity,GetBetaTOF);
	Booker.BookSingleScatter("RigvsBeta_TOF_Z2",600,0,50,100,0.4,1.05,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsGoodTime",GetRigidity,GetBetaTOF);
	
	Booker.BookSingleScatter("RigvsBeta_NaF_data",600,1,50,100,0.7,1.05,"IsPhysTrig&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",GetRigidity,GetBetaRICH);
	Booker.BookSingleScatter("RigvsBeta_NaFD_dataprim",600,1,50,100,0.7,1.05,"IsPhysTrig&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut&IsPrimary",GetRigidity,GetBetaRICH);
	Booker.BookSingleScatter("RigvsBeta_NaF_Z2",600,1,50,100,0.7,1.05,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF",GetRigidity,GetBetaRICH);
	
	Booker.BookSingleScatter("RigvsBeta_Agl_data",600,0,50,100,0.9,1.05,"IsPhysTrig&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",GetRigidity,GetBetaRICH);
	Booker.BookSingleScatter("RigvsBeta_AglD_dataprim",600,0,50,100,0.9,1.05,"IsPhysTrig&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut&IsPrimary",GetRigidity,GetBetaRICH);
	Booker.BookSingleScatter("RigvsBeta_Agl_Z2",600,0,50,100,0.9,1.05,"IsPhysTrig&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl",GetRigidity,GetBetaRICH);


	plottingfunctions.push_back(DrawBetavsRig);

}

void Plotter::BookAcceptanceMatrixAnalysis(){

	BookerMC.BookBinnedScatter("Acceptance Matrix",PRB,"IsProtonMC&IsStandardSel",GetGenMomentum,GetRigidity);

	BookerMC.BookBinnedScatter("Acceptance Matrix v0",PRB,"IsProtonMC&IsStandardSel_v0",GetGenMomentum,GetRigidity);
	void BookHeliumRICHVariablesAnalysis();	BookerMC.BookBinnedScatter("Acceptance Matrix v1",PRB,"IsProtonMC&IsStandardSel_v1",GetGenMomentum,GetRigidity);
	BookerMC.BookBinnedScatter("Acceptance Matrix v2",PRB,"IsProtonMC&IsStandardSel_v2",GetGenMomentum,GetRigidity);
	BookerMC.BookBinnedScatter("Acceptance Matrix v3",PRB,"IsProtonMC&IsStandardSel_v3",GetGenMomentum,GetRigidity);
	BookerMC.BookBinnedScatter("Acceptance Matrix v4",PRB,"IsProtonMC&IsStandardSel_v4",GetGenMomentum,GetRigidity);
	BookerMC.BookBinnedScatter("Acceptance Matrix v5",PRB,"IsProtonMC&IsStandardSel_v5",GetGenMomentum,GetRigidity);
	BookerMC.BookBinnedScatter("Acceptance Matrix v6",PRB,"IsProtonMC&IsStandardSel_v6",GetGenMomentum,GetRigidity);
	BookerMC.BookBinnedScatter("Acceptance Matrix v7",PRB,"IsProtonMC&IsStandardSel_v7",GetGenMomentum,GetRigidity);
	BookerMC.BookBinnedScatter("Acceptance Matrix v8",PRB,"IsProtonMC&IsStandardSel_v8",GetGenMomentum,GetRigidity);

	plottingfunctions.push_back(DrawAcceptanceMatrix);
}


void Plotter::BookMassResoAnalysis(){
	 BookerMC.BookSingleScatter(     "MassvsBetaTOF",100,0.3,0.95,100,0,4.5,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsOnlyFromToF",GetBetaTOF,GetRecMassTOF,GetRecMassTOF);
	 BookerMC.BookSingleScatter(     "MassvsBetaNaF",100,0.8,0.995,100,0,4.5,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",GetBetaRICH,GetRecMassRICH,GetRecMassRICH);
	 BookerMC.BookSingleScatter(     "MassvsBetaAgl",100,0.96,0.998,100,0,4.5,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",GetBetaRICH,GetRecMassRICH,GetRecMassRICH);
	
	plottingfunctions.push_back(DrawMassRes);
}


void Plotter::BookCutVariablesAnalysis(){


	//no L1
	Booker.BookSingleScatter("N Tof Clusters"	,50,0,15,100,0,100,		"IsPositive&IsPrimary&IsBaseline",GetNToFClusters,GetRigidity);
	Booker.BookSingleScatter("Q L tof"		,100,0,4,100,0,100,			"IsPositive&IsPrimary&IsBaseline",GetLtofQ,GetRigidity);
	Booker.BookSingleScatter("Q U tof"		,100,0,4,100,0,100,			"IsPositive&IsPrimary&IsBaseline",GetUtofQ,GetRigidity);
	Booker.BookSingleScatter("tof Coo chi"		,200,0,20,100,0,100,		"IsPositive&IsPrimary&IsBaseline",GetTofChisqcn,GetRigidity);
	Booker.BookSingleScatter("tof Time chi"		,200,0,20,100,0,100,		"IsPositive&IsPrimary&IsBaseline",GetTofChisqtn,GetRigidity);
	Booker.BookSingleScatter("N Tracks"		,100,0,10,100,0,100,			"IsPositive&IsPrimary&IsBaseline",GetNTracks,GetRigidity);
	Booker.BookSingleScatter("On Time"		,100,0,10,100,0,100,			"IsPositive&IsPrimary&IsBaseline",GetTofOnTime,GetRigidity);

	BookerMC.BookSingleScatter("N Tof Clusters MC",50,0,15,100,0,100,		"IsPositive&IsProtonMC&IsBaseline",GetNToFClusters,GetRigidity);
	BookerMC.BookSingleScatter("Q L tof MC",100,0,4,100,0,100,			"IsPositive&IsProtonMC&IsBaseline",GetLtofQ,GetRigidity);
	BookerMC.BookSingleScatter("Q U tof MC",100,0,4,100,0,100,			"IsPositive&IsProtonMC&IsBaseline",GetUtofQ,GetRigidity);
	BookerMC.BookSingleScatter("tof Coo chi MC",200,0,20,100,0,100,		"IsPositive&IsProtonMC&IsBaseline",GetTofChisqcn,GetRigidity);
	BookerMC.BookSingleScatter("tof Time chi MC",200,0,20,100,0,100,		"IsPositive&IsProtonMC&IsBaseline",GetTofChisqtn,GetRigidity);
	BookerMC.BookSingleScatter("N Tracks MC",100,0,10,100,0,100,		"IsPositive&IsProtonMC&IsBaseline",GetNTracks,GetRigidity);
	BookerMC.BookSingleScatter("On Time MC",100,0,10,100,0,100,		"IsPositive&IsProtonMC&IsBaseline",GetTofOnTime,GetRigidity);

	Booker.BookSingleScatter("TofClustersVsQLtof",50,0,15,100,0,4,"IsPositive&IsPrimary&IsBaseline",GetNToFClusters,GetLtofQ);
	Booker.BookSingleScatter("TofClustersVsQUtof",50,0,15,100,0,4,"IsPositive&IsPrimary&IsBaseline",GetNToFClusters,GetUtofQ);
	Booker.BookSingleScatter("TofClustersVsNTrtracks",50,0,15,100,0,10,"IsPositive&IsPrimary&IsBaseline",GetNToFClusters,GetNTracks);
	Booker.BookSingleScatter("TofClustersVsOnTime",50,0,15,100,0,10,"IsPositive&IsPrimary&IsBaseline",GetNToFClusters,GetTofOnTime);
	Booker.BookSingleScatter("OnTimeVsQLtof",50,0,15,100,0,4,"IsPositive&IsPrimary&IsBaseline",GetTofOnTime,GetLtofQ);
	Booker.BookSingleScatter("OnTimeVsQUtof",50,0,15,100,0,4,"IsPositive&IsPrimary&IsBaseline",GetTofOnTime,GetUtofQ);
	Booker.BookSingleScatter("OnTimeVsNTrTracks",50,0,15,100,0,10,"IsPositive&IsPrimary&IsBaseline",GetTofOnTime,GetNTracks);
	Booker.BookSingleScatter("NTrTracksvsQLtof",50,0,15,100,0,4,"IsPositive&IsPrimary&IsBaseline",GetNTracks,GetLtofQ);
	Booker.BookSingleScatter("NTrTracksvsQUtof",50,0,15,100,0,4,"IsPositive&IsPrimary&IsBaseline",GetNTracks,GetUtofQ);

	BookerMC.BookSingleScatter("TofClustersVsQLtof MC",50,0,15,100,0,4,"IsPositive&IsProtonMC&IsBaseline",GetNToFClusters,GetLtofQ);
	BookerMC.BookSingleScatter("TofClustersVsQUtof MC",50,0,15,100,0,4,"IsPositive&IsProtonMC&IsBaseline",GetNToFClusters,GetUtofQ);
	BookerMC.BookSingleScatter("TofClustersVsNTrtracks MC",50,0,15,100,0,10,"IsPositive&IsProtonMC&IsBaseline",GetNToFClusters,GetNTracks);
	BookerMC.BookSingleScatter("TofClustersVsOnTime MC",50,0,15,100,0,10,"IsPositive&IsProtonMC&IsBaseline",GetNToFClusters,GetTofOnTime);
	BookerMC.BookSingleScatter("OnTimeVsQLtof MC",50,0,15,100,0,4,"IsPositive&IsProtonMC&IsBaseline",GetTofOnTime,GetLtofQ);
	BookerMC.BookSingleScatter("OnTimeVsQUtof MC",50,0,15,100,0,4,"IsPositive&IsProtonMC&IsBaseline",GetTofOnTime,GetUtofQ);
	BookerMC.BookSingleScatter("OnTimeVsNTrTracks MC",50,0,15,100,0,10,"IsPositive&IsProtonMC&IsBaseline",GetTofOnTime,GetNTracks);
	BookerMC.BookSingleScatter("NTrTracksvsQLtof MC",50,0,15,100,0,4,"IsPositive&IsProtonMC&IsBaseline",GetNTracks,GetLtofQ);
	BookerMC.BookSingleScatter("NTrTracksvsQUtof MC",50,0,15,100,0,4,"IsPositive&IsProtonMC&IsBaseline",GetNTracks,GetUtofQ);

	//with L1
	Booker.BookSingleScatter("N Tof Clusters L1"	,50,0,15,100,0,100,		"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetNToFClusters,GetRigidity);
	Booker.BookSingleScatter("Q L tof L1"		,100,0,4,100,0,100,			"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetLtofQ,GetRigidity);
	Booker.BookSingleScatter("Q U tof L1"		,100,0,4,100,0,100,			"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetUtofQ,GetRigidity);
	Booker.BookSingleScatter("tof Coo chi L1"		,200,0,20,100,0,100,		"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetTofChisqcn,GetRigidity);
	Booker.BookSingleScatter("tof Time chi L1"		,200,0,20,100,0,100,		"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetTofChisqtn,GetRigidity);
	Booker.BookSingleScatter("N Tracks L1"		,100,0,10,100,0,100,			"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetNTracks,GetRigidity);
	Booker.BookSingleScatter("On Time L1"		,100,0,10,100,0,100,			"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetTofOnTime,GetRigidity);

	BookerMC.BookSingleScatter("N Tof Clusters MC L1",50,0,15,100,0,100,		"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetNToFClusters,GetRigidity);
	BookerMC.BookSingleScatter("Q L tof MC L1",100,0,4,100,0,100,			"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetLtofQ,GetRigidity);
	BookerMC.BookSingleScatter("Q U tof MC L1",100,0,4,100,0,100,			"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetUtofQ,GetRigidity);
	BookerMC.BookSingleScatter("tof Coo chi MC L1",200,0,20,100,0,100,		"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetTofChisqcn,GetRigidity);
	BookerMC.BookSingleScatter("tof Time chi MC L1",200,0,20,100,0,100,		"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetTofChisqtn,GetRigidity);
	BookerMC.BookSingleScatter("N Tracks MC L1",100,0,10,100,0,100,		"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetNTracks,GetRigidity);
	BookerMC.BookSingleScatter("On Time MC L1",100,0,10,100,0,100,		"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetTofOnTime,GetRigidity);

	Booker.BookSingleScatter("TofClustersVsQLtof L1",50,0,15,100,0,4,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetNToFClusters,GetLtofQ);
	Booker.BookSingleScatter("TofClustersVsQUtof L1",50,0,15,100,0,4,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetNToFClusters,GetUtofQ);
	Booker.BookSingleScatter("TofClustersVsNTrtracks L1",50,0,15,100,0,10,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetNToFClusters,GetNTracks);
	Booker.BookSingleScatter("TofClustersVsOnTime L1",50,0,15,100,0,10,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetNToFClusters,GetTofOnTime);
	Booker.BookSingleScatter("OnTimeVsQLtof L1",50,0,15,100,0,4,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetTofOnTime,GetLtofQ);
	Booker.BookSingleScatter("OnTimeVsQUtof L1",50,0,15,100,0,4,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetTofOnTime,GetUtofQ);
	Booker.BookSingleScatter("OnTimeVsNTrTracks L1",50,0,15,100,0,10,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetTofOnTime,GetNTracks);
	Booker.BookSingleScatter("NTrTracksvsQLtof L1",50,0,15,100,0,4,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetNTracks,GetLtofQ);
	Booker.BookSingleScatter("NTrTracksvsQUtof L1",50,0,15,100,0,4,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1",GetNTracks,GetUtofQ);

	BookerMC.BookSingleScatter("TofClustersVsQLtof MC L1",50,0,15,100,0,4,"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetNToFClusters,GetLtofQ);
	BookerMC.BookSingleScatter("TofClustersVsQUtof MC L1",50,0,15,100,0,4,"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetNToFClusters,GetUtofQ);
	BookerMC.BookSingleScatter("TofClustersVsNTrtracks MC L1",50,0,15,100,0,10,"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetNToFClusters,GetNTracks);
	BookerMC.BookSingleScatter("TofClustersVsOnTime MC L1",50,0,15,100,0,10,"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetNToFClusters,GetTofOnTime);
	BookerMC.BookSingleScatter("OnTimeVsQLtof MC L1",50,0,15,100,0,4,"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetTofOnTime,GetLtofQ);
	BookerMC.BookSingleScatter("OnTimeVsQUtof MC L1",50,0,15,100,0,4,"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetTofOnTime,GetUtofQ);
	BookerMC.BookSingleScatter("OnTimeVsNTrTracks MC L1",50,0,15,100,0,10,"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetTofOnTime,GetNTracks);
	BookerMC.BookSingleScatter("NTrTracksvsQLtof MC L1",50,0,15,100,0,4,"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetNTracks,GetLtofQ);
	BookerMC.BookSingleScatter("NTrTracksvsQUtof MC L1",50,0,15,100,0,4,"IsPositive&IsProtonMC&IsBaseline&L1LooseCharge1",GetNTracks,GetUtofQ);


	Booker.BookSingleScatter("Multitracks R"          ,200,0,30,200,0,30,             "IsPositive&IsPrimary&IsBaseline&Is2Tracks",GetRigidity,GetRigiditySecondTrack);
	Booker.BookSingleScatter("Multitracks R MC"       ,200,0,30,200,0,30,             "IsPositive&IsProtonMC&IsBaseline&Is2Tracks",GetRigidity,GetRigiditySecondTrack);


	plottingfunctions.push_back(DrawVariables);
}

void Plotter::DoAllAnalyses(){

	for(int i =0; i<plottingfunctions.size();i++) plottingfunctions[i](finalHistos,finalResults);
}




#endif 
