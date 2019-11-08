#include "fastcut.h"
{
	pp=TProof::Open("");

	gROOT->ProcessLine(".L ~/Deutons/DirectAnalysis/InputFileReader.h");

	TH1::SetDefaultSumw2();

	TChain * chainMC = InputFileReader("/data1/home/fdimicco/Deutons/DirectAnalysis/InputFileLists/1392643918-1395014410/FileListMC2.txt","Compact");
	TChain * chainDT = InputFileReader("/data1/home/fdimicco/Deutons/DirectAnalysis/InputFileLists/1392643918-1395014410/FileListDT2.txt","Compact");

	chainMC->SetProof();
	chainDT->SetProof();

	//cuts

	string beforecutMC =  (IsDownGoing + "&&" + IsPhysTrig + "&&" +  IsCompact_SA + "&&" + IsGoodTrackPattern + "&&" + HasL2 + "&&" + HasL1 +"&&" + qInnerCut+"&&" + qUToFCut+"&&" + qLToFCut).c_str(); 
	string aftercutMC  =  (IsDownGoing + "&&" + IsUnbias i + "&&" +  IsCompact_SA + "&&" + IsGoodTrackPattern + "&&" + HasL2 + "&&" + HasL1 +"&&" + qInnerCut+"&&" + qUToFCut+"&&" + qLToFCut).c_str(); 
	string beforecutDT =  (IsDownGoing + "&&" + IsPhysTrig + "&&" +  IsCompact_SA + "&&" + IsGoodTrackPattern + "&&" + HasL2 + "&&" + HasL1 +"&&" + qInnerCut+"&&" + qUToFCut+"&&" + qLToFCut).c_str(); 
	string aftercutDT  =  (IsDownGoing + "&&" + IsUnbias   + "&&" +  IsCompact_SA + "&&" + IsGoodTrackPattern + "&&" + HasL2 + "&&" + HasL1 +"&&" + qInnerCut+"&&" + qUToFCut+"&&" + qLToFCut).c_str(); 
	
	
	TH1F * hbeforeMC = new TH1F("hbeforeMC","hbeforeMC",30,-1,3);
	TH1F * hafterMC   = new TH1F("hafterMC","hafterMC",30,-1,3) ; 
	TH1F * hbeforeDT = new TH1F("hbeforeDT","hbeforeDT",30,-1,3);
	TH1F * hafterDT   = new TH1F("hafterDT","hafterDT",30,-1,3) ; 

	chainMC->Draw("1>>hbeforeMC",beforecutMC.c_str());
	chainMC->Draw("1>>hafterMC",aftercutMC.c_str());
	chainDT->Draw("1>>hbeforeDT",beforecutDT.c_str());
	chainDT->Draw("1>>hafterDT",aftercutDT.c_str());

	TH1F * MCEff     = (TH1F*) hafterMC->Clone("RatioMC");

	MCEff->Divide(hbeforeMC);	

	MCEff->SetLineColor(2);	

	TH1F * DTEff     = (TH1F*)hafterDT->Clone("RatioDT");

	DTEff->Divide(hbeforeDT);	

	DTEff->SetLineColor(1);	
	MCEff->Draw();	
	DTEff->Draw("same");

}
