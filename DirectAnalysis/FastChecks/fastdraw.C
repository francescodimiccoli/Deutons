#include "fastcut.h"

{
	pp=TProof::Open("");

	gROOT->ProcessLine(".L ./InputFileReader.h");

	TH1::SetDefaultSumw2();

	TChain * chainMC = InputFileReader("/data1/home/fdimicco/Deutons/DirectAnalysis/InputFileLists/1392643918-1395014410/FileListMC2.txt","Compact");
	TChain * chainDT = InputFileReader("/data1/home/fdimicco/Deutons/DirectAnalysis/InputFileLists/1392643918-1395014410/FileListDT2.txt","Compact");

	chainMC->SetProof();
	chainDT->SetProof();

	TH1F * Distr_DT = new TH1F("Distr_DT","Distr_DT",100,0,5);
	TH1F * Distr_MC = new TH1F("Distr_MC","Distr_MC",100,0,5);

	string cutMC =  (IsDownGoing + "&&" + IsPhysTrig + "&&" + IsCleanL1Hit + "&&" + IsGoodTOFStandaloneQ1 + "&&" + IsGoodTrackPattern + "&&" + IsNotL1HitMultiplX + "&&" +HasL1).c_str(); 
	string cutDT =  (IsDownGoing + "&&" + IsPhysTrig + "&&" + IsCleanL1Hit + "&&" + IsGoodTOFStandaloneQ1 + "&&" + IsGoodTrackPattern + "&&" + IsNotL1HitMultiplX + "&&" +HasL1).c_str(); 
	
	chainMC->Draw("(sa_exthit_dl1[0]^2+sa_exthit_dl1[1]^2)^0.5>>Distr_MC",cutMC.c_str());
	chainDT->Draw("(sa_exthit_dl1[0]^2+sa_exthit_dl1[1]^2)^0.5>>Distr_DT",cutDT.c_str());

	
	Distr_DT->Scale(1/Distr_DT->Integral());
	Distr_MC->Scale(1/Distr_MC->Integral());

	Distr_MC->SetLineColor(2);
	Distr_DT->Draw();
	Distr_MC->Draw("same");
	gPad->SetLogy();
}
