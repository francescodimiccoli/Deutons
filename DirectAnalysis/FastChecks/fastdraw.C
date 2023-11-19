#include "fastcut_v7.h"

{
	pp=TProof::Open("");

	gROOT->ProcessLine(".L ./InputFileReader.h");

	TH1::SetDefaultSumw2();

	TChain * chainMC = InputFileReader("../InputFileLists/1475712000-1485043200/FileListMC2.txt_He","Compact");
	TChain * chainDT = InputFileReader("../InputFileLists/1475712000-1485043200/FileListDT0.txt","Compact");

	chainMC->SetProof();
	chainDT->SetProof();

	TH1F * Distr_DT = new TH1F("Distr_DT","Distr_DT",100,0,25);
	TH1F * Distr_MC = new TH1F("Distr_MC","Distr_MC",100,0,25);

	string cutMC =  (IsDownGoing + "&&" + IsPhysTrig + "&&" + IsCleanL1Hit + "&&"  + IsGoodTrackPattern + "&&"  +HasL1+"&&"+HasL2+"&&"+IsKalman+"&&"+trackChicut+"&&"+qInnerCutHe+"&&"+L1LooseCutHe+"&&"+IsSingleTrack).c_str(); 
	string cutDT =  (IsDownGoing + "&&" + IsPhysTrig + "&&" + IsCleanL1Hit + "&&"  + IsGoodTrackPattern + "&&"  +HasL1+"&&"+HasL2+"&&"+IsKalman+"&&"+trackChicut+"&&"+qInnerCutHe+"&&"+L1LooseCutHe+"&&"+IsSingleTrack).c_str(); 
	
	chainMC->Draw("tof_chisqtn>>Distr_MC",cutMC.c_str());
	chainDT->Draw("tof_chisqtn>>Distr_DT",cutDT.c_str());

	
	Distr_DT->Scale(1/Distr_DT->Integral());
	Distr_MC->Scale(1/Distr_MC->Integral());

	Distr_MC->SetLineColor(2);
	Distr_DT->Draw();
	Distr_MC->Draw("same");
	gPad->SetLogy();
}
