#include "fastcut_v7.h"

{
	pp=TProof::Open("");

	gROOT->ProcessLine(".L ./InputFileReader.h");

	TH1::SetDefaultSumw2();

	TChain * chainMC = InputFileReader("../InputFileLists/1475712000-1485043200/FileListMC2.txt_He","Compact");
	TChain * chainDT = InputFileReader("FileListDT.txt","Compact");

	chainMC->SetProof();
	chainDT->SetProof();

	TH1F * Distr_DT = new TH1F("Distr_DT","Distr_DT",60,0,6);
	TH1F * Distr_DT2 = new TH1F("Distr_DT2","Distr_DT2",60,0,6);
	TH1F * Distr_MC = new TH1F("Distr_MC","Distr_MC",300,0,6);

	TH1F * BDTGood= new TH1F("BDTgood","BDTgood",100,-1,1);
	TH1F * BDTBad= new TH1F("BDTbad","BDTbad",100,-1,1);



	string cutDT =  (IsDownGoing + "&&" + IsPhysTrig + "&&" + IsCleanL1Hit + "&&"  + IsGoodTrackPattern + "&&"  +HasL1+"&&"+HasL2+"&&"+IsKalman+"&&"+trackChicut+"&&"+qInnerCut+"&&"+L1LooseCutHe+"&&"+IsSingleTrack+"&&"+qUToFCut+"&&"+qLToFCut+"&&"+IsAgl+"&&"+IsAglWindow).c_str(); 
	string cutDT2 =  (IsDownGoing + "&&" + IsPhysTrig + "&&" + IsCleanL1Hit + "&&"  + IsGoodTrackPattern + "&&"  +HasL1+"&&"+HasL2+"&&"+IsKalman+"&&"+trackChicut+"&&"+qInnerCut+"&&"+L1LooseCutHe+"&&"+IsSingleTrack+"&&"+qUToFCut+"&&"+qLToFCut+"&&"+IsAgl+"&&"+IsAglWindow+"&&"+IsBDTcut).c_str(); 

	string cutDT3 =  (IsDownGoing + "&&" + IsPhysTrig + "&&" + IsCleanL1Hit + "&&"  + IsGoodTrackPattern + "&&"  +HasL1+"&&"+HasL2+"&&"+IsKalman+"&&"+trackChicut+"&&"+qInnerCut+"&&"+L1LooseCutHe+"&&"+IsSingleTrack+"&&"+qUToFCut+"&&"+qLToFCut+"&&"+IsAgl+"&&"+IsAglWindow+"&&"+BDTgood).c_str(); 
	string cutDT4 =  (IsDownGoing + "&&" + IsPhysTrig + "&&" + IsCleanL1Hit + "&&"  + IsGoodTrackPattern + "&&"  +HasL1+"&&"+HasL2+"&&"+IsKalman+"&&"+trackChicut+"&&"+qInnerCut+"&&"+L1LooseCutHe+"&&"+IsSingleTrack+"&&"+qUToFCut+"&&"+qLToFCut+"&&"+IsAgl+"&&"+IsAglWindow+"&&"+BDTbad).c_str(); 
	
	//chainDT->Draw("(trk_rig[3]/rich_beta)*pow(1-pow(rich_beta,2),0.5)>>Distr_DT",cutDT.c_str());
	//chainDT->Draw("(trk_rig[3]/rich_beta)*pow(1-pow(rich_beta,2),0.5)>>Distr_DT2",cutDT2.c_str());
	chainDT->Draw("rich_bdt>>BDTgood",cutDT3.c_str());
	chainDT->Draw("rich_bdt>>BDTbad",cutDT4.c_str());
	
//	Distr_DT->Scale(1/Distr_DT->Integral());
	TCanvas *c1=new TCanvas("Mass");
	Distr_DT->Draw();
	Distr_DT2->Draw("same");
	gPad->SetLogy();
	TCanvas * c2 = new TCanvas("BDT");
	BDTGood->Scale(1/BDTGood->Integral());
	BDTBad->Scale(1/BDTBad->Integral());
	BDTGood->Draw();
	BDTBad->Draw("same");

}
