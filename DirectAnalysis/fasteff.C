{
	pp=TProof::Open("");

	gROOT->ProcessLine(".L ~/Deutons/DirectAnalysis/InputFileReader.h");

	TH1::SetDefaultSumw2();

	TChain * chainMC = InputFileReader("./FileListMC.txt","Compact");
	TChain * chainDT = InputFileReader("./FileListDT.txt","Compact");

	chainMC->SetProof();
	chainDT->SetProof();

	//cuts

	string IsDownGoing = "sa_tof_beta>0&&trk_rig_choutko>0";
	string IsPhysTrig = "((trigpatt & 0x2) != 0) && ((sublvl1&0x3E) !=0)";
	string IsUnbias = "((sublvl1&0x3E) !=0)";
	string IsL1HitNearExtrapol = "(sa_exthit_dl1[0]^2+sa_exthit_dl1[1]^2)^0.5";
	string IsCleanL1Hit = "abs(sa_exthit_ql1-1)<0.4";
	string IsGoodTOFStandaloneQ1 = "sa_tof_qup > 0.8 && sa_tof_qup < 1.3 && sa_tof_qdw > 0.8 && sa_tof_qdw < 1.3 && sa_trd_q > 0.4 &&  sa_trd_q < 1.7";  
	string IsGoodTrackPattern = "((trk_patty&0x2)!=0)&&((trk_patty&0xc)!=0)&&((trk_patty&0x30)!=0)&&((trk_patty&0xc0)!=0)";  
	string HasL1 = "(trk_pattxy&0x1)!=0";

	//L1pickup
	string beforecutMC =  (IsDownGoing + "&&" + IsPhysTrig + "&&" + IsL1HitNearExtrapol + "&&" + IsCleanL1Hit + "&&" + IsGoodTOFStandaloneQ1 + "&&" + IsGoodTrackPattern).c_str(); 
	string aftercutMC  =  (IsDownGoing + "&&" + IsPhysTrig + "&&" + IsL1HitNearExtrapol + "&&" + IsCleanL1Hit + "&&" + IsGoodTOFStandaloneQ1 + "&&" + IsGoodTrackPattern + "&&" + HasL1).c_str(); 
	string beforecutDT =  (IsDownGoing + "&&" + IsPhysTrig + "&&" + IsL1HitNearExtrapol + "&&" + IsCleanL1Hit + "&&" + IsGoodTOFStandaloneQ1 + "&&" + IsGoodTrackPattern).c_str(); 
	string aftercutDT  =  (IsDownGoing + "&&" + IsPhysTrig + "&&" + IsL1HitNearExtrapol + "&&" + IsCleanL1Hit + "&&" + IsGoodTOFStandaloneQ1 + "&&" + IsGoodTrackPattern + "&&" + HasL1).c_str(); 
	
	
	TH1F * hbeforeMC = new TH1F("hbeforeMC","hbeforeMC",30,-1,2);
	TH1F * hafterMC   = new TH1F("hafterMC","hafterMC",30,-1,2) ; 
	TH1F * hbeforeDT = new TH1F("hbeforeDT","hbeforeDT",30,-1,2);
	TH1F * hafterDT   = new TH1F("hafterDT","hafterDT",30,-1,2) ; 

	chainMC->Draw("log10(trk_rig_choutko)>>hbeforeMC",beforecutMC.c_str());
	chainMC->Draw("log10(trk_rig_choutko)>>hafterMC",aftercutMC.c_str());
	chainDT->Draw("log10(trk_rig_choutko)>>hbeforeDT",beforecutDT.c_str());
	chainDT->Draw("log10(trk_rig_choutko)>>hafterDT",aftercutDT.c_str());

	TH1F * MCEff     = (TH1F*) hafterMC->Clone("RatioMC");

	MCEff->Divide(hbeforeMC);	

	MCEff->SetLineColor(2);	

	TH1F * DTEff     = (TH1F*)hafterDT->Clone("RatioDT");

	DTEff->Divide(hbeforeDT);	

	DTEff->SetLineColor(1);	
	MCEff->Draw();	
	DTEff->Draw("same");	


}
