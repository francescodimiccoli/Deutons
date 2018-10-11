{
	pp=TProof::Open("");

	gROOT->ProcessLine(".L ~/Work/Deutons/DirectAnalysis/include/InputFileReader.h");

	TH1::SetDefaultSumw2();

	TChain * chainMC = InputFileReader("./FileListMC.txt","Event");
	TChain * chainDT = InputFileReader("./FileListDT.txt","Event");


	chainMC->SetProof();
	chainDT->SetProof();

	//trigger	
	string beforecutMC = "abs(Tracker.q_inn-1)<0.3 && (Tof.z_nhit>=3)&&(Tof.z_like>0)&&((Tof.clsn[0]+Tof.clsn[2])<=4) && Tof.beta>0&&(trigpatt&0x2)!=0&&Tracker.q_lay[1][0]>0";
	string aftercutMC  = "abs(Tracker.q_inn-1)<0.3 && (Tof.z_nhit>=3)&&(Tof.z_like>0)&&((Tof.clsn[0]+Tof.clsn[2])<=4) && Tof.beta>0&&(trigpatt&0x2)!=0&&Tracker.q_lay[1][0]>0&&(Header.sublvl1&0x1E)!=0";		
	string beforecutDT  = "Header.pres_weight*(100*(abs(Tracker.q_inn-1)<0.3 && (Tof.z_nhit>=3)&&(Tof.z_like>0)&&((Tof.clsn[0]+Tof.clsn[2])<=4) && Tof.beta>0&&Tracker.q_lay[1][0]>0&&((trigpatt&0x2)!=0)&&(Header.sublvl1&0x3E)==0 ))";		
	string aftercutDT   = "Header.pres_weight*((abs(Tracker.q_inn-1)<0.3 && (Tof.z_nhit>=3)&&(Tof.z_like>0)&&((Tof.clsn[0]+Tof.clsn[2])<=4) && Tof.beta>0&&Tracker.q_lay[1][0]>0&&((trigpatt&0x2)!=0)&&(Header.sublvl1&0x3E)!=0))";
	
	//tracking
	//string beforecutMC = "(Tof.z_nhit>=3)&&(Tof.z_like>0)&&((Tof.clsn[0]+Tof.clsn[2])<=4) && Tof.beta>0&&(trigpatt&0x2)!=0&&(Header.sublvl1&0x3E)!=0 && SA.beta>0 && SA.beta_patt==0";
	//string aftercutMC  = "(Tof.z_nhit>=3)&&(Tof.z_like>0)&&((Tof.clsn[0]+Tof.clsn[2])<=4) && Tof.beta>0&&(trigpatt&0x2)!=0&&(Header.sublvl1&0x3E)!=0 && SA.beta>0 && SA.beta_patt==0 && Tracker.rig[1]!=0&&Tracker.q_lay_status[1][1]==0";		

	//minTOF
//	string beforecutMC = "abs(Tracker.q_inn-1)<0.3 && (trigpatt&0x2)!=0 &&(Header.sublvl1&0x1E)!=0 && Tracker.q_lay[1][0] > 0 && Tracker.rig[1]>0&& Tof.beta>0";
//	string aftercutMC  = "abs(Tracker.q_inn-1)<0.3 && (trigpatt&0x2)!=0 &&(Header.sublvl1&0x1E)!=0 && Tracker.q_lay[1][0] > 0 && Tracker.rig[1]>0&& Tof.beta>0&& (Tof.z_nhit>=3)&&(Tof.z_like>0)&&((Tof.clsn[0]+Tof.clsn[2])<=4)";		
//	string beforecutDT  = "Header.pres_weight*(abs(Tracker.q_inn-1)<0.3 && (trigpatt&0x2)!=0 &&(Header.sublvl1&0x1E)!=0 && Tracker.q_lay[1][0] > 0 && Tracker.rig[1]>0&& Tof.beta>0)";		
//	string aftercutDT   = "Header.pres_weight*(abs(Tracker.q_inn-1)<0.3 && (trigpatt&0x2)!=0 &&(Header.sublvl1&0x1E)!=0 && Tracker.q_lay[1][0] > 0 && Tracker.rig[1]>0&& Tof.beta>0&& (Tof.z_nhit>=3)&&(Tof.z_like>0)&&((Tof.clsn[0]+Tof.clsn[2])<=4))";
	
	TH1F * hbeforeMC = new TH1F("hbeforeMC","hbeforeMC",30,-1,2);
	TH1F * hafterMC   = new TH1F("hafterMC","hafterMC",30,-1,2) ; 
	TH1F * hbeforeDT = new TH1F("hbeforeDT","hbeforeDT",30,-1,2);
	TH1F * hafterDT   = new TH1F("hafterDT","hafterDT",30,-1,2) ; 

	chainMC->Draw("log10(rig[1])>>hbeforeMC",beforecutMC.c_str());
	chainMC->Draw("log10(rig[1])>>hafterMC",aftercutMC.c_str());
	chainDT->Draw("log10(rig[1])>>hbeforeDT",beforecutDT.c_str());
	chainDT->Draw("log10(rig[1])>>hafterDT",aftercutDT.c_str());

	TH1F * MCEff     = hafterMC->Clone("RatioMC");

	MCEff->Divide(hbeforeMC);	

	MCEff->SetLineColor(2);	

	TH1F * DTEff     = hafterDT->Clone("RatioDT");

	//for trigger
	hbeforeDT->Add(hafterDT);

	DTEff->Divide(hbeforeDT);	

	DTEff->SetLineColor(1);	
	DTEff->Draw();	
	MCEff->Draw("same");	


}
