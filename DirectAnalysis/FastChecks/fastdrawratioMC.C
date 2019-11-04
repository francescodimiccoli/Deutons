{
	pp=TProof::Open("");

	gROOT->ProcessLine(".L ~/Work/Deutons/DirectAnalysis/include/InputFileReader.h");

	TH1::SetDefaultSumw2();

	TChain * chainMC = InputFileReader("./FileListMC.txt","RTI");
	TChain * chainDT = InputFileReader("./FileListDT.txt","RTI");


	chainMC->SetProof();
	chainDT->SetProof();

	//trigger	
	string beforecutMC  = "abs(Tracker.q_inn-1)<0.3 && (trigpatt&0x2)!=0 &&(Header.sublvl1&0x1E)!=0 && Tracker.q_lay[1][0] > 0 && Tracker.rig[1]>0&& Tof.beta>0&& (Tof.z_nhit>=3)&&(Tof.z_like>0)&&((Tof.clsn[0]+Tof.clsn[2])<=4 && Header.ntrtrack==1 )";
	string beforecutDT  = "abs(Tracker.q_inn-1)<0.3 && (trigpatt&0x2)!=0 &&(Header.sublvl1&0x1E)!=0 && Tracker.q_lay[1][0] > 0 && Tracker.rig[1]>0&& Tof.beta>0&& (Tof.z_nhit>=3)&&(Tof.z_like>0)&&((Tof.clsn[0]+Tof.clsn[2])<=4 && Header.ntrtrack==1 )";		
	beforecutDT = "Header.pres_weight*(" + beforecutDT + "&&Tracker.rig[1]>1.2*Tracker.stoermer_cutoff[0])";

	TH2F * hbeforeMC = new TH2F("hbeforeMC","hbeforeMC",100,0,40,100,0,40);
	TH2F * hbeforeDT = new TH2F("hbeforeDT","hbeforeDT",100,0,40,100,0,40);


//	chainMC->Draw("(Tof.q_lay[0] + Tof.q_lay[1] ) / 2.0 : (Tof.q_lay[2] + Tof.q_lay[3] ) / 2.0>>hbeforeMC",beforecutMC.c_str());
//	chainDT->Draw("(Tof.q_lay[0] + Tof.q_lay[1] ) / 2.0 : (Tof.q_lay[2] + Tof.q_lay[3] ) / 2.0>>hbeforeDT",beforecutDT.c_str());
//	chainMC->Draw("Tof.beta : Tracker.rig[1]>>hbeforeMC",beforecutMC.c_str());
//	chainDT->Draw("Tof.beta : Tracker.rig[1]>>hbeforeDT",beforecutDT.c_str());
	chainMC->Draw("cf[0][2][1] :cf[1][2][1] >>hbeforeMC",beforecutMC.c_str());
	chainDT->Draw("cf[0][2][1] :cf[1][2][1] >>hbeforeDT",beforecutDT.c_str());

	hbeforeMC->Scale(1/hbeforeMC->Integral());
	hbeforeDT->Scale(1/hbeforeDT->Integral());

	TH2F * Ratio = hbeforeDT->Clone("Ratio");
	Ratio -> Add(hbeforeMC,-1);
	Ratio -> Divide(hbeforeMC);

	TCanvas * c1 = new TCanvas("Plot MC");
	c1->cd();
	c1->SetLogz();
	hbeforeMC->Draw("colz");

	TCanvas * c2 = new TCanvas("Plot DT");
	c2->cd();
	c2->SetLogz();
	hbeforeDT->Draw("colz");
	
	TCanvas * c3 = new TCanvas("Ratio DT/MC");
	c3->cd();
	Ratio->GetZaxis()->SetRangeUser(-2,2);
	Ratio->Draw("colz");



}
