{
	pp=TProof::Open("");

	gROOT->ProcessLine(".L ./InputFileReader.h");

	TH1::SetDefaultSumw2();

	TChain * chainMC = InputFileReader("./FileListMC.txt","Event");
	TChain * chainDT = InputFileReader("./FileListDT.txt","Compact");

	chainMC->SetProof();
	chainDT->SetProof();


	string beforecutDT = "trk_ql1>0.7 && trk_ql1<1.6 && ((trigpatt & 0x2) != 0) && ((sublvl1&0x3E) !=0) && trk_rig[1]>0&&trk_chisqn[0]<100&&tof_beta>0&&trk_qinn>0.07&&trk_qinn<13&&tof_beta<2.85&&fabs(tof_qdw-1)<0.3&&abs(tof_qup-1)<0.3&&trk_rig[1]<1.5";
	
	//trigger	
/*	TH2F * hbeforeDT1 = new TH2F("hbeforeDT1","hbeforeDT1",500,0,4,500,-4,4);
	TH2F * hbeforeDT2 = new TH2F("hbeforeDT2","hbeforeDT2",500,0,4,500,-4,4);
	TH2F * hbeforeDT3 = new TH2F("hbeforeDT3","hbeforeDT3",500,0,4,500,-4,4);
*/

	TH2F * hbeforeDT1 = new TH2F("hbeforeDT1","hbeforeDT1",500,0.,1.1,500,-4,4);
	TH2F * hbeforeDT2 = new TH2F("hbeforeDT2","hbeforeDT2",500,0.,1.1,500,-4,4);
	TH2F * hbeforeDT3 = new TH2F("hbeforeDT3","hbeforeDT3",500,0.,1.1,500,-4,4);

//	TH3F * hbeforeDT1 = new TH3F("hbeforeDT1","hbeforeDT1",100,0.4,1.1,100,0,3,100,0.2,1.8);

/*	chainDT->Draw("trk_qinn :trk_rig[1]/tof_beta*(1-tof_beta^2)^0.5>>hbeforeDT1",beforecutDT.c_str());
	chainDT->Draw("tof_qdw  :trk_rig[1]/tof_beta*(1-tof_beta^2)^0.5>>hbeforeDT2",beforecutDT.c_str());
	chainDT->Draw("tof_qup  :trk_rig[1]/tof_beta*(1-tof_beta^2)^0.5>>hbeforeDT3",beforecutDT.c_str());
*/
	chainDT->Draw("trk_qinn :tof_beta>>hbeforeDT1",beforecutDT.c_str());
	chainDT->Draw("tof_qdw  :tof_beta>>hbeforeDT2",beforecutDT.c_str());
	chainDT->Draw("tof_qup  :tof_beta>>hbeforeDT3",beforecutDT.c_str());



//	chainDT->Draw("tof_qup  :trk_rig[1]:tof_beta>>hbeforeDT1",beforecutDT.c_str());



	TCanvas * c1 = new TCanvas("quinn");
	c1->cd();
	hbeforeDT1->Draw("col");;
	TCanvas * c2 = new TCanvas("qdown");
	c2->cd();
	hbeforeDT2->Draw("col");;
	TCanvas * c3 = new TCanvas("qup");
	c3->cd();
	hbeforeDT3->Draw("col");;



}
