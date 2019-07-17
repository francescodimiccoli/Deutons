{
	pp=TProof::Open("");

	gROOT->ProcessLine(".L ./InputFileReader.h");

	TH1::SetDefaultSumw2();

	TChain * chainMC = InputFileReader("./FileListMC.txt","Event");
	TChain * chainDT = InputFileReader("./FileListDT.txt","Compact");

	chainMC->SetProof();
	chainDT->SetProof();


	string beforecutDT = "trk_q_lay[0]>0.7 && trk_q_lay[0]<1.6 && ((trigpatt & 0x2) != 0) && ((sublvl1&0x3E) !=0) && trk_rig_kalman[1]>0&&trk_chisqn[0]<10&&tof_beta>0.5&&trk_qinn>0.7&&trk_qinn<1.3&&tof_beta<0.78&&fabs(tof_qdw-1)<0.3&&abs(tof_qup-1)<0.3";
	string beforecutDTfragm = "trk_q_lay[0]>1.7 && trk_q_lay[0]<2.6 && ((trigpatt & 0x2) != 0) && ((sublvl1&0x3E) !=0) && trk_rig_kalman[1]>0&&trk_chisqn[0]<10&&tof_beta>0.5&&trk_qinn>0.7&&trk_qinn<1.3&&tof_beta<0.78&&fabs(tof_qdw-1)<0.3&&abs(tof_qup-1)<0.3";

	//string beforecutDT = "trk_q_lay[0]>0.7 && trk_q_lay[0]<1.6 && ((trigpatt & 0x2) != 0) && ((sublvl1&0x3E) !=0) && trk_chisqn[0]<10&&trk_qinn>0.7&&trk_qinn<1.3&&fabs(tof_qdw-1)<0.3&&abs(tof_qup-1)<0.3";
	
/*
 *fragm charge
	TH2F * hbeforeDT1 = new TH2F("hbeforeDT1","hbeforeDT1",500,0.,4,500,0,2);
	TH2F * hbeforeDT2 = new TH2F("hbeforeDT2","hbeforeDT2",500,0.,4,500,0,2);
	TH2F * hbeforeDT3 = new TH2F("hbeforeDT3","hbeforeDT3",500,0.,4,500,0,2);
	chainDT->Draw("tof_qup  :trk_rig_kalman[1]/tof_beta*(1-tof_beta^2)^0.5>>hbeforeDT1",beforecutDT.c_str());
	chainDT->Draw("trk_qinn  :trk_rig_kalman[1]/tof_beta*(1-tof_beta^2)^0.5>>hbeforeDT2",beforecutDT.c_str());
	chainDT->Draw("trk_qinn  :trk_rig_kalman[1]/tof_beta*(1-tof_beta^2)^0.5>>hbeforeDT3",beforecutDT.c_str());

	TCanvas * c1 = new TCanvas("quinn");
	c1->cd();
	hbeforeDT1->Draw("col");;
	TCanvas * c2 = new TCanvas("qdown");
	c2->cd();
	hbeforeDT2->Draw("col");;
	TCanvas * c3 = new TCanvas("qup");
	c3->cd();
	hbeforeDT3->Draw("col");;


*/

	//fragm mass
//	TH1F * hbeforeDT1 = new TH1F("hbeforeDT1","hbeforeDT1",500,0.,4);
//	TH1F * hbeforeDT2 = new TH1F("hbeforeDT2","hbeforeDT2",500,0.,4);
	//chainDT->Draw("trk_rig_kalman[1]/tof_beta*(1-tof_beta^2)^0.5>>hbeforeDT2",beforecutDTfragm.c_str());
//	chainDT->Draw("trk_rig_kalman[1]/tof_beta*(1-tof_beta^2)^0.5>>hbeforeDT1",beforecutDT.c_str());
	

	//fragm rigvsbeta
//	TH2F * hbeforeDT1 = new TH2F("hbeforeDT1","hbeforeDT1",100,0.,4,100,0,4);
//	chainDT->Draw("trk_rig_kalman[1]:trk_rig_choutko[1]>>hbeforeDT1",beforecutDT.c_str());
	
	//fragm rig
	//TH1F * hbeforeDT1 = new TH1F("hbeforeDT1","hbeforeDT1",500,0.,4);
	//TH1F * hbeforeDT2 = new TH1F("hbeforeDT2","hbeforeDT2",500,0.,4);
	//chainDT->Draw("trk_rig_kalman[1]>>hbeforeDT2",beforecutDTfragm.c_str());
	//chainDT->Draw("trk_rig_kalman[1]>>hbeforeDT1",beforecutDT.c_str());

	//fragm qL1vsmass
        beforecutDT = "((trigpatt & 0x2) != 0) && ((sublvl1&0x3E) !=0) && trk_rig_kalman[1]>0&&trk_chisqn[0]<10&&tof_beta>0.5&&trk_qinn>0.7&&trk_qinn<1.3&&tof_beta<0.78&&fabs(tof_qdw-1)<0.3&&abs(tof_qup-1)<0.3&&trk_rig_kalman[1]/tof_beta*(1-tof_beta^2)^0.5<1.25";
	
	TH2F * hbeforeDT1 = new TH2F("hbeforeDT1","hbeforeDT1",200,0.,2,200,0,4);
	chainDT->Draw("trk_q_lay[0]:trk_rig_kalman[1]/tof_beta*(1-tof_beta^2)^0.5>>hbeforeDT1",beforecutDT.c_str());
		
	TCanvas * c1 = new TCanvas("mass");
//	hbeforeDT1->Draw("col");
	//hbeforeDT2->Draw("same");

}
