{
	pp=TProof::Open("");

	gROOT->ProcessLine(".L ./InputFileReader.h");

	TH1::SetDefaultSumw2();

	TChain * chainMC = InputFileReader("./FileListMC.txt","Event");
	TChain * chainDT = InputFileReader("./FileListDT.txt","Compact");
	TChain * chainRTI= InputFileReader("./FileListDT.txt","RTI");

	chainMC->SetProof();
	chainDT->SetProof();

	
	//rich efficiencies from data:
	cout<<"************************************* RICH EFFICIENCY FROM DATA ****************************************"<<endl;

	//naf
	TH1F * hbeforeDTNaF = new TH1F("hbeforeDTNaF","hbeforeDTNaF",100,0.5,5);
	TH1F * hafterDTNaF   = new TH1F("hafterDTNaF","hafterDTNaF",100,0.5,5) ; 

	string beforecutDTNaF = "tof_beta>0.3 && (trigpatt&0x2)!=0&&(sublvl1&0x3E)!=0 && trk_rig[1]>0 &&fabs(trk_ql1-1)>0.99 &&fabs(trk_ql2-1)>0.99 && fabs(trk_qinn-1)<0.4";
	string aftercutDTNaF  = "tof_beta>0.3 && (trigpatt&0x2)!=0&&(sublvl1&0x3E)!=0 && trk_rig[1]>0 &&fabs(trk_ql1-1)>0.99 &&fabs(trk_ql2-1)>0.99 &&  fabs(trk_qinn-1)<0.4 && rich_select==1 && rich_bdt>0.26 ";

	chainDT->Draw("(trk_rig[1]^2+0.938^2)^0.5 - 0.938>>hbeforeDTNaF",beforecutDTNaF.c_str());
	chainDT->Draw("(trk_rig[1]^2+0.938^2)^0.5 - 0.938>>hafterDTNaF" ,aftercutDTNaF.c_str());

	TCanvas * r = new TCanvas("rich naf eff");
	TH1F * EffNaF = (TH1F*) hafterDTNaF->Clone();
	EffNaF->Divide(hbeforeDTNaF);
	EffNaF->Draw();

	//agl
	TH1F * hbeforeDTAgl = new TH1F("hbeforeDTAgl","hbeforeDTAgl",100,2.5,10);
	TH1F * hafterDTAgl   = new TH1F("hafterDTAgl","hafterDTAgl",100,2.5,10) ; 

	string beforecutDTAgl = "tof_beta>0.3 && (trigpatt&0x2)!=0&&(sublvl1&0x3E)!=0 && trk_rig[1]>0 &&fabs(trk_ql1-1)>0.99 &&fabs(trk_ql2-1)>0.99 && fabs(trk_qinn-1)<0.4";
	string aftercutDTAgl  = "tof_beta>0.3 && (trigpatt&0x2)!=0&&(sublvl1&0x3E)!=0 && trk_rig[1]>0 &&fabs(trk_ql1-1)>0.99 &&fabs(trk_ql2-1)>0.99 &&  fabs(trk_qinn-1)<0.4 && rich_select==2 && rich_bdt>0.25 ";

	chainDT->Draw("(trk_rig[1]^2+0.938^2)^0.5 - 0.938>>hbeforeDTAgl",beforecutDTAgl.c_str());
	chainDT->Draw("(trk_rig[1]^2+0.938^2)^0.5 - 0.938>>hafterDTAgl" ,aftercutDTAgl.c_str());

	TCanvas * l = new TCanvas("rich agl eff");
	TH1F * EffAgl = (TH1F*) hafterDTAgl->Clone();
	EffAgl->Divide(hbeforeDTAgl);
	EffAgl->Draw();


	cout<<"************************************* ACCEPTANCE  ****************************************"<<endl;


	//tof
	string beforecutMC = "";
	string aftercutMC  = "Tof.beta>0.3&&(trigpatt&0x2)!=0&&(Header.sublvl1&0x3E)!=0 && Tracker.rig[1][1]!=0&&fabs(Tracker.q_lay[1][0]-1)<0.99&&Tracker.q_clu_status[1][0]==0&&fabs(Tracker.q_inn[0]-1)<4";

	TH1F * hbeforeMC = new TH1F("hbeforeMC","hbeforeMC",100,0.05,2);
	TH1F * hafterMC   = new TH1F("hafterMC","hafterMC",100,0.05,2) ; 

	chainMC->Draw("(rig[1]^2+0.938^2)^0.5 - 0.938>>hafterMC",aftercutMC.c_str());


	float normalization = chainMC->GetEntries()*(pow(0.0092,-1));
	float range = log(100)-log(0.5);

	for(int i=0;i<hbeforeMC->GetNbinsX();i++){
		float rigmin=pow((pow(hbeforeMC->GetBinLowEdge(i)+0.938,2)-0.938*0.938),0.5)  ;
		float rigmax=pow((pow(hbeforeMC->GetBinLowEdge(i+1)+0.938,2)-0.938*0.938),0.5)  ;
		hbeforeMC->SetBinContent(i+1,normalization*(log(rigmax)-log(rigmin))/range);
	}
	TH1F * MCEff     = (TH1F*) hafterMC->Clone("RatioMC");
	MCEff->Divide(hbeforeMC);	
	MCEff->Scale(30);
	MCEff->Rebin(3);

	//naf
	string beforecutMCNaF = "";
	string aftercutMCNaF  = "Tof.beta>0.3&&(trigpatt&0x2)!=0&&(Header.sublvl1&0x3E)!=0 && Tracker.rig[1][1]!=0&&fabs(Tracker.q_lay[1][0]-1)<0.99&&Tracker.q_clu_status[1][0]==0&&fabs(Tracker.q_inn[0]-1)<0.4";

	TH1F * hbeforeMCNaF = new TH1F("hbeforeMCNaF","hbeforeMCNaF",100,0.5,5);
	TH1F * hafterMCNaF   = new TH1F("hafterMCNaF","hafterMCNaF",100,0.5,5) ; 

	chainMC->Draw("(rig[1]^2+0.938^2)^0.5 - 0.938>>hafterMCNaF",aftercutMCNaF.c_str());


	normalization = chainMC->GetEntries()*(pow(0.0092,-1));
	range = log(100)-log(0.5);

	for(int i=0;i<hbeforeMC->GetNbinsX();i++){
		float rigmin=pow((pow(hbeforeMCNaF->GetBinLowEdge(i)+0.938,2)-0.938*0.938),0.5)  ;
		float rigmax=pow((pow(hbeforeMCNaF->GetBinLowEdge(i+1)+0.938,2)-0.938*0.938),0.5)  ;
		hbeforeMCNaF->SetBinContent(i+1,normalization*(log(rigmax)-log(rigmin))/range);
	}
	TH1F * MCNaFEff     = (TH1F*) hafterMCNaF->Clone("RatioMCNaF");
	MCNaFEff->Divide(hbeforeMCNaF);	
	MCNaFEff->Scale(30);
	MCNaFEff->Multiply(EffNaF);
	MCNaFEff->Rebin(3);


	//agl
	string beforecutMCAgl = "";
	string aftercutMCAgl  = "Tof.beta>0.3&&(trigpatt&0x2)!=0&&(Header.sublvl1&0x3E)!=0 && Tracker.rig[1][1]!=0&&fabs(Tracker.q_lay[1][0]-1)<0.99&&Tracker.q_clu_status[1][0]==0&&fabs(Tracker.q_inn[0]-1)<0.4";

	TH1F * hbeforeMCAgl = new TH1F("hbeforeMCAgl","hbeforeMCAgl",100,2.5,10);
	TH1F * hafterMCAgl   = new TH1F("hafterMCAgl","hafterMCAgl",100,2.5,10) ; 

	chainMC->Draw("(rig[1]^2+0.938^2)^0.5 - 0.938>>hafterMCAgl",aftercutMCAgl.c_str());


	normalization = chainMC->GetEntries()*(pow(0.0092,-1));
	range = log(100)-log(0.5);

	for(int i=0;i<hbeforeMC->GetNbinsX();i++){
		float rigmin=pow((pow(hbeforeMCAgl->GetBinLowEdge(i)+0.938,2)-0.938*0.938),0.5)  ;
		float rigmax=pow((pow(hbeforeMCAgl->GetBinLowEdge(i+1)+0.938,2)-0.938*0.938),0.5)  ;
		hbeforeMCAgl->SetBinContent(i+1,normalization*(log(rigmax)-log(rigmin))/range);
	}
	TH1F * MCAglEff     = (TH1F*) hafterMCAgl->Clone("RatioMCAgl");
	MCAglEff->Divide(hbeforeMCAgl);	
	MCAglEff->Scale(30);
	MCAglEff->Multiply(EffAgl);
	MCAglEff->Rebin(3);


	//draw
	TH2F * Frame = new TH2F("Frame","Frame",1000,0,15,1000,1e-4,1);
	TCanvas *f=new TCanvas("Acceptance");
	gPad->SetLogy();
	gPad->SetLogx();
	Frame->Draw();
	MCEff->SetLineColor(2);
	MCNaFEff->SetLineColor(3);
	MCAglEff->SetLineColor(4);
	MCEff->Draw("same");
	MCNaFEff->Draw("same");
	MCAglEff->Draw("same");

	cout<<"********************************* EXPOSURE TIME ********************************"<<endl;
	TH2F * ExposurevsCut = new TH2F("ExposurevsCut","ExposurevsCut",500,0,50,500,0,1.5);
	chainRTI->Draw("lf:cf[0][2][1]>>ExposurevsCut");
	TH1F * ExposureTime = new TH1F("ExposureTime","ExposureTime",500,0,50);


	for(int i=0;i<ExposureTime->GetNbinsX();i++){
		float R=ExposureTime->GetBinLowEdge(i);
		float time = 0;
			for(int j=0;j<ExposurevsCut->GetNbinsX();j++){
				float Rcut = 1.2*ExposurevsCut->GetXaxis()->GetBinLowEdge(j);
				if(Rcut<R)
					for(int k=0;k<ExposurevsCut->GetNbinsY();k++)
						time+=ExposurevsCut->GetBinContent(j,k)*ExposurevsCut->GetYaxis()->GetBinCenter(k);
			}
		ExposureTime->SetBinContent(i,time);
	}



	TCanvas * expo = new TCanvas("expo");
	ExposureTime->Draw();
}
