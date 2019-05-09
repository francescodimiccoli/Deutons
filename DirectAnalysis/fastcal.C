#include "InputFileReader.h"

TH1D * GetMeans(TH2F * h1){
	string nameh1 = h1->GetName();
	nameh1 = nameh1 + "_1";
	TH1D * h1_1 = (TH1D*)gDirectory->Get(nameh1.c_str());	
	return h1_1;

}

TSpline3 * ExtractSpline(TH1D * Cal){

	std::vector<double> x;
	std::vector<double> y;
	int j=0;
        for(int i=0;i<Cal->GetNbinsX();i++){
             if(Cal->GetBinContent(i+1)!=0){
			x.push_back(Cal->GetBinCenter(i+1));
                	y.push_back(Cal->GetBinContent(i+1));
			cout<<x[j]<<" "<<y[j]<<endl;
			j++;
		}
        }
	double X[x.size()];
	double Y[y.size()];
	
	for(int i=0;i<x.size();i++){
		X[i] = x[i];
		Y[i] = y[i];
	}

        TSpline3 * Calibration = new TSpline3(Cal->GetName(),X,Y,x.size());
        Calibration->SetName(Cal->GetName());
        return Calibration;

}


void DrawSlices(TCanvas * c, TH2F * hMC, TH2F *hDT);

int fastcal(){

	TProof * pp=TProof::Open("");

	//gROOT->ProcessLine(".L ./InputFileReader.h");


	TH1::SetDefaultSumw2();

	TChain * chainMC = InputFileReader("./FileListMC.txt","Event");
	TChain * chainDT = InputFileReader("./FileListDT.txt","Compact");


	chainMC->SetProof();
	chainDT->SetProof();

	//cuts
	string beforecutMC = "q_lay[1][0]>0 && q_lay[1][0]<1.8 && ((trigpatt & 0x2) != 0) && ((sublvl1&0x3E) !=0) && rig[1][1]>0&&chisqn[1][0]<100&&beta>0";
	string beforecutDT = "trk_ql1>0 && trk_ql1< 1.8 && ((trigpatt & 0x2) != 0) && ((sublvl1&0x3E) !=0) && trk_rig[1]>0&&trk_chisqn[0]<100&&tof_beta>0";

	//drawings
	TH2F * hMCqinn = new TH2F("hMCqinn","hMCqinn",100,0.4,1.1,100,0,3);
	TH2F * hDTqinn = new TH2F("hDTqinn","hDTqinn",100,0.4,1.1,100,0,3);
	chainMC->Draw("q_inn[0]:beta>>hMCqinn",beforecutMC.c_str());
	chainDT->Draw("trk_qinn:tof_beta>>hDTqinn",beforecutDT.c_str());
	TH2F * hMCqutof = new TH2F("hMCqutof","hMCqutof",100,0.4,1.1,100,0.6,1.7);
	TH2F * hDTqutof = new TH2F("hDTqutof","hDTqutof",100,0.4,1.1,100,0.6,1.7);
	chainMC->Draw("(Tof.q_lay[0]+Tof.q_lay[1])/2:beta>>hMCqutof",beforecutMC.c_str());
	chainDT->Draw("tof_qup:tof_beta>>hDTqutof",beforecutDT.c_str());
	TH2F * hMCqltof = new TH2F("hMCqltof","hMCqltof",100,0.4,1.1,100,0.6,1.7);
	TH2F * hDTqltof = new TH2F("hDTqltof","hDTqltof",100,0.4,1.1,100,0.6,1.7);
	chainMC->Draw("(Tof.q_lay[2]+Tof.q_lay[3])/2:beta>>hMCqltof",beforecutMC.c_str());
	chainDT->Draw("tof_qdw:tof_beta>>hDTqltof",beforecutDT.c_str());
	TH2F * hMCql1 = new TH2F("hMCql1","hMCql1",100,0.4,1.1,100,0.6,1.7);
	TH2F * hDTql1 = new TH2F("hDTql1","hDTql1",100,0.4,1.1,100,0.6,1.7);
	chainMC->Draw("q_lay[1][0]:beta>>hMCql1",beforecutMC.c_str());
	chainDT->Draw("trk_ql1:tof_beta>>hDTql1",beforecutDT.c_str());






	//plottings
	TCanvas * c1 = new TCanvas("qinn_MC");
	c1->cd();
	hMCqinn->Draw("col");
	TCanvas * c2 = new TCanvas("qinn_DT");
	c2->cd();
	hDTqinn->Draw("col");

	TCanvas * c1_ = new TCanvas("qutof_MC");
	c1_->cd();
	hMCqutof->Draw("col");
	TCanvas * c2_ = new TCanvas("qutof_DT");
	c2_->cd();
	hDTqutof->Draw("col");

	TCanvas * c1__ = new TCanvas("qltof_MC");
	c1__->cd();
	hMCqutof->Draw("col");
	TCanvas * c2__ = new TCanvas("qltof_DT");
	c2__->cd();
	hDTqltof->Draw("col");

	TCanvas * c1___ = new TCanvas("ql1_MC");
	c1___->cd();
	hMCql1->Draw("col");
	TCanvas * c2___ = new TCanvas("ql1_DT");
	c2___->cd();
	hDTql1->Draw("col");



	TFile * Out = TFile::Open("ChargeCalibrations.root","RECREATE");
	TCanvas * c3 = new TCanvas("qinn Means");
	DrawSlices(c3,hMCqinn,hDTqinn);
	TCanvas * c4 = new TCanvas("qutof Means");
	DrawSlices(c4,hMCqutof,hDTqutof);
	TCanvas * c5 = new TCanvas("qltof Means");
	DrawSlices(c5,hMCqltof,hDTqltof);
	TCanvas * c6 = new TCanvas("ql1 Means");
	DrawSlices(c6,hMCql1,hDTql1);

	Out->Save();

	//test

	TFile * Out2 = TFile::Open("ChargeCalibrations.root","READ");
	TSpline3 * qinn_cal = (TSpline3*) Out2->Get("hDTqinn_1");
	
	TH2F * hMCqinn_corr = new TH2F("hMCqinn_corr","hMCqinn_corr",100,0.4,1.1,100,0,3);

	float qinn=0;
	float beta=0;
	
	chainMC->SetBranchAddress("q_inn[0]",&qinn);		
	chainMC->SetBranchAddress("beta",&beta);		

	for(int i=0; i< chainMC->GetEntries();i++) {
		chainMC->GetEvent(i);
		cout<<qinn<<" "<<beta<<" "<<qinn_cal->Eval(0.7)<<endl;
	}



	return 0;
}


void DrawSlices(TCanvas * c, TH2F * hMC, TH2F *hDT){
	//slices
	hMC->FitSlicesY(0,10,hMC->GetNbinsX(),0,"G4"); 	
	hDT->FitSlicesY(0,10,hMC->GetNbinsX(),0,"G4"); 	

	TH1D * MeansMC = GetMeans(hMC);
	TH1D * MeansDT = GetMeans(hDT);

	
	c->cd(1);
	MeansMC -> SetMarkerStyle(8);
	MeansMC -> SetMarkerColor(2);
	MeansDT -> SetMarkerStyle(8);
	MeansDT -> SetMarkerColor(1);
	MeansMC->Draw();
	MeansDT->Draw("same");

	TH1D * cal = (TH1D*) MeansDT->Clone();
	cal->Divide(MeansMC);
	cal->Write();

	TSpline3 * calibration = ExtractSpline(cal);
	calibration->Write();
}


