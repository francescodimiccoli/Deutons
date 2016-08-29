#include "TH2.h"
#include "TH3.h"
#include <TVector3.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cstring>
#include <vector>
#include "TMath.h"
#include "TCanvas.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
#include <stdarg.h>
#include <TSpline.h>
#include "TFractionFitter.h"
#include "THStack.h"
#include "TNtuple.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF2.h"
#include <math.h>
#include <cstring>
#include "Functions_auto.cpp"
#include <TMVA/Factory.h>
#include <TMVA/Reader.h>
#include <TMVA/Tools.h>

using namespace std;

float NAnticluster=0;
float Clusterinutili=0;
float ClusterTOFusati=0;
float ClusterTOFtotali=0;
float NTRDSegments=0;
float Massa=0;
float Beta=0;
float Betacorr=0;
float DiffR=0;
float ProbQ=0;
float layernonusati=0;
float residuoX=0;
float residuoY=0;
float fuoriX=0;
float fuoriY=0;
float CaricaTOF=0;
float CaricaTRD=0;
float CaricaTrack=0;
float EdepTOF=0;
float EdepTRD=0;
float EdepTrack=0;
float EdepTOFteo=0;
float R=0;
float probP=0;
double geomag[12]= {0,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.3};
int zone=0;
float R_corr=0;
float probD=0;
float Chisquare=0;
float Momentogen=0;
float Betagen=0;
float clusterTOFfuori=0;
float clusterTrackfuori=0;
float clusterTRDfuori=0;
float DistTrack=0;
float Richtotused=0;
float RichPhEl=0;
float DistTrack_p=0;
float ThetaTrack=0;
float PhiProjTrack=0;
float BetaRICH=0;
float PhiTrack=0;
float DistD=0;
float EdepTOFU=0;
float EdepTOFD=0;
float ThetaTOF=0;
float PhiTOF=0;
float PhiTRD=0;
float ThetaTRD=0;
float Massagen=0;
float LDiscr=0;
float LDiscriminant=0;
float Latitude=0;
float RminTOF=0; float RminTrack=0; float RminTRD=0;
float P_ID=0;
float TOF_Up_Down=0;
float EdepL1=0;
float DiffTrackEdep=0;
float Track_Up_Down=0;
float BDTresponse=0;
float S=0;
float B=0;
float IsCharge1=0;
float Rcutoff=0;
float MassLimit=0.5;
float Cutmask=0;
float MC_type=0;
TMVA::Reader *reader;

void BDTreader();


float ReturnMass_Gen( float MC_type)
{
   float Mass_gen=0;
   if ( ( ( (int) MC_type) &0xFF    ) >0)      Mass_gen = 0.938;
   if ( ( ( (int) MC_type) &0xFF00  ) >0)      Mass_gen = 1.875;
   if ( ( ( (int) MC_type) &0xFF0000) >0)      Mass_gen = 3.725;
   return Mass_gen;
}


int main()
{
	cout<<"*********************** CALIB. READING *********************"<<endl;

        string nomecal=("/storage/gpfs_ams/ams/users/fdimicco/Deutons/CodesforAnalysis/CALIBRAZIONI/2012_01.root");
        TFile *calib = TFile::Open(nomecal.c_str());
        cout<<"calibrazione: "<<calib<<endl;
        Rig = (TSpline3 *) calib->Get("Fit Results/Splines/Rig");
        beta = (TSpline3 *) calib->Get("Fit Results/Splines/beta");
        betaNaF = (TF1 *) calib->Get("Fit Results/Splines/SigmaInvBetaNaF_spl");
        betaAgl = (TF1 *) calib->Get("Fit Results/Splines/SigmaInvBetaAgl_spl");
        eL1 = (TSpline3 *) calib->Get("Fit Results/Splines/eL1");
        etofu =  (TSpline3 *) calib->Get("Fit Results/Splines/etofu");
        etrack =  (TSpline3 *) calib->Get("Fit Results/Splines/etrack");
        etofd =  (TSpline3 *) calib->Get("Fit Results/Splines/etofd");
        EdepL1beta =  (TSpline3 *) calib->Get("Fit Results/Splines/EdepL1beta");
        EdepTOFbeta =  (TSpline3 *) calib->Get("Fit Results/Splines/EdepTOFbeta");
        EdepTrackbeta =  (TSpline3 *) calib->Get("Fit Results/Splines/EdepTrackbeta");
        EdepTOFDbeta =  (TSpline3 *) calib->Get("Fit Results/Splines/EdepTOFDbeta");
        Corr_L1 =  (TSpline3 *) calib->Get("Fit Results/Splines/Corr_L1");
        Corr_TOFU =  (TSpline3 *) calib->Get("Fit Results/Splines/Corr_TOFU");
        Corr_Track =  (TSpline3 *) calib->Get("Fit Results/Splines/Corr_Track");
        Corr_TOFD =  (TSpline3 *) calib->Get("Fit Results/Splines/Corr_TOFD");
        cout<<Rig<<" "<<beta<<" "<<" "<<betaNaF<<" "<<betaAgl<<" "<<eL1<<" "<<etofu<<" "<<etrack<<" "<<etofd<<" "<<EdepL1beta<<" "<<EdepTOFbeta<<" "<<EdepTrackbeta<<" "<<EdepTOFDbeta<<" "<<Corr_L1<<" "<<Corr_TOFU<<" "<<Corr_Track<<" "<<Corr_TOFD<<endl;
        cout<<"******************************"<<endl;

	TFile *file1 =TFile::Open("../Risultati/QualityTraining/RisultatiMC.root");
	TFile *file2 =TFile::Open("../Risultati/QualityTraining/RisultatiDATI.root");
	TNtuple *ntupla1=(TNtuple*)file1->Get("grandezzequal");
	TNtuple *ntupla2=(TNtuple*)file2->Get("grandezzequal");
	//grandezzequal

	ntupla1->SetBranchAddress("NAnticluster",&NAnticluster);
	ntupla1->SetBranchAddress("NTRDSegments",&NTRDSegments);
	ntupla1->SetBranchAddress("Velocity",&Beta);
	ntupla1->SetBranchAddress("R",&R);
	ntupla1->SetBranchAddress("MC_type",&MC_type);
	ntupla1->SetBranchAddress("ClusterTOFusati",&ClusterTOFusati);
	ntupla1->SetBranchAddress("ClusterTOFtotali",&ClusterTOFtotali);
	ntupla1->SetBranchAddress("Clusterinutili",&Clusterinutili);
	ntupla1->SetBranchAddress("Cutmask",&Cutmask);
	ntupla1->SetBranchAddress("DiffR",&DiffR);
	ntupla1->SetBranchAddress("ProbQ",&ProbQ);
	ntupla1->SetBranchAddress("layernonusati",&layernonusati);
	ntupla1->SetBranchAddress("residuoX",&residuoX);
	ntupla1->SetBranchAddress("residuoY",&residuoY);
	ntupla1->SetBranchAddress("fuoriX",&fuoriX);
	ntupla1->SetBranchAddress("fuoriY",&fuoriY);
	ntupla1->SetBranchAddress("Chisquare",&Chisquare);
	ntupla1->SetBranchAddress("EdepTOFU",&EdepTOFU);	
	ntupla1->SetBranchAddress("EdepTOFD",&EdepTOFD);
	ntupla1->SetBranchAddress("DistD",&DistD);
	ntupla1->SetBranchAddress("EdepL1",&EdepL1);
	ntupla1->SetBranchAddress("Track_Up_Down",&Track_Up_Down);
	ntupla1->SetBranchAddress("DiffTrackEdep",&DiffTrackEdep);
	ntupla1->SetBranchAddress("Richtotused",&Richtotused);
        ntupla1->SetBranchAddress("RichPhEl",&RichPhEl);
	ntupla1->SetBranchAddress("Momentogen",&Momentogen);
	ntupla1->SetBranchAddress("BDT_response",&BDTresponse);
	ntupla1->SetBranchAddress("IsCharge1",&IsCharge1);


	ntupla2->SetBranchAddress("NAnticluster",&NAnticluster);
	ntupla2->SetBranchAddress("NTRDSegments",&NTRDSegments);
	ntupla2->SetBranchAddress("Velocity",&Beta);
	ntupla2->SetBranchAddress("R",&R);
	ntupla2->SetBranchAddress("Rcutoff",&Rcutoff);
	ntupla2->SetBranchAddress("Latitude",&Latitude);
	ntupla2->SetBranchAddress("Cutmask",&Cutmask);
	ntupla2->SetBranchAddress("ClusterTOFusati",&ClusterTOFusati);
	ntupla2->SetBranchAddress("ClusterTOFtotali",&ClusterTOFtotali);
	ntupla2->SetBranchAddress("Clusterinutili",&Clusterinutili);
	ntupla2->SetBranchAddress("DiffR",&DiffR);
	ntupla2->SetBranchAddress("ProbQ",&ProbQ);
	ntupla2->SetBranchAddress("layernonusati",&layernonusati);
	ntupla2->SetBranchAddress("residuoX",&residuoX);
	ntupla2->SetBranchAddress("residuoY",&residuoY);
	ntupla2->SetBranchAddress("fuoriX",&fuoriX);
	ntupla2->SetBranchAddress("fuoriY",&fuoriY);
	ntupla2->SetBranchAddress("Chisquare",&Chisquare);
	ntupla2->SetBranchAddress("EdepTOFU",&EdepTOFU);  
	ntupla2->SetBranchAddress("EdepTOFD",&EdepTOFD);
	ntupla2->SetBranchAddress("DistD",&DistD);
	ntupla2->SetBranchAddress("EdepL1",&EdepL1);
	ntupla2->SetBranchAddress("Track_Up_Down",&Track_Up_Down);
	ntupla2->SetBranchAddress("DiffTrackEdep",&DiffTrackEdep);
	ntupla2->SetBranchAddress("Richtotused",&Richtotused);
        ntupla2->SetBranchAddress("RichPhEl",&RichPhEl);
	ntupla2->SetBranchAddress("BDT_response",&BDTresponse);
	ntupla2->SetBranchAddress("IsCharge1",&IsCharge1);


	TCanvas *f=new TCanvas("R vs Rgen (BAD)");
	TCanvas *d_tris=new TCanvas("Sigma distance R vs Beta");
	TCanvas *c1[9];
	TH1F *grafico1[9]; 
	TH1F *grafico2[9]; 
	TH1F *grafico3[9];
	TH1F *graficoz[9][11];
	TH1F *RisoluzioniBeta[30];
	TH1F *RisoluzioniR[24];
	TH2F *RvsRgen_bad =new TH2F("","",500,0,5,500,0,5);
	TH2F *RvsRgen_sigma =new TH2F("","",500,0,5,500,0,5);
	TH2F *sigmagen_bad =new TH2F("","",500,0,30,500,0,30);
	for(int i=0;i<RvsRgen_sigma->GetNbinsX();i++)
                        for(int j=0;j<RvsRgen_sigma->GetNbinsY();j++){
                                R=RvsRgen_sigma->GetXaxis()->GetBinCenter(i);
				Momentogen=RvsRgen_sigma->GetYaxis()->GetBinCenter(j);
				if(R<Momentogen-1*pow(Momentogen,2)*Rig->Eval(Momentogen)||R>Momentogen+1*pow(Momentogen,2)*Rig->Eval(Momentogen))
                                                        RvsRgen_sigma->SetBinContent(i,j,1);
				if(R<Momentogen-6*pow(Momentogen,2)*Rig->Eval(Momentogen)||R>Momentogen+6*pow(Momentogen,2)*Rig->Eval(Momentogen))
                                                        RvsRgen_sigma->SetBinContent(i,j,4);

				}

	string Variables[9]={"N. Anti-clusters","Unused TOF Clusters","|Rup-Rdown|:R","Unused Tracker layers","Tracker: Y Hits without X","Track Chi^2","RICH Hits: tot - used","Rich Photoelectrons","|E. dep.(tot)-E.dep.(track)|"};
	grafico1[0]=new TH1F("Anticl",Variables[0].c_str(),5,0,10);
	grafico1[1]=new TH1F("Un. TOF",Variables[1].c_str(),20,0,20);
	grafico1[2]=new TH1F("R Diff.",Variables[2].c_str(),50,0,1);
	grafico1[3]=new TH1F("Un. layers",Variables[3].c_str(),10,0,10);
	grafico1[4]=new TH1F("fuori X",Variables[4].c_str(),10,0,10);
	grafico1[5]=new TH1F("Track Chi^2",Variables[5].c_str(),1000,0,25);
	grafico1[6]=new TH1F("RICH Hits: tot - used",Variables[6].c_str(),100,0,30);
        grafico1[7]=new TH1F("Rich Photoelectrons",Variables[7].c_str(),300,0,30);
        grafico1[8]=new TH1F("Track Edep: Tot - Track",Variables[8].c_str(),300,0,30);

	grafico2[0]=new TH1F("Anticl-good",Variables[0].c_str(),5,0,10);
	grafico2[1]=new TH1F("Un. TOF-good",Variables[1].c_str(),20,0,20);
	grafico2[2]=new TH1F("R Diff.-good",Variables[2].c_str(),50,0,1);
	grafico2[3]=new TH1F("Un. layers-good",Variables[3].c_str(),10,0,10);
	grafico2[4]=new TH1F("fuori X-good",Variables[4].c_str(),10,0,10);
	grafico2[5]=new TH1F("Track Chi^2-good",Variables[5].c_str(),1000,0,25);
	grafico2[6]=new TH1F("RICH Hits: tot - used-good",Variables[6].c_str(),100,0,30);
        grafico2[7]=new TH1F("Rich Photoelectrons-good",Variables[7].c_str(),300,0,30);	
	grafico2[8]=new TH1F("Track Edep: Tot - Track-good",Variables[8].c_str(),300,0,30);

	grafico3[0]=new TH1F("Anticl-He",Variables[0].c_str(),5,0,10);
	grafico3[1]=new TH1F("Un. TOF-He",Variables[1].c_str(),20,0,20);
	grafico3[2]=new TH1F("R Diff.-He",Variables[2].c_str(),50,0,1);
	grafico3[3]=new TH1F("Un. layers-He",Variables[3].c_str(),10,0,10);
	grafico3[4]=new TH1F("fuori X-He",Variables[4].c_str(),10,0,10);
	grafico3[5]=new TH1F("Track Chi^2-He",Variables[5].c_str(),1000,0,25);
	grafico3[6]=new TH1F("RICH Hits: tot - used-He",Variables[6].c_str(),100,0,30);
        grafico3[7]=new TH1F("Rich Photoelectrons-He",Variables[7].c_str(),300,0,30);
	grafico3[8]=new TH1F("Track Edep: Tot - Track-good-He",Variables[8].c_str(),300,0,30);


	for(int z=0;z<11;z++){
		graficoz[0][z]=new TH1F(("Anticl"+to_string(z)).c_str(),(Variables[0]+"_"+to_string(z)).c_str(),5,0,10);
		graficoz[1][z]=new TH1F(("Un. TOF"+to_string(z)).c_str(),(Variables[1]+"_"+to_string(z)).c_str(),20,0,20);
		graficoz[2][z]=new TH1F(("R Diff."+to_string(z)).c_str(),(Variables[2]+"_"+to_string(z)).c_str(),50,0,1);
		graficoz[3][z]=new TH1F(("Un. layers"+to_string(z)).c_str(),(Variables[3]+"_"+to_string(z)).c_str(),10,0,10);
		graficoz[4][z]=new TH1F(("fuori X"+to_string(z)).c_str(),(Variables[4]+"_"+to_string(z)).c_str(),10,0,10);
		graficoz[5][z]=new TH1F(("Track Chi^2"+to_string(z)).c_str(),(Variables[5]+"_"+to_string(z)).c_str(),50,0,25);
		graficoz[6][z]=new TH1F(("RICH Hits: tot - used"+to_string(z)).c_str(),(Variables[6]+"_"+to_string(z)).c_str(),100,0,30);
		graficoz[7][z]=new TH1F(("Rich Photoelectrons"+to_string(z)).c_str(),(Variables[7]+"_"+to_string(z)).c_str(),300,0,30);	    
		graficoz[8][z]=new TH1F(("Track Edep: Tot - Track"+to_string(z)).c_str(),(Variables[8]+"_"+to_string(z)).c_str(),300,0,30);
	}


	cout<<"**************************** BETA BINS NaF***********************************"<<endl;
	float B=0.4;
        float B1=0;
        float B2=0;
        float E=0.1;
        int binnum=0;
        float Betabins[18]={0.4};
        float Betacent[18]={0};
        float Ekincent[18]={0};
        float a=(log(4.025)-log(0.666))/18;
        float E2=exp(log(0.666)+1.5*a);
        binnum=0;
        while(B1<0.98){
                E=exp(log(0.666)+binnum*a);
                E2=exp(log(0.666)+(binnum+0.5)*a);
                B1=sqrt(1-1/(pow(E+1,2)));
                B2=sqrt(1-1/(pow(E2+1,2)));
                Betabins[binnum]=B1;
                Betacent[binnum]=B2;
                Ekincent[binnum]=1/pow(1-pow(B1,2),0.5)-1;
                binnum++;
        }
        string TitoliNaF[18];
        for(int i=0;i<18;i++){
                ostringstream ss;
                ss<<Betabins[i];
                TitoliNaF[i]= ss.str();
                cout<<TitoliNaF[i]<<", ";
        }


	float Rcut[17]={0.6,0.6,0.67,0.68,0.7,0.72,0.75,0.8,0.85,1,1.1,1.17,1.22,1.4,1.6,1.7,1.8};
	TF1 *RBeta = new TF1("f1","pow((pow(0.938,2)*(pow(x,2)/(1-pow(x,2)))),0.5)",0.1,0.999999999999999999999999999999);
	TF1 *RBeta_P = new TF1("f1","pow((pow(1.875,2)*(pow(x,2)/(1-pow(x,2)))),0.5)",0.1,0.999999999999999999999999999999);

	cout<<"************************************** TRAINING *****************************************************"<<endl;
	int avanzamento=0;
	bool Hecut=false;
	for(int i=0; i<ntupla1->GetEntries()/2;i++) {
		 ntupla1->GetEvent(i);
		if(100*(i/(float)(ntupla1->GetEntries()))>avanzamento) {cout<<avanzamento<<endl;avanzamento++;}
		Massagen= ReturnMass_Gen(MC_type);
		Massa=pow(fabs(pow(fabs(R)*pow((1-pow(Beta,2)),0.5)/Beta,2)),0.5);
		float EdepTOFud=(EdepTOFU+EdepTOFD)/2;
		//Clusterinutili=ClusterTOFtotali-ClusterTOFusati;
		TOF_Up_Down=fabs(EdepTOFD-EdepTOFU);
		Betagen=(pow(pow(Momentogen/Massagen,2)/(1+pow(Momentogen/Massagen,2)),0.5));
		
		if(((int)Cutmask>>11)==512&&Beta<0.97&&Beta>0.8&&(1/Massa)<MassLimit&&R>=0){ 	
			if(Massagen<1&&Massagen>0.5)
			{
				grafico1[0]->Fill(NAnticluster);
				grafico1[1]->Fill(Clusterinutili);
				grafico1[2]->Fill(DiffR);
				grafico1[3]->Fill(layernonusati);
				grafico1[4]->Fill(fuoriX);
				grafico1[5]->Fill(Chisquare);
				grafico1[6]->Fill(Richtotused);
				grafico1[7]->Fill(RichPhEl);
				grafico1[8]->Fill(DiffTrackEdep);
			}
			if(Massagen<2&&Massagen>1.5)
			{
				grafico2[0]->Fill(NAnticluster);
				grafico2[1]->Fill(Clusterinutili);
				grafico2[2]->Fill(DiffR);
				grafico2[3]->Fill(layernonusati);
				grafico2[4]->Fill(fuoriX);
				grafico2[5]->Fill(Chisquare);
				grafico2[6]->Fill(Richtotused);
				grafico2[7]->Fill(RichPhEl);
				grafico2[8]->Fill(DiffTrackEdep);
			}
			if(Massagen<4.5&&Massagen>2.5&&IsCharge1==1)
			{
				grafico3[0]->Fill(NAnticluster);
				grafico3[1]->Fill(Clusterinutili);
				grafico3[2]->Fill(DiffR);
				grafico3[3]->Fill(layernonusati);
				grafico3[4]->Fill(fuoriX);
				grafico3[5]->Fill(Chisquare);
				grafico3[6]->Fill(Richtotused);
				grafico3[7]->Fill(RichPhEl);
				grafico3[8]->Fill(DiffTrackEdep);
			}

		}

	}

	/////////////////////////////////////////////////
	double X[9][1000];
	double Y[9][1000];
	double Y_2[9][1000];
	TSpline3 *Bkgnd[9]; 
	TSpline3 *Signal[9];

	// NORMALIZZAZIONE IN AREA

	for(int m=0;m<9;m++){
		int f=0;
		float area1=0;
		float area2=0;
		float area3=0;
		int entries1=grafico1[m]->GetEntries();
		int entries2=grafico2[m]->GetEntries();
		int entries3=grafico3[m]->GetEntries();
		cout<<"Tot entries 1: "<<entries1<<endl;
		cout<<"Tot entries 2: "<<entries2<<endl;
		cout<<"Tot entries 3: "<<entries3<<endl;
		cout<<"Bins 1: "<<grafico1[m]->GetNbinsX()<<endl;
		cout<<"Bins 2: "<<grafico2[m]->GetNbinsX()<<endl;
		cout<<"Bins 3: "<<grafico3[m]->GetNbinsX()<<endl;
		for(int i=0;i<=grafico1[m]->GetNbinsX()+1;i++){
			grafico1[m]->SetBinContent(i,grafico1[m]->GetBinContent(i)/(double)entries1);
			area1 = area1 + grafico1[m]->GetBinContent(i);
		}


		for(int i=0;i<=grafico2[m]->GetNbinsX()+1;i++){
			grafico2[m]->SetBinContent(i,grafico2[m]->GetBinContent(i)/(double)entries2);
			area2 = area2 + grafico2[m]->GetBinContent(i);
		}
		for(int i=0;i<=grafico3[m]->GetNbinsX()+1;i++){
			grafico3[m]->SetBinContent(i,grafico3[m]->GetBinContent(i)/(double)entries3);
			area3 = area3 + grafico3[m]->GetBinContent(i);
		}

		cout<<"Area1: "<<entries1<<endl;
		cout<<"Area2: "<<entries2<<endl;
		cout<<"Area3: "<<entries3<<endl;

		////////////////////////////////////////
		///////////// DEFINING PDF ////////////////
		int g=0;
		for(int i=0;i<=grafico1[m]->GetNbinsX();i++){
			if(grafico1[m]->GetBinContent(i)>0){
				X[m][g]=grafico1[m]->GetBinLowEdge(i);
				Y[m][g]=grafico1[m]->GetBinContent(i);
				Y_2[m][g]=grafico2[m]->GetBinContent(i);
				g++;
			}
		}
		cout<<endl;
		Bkgnd[m]= new TSpline3("Cubic Spline",X[m],Y[m],g);
		Signal[m]=new TSpline3("Cubic Spline",X[m],Y_2[m],g);
		///////////////////////////////////////////

		cout<<m<<endl;

		c1[m] = new TCanvas(Variables[m].c_str());
		c1[m]->Divide(2,1);
		c1[m]->cd(1);
		gPad->SetGridx();
		gPad->SetGridy();
		grafico1[m]->SetBarOffset(0.25);
		/*grafico3[m]->SetBarOffset(1.5);*/
		grafico1[m]->SetBarWidth(0.5);
                grafico2[m]->SetBarWidth(0.5);
		grafico2[m]->GetXaxis()->SetTitle(Variables[m].c_str());
		grafico2[m]->Draw("B");
		grafico1[m]->Draw("B,SAME");
		//grafico3[m]->Draw("SAME");
		grafico1[m]->SetFillColor(2);
		grafico2[m]->SetFillColor(4);
		grafico3[m]->SetFillColor(3);
		grafico3[m]->SetFillStyle(3001);
		Bkgnd[m]->SetLineWidth(2);
		Bkgnd[m]->SetLineColor(2);
		//Bkgnd[m]->Draw("SAME");
		Signal[m]->SetLineWidth(2);
		Signal[m]->SetLineColor(4);
		//Signal[m]->Draw("SAME");
	}
	cout<<"************************************** TEST *****************************************************"<<endl;
	/////////////LKLHD CALCULATION//////////
	TCanvas *c2=new TCanvas("Likelihood response");
	cout<<"Likelihood calculation..."<<endl;
	TH1F *Sl =new TH1F("Signal","Signal",50,0,7);
	TH1F *Bd =new TH1F("Bkgnd","Bkgnd",50,0,7);
	TH1F *He =new TH1F("Bkgnd","Bkgnd",50,0,7);
	TH1F *LKL_Data=new TH1F("Lk_Data","Signal",50,0,7);
	TH1F *LKL_Dataz[11];
	for(int z=0;z<11;z++) LKL_Dataz[z]=new TH1F(("Lk_Data_"+to_string(z)).c_str(),("Lk_Data_"+to_string(z)).c_str(),50,0,7);
	TH2F *LkvsDist_P=new TH2F("LkvsDist_P","",500,0,1,500,-1,1);
	TH2F *LkvsDist_D=new TH2F("LkvsDist_D","",500,0,1,500,-1,1);
	
	float var[9]={0,0,0,0,0,0,0,0,0};
	BDTreader();
	avanzamento=0;
	int qu=0;
	for(int l=0; l<ntupla1->GetEntries()/3;l++) {
		 ntupla1->GetEvent(l);
		if(100*(l/(float)(ntupla1->GetEntries()))>avanzamento) {cout<<avanzamento<<endl;avanzamento++;}
		Massagen= ReturnMass_Gen(MC_type);
		Massa=pow(fabs(pow(fabs(R)*pow((1-pow(Beta,2)),0.5)/Beta,2)),0.5);
		Betagen=(pow(pow(Momentogen/Massagen,2)/(1+pow(Momentogen/Massagen,2)),0.5));
		//Clusterinutili=fabs(ClusterTOFtotali-ClusterTOFusati);
		
		TOF_Up_Down=fabs(EdepTOFD-EdepTOFU);
		var[0]=NAnticluster;
		var[1]=Clusterinutili;
		var[2]=DiffR;
		var[3]=layernonusati;
		var[4]=fuoriX;
		var[5]=Chisquare;
		var[6]=Richtotused;
		var[7]=RichPhEl;
		var[8]=DiffTrackEdep;
		double L_true=1;
		double L_false=1;
		for(int m=0;m<8;m++)
		{
			if(m!=10){
			L_false=L_false*Bkgnd[m]->Eval(var[m]);
			L_true=L_true*Signal[m]->Eval(var[m]);
			}
		}

		float L_Discr=-log(1-(L_true/(L_false+L_true)));
		if(((int)Cutmask>>11)==512&&Beta<0.97&&Beta>0.8&&(1/Massa)<MassLimit&&R>=0&&IsCharge1==1)
		{
			BDTresponse=reader->EvaluateMVA("BDTmethod");
			if(Massagen<1&&Massagen>0.5){	    
				qu++;
				Bd->Fill(L_Discr);
				LkvsDist_P->Fill(L_Discr,DistD);
			}
			if(Massagen<2&&Massagen>1.5){
				Sl->Fill(L_Discr);
				LkvsDist_D->Fill(L_Discr,DistD);
			}
			if(Massagen<4&&Massagen>2.5){
				He->Fill(L_Discr);
			}

		}
		
	}


	avanzamento=0;
	for(int l=0; l<ntupla2->GetEntries()/3;l++) {
		 ntupla2->GetEvent(l);
		if(100*(l/(float)(ntupla2->GetEntries()))>avanzamento) {cout<<avanzamento<<endl;avanzamento++;}
		Massa=pow(fabs(pow(fabs(R)*pow((1-pow(Beta,2)),0.5)/Beta,2)),0.5);
		//Clusterinutili=fabs(ClusterTOFtotali-ClusterTOFusati);
		
		TOF_Up_Down=fabs(EdepTOFD-EdepTOFU);
		Betagen=(pow(pow(Momentogen/Massagen,2)/(1+pow(Momentogen/Massagen,2)),0.5));
		zone = 0;

		for(int i=0; i<11; i++) {
                        double geo= geomag[i]  ;
                        double geo2=geomag[i+1];
                        if(Latitude>geo && Latitude<=geo2)
                                zone=i;
                }
	

	
		if(((int)Cutmask>>11)==512&&IsCharge1==1&&R>1.2*Rcutoff&&R>30){
                        graficoz[0][zone]->Fill(NAnticluster);
                        graficoz[1][zone]->Fill(Clusterinutili);
                        graficoz[2][zone]->Fill(DiffR);
                        graficoz[3][zone]->Fill(layernonusati);
                        graficoz[4][zone]->Fill(fuoriX);
                        graficoz[5][zone]->Fill(Chisquare);
                        graficoz[6][zone]->Fill(Richtotused);
                        graficoz[7][zone]->Fill(RichPhEl);
                        graficoz[8][zone]->Fill(DiffTrackEdep);
                }

		var[0]=NAnticluster;
		var[1]=Clusterinutili;
		var[2]=DiffR;
		var[3]=layernonusati;
		var[4]=fuoriX;
		var[5]=Chisquare;
		var[6]=Richtotused;
		var[7]=RichPhEl;
		var[8]=DiffTrackEdep;
			

		double L_true=1;
		double L_false=1;
		for(int m=0;m<8;m++)
		{
			if(m!=10){
			L_false=L_false*Bkgnd[m]->Eval(var[m]);
			L_true=L_true*Signal[m]->Eval(var[m]);
			}
		}
		float L_Discr=-log(1-(L_true/(L_false+L_true)));
		if(((int)Cutmask>>11)==512&&IsCharge1==1&&R>1.2*Rcutoff&&R>30){
			LKL_Dataz[zone]->Fill(L_Discr);
		}
		if(((int)Cutmask>>11)==512&&Beta<0.97&&Beta>0.8&&(1/Massa)<MassLimit&&R>=0&&IsCharge1==1)
		{
			LKL_Data->Fill(L_Discr);
		}
	}
	
	//////////// LAT DEPENDENCE ///////////////
	for(int m=0;m<9;m++){
		for(int z=1;z<11;z++){
			graficoz[m][z]->Scale(1/graficoz[m][z]->GetEntries());
		}
	}

	for(int m=0;m<9;m++){
		c1[m]->cd(2);
		gPad->SetGridx();
		gPad->SetGridy();
		graficoz[m][1]->SetBarOffset(0.25);
		/*grafico3[m]->SetBarOffset(1.5);*/
		graficoz[m][1]->SetBarWidth(0.5);
		graficoz[m][1]->SetBarWidth(0.5);
		graficoz[m][1]->GetXaxis()->SetTitle(Variables[m].c_str());
		graficoz[m][1]->SetLineColor(1);
		graficoz[m][1]->Draw("L");
		for(int z=2;z<11;z++){
			graficoz[m][z]->SetBarWidth(0.5);
			graficoz[m][z]->SetBarWidth(0.5);
			graficoz[m][z]->GetXaxis()->SetTitle(Variables[m].c_str());
			graficoz[m][z]->SetLineColor(z-1);
			graficoz[m][z]->Draw("L,SAME");
		}
	}


	/////////////// DISTR. TOT E BAD //////////////////////
	f->cd();
	gPad->SetLogz();
	RvsRgen_bad->GetXaxis()->SetRangeUser(0.55,4.5);
	RvsRgen_bad->GetYaxis()->SetRangeUser(0,4.5);
	RvsRgen_bad->GetXaxis()->SetTitle("R [GV]");
	RvsRgen_bad->GetYaxis()->SetTitle("R gen [GV]");
	RvsRgen_bad->Draw();
	RvsRgen_sigma->SetLineColor(2);
	RvsRgen_sigma->Draw("same,CONT2");
	d_tris->cd();
	gPad->SetLogz();
        sigmagen_bad->GetXaxis()->SetRangeUser(0.55,4.5);
        sigmagen_bad->GetYaxis()->SetRangeUser(0,4.5);
        sigmagen_bad->GetXaxis()->SetTitle("(R-Rgen)/sigma");
        sigmagen_bad->GetYaxis()->SetTitle("(Beta-Betagen)/sigma");
        sigmagen_bad->Draw();

	string Titolo;


	/////////////// LK VS DIST ////////////////////////////
	TCanvas *c7_bis=new TCanvas("Likelihood vs Distance");
	c7_bis->cd();
	gPad->SetLogy();
	LkvsDist_P->SetMarkerStyle(7);
	LkvsDist_P->SetMarkerColor(2);
	LkvsDist_D->Draw("col");
	LkvsDist_P->Draw("same");


	/////////////// LAT LK ////////////////////////////////////	
	TCanvas *c7_tris=new TCanvas("Likelihood LAT dep.");
	c7_tris->cd();
	gPad->SetLogy();
	LKL_Dataz[1]->GetXaxis()->SetTitle("Likelihood Distribution");
	LKL_Dataz[1]->SetTitle("Likelihood Distribution");
	LKL_Dataz[1]->SetLineColor(1);
	LKL_Dataz[1]->Scale(1/LKL_Dataz[1]->GetEntries());
	LKL_Dataz[1]->Draw("L");
	for(int z=2;z<11;z++){
		LKL_Dataz[z]->SetLineColor(z-1);
		LKL_Dataz[z]->Scale(1/LKL_Dataz[z]->GetEntries());
		LKL_Dataz[z]->Draw("L,SAME");
	}



	/////////////// TEMPLATE FIT //////////////////////////
	c2->cd();
	TObjArray *mc = new TObjArray(2);
	mc->Add(Sl);
	mc->Add(Bd);
	TFractionFitter* fit = new TFractionFitter(LKL_Data, mc);
	int status = fit->Fit();
	std::cout << "fit status: " << status << std::endl;
	float pesoP=0;
	float pesoD=0;
	if (status == 0) {                       // check on fit status
		double w1,e1,w2,e2,w3,e3=1;
		fit->GetResult(0,w1,e1);
		fit->GetResult(1,w2,e2);
		//fit->GetResult(2,w3,e3);
		TH1F * result = (TH1F*) fit->GetPlot();
		float itot= result->Integral();
		float i1 = Sl->Integral();
		float i2 = Bd->Integral();
		float i3 = He->Integral();
		cout<<itot<<endl;
		cout<<w1<<" "<<w2<<" "<<w3<<endl;
		pesoP=w2/i2*itot;
		pesoD=w1/i1*itot;
		for(int i=0; i<Sl->GetNbinsX();i++)
			Sl->SetBinContent(i,w1*Sl->GetBinContent(i)/i1*itot);
		for(int i=0; i<Bd->GetNbinsX();i++)
			Bd->SetBinContent(i,w2*Bd->GetBinContent(i)/i2*itot);
		//for(int i=0; i<He->GetNbinsX();i++)
		//	He->SetBinContent(i,w3*He->GetBinContent(i)/i3*itot);
		Sl->SetFillColor(4);
		He->SetFillColor(3);
		Bd->SetFillColor(2);
		//Bd->Draw();
		//Sl->Draw("same");
		//He->Draw("same");
		THStack *hs= new THStack("hs","Likelihood response (stacked)");
		Sl->SetFillStyle(3001);
        	Bd->SetFillStyle(3001);
		hs->Add(Sl);
		hs->Add(Bd);
		//hs->Add(He);        
		LKL_Data->GetXaxis()->SetTitle("Likelihood response");
		LKL_Data->SetMarkerStyle(8);
		//LKL_Data->Draw("ep");
		//hs->Draw("same");
		result->SetLineColor(5);
		//result->Draw("SAME");
	}
	if (status != 0) {
        Sl->SetFillColor(4);
        He->SetFillColor(3);
        Bd->SetFillColor(2);
        Sl->SetFillStyle(3001);
        He->SetFillStyle(3001);
        Bd->SetFillStyle(3001);
        //Bd->Draw();
        //Sl->Draw("same");
        //He->Draw("same");
        LKL_Data->SetMarkerStyle(8);
        //LKL_Data->Draw("epsame");
        }
	Bd->Scale(1/Bd->Integral());
	Sl->Scale(1/Sl->Integral());
	Bd->Draw();
	Sl->Draw("same");
	//////////////////// CUT OPTIMIZATION /////////

	TGraph *RPvsEff = new TGraph();
	RPvsEff->SetTitle("Likelihood (new Variables)");
	TGraph *EffvsCut = new TGraph();
	EffvsCut->SetTitle("Likelihood (new Variables)");
	TGraph *Ott = new TGraph();
	Ott->SetTitle("Likelihood (new Variables)");
	float efficienza1=0;
	float efficienza2=0;
	for(int j=0;j<=Sl->GetNbinsX()+1;j++){
		double taglio=Sl->GetBinCenter(j);
		efficienza1=0;
		efficienza2=0;
		efficienza1=0;
		for(int i=Sl->GetNbinsX(); i>=0;i--)
			if(Sl->GetBinCenter(i)>=taglio) {efficienza1=efficienza1+Sl->GetBinContent(i);}
		efficienza1=efficienza1/Sl->Integral();
		for(int i=Bd->GetNbinsX(); i>=0;i--)
			if(Bd->GetBinCenter(i)>=taglio) {efficienza2=efficienza2+Bd->GetBinContent(i);}
		efficienza2=1-efficienza2/Bd->Integral();
		RPvsEff->SetPoint(j,efficienza1,efficienza2);
		EffvsCut->SetPoint(j,taglio,efficienza1);
		S=efficienza1*Sl->Integral(); B=182.601*(1-efficienza2)*Bd->Integral();
		if(S+B>0) Ott->SetPoint(j,taglio,S/pow(B,0.5));
	}
	RPvsEff->GetXaxis()->SetTitle("Efficiency");
	RPvsEff->GetYaxis()->SetTitle("Background rejection");
	RPvsEff->SetLineColor(2);
	RPvsEff->SetLineWidth(2);
	EffvsCut->GetXaxis()->SetTitle("Cut");
	EffvsCut->GetYaxis()->SetTitle("Efficiency");
	EffvsCut->SetLineColor(2);
	EffvsCut->SetLineWidth(2);
	Ott->GetXaxis()->SetTitle("Cut");
	Ott->GetYaxis()->SetTitle("Significance (S/sqrt(S+B))");
	Ott->SetLineColor(2);
	Ott->SetLineWidth(2);
	RPvsEff->GetXaxis()->SetRangeUser(-1,1);
	RPvsEff->GetYaxis()->SetRangeUser(0.6,1);
	EffvsCut->GetXaxis()->SetRangeUser(-1,1);
	Ott->GetXaxis()->SetRangeUser(-1,1);
	Ott->GetYaxis()->SetRangeUser(0.1,200);
	TCanvas *c3=new TCanvas("R.P. vs Eff");
	c3->cd();
	c3->SetGridx();
	c3->SetGridy();
	RPvsEff->Draw("AL");
	TCanvas *c5=new TCanvas("Eff vs Cut");
	c5->cd();
	c5->SetGridx();
	c5->SetGridy();
	EffvsCut->Draw("AL");
	TCanvas *c6=new TCanvas("Cut Optimization");
	c6->cd();
	c6->SetGridx();
	c6->SetGridy();
	c6->SetLogy();
	Ott->Draw("AL");
	//////////////////////////////////////////////

	////////////////////////////////
	/////////OUTPUT////////////
	string nomefile="./QualityVariables_NaF.root";
	TFile *f_out=new TFile(nomefile.c_str(), "RECREATE");
	f_out->mkdir("Bad Events Study");
	f_out->mkdir("Variables");
	f_out->mkdir("Distributions");
	f_out->mkdir("Performances");
	f_out->mkdir("Splines");
	f_out->cd("RvsBeta");
	f_out->cd("Variables");
	for(int m=0;m<9;m++) {c1[m]->Write();}
	f_out->cd("Bad Events Study");
	f->Write();
        d_tris->Write();
	f_out->cd("Distributions");
	c2->Write();
	c7_tris->Write();
	f_out->cd("Performances");
	c3->Write();
	c5->Write();
	c6->Write();
	c7_bis->Write();
	f_out->cd("Splines");
	string nome;
	string Spline="Spline NaF: ";
	string signal="_SGNL";
	string bkgnd="_BKGND";
	for(int m=0;m<9;m++) {
		nome=Spline+Variables[m]+bkgnd;	 
		Bkgnd[m]->Write(nome.c_str());
		nome=Spline+Variables[m]+signal;
		Signal[m]->Write(nome.c_str());
	}
	f_out->Write();
	f_out->Close();
	////////////////////////////
	return 1;
}


void BDTreader()
{
        TMVA::Tools::Instance();
        reader = new TMVA::Reader( "V:Color:!Silent" );
        reader->AddVariable("NAnticluster", &NAnticluster);
        reader->AddVariable("Chisquare",&Chisquare);
        reader->AddVariable("layernonusati",&layernonusati);
        reader->AddVariable("NTofUsed := NTofClusters - NTofClustersusati",&Clusterinutili);
        reader->AddVariable("diffR := TMath::Abs(Rup-Rdown)/R",&DiffR);
        reader->AddVariable("TOF_Up_Down := TMath::Abs(Endep[2]+Endep[3]-Endep[0]-Endep[1])", &TOF_Up_Down);

        reader->BookMVA("BDTmethod", "../TMVA/QualityBDT_BDT.weights.xml");
}

