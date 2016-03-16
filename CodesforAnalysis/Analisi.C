#include "TH2.h"
#include "TH3.h"
#include <TVector3.h>
#include <fstream>
#include <sstream>
#include <vector>
#include "TCanvas.h"
#include "TLegend.h"
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
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include <cstring>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <math.h>

#include "Functions_auto.h"
#include "Parte2/Definitions.h"
#include "Parte2/FitError.h"
#include "Parte2/MCpreeff.h"
#include "Parte2/MCUnbiaseff.h"
#include "Parte2/Hecut.h"
#include "Parte2/SlidesforPlot.h"
#include "Parte2/MCQualeff.h"
#include "Parte2/MCTrackeff.h"
#include "Parte2/MCpreSeleff.h"
#include "Parte2/MCpreCheck.h"
#include "Parte2/MigrationMatrix.h"
#include "Parte2/DATAUnbiaseff.h"
#include "Parte2/CorrelazionePreselezioni.h"
#include "Parte2/DATApreSeleff.h"
#include "Parte2/DATAQualeff.h"
#include "Parte2/DVSMCpreSeleff.h"
#include "Parte2/DVSMCQualeff2.h"
#include "Parte2/DVSMCTrackeff.h"
#include "Parte2/AcceptanceP.h"
#include "Parte2/AcceptanceD.h"
#include "Parte2/CorrLAT.h"
#include "Parte2/ProtonFlux.h"
#include "Parte2/Deutons.h"
#include "Parte2/DeutonsDist.h"
#include "Parte2/MCMC.h"
#include "Parte2/DeutonsFlux.h"
#include "Parte2/DeutonsFlux_Dist.h"
#include "Parte2/Cuts.h"
#include "FillIstogram.h"
using namespace std;


int main(int argc, char * argv[])
{
	cout<<"Month _ Indx _ Frac "<<endl;
	//string percorso="/home/francesco/PhD/LocalCNAF/";
	string percorso="/storage/gpfs_ams/ams/users/fdimicco/Deutons";
	string frac=argv[3];
	FRAC=atoi(argv[3]);
	INDX=atoi(argv[2]);
	string mese=argv[1];
	cout<<"********************************************** R BINS ******************************************************************************"<<endl;
	for(int i=0;i<44;i++)
	{
		float temp=i+14;
		bin[i]=0.1*pow(10,temp/(9.5*2));
		if(i<43) {R_cent[i]=0.1*pow(10,(temp+0.5)/(9.5*2));
			encindeut[i]=pow(((1+pow((R_cent[i]/1.875),2))),0.5)-1;
			encinprot[i]=pow(((1+pow((R_cent[i]/0.938),2))),0.5)-1;
		}
		cout<<bin[i]<<endl;
	}
	for(int i=0;i<43;i++) {
		deltaencinprot[i]=(pow(((1+pow((bin[i+1]/0.938),2))),0.5)-1)-(pow(((1+pow((bin[i]/0.938),2))),0.5)-1);
		deltaencindeut[i]=(pow(((1+pow((bin[i+1]/1.875),2))),0.5)-1)-(pow(((1+pow((bin[i]/1.875),2))),0.5)-1);
	}

	cout<<"**************************** BETA BINS TOF***********************************"<<endl;
	float B=0.4;
	float B1=0;
	float B2=0;
	float E=0.1;
	int binnum=0;
	float a=(log(1)-log(0.1))/18;
	float E2=exp(log(0.1)+1.5*a);
	while(B1<0.85){
		E=exp(log(0.1)+binnum*a);
		E2=exp(log(0.1)+(binnum+0.5)*a);
		B1=sqrt(1-1/(pow(E+1,2)));
		B2=sqrt(1-1/(pow(E2+1,2)));
		Betabins[binnum]=B1;
		BetabinsR_P[binnum]=0.938*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		BetabinsR_D[binnum]=1.875*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		Betacent[binnum]=B2;
		Ekincent[binnum]=1/pow(1-pow(B2,2),0.5)-1;
		binnum++;
	}
	for(int i=0;i<18;i++) deltaencinTOF[i]=(1/pow(1-pow(Betabins[i+1],2),0.5)-1)-(1/pow(1-pow(Betabins[i],2),0.5)-1);
	for(int i=0;i<18;i++) cout<<Betabins[i]<<" "<<BetabinsR_D[i]<<" "<<BetabinsR_P[i]<<endl;
	string TitoliTOF[18];
	for(int i=0;i<18;i++){
		ostringstream ss;
		ss<<Betabins[i];
		TitoliTOF[i]= ss.str();
		cout<<TitoliTOF[i]<<", ";
	}
	cout<<endl;

	cout<<"**************************** BETA BINS NaF***********************************"<<endl;
	a=(log(4.025)-log(0.666))/18;
	E2=exp(log(0.666)+1.5*a);
	binnum=0;
	while(B1<0.98){
		E=exp(log(0.666)+binnum*a);
		E2=exp(log(0.666)+(binnum+0.5)*a);
		B1=sqrt(1-1/(pow(E+1,2)));
		B2=sqrt(1-1/(pow(E2+1,2)));
		BetabinsNaF[binnum]=B1;
		BetabinsNaFR_P[binnum]=0.938*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		BetabinsNaFR_D[binnum]=1.875*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		BetacentNaF[binnum]=B2;
		EkincentNaF[binnum]=1/pow(1-pow(B1,2),0.5)-1;
		binnum++;
	}
	for(int i=0;i<18;i++) deltaencinNaF[i]=(1/pow(1-pow(BetabinsNaF[i+1],2),0.5)-1)-(1/pow(1-pow(BetabinsNaF[i],2),0.5)-1);
	for(int i=0;i<18;i++) cout<<BetabinsNaF[i]<<" "<<BetabinsNaFR_D[i]<<" "<<EkincentNaF[i]<<" "<<deltaencinNaF[i]<<endl;
	string TitoliNaF[18];
	for(int i=0;i<18;i++){
		ostringstream ss;
		ss<<BetabinsNaF[i];
		TitoliNaF[i]= ss.str();
		cout<<TitoliNaF[i]<<", ";
	}
	cout<<endl;
	cout<<"**************************** BETA BINS Agl***********************************"<<endl;
	a=(log(9.01)-log(2.57))/18;
	E2=exp(log(2.57)+1.5*a);
	binnum=0;
	while(B1<0.995){
		E=exp(log(2.57)+binnum*a);
		E2=exp(log(2.57)+(binnum+0.5)*a);
		B1=sqrt(1-1/(pow(E+1,2)));
		B2=sqrt(1-1/(pow(E2+1,2)));
		BetabinsAgl[binnum]=B1;
		BetabinsAglR_P[binnum]=0.938*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		BetabinsAglR_D[binnum]=1.875*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		BetacentAgl[binnum]=B2;
		EkincentAgl[binnum]=1/pow(1-pow(B2,2),0.5)-1;
		binnum++;
	}
	for(int i=0;i<18;i++) deltaencinAgl[i]=(1/pow(1-pow(BetabinsAgl[i+1],2),0.5)-1)-(1/pow(1-pow(BetabinsAgl[i],2),0.5)-1);
	for(int i=0;i<18;i++) cout<<BetabinsAgl[i]<<" "<<BetabinsAglR_D[i]<<" "<<EkincentAgl[i]<<" "<<deltaencinAgl[i]<<endl;
	string TitoliAgl[18];
	for(int i=0;i<18;i++){
		ostringstream ss;
		ss<<BetabinsAgl[i];
		TitoliAgl[i]= ss.str();
		cout<<TitoliAgl[i]<<", ";
	}
	cout<<endl;

	////////// BINNAGGIO IN BETA
	/*for(int i=0;i<19;i++){
	  BetaD[i]=Betabins[i];
	  BetaP[i]=Betabins[i];
	  BetaNaFD[i]=BetabinsNaF[i];
	  BetaNaFP[i]=BetabinsNaF[i];
	  BetaAglD[i]=BetabinsAgl[i];
	  BetaAglP[i]=BetabinsAgl[i];
	  }*/	
	////////////////////////////
	/////////// BINNAGGIO IN RIGIDITA'
	for(int i=0;i<20;i++){
		BetaD[i]=BetabinsR_D[i];
		BetaP[i]=BetabinsR_P[i];
		BetaNaFD[i]=BetabinsNaFR_D[i];
		BetaNaFP[i]=BetabinsNaFR_P[i];
		BetaAglD[i]=BetabinsAglR_D[i];
		BetaAglP[i]=BetabinsAglR_P[i];
		cout<<"TOF: "<<BetaD[i]<<" "<<BetaP[i]<<" NaF "<<BetaNaFD[i]<<" "<<BetaNaFP[i]<<" Agl "<<BetaAglD[i]<<" "<<BetaAglP[i]<<endl;
	}
	////////////////////////////

	cout<<"************************************* ISTOGRAM FILLING **************************************************************************"<<endl;
	FillIstogram(INDX,frac,mese);
	string nomefile=percorso + "/Risultati/risultati/"+mese+"_"+frac+"_P1.root";
	TFile *file1 =TFile::Open(nomefile.c_str());
	if(!file1){
		nomefile=percorso + "/Risultati/"+mese+"/"+mese+"_"+frac+"_P1.root";
		file1 =TFile::Open(nomefile.c_str());
	}

	cout<<"********************************** ANALYSIS *************************************************************************************"<<endl;
	MCpreeff(file1);
	MCUnbiaseff(file1);
	Hecut(file1);
	SlidesforPlot(file1);
	MCQualeff(file1);
	MCTrackeff(file1);
	MCpreSeleff(file1);
	Correlazione_Preselezioni(file1);
	MCpreCheck();
	MigrationMatrix();
	DATAUnbiaseff(file1);
	DATApreSeleff(file1);
	DVSMCTrackeff(file1);
	DATAQualeff(file1);
	DVSMCpreSeleff(file1);
	DVSMCQualeff2(file1);
	AcceptanceP(file1);
	AcceptanceD(file1);
	CorrLAT(file1);
	ProtonFlux(file1);
	if(frac=="tot")DeutonsTemplFits(file1);
	if(frac=="tot")DeutonsTemplFits_Dist(file1);
	DeutonFlux(file1);
	DeutonFlux_Dist(file1);
	cout<<"************************************** OUTPUT **************************************************************"<<endl;
	nomefile=percorso + "/CodesforAnalysis/Final_plots/"+mese+".root";
	TFile *f_out=new TFile(nomefile.c_str(), "RECREATE");
	string nome;
	f_out->mkdir("Slides common plots");
	f_out->mkdir("MC Results");
	f_out->mkdir("MC Results/He cut");
	f_out->mkdir("MC Results/Preselections");
	f_out->mkdir("MC Results/Quality");
	f_out->mkdir("Data-driven Results");
	f_out->mkdir("Data-driven Results/Latitude Effect");
	f_out->mkdir("Data-driven Results/Data vs MC");
	f_out->mkdir("Data-driven Results/Data vs MC/Preselections");
	f_out->mkdir("Data-driven Results/Data vs MC/Quality");
	f_out->mkdir("Acceptance");
	f_out->mkdir("P Flux");
	f_out->mkdir("Mass Template Fits");
	f_out->mkdir("Dist Template Fits");
	f_out->mkdir("D Flux (Mass fit)");
	f_out->mkdir("D Flux (Dist. fit)");
	for(int l=0;l<11;l++) {
		nome="Mass Template Fits/Mass TOF Template Fits/Geo.Zones/Zona "+numero[l];
		f_out->mkdir(nome.c_str());
	}
	f_out->mkdir("Mass Template Fits/Mass TOF Template Fits/Primaries/D.bins");
	f_out->mkdir("Mass Template Fits/Mass TOF Template Fits/Primaries/P.bins");
	f_out->mkdir("Mass Template Fits/Mass TOF Template Fits/Geo.Zones");
	for(int l=0;l<11;l++) {
		nome="Mass Template Fits/Mass NaF Template Fits/Geo.Zones/Zona "+numero[l];
		f_out->mkdir(nome.c_str());
	}
	f_out->mkdir("Mass Template Fits/Mass NaF Template Fits/Primaries/D.bins");
	f_out->mkdir("Mass Template Fits/Mass NaF Template Fits/Primaries/P.bins");
	f_out->mkdir("Mass Template Fits/Mass NaF Template Fits/Geo.Zones");

	for(int l=0;l<11;l++) {
		nome="Mass Template Fits/Mass Agl Template Fits/Geo.Zones/Zona "+numero[l];
		f_out->mkdir(nome.c_str());
	}
	f_out->mkdir("Mass Template Fits/Mass Agl Template Fits/Primaries/D.bins");
	f_out->mkdir("Mass Template Fits/Mass Agl Template Fits/Primaries/P.bins");
	f_out->mkdir("Mass Template Fits/Mass Agl Template Fits/Geo.Zones");

	for(int l=0;l<11;l++) {
                nome="Dist Template Fits/Dist TOF Template Fits/Geo.Zones/Zona "+numero[l];
                f_out->mkdir(nome.c_str());
        }
        f_out->mkdir("Dist Template Fits/Dist TOF Template Fits/Primaries/D.bins");
	f_out->mkdir("Dist Template Fits/Dist TOF Template Fits/Primaries/P.bins");
        f_out->mkdir("Dist Template Fits/Dist TOF Template Fits/Geo.Zones");

	for(int l=0;l<11;l++) {
                nome="Dist Template Fits/Dist NaF Template Fits/Geo.Zones/Zona "+numero[l];
                f_out->mkdir(nome.c_str());
        }
        f_out->mkdir("Dist Template Fits/Dist NaF Template Fits/Primaries/D.bins");
	f_out->mkdir("Dist Template Fits/Dist NaF Template Fits/Primaries/P.bins");
        f_out->mkdir("Dist Template Fits/Dist NaF Template Fits/Geo.Zones");
	for(int l=0;l<11;l++) {
                nome="Dist Template Fits/Dist Agl Template Fits/Geo.Zones/Zona "+numero[l];
                f_out->mkdir(nome.c_str());
        }
        f_out->mkdir("Dist Template Fits/Dist Agl Template Fits/Primaries/D.bins");
	f_out->mkdir("Dist Template Fits/Dist Agl Template Fits/Primaries/P.bins");	
	f_out->mkdir("Dist Template Fits/Dist Agl Template Fits/Geo.Zones");

	f_out->mkdir("Export");
	f_out->mkdir("Export/Hecut");
	f_out->mkdir("Export/MCMC");
	f_out->mkdir("Export/Qual.Sel.Eff.");
	f_out->mkdir("Export/Acceptance");
	f_out->mkdir("Export/Slidesplot");
	
	f_out->cd("Slides common plots");
	p1->Write();
	p2->Write();
	p3->Write();
	p4->Write();
	p5->Write();
	p6->Write();
	p7->Write();
	p8->Write();
	p9->Write();
	p10->Write();
        p11->Write();
        p12->Write();
	p13->Write();
        p14->Write();
        p15->Write();
	p10Q->Write();
        p11Q->Write();
        p12Q->Write();
	p13Q->Write();
        p14Q->Write();
        p15Q->Write();
	p16->Write();
	p17->Write();
	p18->Write();
	p19->Write();
	p20->Write();
	p21->Write();
	p22->Write();
	p23->Write();
	p24->Write();
	f_out->cd("MC Results");
	f_out->cd("MC Results/He cut");
        c36->Write();
        c36_bis->Write();
	c37->Write();
	f_out->cd("MC Results/Preselections");
	c4->Write();
	c4_bis->Write();
	c11->Write();
	c_7->Write();
	c7->Write();
	c8->Write();
	c10->Write();
	c13->Write();
	for(int S=0;S<3;S++) c9[S]->Write();
	c27->Write();
	f_out->cd("MC Results/Quality");
	c5->Write();
	c5_bis->Write();
	c6->Write();
	c6_bis->Write();
	f_out->cd("Data-driven Results");
	c12->Write();
	c28->Write();
	f_out->cd("Data-driven Results/Latitude Effect");
	for(int S=0;S<3;S++) c14[S]->Write();
	c15->Write();
	c16->Write();
	c26->Write();
	c26_bis->Write();
	f_out->cd("Data-driven Results/Data vs MC");
	f_out->cd("Data-driven Results/Data vs MC/Preselections");
	for(int S=0;S<3;S++) c17[S]->Write();
	f_out->cd("Data-driven Results/Data vs MC/Quality");
	c20->Write();
	c21->Write();
	f_out->cd("Acceptance");
	c22->Write();
	c31->Write();
	c31_tris->Write();
	c31_bis->Write();
	f_out->cd("P Flux");
	c23->Write();
	c24->Write();
	c25->Write();
	if(frac=="tot"){
		for(int l=1;l<11;l++) {
			nome="Mass Template Fits/Mass TOF Template Fits/Geo.Zones/Zona "+numero[l];
			f_out->cd(nome.c_str());
			for(int m=0;m<18;m++) c30[l][m]->Write();
		}
		f_out->cd("Mass Template Fits/Mass TOF Template Fits/Primaries/D.bins");
		for(int m=0;m<18;m++) c30[11][m]->Write();
		f_out->cd("Mass Template Fits/Mass TOF Template Fits/Primaries/P.bins");
		for(int m=0;m<18;m++) c30[0][m]->Write();
		for(int l=1;l<11;l++) {
			nome="Mass Template Fits/Mass NaF Template Fits/Geo.Zones/Zona "+numero[l];
			f_out->cd(nome.c_str());
			for(int m=0;m<18;m++) c30_bis[l][m]->Write();
		}
		f_out->cd("Mass Template Fits/Mass NaF Template Fits/Primaries/D.bins");
		for(int m=0;m<18;m++) c30_bis[11][m]->Write();
		f_out->cd("Mass Template Fits/Mass NaF Template Fits/Primaries/P.bins");
		for(int m=0;m<18;m++) c30_bis[0][m]->Write();
		for(int l=1;l<11;l++) {
			nome="Mass Template Fits/Mass Agl Template Fits/Geo.Zones/Zona "+numero[l];
			f_out->cd(nome.c_str());
			for(int m=0;m<18;m++) c30_tris[l][m]->Write();
		}
		f_out->cd("Mass Template Fits/Mass Agl Template Fits/Primaries/D.bins");
		for(int m=0;m<18;m++) c30_tris[11][m]->Write();
		f_out->cd("Mass Template Fits/Mass Agl Template Fits/Primaries/P.bins");
		for(int m=0;m<18;m++) c30_tris[0][m]->Write();
		for(int l=1;l<11;l++) {
			nome="Dist Template Fits/Dist TOF Template Fits/Geo.Zones/Zona "+numero[l];
			f_out->cd(nome.c_str());
			for(int m=0;m<18;m++) c40[l][m]->Write();
		}
		f_out->cd("Dist Template Fits/Dist TOF Template Fits/Primaries/D.bins");
		for(int m=0;m<18;m++) c40[11][m]->Write();
		f_out->cd("Dist Template Fits/Dist TOF Template Fits/Primaries/P.bins");
		for(int m=0;m<18;m++) c40[0][m]->Write();
		for(int l=1;l<11;l++) {
			nome="Dist Template Fits/Dist NaF Template Fits/Geo.Zones/Zona "+numero[l];
			f_out->cd(nome.c_str());
			for(int m=0;m<18;m++) c40_bis[l][m]->Write();
		}
		f_out->cd("Dist Template Fits/Dist NaF Template Fits/Primaries/D.bins");
		for(int m=0;m<18;m++) c40_bis[11][m]->Write();
		f_out->cd("Dist Template Fits/Dist NaF Template Fits/Primaries/P.bins");
		for(int m=0;m<18;m++) c40_bis[0][m]->Write();
		for(int l=1;l<11;l++) {
			nome="Dist Template Fits/Dist Agl Template Fits/Geo.Zones/Zona "+numero[l];
			f_out->cd(nome.c_str());
			for(int m=0;m<18;m++) c40_tris[l][m]->Write();
		}
		f_out->cd("Dist Template Fits/Dist Agl Template Fits/Primaries/D.bins");
		for(int m=0;m<18;m++) c40_tris[11][m]->Write();
		f_out->cd("Dist Template Fits/Dist Agl Template Fits/Primaries/P.bins");
		for(int m=0;m<18;m++) c40_tris[0][m]->Write();

	}
	f_out->cd("D Flux (Mass fit)");
	c33->Write();
	c32->Write();
	c34->Write();
	c35->Write();
	f_out->cd("D Flux (Dist. fit)");
        c43->Write();
        c42->Write();
        c44->Write();
        c45->Write();
	f_out->cd("Export");
	MigrMatrix->Write();
	EffTrackerMCP_R_TH1F->Write();
	EffTrackerMCP_TH1F->Write();
	EffTrackerMCD_R_TH2F->Write();
	EffTrackerMCD_TH2F->Write();
	EffTOF_MCP_R_TH1F->Write();
	EffTOF_MCP_TH1F->Write();
	EffTOF_MCD_R_TH2F->Write();
	EffTOF_MCD_TH2F->Write();
	EffMCLikP_TH1F->Write();
	EffMCLikD_TH2F->Write();
	EffMCLikP_Beta_TH1F->Write();
	EffMCLikD_Beta_TH2F->Write();
	EffPreSelMCP_R_TH2F->Write();
	EffPreSelMCP_TH2F->Write();
	EffPreSelMCD_R_TH3F->Write();
	EffPreSelMCD_TH3F->Write();
	for(int S=0;S<3;S++) CorrLAT_pre[S]->Write();
	for(int S=0;S<3;S++) CorrLATpre_spl[S]->Write();
	for(int S=0;S<3;S++) CorrLATpre_Spl[S]->Write();
	CorrLAT_Lik->Write();
	CorrLAT_Dist->Write();
	CorrLAT_Dist_spl->Write();
	CorrLAT_Dist_Spl->Write();
	CorrLAT_Lik_spl->Write();
	CorrLAT_Lik_Spl->Write();
	CorrLAT_LikNaF->Write();
        CorrLAT_DistNaF->Write();
        CorrLAT_DistNaF_spl->Write();
        CorrLAT_DistNaF_Spl->Write();
        CorrLAT_LikNaF_spl->Write();
        CorrLAT_LikNaF_Spl->Write();
	CorrLAT_LikAgl->Write();
        CorrLAT_DistAgl->Write();
        CorrLAT_DistAgl_spl->Write();
        CorrLAT_DistAgl_Spl->Write();
        CorrLAT_LikAgl_spl->Write();
        CorrLAT_LikAgl_Spl->Write();
	LikDVSMC_P_Graph->Write();
	DistDVSMC_P_Graph->Write();
	LikDVSMC_P_graph->Write();
	DistDVSMC_P_graph->Write();
	for(int S=0;S<3;S++) PreDVSMC_P_Graph[S]->Write();
	for(int S=0;S<3;S++) PreDVSMC_P[S]->Write();
	CorrLAT_tot_spl->Write();
	CorrLAT_pre_spl->Write();
	CorrLAT_totM2->Write();
	CorrLAT_preM2->Write();
	for(int i=0;i<11;i++) P_Fluxgeo[i]->Write();
	PFlux->Write();
	P_Flux->Write();
	P_Flux_pre->Write();
	DFluxgeoTOF->Write();
	f_out->cd("Export/Hecut");
	HecutMC_P->Write();
        HecutMC_He->Write();
        Hecut_D->Write();
        HecutMC_P1->Write();
        HecutMC_P2->Write();
        HecutMC_He1->Write();
        HecutMC_He2->Write();
	f_out->cd("Export/MCMC");
	MCMCPTemplatesTOF->Write();
        MCMCDTemplatesTOF->Write();
        MCMCHeTemplatesTOF->Write();
        MCMCPTemplatesNaF->Write();
        MCMCDTemplatesNaF->Write();
        MCMCHeTemplatesNaF->Write();
        MCMCPTemplatesAgl->Write();
        MCMCDTemplatesAgl->Write();
        MCMCHeTemplatesAgl->Write();
        MCMCDataTOF->Write();
        MCMCDataNaF->Write();
        MCMCDataAgl->Write();
	f_out->cd("Export/Qual.Sel.Eff.");
	EffMCDistP_TH1F->Write();
	EffMCDistD_TH2F->Write();
	EffMCDistP_Beta_TH1F->Write();
	EffMCDistD_Beta_TH2F->Write();
	EffMCDistP_BetaNaF_TH1F->Write();
	EffMCDistD_BetaNaF_TH2F->Write();
	EffMCDistP_BetaAgl_TH1F->Write();
	EffMCDistD_BetaAgl_TH2F->Write();	
	f_out->cd("Export/Acceptance");
	AcceptPzone->Write();
        AcceptDzoneTOF->Write();
        AcceptDzoneNaF->Write();
        AcceptDzoneAgl->Write();
	AcceptPTOF->Write();
	AcceptPNaF->Write();
	AcceptPAgl->Write();
	AcceptDTOF->Write();
	AcceptDNaF->Write();
	AcceptDAgl->Write();
	f_out->cd("Export/Slidesplot");
	SlidesforPlot_Write();	
	f_out->Write();
        f_out->Close();
	
	return 1;
}
