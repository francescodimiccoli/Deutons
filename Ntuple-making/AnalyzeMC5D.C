#include "TROOT.h"
#include "TFile.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include <TVector3.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cstring>
#include <vector>
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <TSpline.h>
//#include "Functions_auto.h"
#include "Selections5D.h"

using namespace std;
void aggiungiantupla (TTree *albero,int i,TNtuple *ntupla,int P_ID);
void Grandezzesep (TTree *albero,int i,TNtuple *ntupla);
void Grandezzesepd (TTree *albero,int i,TNtuple *ntupla);
void Grandezzequal (TTree *albero,int i,TNtuple *ntupla);
void GrandezzequalRICH (TTree *albero,int i,TNtuple *ntupla);
void Trigg (TTree *albero,int i,TNtuple *ntupla);
int AssegnaCutmask(int CT,float MG);
int Ev_Num;
int Trig_Num;
double Trig;
double TotalTrig=0;
double Totalevents=0;
double totaltrig=0;
double totaltrig2=0;
double response[44][44];
double norm[44];
int Number=0;
int scelta=0;
int efficienzagenbeta_P[18]={0};
int efficienzagenbeta_D[18][6]={{0}};
int preselectedbeta_P[18]={0};
int preselectedbeta_D[18][6]={{0}};
TH1F * selezioni_PHLMC[10];
TH1F * selected_PHLMC[10];
TH2F * selezioni_DHLMC[10];
TH2F * selected_DHLMC[10];
TH1F * PrescaledMC = new TH1F("PrescaledMC","PrescaledMC",43,0,43);
TH1F * UnbiasMC= new TH1F("UMC","UMC",43,0,43);
TH1F * UPreselectedMC= new TH1F("PreselectedMC","PreselectedMC",43,0,43);

int main(int argc, char * argv[]){

	gROOT->ProcessLine("#include <vector>");
	///////////////////////////// CALIBR.
	int control=0;
	string calib=argv[2];
	string nomecal=("/storage/gpfs_ams/ams/users/fdimicco/CodesforAnalysis/CALIBRAZIONI/"+calib+".root");
	TFile *_file2 = TFile::Open(nomecal.c_str());
	cout<<"calibrazione: "<<_file2<<endl;
	if(!_file2) {nomecal=("/storage/gpfs_ams/ams/users/fdimicco/CodesforAnalysis/2011_07.root");
	_file2 = TFile::Open(nomecal.c_str());
	control=1;
	}
	Rig = (TSpline3 *) _file2->Get("Fit Results/Splines/Rig");
        beta = (TSpline3 *) _file2->Get("Fit Results/Splines/beta");
        betaNaF = (TF1 *) _file2->Get("Fit Results/Splines/SigmaInvBetaNaF_spl");
	betaAgl = (TF1 *) _file2->Get("Fit Results/Splines/SigmaInvBetaAgl_spl");
	eL1 = (TSpline3 *) _file2->Get("Fit Results/Splines/eL1");
        etofu =  (TSpline3 *) _file2->Get("Fit Results/Splines/etofu");
        etrack =  (TSpline3 *) _file2->Get("Fit Results/Splines/etrack");
        etofd =  (TSpline3 *) _file2->Get("Fit Results/Splines/etofd");
        EdepL1beta =  (TSpline3 *) _file2->Get("Fit Results/Splines/EdepL1beta");
        EdepTOFbeta =  (TSpline3 *) _file2->Get("Fit Results/Splines/EdepTOFbeta");
        EdepTrackbeta =  (TSpline3 *) _file2->Get("Fit Results/Splines/EdepTrackbeta");
        EdepTOFDbeta =  (TSpline3 *) _file2->Get("Fit Results/Splines/EdepTOFDbeta");
        Corr_L1 =  (TSpline3 *) _file2->Get("Fit Results/Splines/Corr_L1");
        Corr_TOFU =  (TSpline3 *) _file2->Get("Fit Results/Splines/Corr_TOFU");
        Corr_Track =  (TSpline3 *) _file2->Get("Fit Results/Splines/Corr_Track");
        Corr_TOFD =  (TSpline3 *) _file2->Get("Fit Results/Splines/Corr_TOFD");
	cout<<Rig<<" "<<beta<<" "<<" "<<betaNaF<<" "<<betaAgl<<" "<<eL1<<" "<<etofu<<" "<<etrack<<" "<<etofd<<" "<<EdepL1beta<<" "<<EdepTOFbeta<<" "<<EdepTrackbeta<<" "<<EdepTOFDbeta<<" "<<Corr_L1<<" "<<Corr_TOFU<<" "<<Corr_Track<<" "<<Corr_TOFD<<endl;
	
	string tagli[10]={"Trigger","3of4 TOF","TRD Segments","Rigidity exists","Chi^2 R","Matching TOF","Matching TRD","In TRD Accept.","1 Particle","1 Tr. Track"};
	string nome;

	//////////////////////// DICHIARAZIONE IST. //////////////////////
	for(int j=0;j<10;j++){
		nome="SelezioniHL_MC"+tagli[j];
		selezioni_PHLMC[j]=new TH1F(nome.c_str(),nome.c_str(),43,0,43);
		nome="SelectedHL_MC"+tagli[j];
		selected_PHLMC[j]=new TH1F(nome.c_str(),nome.c_str(),43,0,43);
		nome="SelezioniHL_MCD"+tagli[j];
		selezioni_DHLMC[j]=new TH2F(nome.c_str(),nome.c_str(),18,0,18,6,0,6);
		nome="SelectedHL_MCD"+tagli[j];
		selected_DHLMC[j]=new TH2F(nome.c_str(),nome.c_str(),18,0,18,6,0,6);
	}

	/////////////////////////////////////////////////////////////////

	string Variables[9]={"N. Anti-clusters","Unused TOF Clusters","|Rup-Rdown|:R","Unused Tracker layers","Tracker: Y Hits without X","Track Chi^2","|E.dep(lower TOF) - E.dep(upper TOF)|","|E.dep(layer 2)-E.dep(layer 1)|","|E. dep.(tot)-E.dep.(track)|"};
	for(int u2=0;u2<9;u2++) {
		nome="Splines/Spline: "+Variables[u2]+"_SGNL"; Signal[u2]=(TSpline3 *) _file1->Get(nome.c_str());
		nome="Splines/Spline: "+Variables[u2]+"_BKGND"; Bkgnd[u2]=(TSpline3 *) _file1->Get(nome.c_str());
	}
	cout<<_file1<<endl;
	string VariablesRICH[9]={"N. Anti-clusters","Unused TOF Clusters","|Rup-Rdown|:R","Unused Tracker layers","Tracker: Y Hits without X","Track Chi^2","RICH Hits: tot - used","Rich Photoelectrons","|E. dep.(tot)-E.dep.(track)|"};
	for(int u2=0;u2<9;u2++) {
                nome="Splines/Spline NaF: "+VariablesRICH[u2]+"_SGNL"; SignalNaF[u2]=(TSpline3 *) _file3->Get(nome.c_str());
                nome="Splines/Spline NaF: "+VariablesRICH[u2]+"_BKGND"; BkgndNaF[u2]=(TSpline3 *) _file3->Get(nome.c_str());
	}
        cout<<_file3<<endl;
	for(int u2=0;u2<9;u2++) {
                nome="Splines/Spline Agl: "+VariablesRICH[u2]+"_SGNL"; SignalAgl[u2]=(TSpline3 *) _file3b->Get(nome.c_str());
                nome="Splines/Spline Agl: "+VariablesRICH[u2]+"_BKGND"; BkgndAgl[u2]=(TSpline3 *) _file3b->Get(nome.c_str());
        }
        cout<<_file3b<<endl;

	for(int qs=0;qs<9;qs++) cout<<Signal[qs]<<" ";
	cout<<endl;
	for(int qs=0;qs<9;qs++) cout<<BkgndNaF[qs]<<" ";
        cout<<endl;
	for(int qs=0;qs<9;qs++) cout<<SignalAgl[qs]<<" ";
        cout<<endl;
	
	cout<<"Vuoi Ntuple? (1-Sì;2-No)"<<endl;
	//cin>>scelta;
	scelta=1;
	cout<<"**************************** R BINS ***********************************"<<endl;
	for(int i=0;i<44;i++)
	{
		float temp=i+14;
		bin[i]=0.1*pow(10,temp/(9.5*2));
		//************** bin DAV
		/*        bin[i]=exp(E);
			  E=E+a;*/
		cout<<bin[i]<<endl;

	}

	cout<<"**************************** BETA BINS ***********************************"<<endl;
	float B=0.4;
	float B1=0;
	float B2=0;
	float E=0.1;
	int binnum=1;
	float a=(log(0.9)-log(0.1))/18;
	float E2=exp(log(0.1)+1.5*a);
	float Betabins[18]={0};
	float Betacent[18]={0};
	while(B<0.825){
		B=B+2*(pow(B,2)*beta->Eval(B));
		E=exp(log(0.1)+binnum*a);
		E2=exp(log(0.1)+(binnum+0.5)*a);
		B1=sqrt(1-1/(pow(E+1,2)));
		B2=sqrt(1-1/(pow(E2+1,2)));
		cout<<B<<" "<<binnum<<" "<<B1<<" "<<B2<<endl;
		Betabins[binnum-1]=B1;
		Betacent[binnum-1]=B2;
		binnum++;
	}
	
	float encinTOF[18];
	for(int i=0;i<18;i++) encinTOF[i]=1/pow(1-pow(Betabins[i],2),0.5)-1;
	

	TFile *file =TFile::Open("/storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/sommaMC.root");
	TTree *geo_stuff = (TTree *)file->Get("parametri_geo");
	string ARGV(argv[1]);
	string indirizzo_out="/storage/gpfs_ams/ams/users/fdimicco/Risultati/risultati/RisultatiMC_"+ARGV+".root";
	TFile * File = new TFile(indirizzo_out.c_str(), "RECREATE");
	TNtuple *grandezzequal = new TNtuple("grandezzequal","grandezzequal","Beta:Massa_gen:R:NAnticluster:Clusterinutili:DiffR:fuoriX:EdepL1:layernonusati:Chisquare:EdepTOFU:EdepTOFD:Cutmask:Momentogen:DistD");
	TNtuple *grandezzequalRICH = new TNtuple("grandezzequalRICH","grandezzequalRICH","BetaRICH_new:Massa_gen:R:NAnticluster:Clusterinutili:DiffR:fuoriX:EdepL1:layernonusati:Chisquare:Richtotused:RichPhEl:Cutmask:Momentogen:DistD");
	TNtuple *grandezzesepd = new TNtuple("grandezzesepd","grandezzesepd","R:Beta:EdepL1:Massagen:Cutmask:X:YTOFU:YTrack:YTOFD:Momentogen:BetaRICH_new:LDiscriminant:BDT_response:Dist5D:Dist5D_P");
	TNtuple * pre = new TNtuple("pre","distr for giov","R:Beta:EdepL1:EdepTOFU:EdepTrack:EdepTOFD:EdepECAL:Massagen:Momentogen:Betagen:Dist5D:Dist5D_P:BetaRICH_new:Cutmask:LDiscriminant");
	TNtuple * trig = new TNtuple("trig","trig","Massagen:Momento_gen:Ev_Num:Trig_Num:R_pre:Beta_pre:Cutmask:EdepL1:EdepTOFU:EdepTOFD:EdepTrack:BetaRICH:EdepECAL:Unbias");
	
	BDTreader();
	geo_stuff->SetBranchAddress("Momento_gen",&Momento_gen);
	geo_stuff->SetBranchAddress("Beta_gen",&Beta_gen);
	geo_stuff->SetBranchAddress("Massa_gen",&Massa_gen);
	geo_stuff->SetBranchAddress("Ev_Num",&Ev_Num);
	geo_stuff->SetBranchAddress("Trig_Num",&Trig_Num);
	geo_stuff->SetBranchAddress("Pres_Unbias",&Pres_Unbias);
	geo_stuff->SetBranchAddress("Preselected",&Preselected);
	geo_stuff->SetBranchAddress("Beta_pre",&Beta_pre);
	geo_stuff->SetBranchAddress("R_pre",&R_pre);
	geo_stuff->SetBranchAddress("CUTMASK",&CUTMASK);
	geo_stuff->SetBranchAddress("trtrack_edep",&trtrack_edep);
	geo_stuff->SetBranchAddress("trtot_edep",&trtot_edep);
	geo_stuff->SetBranchAddress("Endep",&Endep);
	geo_stuff->SetBranchAddress("BetaRICH_new",&BetaRICH_new);
	geo_stuff->SetBranchAddress("RICHmask_new",&RICHmask_new);
	geo_stuff->SetBranchAddress("EdepECAL",&EdepECAL);	
	geo_stuff->SetBranchAddress("Unbias",&Unbias);
	geo_stuff->SetBranchAddress("Beta_gen",&Beta_gen);
	geo_stuff->SetBranchAddress("Massa_gen",&Massa_gen);
	geo_stuff->SetBranchAddress("layernonusati",&layernonusati);
	geo_stuff->SetBranchAddress("NAnticluster",&NAnticluster);
	geo_stuff->SetBranchAddress("NTRDSegments",&NTRDSegments);
	geo_stuff->SetBranchAddress("NTofClusters",&NTofClusters);
	geo_stuff->SetBranchAddress("NTofClustersusati",&NTofClustersusati);
	geo_stuff->SetBranchAddress("Rup",&Rup);
	geo_stuff->SetBranchAddress("Rdown",&Rdown);
	geo_stuff->SetBranchAddress("R",&R);
	geo_stuff->SetBranchAddress("Chisquare",&Chisquare);
	geo_stuff->SetBranchAddress("ResiduiX",&ResiduiX);
	geo_stuff->SetBranchAddress("ResiduiY",&ResiduiY);
	geo_stuff->SetBranchAddress("chiq",&chiq);
	geo_stuff->SetBranchAddress("R_",&R_);
	geo_stuff->SetBranchAddress("Beta",&Beta);
	geo_stuff->SetBranchAddress("Betacorr",&Betacorr);
	geo_stuff->SetBranchAddress("NTrackHits",&NTrackHits);
	geo_stuff->SetBranchAddress("clustertrack",&clusterTrack);
	geo_stuff->SetBranchAddress("clustertottrack",&clustertotTrack);
	geo_stuff->SetBranchAddress("zonageo",&zonageo);
	geo_stuff->SetBranchAddress("Massa",&Massa);
	geo_stuff->SetBranchAddress("Richtotused",&Richtotused);
	geo_stuff->SetBranchAddress("RichPhEl",&RichPhEl);

	int giov=0;
	int nprotoni=0;
	int events =geo_stuff->GetEntries();
	int INDX=atoi(argv[1]);
	cout<<endl;
	cout<<INDX<<endl;
	cout<<endl;
	float entries=0;
	cout<<"Eventi: "<<events<<endl;
	//entries=entries/5;
	for(int i=(events/100)*INDX;i<(events/100)*(INDX+1);i++)
	{
		if(i%1300==0) cout<<i/(float)events*100<<"%"<<endl;
		int k = geo_stuff->GetEvent(i);	
		if(Massa_gen>1.8569&&Massa_gen<1.8571) Massa_gen=1.8570; 
		if(Trig_Num<Trig) TotalTrig=TotalTrig+(double)Trig; 	
		totaltrig=TotalTrig+Trig;
		Trig=Trig_Num;
		R_corr=R;
		Particle_ID=0;
		Cutmask=CUTMASK;
		Cutmask=CUTMASK|(1<<10);
		Cutmask = Cutmask|(RICHmask_new<<11);
		if(!(((Cutmask&187)==187))) continue;
		entries++;
		if(Unbias==1) continue;
		if (Quality(geo_stuff,i)){
			giov++;
			if(scelta==1) aggiungiantupla(geo_stuff,i,pre,0);
			if(control==1) continue;
			//////////////////// CORR EDEP /////////////////////
			EdepTOFU=((EdepTOFU)*Corr_TOFU->Eval(Beta));
			EdepTOFD=((EdepTOFD)*Corr_TOFD->Eval(Beta));
			//E_depTRD=((E_depTRD)*Corr_TRD->Eval(Beta));
			EdepTrack=((EdepTrack)*Corr_Track->Eval(Beta));
			///////////////////////////////////////////////////// 
			//////////////// MATRICE DI RISPOSTA ///////////////
			if(Massa_gen<1&&Massa_gen>0.5){
				for(int I=0;I<44;I++)
					if(fabs(Momento_gen)<bin[I+1]&&fabs(Momento_gen)>bin[I])
						for(int J=0;J<44;J++)
							if(fabs(R_corr)<bin[J+1]&&fabs(R_corr)>bin[J])  response[J][I]++;
				for(int I=0;I<44;I++)
					if(fabs(Momento_gen)<bin[I+1]&&fabs(Momento_gen)>bin[I]) norm[I]++;
			}			
			////////////////////////////////////////////////////
			Protoni(geo_stuff,i);
			if (Deutoni(geo_stuff,i)) if(scelta==1) Grandezzesepd(geo_stuff,i,grandezzesepd); 
		} 
		if(control==1) continue;
		if(scelta==1) Grandezzequal(geo_stuff,i,grandezzequal);
		if(scelta==1) GrandezzequalRICH(geo_stuff,i,grandezzequalRICH);													
	}

	//ANALISI SET COMPLETO
	cout<<"Trigger tot: "<<totaltrig<<endl;
	TotalTrig=0;
	Trig=0;
	int z;
	for(z=(events/100)*INDX;z<(events/100)*(INDX+1);z++)
	{
		int k = geo_stuff->GetEvent(z);
		if(Massa_gen>1.8569&&Massa_gen<1.8571) Massa_gen=1.8570;
		if(Trig_Num<Trig)TotalTrig=TotalTrig+(double)Trig;
		totaltrig2=TotalTrig+Trig;
		Trig=Trig_Num;
		if(z%100000==0) cout<<z/(float)events*100<<" "<<totaltrig2<<endl;
		
		Cutmask=CUTMASK;
                Cutmask=CUTMASK|(1<<10);
                Cutmask = Cutmask|(RICHmask_new<<11);
                //Cutmask = AssegnaCutmask(Cutmask,Massa_gen);
		
		EdepTrack=0;
	        EdepTOFU=((*Endep)[0]+(*Endep)[1])/2;
	        EdepTOFD=((*Endep)[2]+(*Endep)[3])/2;
       	 	EdepTOFU=((EdepTOFU)*Corr_TOFU->Eval(Beta));
                EdepTOFD=((EdepTOFD)*Corr_TOFD->Eval(Beta));
                //E_depTRD=((E_depTRD)*Corr_TRD->Eval(Beta));
		for(int layer=1;layer<8;layer++) EdepTrack+=(*trtot_edep)[layer]; EdepTrack=EdepTrack/7;
		EdepTrack=((EdepTrack)*Corr_Track->Eval(Beta));	
	
		if(scelta==1) Trigg(geo_stuff,z,trig);
		
	}




	cout<<"Eventi Tot: "<<z<<endl;
	cout<<"Preselezionate Tot: "<<entries<<endl;
	cout<<endl;
	cout<<"selezioni di qualità: "<<giov<<endl;
	cout<<"N. Protoni: "<<nprotoni<<endl;
	cout<<endl;
	cout<<"-----------------------"<<endl;
	cout<<"-----EFF. SELEZIONI (risp a PRE)----"<<endl;
	cout<<"Taglio anticluster: "<<f/(float)entries*100<<endl;
	cout<<"Taglio segmenti TRD: "<<g/(float)entries*100<<endl;
	cout<<"Taglio TOF Clusters inutilizzati: "<<h/(float)entries*100<<endl;
	cout<<"Taglio Confronto Rup/Rdown: "<<l/(float)entries*100<<endl;
	cout<<"Taglio Prob. Q: "<<m/(float)entries*100<<endl;
	cout<<"Taglio Fit multipli: "<<n/(float)entries*100<<endl;
	cout<<"Taglio HitX: "<<o/(float)entries*100<<endl;
	cout<<"Taglio Chi-Quadro: "<<r/(float)entries*100<<endl;
	cout<<"Taglio carica TOF: "<<t/(float)entries*100<<endl;
	cout<<"Taglio carica Track: "<<u/(float)entries*100<<endl;
	cout<<"Taglio carica TRD: "<<v/(float)entries*100<<endl;
	cout<<"Taglio layer non usati: "<<s/(float)entries*100<<endl;
	cout<<endl;
	cout<<"-----SELEZIONI DI REIEZIONE----"<<endl;
	cout<<"Cluster TOF: "<<b/(float)giov*100<<endl;
	cout<<endl;
	cout<<"Cluster Track: "<<d/(float)giov*100<<endl;
	cout<<endl;
	cout<<"Cluster TRD: "<<e/(float)giov*100<<endl;
	cout<<endl;

	cout<<"-----ANALISI EFFICIENZA (PROTONI)------------- "<<endl;
	cout<<endl;
	for(int i=0; i<43;i++)
		if(efficienzagen[i]!=0)
			cout<<"efficienza da"<<bin[i]<<" a " << bin[i+1]<<": "<<protoni[i][0]<<" nel bin " <<i<<" "<<( protoni[i][0]/(double)efficienzagen[i]*100)<<"%"<< " rispetto a preselezione:  "<<protoni[i][0]/(double)preselezionate[i][0]*100<<"%"<<endl;

		else cout<<"efficienza da "<<bin[i]<<" a "<< bin[i+1]<<": nessun evento generato nel bin"<<endl;
	cout<<endl;   
	cout<<"-----ANALISI EFFICIENZA (DEUTONI)------------- "<<endl;
	cout<<endl;
	for(int i=0; i<43;i++)
		if(efficienzagen_D[i]!=0)
			cout<<"efficienza da"<<bin[i]<<" a " << bin[i+1]<<": "<<deutoni[i][0]<<" nel bin " <<i<<" "<<( deutoni[i][0]/(double)efficienzagen_D[i][0]*100)<<"%"<<endl;

		else cout<<"efficienza da "<<bin[i]<<" a "<< bin[i+1]<<": nessun evento generato nel bin"<<endl;
	cout<<endl;       
	cout<<"----ANALISI FONDO----"<<endl;
	cout<<endl;
	for(int i=0; i<43;i++) {
		if(efficienzagen[i]!=0&&background[i]>0)
			cout<<"Fondo da"<<bin[i]<<" a " << bin[i+1]<<": "<<background[i]<<" nel bin " <<i<<" su "<< efficienzagen[i]<<" : "<<(background[i]/(double)efficienzagen[i]*100)<<"% ---> R.P. : "<< protoni[i][0]/background[i]<<endl;

		if(efficienzagen[i]!=0&&background[i]==0)
			cout<<"Fondo da"<<bin[i]<<" a " << bin[i+1]<<": "<<background[i]<<" nel bin " <<i<<" su "<<efficienzagen[i]<<" : "<<"Minore di: "<<(1/(double)efficienzagen[i]*100)<<"% ---> R.P. Maggiore di: "<<protoni[i][0]<<endl;

		if(efficienzagen[i]==0) cout<<"Fondo da "<<bin[i]<<" a "<< bin[i+1]<<": nessun evento generato nel bin"<<endl;
	}

	cout<<"------------------- MATRICE DI RISPOSTA -------------"<<endl;
	for(int j=0; j<43;j++){
		for(int i=0; i<43;i++)
			cout<<response[j][i]<<" ";
		cout<<endl;}

		cout<<"------------------- NORM -------------"<<endl;
		for(int i=0; i<43;i++)
			cout<<norm[i]<<" ";
		cout<<endl;

		/*
		TH1F * preselezionate0=new TH1F("preselezionate0","preselezionate0",43,0,43);
		for(int j=0; j<43;j++) preselezionate0->SetBinContent(j+1,preselezionate[j][0]);

		TH2F * preselezionateD=new TH2F("preselezionateD","preselezionateD",43,0,43,6,0,6);
		for(int j=0; j<43;j++) for(int f=0;f<6;f++) preselezionateD->SetBinContent(j+1,f+1,preselezionate_D[j][f]);

		TH1F * Preselectedbeta_P=new TH1F("preselectedbeta_P","preselectedbeta_P",18,0,18);
		for(int j=0; j<18;j++) Preselectedbeta_P->SetBinContent(j+1,preselectedbeta_P[j]);

		TH2F * Preselectedbeta_D=new TH2F("preselectedbeta_D","preselectedbeta_D",18,0,18,6,0,6);
		for(int j=0; j<18;j++) for(int f=0;f<6;f++) Preselectedbeta_D->SetBinContent(j+1,f+1,preselectedbeta_D[j][f]);

		TH1F * tempi0=new TH1F("tempi0","tempi0",43,0,43);
		for(int j=0;j<43;j++) tempi0->SetBinContent(j+1,efficienzagen[j]);

		TH2F * tempi0_D=new TH2F("tempi0_D","tempi0_D",43,0,43,6,0,6);
		for(int j=0;j<43;j++) for(int f=0;f<6;f++) tempi0_D->SetBinContent(j+1,f+1,efficienzagen_D[j][f]);

		TH1F * Efficienzagenbeta_P=new TH1F("efficienzagenbeta_P","efficienzagenbeta_P",18,0,18);
		for(int j=0; j<18;j++) Efficienzagenbeta_P->SetBinContent(j+1,efficienzagenbeta_P[j]);

		TH2F * Efficienzagenbeta_D=new TH2F("efficienzagenbeta_D","efficienzagenbeta_D",18,0,18,6,0,6);
		for(int j=0; j<18;j++) for(int f=0;f<6;f++) Efficienzagenbeta_D->SetBinContent(j+1,f+1,efficienzagenbeta_D[j][f]);

		preselezionate0->Write();
		preselezionateD->Write();
		Preselectedbeta_P->Write();
		Preselectedbeta_D->Write();
		tempi0->Write();
		tempi0_D->Write();
		tempi0_D->Write();
		Efficienzagenbeta_P->Write();
		Efficienzagenbeta_D->Write();
		for(int S=0;S<10;S++)  selezioni_PHLMC[S]->Write();
		for(int S=0;S<10;S++)  selected_PHLMC[S]->Write();
		for(int S=0;S<10;S++)  selezioni_DHLMC[S]->Write();
		for(int S=0;S<10;S++)  selected_DHLMC[S]->Write();
		*/
		
		if(scelta==1) File->Write();
		File->Close();
		return 1;
}

int AssegnaCutmask(int CT,float MG)
{
	int Cutmask=CT;
	if(MG<1) Cutmask=Cutmask | (1<<22);
	if(MG>1&&MG<1.85705) Cutmask=Cutmask | (1<<23);
	if(MG>1.85705&&MG<1.85715) Cutmask=Cutmask | (1<<24);
	if(MG>1.85715&&MG<1.85725) Cutmask=Cutmask | (1<<25);
	if(MG>1.85725&&MG<1.85735) Cutmask=Cutmask | (1<<26);
	if(MG>1.85735&&MG<1.85745) Cutmask=Cutmask | (1<<27);
	if(MG>1.85745&&MG<1.85755) Cutmask=Cutmask | (1<<28);
	if(MG>3&&MG<4) Cutmask=Cutmask | (1<<29);
	return Cutmask;
}

void Trigg (TTree *albero,int i,TNtuple *ntupla)
{
	int k = albero->GetEvent(i);
	ntupla->Fill(Massa_gen,Momento_gen,Ev_Num,Trig_Num,R_pre,Beta_pre,Cutmask,(*trtrack_edep)[0],EdepTOFU,EdepTOFD,EdepTrack,BetaRICH_new,EdepECAL,Unbias);
		
}

void aggiungiantupla (TTree *albero,int i,TNtuple *ntupla,int P_ID)
{
	int k = albero->GetEvent(i);
	ntupla->Fill(R,Beta,(*trtrack_edep)[0],EdepTOFU,EdepTrack,EdepTOFD,EdepECAL,Massa_gen,Momento_gen,Beta_gen,Dist5D,Dist5D_P,BetaRICH_new,Cutmask,LDiscriminant);

}

void Grandezzequal (TTree *albero,int i,TNtuple *ntupla)
{
	int k = albero->GetEvent(i);
	ntupla->Fill(Beta,Massa_gen,R,NAnticluster,NTofClusters-NTofClustersusati,fabs(Rup-Rdown)/R,fuoriX,(*trtrack_edep)[0],layernonusati,Chisquare,EdepTOFU,EdepTOFD,Cutmask,Momento_gen,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
}

void GrandezzequalRICH (TTree *albero,int i,TNtuple *ntupla)
{
        int k = albero->GetEvent(i);
        ntupla->Fill(BetaRICH_new,Massa_gen,R,NAnticluster,NTofClusters-NTofClustersusati,fabs(Rup-Rdown)/R,fuoriX,(*trtrack_edep)[0],layernonusati,Chisquare,Richtotused,RichPhEl,Cutmask,Momento_gen,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
}

void Grandezzesepd (TTree *albero,int i,TNtuple *ntupla)
{

	int k = albero->GetEvent(i);
	ntupla->Fill(R,Beta,(*trtrack_edep)[0],Massa_gen,Cutmask,PSCALTOF2,PSCALTOF3,PSCALTrack3,PSCALTRD3,Momento_gen,BetaRICH_new,LDiscriminant,BDT_response,Dist5D,Dist5D_P);
}


