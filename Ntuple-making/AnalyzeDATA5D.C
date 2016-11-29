#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "Selections5D.h"
#include "Commonglobals.cpp"
#include "../include/binning.h"

void Controllodist (TTree *albero,int i,TNtuple *ntupla);
void Controllodistd (TTree *albero,int i,TNtuple *ntupla);
void Selez (TTree *albero,int i,TNtuple *ntupla);
int giov=0;
int nprotoni=0;
int entriestot=0;
int events=0;
int INDX;
float BetanS=0;
float Encinp;
float Encind;
int PhysBPatt=0;
float esposizionep[18]= {0};
float esposizioned[18]= {0};
float esposizionepgeo[18][11]= {{0}};
float esposizionedgeo[18][11]= {{0}};
float esposizionepNaF[18]= {0};
float esposizionedNaF[18]= {0};
float esposizionepgeoNaF[18][11]= {{0}};
float esposizionedgeoNaF[18][11]= {{0}};
float esposizionepAgl[18]= {0};
float esposizionedAgl[18]= {0};
float esposizionepgeoAgl[18][11]= {{0}};
float esposizionedgeoAgl[18][11]= {{0}};
string numero[11]= {"0","1","2","3","4","5","6","7","8","9","10"};
float Rcut[11]= {18,18,16,14,12,10,8,6,4,2,1};
TH2F * unbias = new TH2F("Unbias","Unbias",43,0,43,11,0,11);
TH2F * UPreselected = new TH2F("Upreselected","Upreselected",43,0,43,11,0,11);
TH1F * UnbiasHL = new TH1F("UnbiasHL","UnbiasHL",43,0,43);
TH1F * UPreselectedHL = new TH1F("UpreselectedHL","Upreselected",43,0,43);
TH2F * selezioni_P[10];
TH2F * selected_P[10];
TH1F * selezioni_PHL[10];
TH1F * selected_PHL[10];
TH1F * selezioni_DHL[10];
TH1F * selected_DHL[10];


float mcweight=1;

int main(int argc, char * argv[])
{

	if (argc<3) return 0;
	string nome;
	string tagli[10]= {"Trigger","3of4 TOF","TRD Segments","Rigidity exists","Chi^2 R","Matching TOF","Matching TRD","In TRD Accept.","1 Particle","1 Tr. Track"};
	/////////// CALIBR.
	int control=0;
	string calib=argv[2];
	string nomecal=("/storage/gpfs_ams/ams/users/fdimicco/Deutons/CodesforAnalysis/CALIBRAZIONI/"+calib+".root");
	TFile *_file2 = TFile::Open(nomecal.c_str());
	cout<<"calibrazione: "<<_file2<<endl;
	if(!_file2) {
		nomecal=("/storage/gpfs_ams/ams/users/fdimicco/Deutons/CodesforAnalysis/2011_07.root");
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
	cout<<Rig<<" "<<beta<<" "<<betaNaF<<" "<<betaAgl<<" "<<eL1<<" "<<etofu<<" "<<etrack<<" "<<etofd<<" "<<EdepL1beta<<" "<<EdepTOFbeta<<" "<<EdepTrackbeta<<" "<<EdepTOFDbeta<<" "<<Corr_L1<<" "<<Corr_TOFU<<" "<<Corr_Track<<" "<<Corr_TOFD<<endl;

	//////////////////////// DICHIARAZIONE IST. //////////////////////
	for(int j=0; j<10; j++) {
		nome="Selezioni"+tagli[j];
		selezioni_P[j]=new TH2F(nome.c_str(),nome.c_str(),43,0,43,11,0,11);
		nome="Selected"+tagli[j];
		selected_P[j]=new TH2F(nome.c_str(),nome.c_str(),43,0,43,11,0,11);
		nome="SelezioniHL"+tagli[j];
		selezioni_PHL[j]=new TH1F(nome.c_str(),nome.c_str(),43,0,43);
		nome="SelectedHL"+tagli[j];
		selected_PHL[j]=new TH1F(nome.c_str(),nome.c_str(),43,0,43);
		nome="SelezioniHL_D"+tagli[j];
		selezioni_DHL[j]=new TH1F(nome.c_str(),nome.c_str(),43,0,43);
		nome="SelectedHL_D"+tagli[j];
		selected_DHL[j]=new TH1F(nome.c_str(),nome.c_str(),43,0,43);
		nome="SelezioniHL_MC"+tagli[j];
	}
	//////////////////////////////////////////////////////////////////
	string Variables[9]= {"N. Anti-clusters","Unused TOF Clusters","|Rup-Rdown|:R","Unused Tracker layers","Tracker: Y Hits without X","Track Chi^2","|E.dep(lower TOF) - E.dep(upper TOF)|","|E.dep(layer 2)-E.dep(layer 1)|","|E. dep.(tot)-E.dep.(track)|"};
	for(int u2=0; u2<9; u2++) {
		nome="Splines/Spline: "+Variables[u2]+"_SGNL";
		Signal[u2]=(TSpline3 *) _file1->Get(nome.c_str());
		nome="Splines/Spline: "+Variables[u2]+"_BKGND";
		Bkgnd[u2]=(TSpline3 *) _file1->Get(nome.c_str());
	}
	cout<<_file1<<endl;
	string VariablesRICH[9]= {"N. Anti-clusters","Unused TOF Clusters","|Rup-Rdown|:R","Unused Tracker layers","Tracker: Y Hits without X","Track Chi^2","RICH Hits: tot - used","Rich Photoelectrons","|E. dep.(tot)-E.dep.(track)|"};
	for(int u2=0; u2<9; u2++) {
		nome="Splines/Spline NaF: "+VariablesRICH[u2]+"_SGNL";
		SignalNaF[u2]=(TSpline3 *) _file3->Get(nome.c_str());
		nome="Splines/Spline NaF: "+VariablesRICH[u2]+"_BKGND";
		BkgndNaF[u2]=(TSpline3 *) _file3->Get(nome.c_str());
	}
	cout<<_file3<<endl;
	for(int u2=0; u2<9; u2++) {
		nome="Splines/Spline Agl: "+VariablesRICH[u2]+"_SGNL";
		SignalAgl[u2]=(TSpline3 *) _file3b->Get(nome.c_str());
		nome="Splines/Spline Agl: "+VariablesRICH[u2]+"_BKGND";
		BkgndAgl[u2]=(TSpline3 *) _file3b->Get(nome.c_str());
	}
	cout<<_file3b<<endl;

	for(int qs=0; qs<9; qs++) cout<<Signal[qs]<<" ";
	cout<<endl;
	for(int qs=0; qs<9; qs++) cout<<BkgndNaF[qs]<<" ";
	cout<<endl;
	for(int qs=0; qs<9; qs++) cout<<SignalAgl[qs]<<" ";
	cout<<endl;


	cout<<"Vuoi Ntuple? (1-Sì;2-No)"<<endl;
	//cin>>scelta;
	scelta=1;
	cout<<"**************************** R BINS ***********************************"<<endl;
	Particle proton(0.9382720813, 1, 1);
	Binning PRB(proton);

	PRB.setBinsFromRigidity(43, 0.5, 100);
	PRB.Print();



	cout<<"**************************** BETA BINS ***********************************"<<endl;

	Particle deuton(1.8756129   , 1, 2);

	Binning ToFDB(deuton);
	Binning ToFPB(proton);
	
	float ekmin=0.1, ekmax=1;
   	ToFDB.setBinsFromEkPerMass (18, ekmin, ekmax);
	ToFPB.setBinsFromEkPerMass (18, ekmin, ekmax);
	
	ToFDB.Print();
	

	cout<<"**************************** BETA BINS NaF***********************************"<<endl;

	Binning NaFDB(deuton);
	Binning NaFPB(proton);        

        ekmin=0.666; ekmax=4.025;
        NaFDB.setBinsFromEkPerMass (18, ekmin, ekmax);
        NaFPB.setBinsFromEkPerMass (18, ekmin, ekmax);
	
	NaFDB.Print();


	cout<<"**************************** BETA BINS Agl***********************************"<<endl;

	Binning AglDB(deuton);
	Binning AglPB(proton);

        ekmin=2.57; ekmax=9.01;
        AglDB.setBinsFromEkPerMass (18, ekmin, ekmax);
	AglPB.setBinsFromEkPerMass (18, ekmin, ekmax);

        AglDB.Print();






	string ARGV(argv[1]);

	string indirizzo_in="/storage/gpfs_ams/ams/users/fdimicco/MAIN/sommadati/sommadati"+ARGV+".root";
	TFile *file =TFile::Open(indirizzo_in.c_str());
	TTree *geo_stuff = (TTree *)file->Get("parametri_geo");
	std::cout<<geo_stuff<<endl;
	bool check=geo_stuff?true:false;

	if(check==false) std::cout<<"Skipping file!"<<endl;
	string indirizzo_out="/storage/gpfs_ams/ams/users/fdimicco/Deutons/Risultati/"+calib+"/RisultatiDATI_"+ARGV+".root";

	TFile * File = new TFile(indirizzo_out.c_str(), "RECREATE");

	TNtuple *grandezzequal = new TNtuple("grandezzequal","grandezzequal","Velocity:Rcutoff:R:NAnticluster:Clusterinutili:DiffR:fuoriX:layernonusati:Chisquare:Richtotused:RichPhEl:Cutmask:Latitude:EdepTRD:IsCharge1");
	TNtuple *grandezzesepd = new TNtuple("grandezzesepd","grandezzesepd","R:Beta:EdepL1:Cutmask:Latitude:PhysBPatt:EdepTOFU:EdepTrack:EdepTOFD:Rcutoff:BetaRICH_new:LDiscriminant:Dist5D:Dist5D_P");
	TNtuple * pre = new TNtuple("Pre","distr for qual","R:Beta:EdepL1:EdepTOFU:EdepTOFD:EdepTrack:EdepECAL:Rcutoff:Latitude:BetaRICH_new:Cutmask:BetanS:BetaR");
	TNtuple * trig = new TNtuple("trig","trig","Seconds:Latitude:Rcutoff:R_L1:R_pre:Beta_pre:Cutmask:EdepL1:EdepTOFU:EdepTOFD:EdepTrack:BetaRICH:EdepECAL:PhysBPatt:Livetime");


	BDTreader();
	if(check) {
		geo_stuff->SetBranchAddress("U_time",&U_time);
		geo_stuff->SetBranchAddress("Latitude",&Latitude);
		geo_stuff->SetBranchAddress("zonageo",&zonageo);
		geo_stuff->SetBranchAddress("Rcutoff",&Rcutoff);
		geo_stuff->SetBranchAddress("Livetime",&Livetime);
		geo_stuff->SetBranchAddress("PhysBPatt",&PhysBPatt);
		geo_stuff->SetBranchAddress("Pres_Unbias",&Pres_Unbias);
		geo_stuff->SetBranchAddress("Preselected",&Preselected);
		geo_stuff->SetBranchAddress("R_pre",&R_pre);
		geo_stuff->SetBranchAddress("Beta_pre",&Beta_pre);
		geo_stuff->SetBranchAddress("CUTMASK",&CUTMASK);
		geo_stuff->SetBranchAddress("trtrack_edep",&trtrack_edep);
		geo_stuff->SetBranchAddress("trtot_edep",&trtot_edep);
		geo_stuff->SetBranchAddress("Endep",&Endep);
		geo_stuff->SetBranchAddress("BetaRICH_new",&BetaRICH_new);
		geo_stuff->SetBranchAddress("RICHmask_new",&RICHmask_new);
		geo_stuff->SetBranchAddress("EdepECAL",&EdepECAL);
		geo_stuff->SetBranchAddress("layernonusati",&layernonusati);
		geo_stuff->SetBranchAddress("NAnticluster",&NAnticluster);
		geo_stuff->SetBranchAddress("NTofClusters",&NTofClusters);
		geo_stuff->SetBranchAddress("NTofClustersusati",&NTofClustersusati);
		geo_stuff->SetBranchAddress("Rup",&Rup);
		geo_stuff->SetBranchAddress("Rdown",&Rdown);
		geo_stuff->SetBranchAddress("R",&R);
		geo_stuff->SetBranchAddress("chiq",&chiq);
		geo_stuff->SetBranchAddress("R_",&R_);
		geo_stuff->SetBranchAddress("Chisquare",&Chisquare);
		geo_stuff->SetBranchAddress("ResiduiX",&ResiduiX);
		geo_stuff->SetBranchAddress("ResiduiY",&ResiduiY);
		geo_stuff->SetBranchAddress("Beta",&Beta);
		geo_stuff->SetBranchAddress("BetaR",&BetaR);
		geo_stuff->SetBranchAddress("NTrackHits",&NTrackHits);
		geo_stuff->SetBranchAddress("zonageo",&zonageo);
		geo_stuff->SetBranchAddress("Massa",&Massa);
		geo_stuff->SetBranchAddress("Richtotused",&Richtotused);
		geo_stuff->SetBranchAddress("RichPhEl",&RichPhEl);
		geo_stuff->SetBranchAddress("R_L1",&R_L1);	
	}

	giov=0;
	nprotoni=0;
	if(check) events =geo_stuff->GetEntries();
	INDX=atoi(argv[1]);
	cout<<endl;
	cout<<INDX<<endl;
	cout<<endl;
	cout<<"Eventi: "<<events<<endl;
	int entries=0;
	if(entries>entriestot) entries=entriestot;

	for(int i=0; i<events; i++) {
		if(!check) break;
		if(i%1300==0) cout<<i/(float)events*100<<"%"<<endl;
		geo_stuff->GetEvent(i);
		R_corr=R;
		if(i==0) { tbeg=U_time; cout <<"Tempo Iniziale: "<<tbeg<<endl;}
		if(i==events-1) {tend=U_time; cout <<"Tempo Finale: "<<tend<<endl;}
		for(int I=0; I<43; I++) if(R<bin[I+1]&&R>bin[I]) preselezionate[I][zonageo]++;
		Particle_ID=0;

		Cutmask=CUTMASK;
		Cutmask=CUTMASK|(1<<10);
		//Temp. rich bug fixing
		if(((Cutmask>>11)==0||(Cutmask>>11)==512)&&BetaRICH_new==-1) RICHmask_new=1;
		Cutmask = Cutmask|(RICHmask_new<<11);
		if(!(((Cutmask&187)==187))) continue;
		entries++;
		if (Quality(geo_stuff,i)) {
			giov++;
			if(scelta==1) aggiungiantupla(geo_stuff,i,pre);
			if(control!=1) {
				Protoni(geo_stuff,i);
				if (Deutoni(geo_stuff,i)) if(scelta==1) Grandezzesepd(geo_stuff,i,grandezzesepd);
			}
		}
		if(scelta==1) Grandezzequal(geo_stuff,i,grandezzequal);
	}

	cout<<"------Calcolo Live Time-----"<<endl;
	zona=0;
	for(int j=0; j<43; j++) for(int i=0; i<11; i++) tempobingeo[j][i]=0;
	int contriniz=0;
	int zona1=0;
	int zona2=0;
	int contaeventi=0;
	int z=0;
	int Tempo=0;
	seconds=0;
	for(z=0; z<events; z++) {
		if(!check) break;
		geo_stuff->GetEvent(z);
		contaeventi++;
		R=R_pre;

		if(z==0) { tbeg=U_time; cout <<"Tempo Iniziale: "<<tbeg<<endl;}	
		seconds = U_time - tbeg;
		
		if(z%1000000==0) cout<<z/(float)events*100<<"%"<<endl;
	
		for(int i=0; i<11; i++) {
			double geo= geomag[i]  ;
			double geo2=geomag[i+1];
			if(Latitude>geo && Latitude<=geo2)
				zona=i;

			else tempozona[i]=U_time;
		}
		if(U_time!=tempozona[zona]) {
			tempozona[zona]=U_time;
			Time[zona]=Time[zona]+Livetime;
			contasecondi[zona]++;
		}


		Cutmask=CUTMASK;
		Cutmask=CUTMASK|(1<<10);
		Cutmask = Cutmask|(RICHmask_new<<11);

		EdepTrack=0;
		EdepTOFU=((*Endep)[0]+(*Endep)[1])/2;
		EdepTOFD=((*Endep)[2]+(*Endep)[3])/2;

		for(int layer=1; layer<8; layer++) EdepTrack+=(*trtot_edep)[layer];
		EdepTrack=EdepTrack/7;
		if(scelta==1) Trigg(geo_stuff,z,trig);
	}
	/////////////////////////////////////////////////////

	for(int i=0; i<11; i++) {
		for(int j=0; j<43; j++) std::cout<<tempobingeo[j][i]<<" ";
		std::cout<<endl;
	}


	if(scelta==1) File->Write();
	if(scelta==1) File->Close();

	return 1;
}

void Trigg (TTree *albero,int i,TNtuple *ntupla)
{
	albero->GetEvent(i);
	ntupla->Fill(seconds,Latitude,Rcutoff,R_L1,R_pre,Beta_pre,Cutmask,(*trtrack_edep)[0],EdepTOFU,EdepTOFD,EdepTrack,BetaRICH_new,EdepECAL,PhysBPatt,Livetime);
}


void aggiungiantupla (TTree *albero,int i,TNtuple *ntupla)
{
	albero->GetEvent(i);
	ntupla->Fill(R,Beta,(*trtrack_edep)[0],EdepTOFU,EdepTOFD,EdepTrack,EdepECAL,Rcutoff,Latitude,BetaRICH_new,Cutmask,Beta_pre,BetaR);
}

void Grandezzequal (TTree *albero,int i,TNtuple *ntupla)
{
	albero->GetEvent(i);
	ntupla->Fill(Velocity,Rcutoff,R,NAnticluster,NTofClusters-NTofClustersusati,fabs(Rup-Rdown)/R,fuoriX,layernonusati,Chisquare,Richtotused,RichPhEl,Cutmask,Latitude,E_depTRD,IsCharge1);
}

void Grandezzesepd (TTree *albero,int i,TNtuple *ntupla)
{
	albero->GetEvent(i);
	ntupla->Fill(R,Beta,(*trtrack_edep)[0],Cutmask,Latitude,PhysBPatt,EdepTOFU,EdepTrack,EdepTOFD,Rcutoff,BetaRICH_new,LDiscriminant,Dist5D,Dist5D_P);
}


