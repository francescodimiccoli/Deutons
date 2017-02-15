#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "../include/binning.h"
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include <TFile.h>
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include <TMVA/Reader.h>
#include <TMVA/Tools.h>


#include "Commonglobals.cpp"
#include "Variables.cpp"
#include "Discriminants5D.h"

#include "reweight.h"
#include "histUtils.h"
#include "../include/binning.h"



bool ReadCalibration(string month);
bool ReadPdfForLikelihood();
int AssignMC_type(float Massa_gen);
float mcweight=0;
bool IsMC(TTree * tree);
Reweighter ReweightInitializer();
void CalibrateEdep(Variables *vars);

void ProcessEvent(Variables *vars,bool isMC,Reweighter reweighter);


int main(int argc, char * argv[])
{
	string month = argv[1];
	if(!ReadCalibration(month)) return 0;
	if(!ReadPdfForLikelihood()) return 0;

	cout<<"**************************** R BINS ***********************************"<<endl;
        Particle proton(0.9382720813, 1, 1);
        Binning PRB(proton);

        PRB.setBinsFromRigidity(43, 0.5, 100);
        PRB.Print();

        cout<<"**************************** BETA BINS TOF***********************************"<<endl;

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

	
	string INPUT(argv[2]);
	TFile *file =TFile::Open(INPUT.c_str());
	TTree *tree = (TTree *)file->Get("parametri_geo");
	bool check=tree?true:false;
	if(!check) {cout<<"Corrupted input file"<<endl;return 0;}
	bool isMC = IsMC(tree);


 	 Variables * vars = new Variables;
      	 vars->ReadBranches(tree);

	 TNtuple * grandezzesepd;
	 TNtuple * trig;
	 TNtuple * Q;

	 if(isMC) cout<<"It is a MC file"<<endl; 
	 if(isMC){
	 	grandezzesepd=new TNtuple("grandezzesepd","grandezzesepd","R:Beta:EdepL1:MC_type:Cutmask:PhysBPatt:EdepTOF:EdepTrack:EdepTOFD:Momentogen:BetaRICH_new:LDiscriminant:mcweight:Dist5D:Dist5D_P");
		trig=new TNtuple("trig","trig","MC_type:Momento_gen:Ev_Num:R_L1:R_pre:Beta_pre:Cutmask:EdepL1:EdepTOFU:EdepTOFD:EdepTrack:BetaRICH:EdepECAL:PhysBPatt:mcweight");
		Q = new TNtuple("Q","Q","R:Beta:MC_type:Cutmask:BetaRICH_new:Dist5D:Dist5D_P:LDiscriminant:Rmin:qL1:qInner:qUtof:qLtof:Momentogen");
	 }
	 else{
		grandezzesepd = new TNtuple("grandezzesepd","grandezzesepd","R:Beta:EdepL1:Cutmask:Latitude:PhysBPatt:EdepTOFU:EdepTrack:EdepTOFD:Rcutoff:BetaRICH_new:LDiscriminant:Dist5D:Dist5D_P");
		trig=new TNtuple("trig","trig","Seconds:Latitude:Rcutoff:R_L1:R_pre:Beta_pre:Cutmask:EdepL1:EdepTOFU:EdepTOFD:EdepTrack:BetaRICH:EdepECAL:PhysBPatt:Livetime");
		Q = new TNtuple("Q","Q","R:Beta:Cutmask:BetaRICH_new:Dist5D:Dist5D_P:LDiscriminant:Rmin:qL1:qInner:qUtof:qLtof"); 
	}

	Reweighter reweighter;
	if(isMC) reweighter=ReweightInitializer(); 	
	
	for(int i=0; i<100; i++){
		tree->GetEvent(i);
		vars->Update();
		ProcessEvent(vars,isMC,reweighter);	
	}	
	return 0;
}


void ProcessEvent(Variables *vars,bool isMC,Reweighter reweighter){
	vars->joinCutmask=vars->CUTMASK;
	vars->joinCutmask=vars->CUTMASK|(1<<10);
	vars->joinCutmask = vars->joinCutmask|(vars->RICHmask_new<<11);	
	if(isMC){
		vars->mcweight=reweighter.getWeight(fabs(vars->Momento_gen));
                        if(vars->Momento_gen<1) vars->mcweight=1;
			CalibrateEdep(vars);	
		
	}
	Likelihood(vars);	
	vars->PrintCurrentState();
	return;		

}




bool IsMC(TTree * tree){
	
	Variables * Vars = new Variables;
        Vars->ReadBranches(tree);
	tree->GetEvent(10);
	
	if (Vars->Massa_gen>0) return true;
	else return false;
}

bool ReadCalibration(string month){

	cout<<"****************** CALIB. READING **************************"<<endl;	
	string nomecal=("/storage/gpfs_ams/ams/users/fdimicco/Deutons/CodesforAnalysis/CALIBRAZIONI/"+month+".root");
        TFile *_file2 = TFile::Open(nomecal.c_str());
        if(_file2) cout<<"calibration file found: "<<(month + ".root").c_str()<<endl;
	else { cout<<"calibration file not found"<<endl; return false;}
	
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
	
	return true;

}


bool ReadPdfForLikelihood(){
	cout<<"*************** PDF for LIKELIHOOD READING ********************"<<endl;
	string nome;
	TFile *_file1 = TFile::Open("/storage/gpfs_ams/ams/users/fdimicco/Deutons/Ntuple-making/QualityVariables.root");
	TFile *_file3 = TFile::Open("/storage/gpfs_ams/ams/users/fdimicco/Deutons/Ntuple-making/QualityVariables_NaF.root");
	TFile *_file3b = TFile::Open("/storage/gpfs_ams/ams/users/fdimicco/Deutons/Ntuple-making/QualityVariables_Agl.root");
	if(!_file1||!_file3||!_file3b) {cout<<"PDF for likelihood files not found"<<endl; return false;}
	else cout<<"PDF for likelihood files found"<<endl;

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
        for(int qs=0; qs<9; qs++) cout<<SignalNaF[qs]<<" ";
        cout<<endl;
        for(int qs=0; qs<9; qs++) cout<<SignalAgl[qs]<<" ";
        cout<<endl;

	return true;	
}


Reweighter ReweightInitializer(){
	Histogram   mcFlux = makeLogUniform(500, 0.5, 100);
	Histogram dataFlux = loadGalpropFile("/storage/gpfs_ams/ams/users/fdimicco/Deutons/include/CRDB_ProtonsAMS_R.galprop");
	dataFlux.multiply( mcFlux.at(1.05) / dataFlux.getContent()[0] );
	Reweighter reweighter(mcFlux, dataFlux);
	return reweighter;
}


void CalibrateEdep(Variables *vars){
	
	float Velocity = vars->Beta;
	if(vars->BetaRICH_new>0) Velocity=vars->BetaRICH_new;
	Velocity=fabs(Velocity);	
	
	if(Velocity>0&&Velocity<1){
		vars->EdepTOFU =((vars->EdepTOFU)*Corr_TOFU->Eval(Velocity));
		vars->EdepTOFD =((vars->EdepTOFD)*Corr_TOFD->Eval(Velocity));
		vars->EdepTrack=((vars->EdepTrack)*Corr_Track->Eval(Velocity));
	}

	return;
}


int AssignMC_type(float Massa_gen)
{
        int MC_type=0;
        int cursor=0;
        if(Massa_gen<1) {
                if(Massa_gen<0.93805) MC_type = MC_type|(1<<(cursor+0));
                if(Massa_gen<0.939  ) MC_type = MC_type|(1<<(cursor+1));
                MC_type = MC_type|(1<<(cursor+2));
        }
        cursor=8;
        if(Massa_gen>1&&Massa_gen<2) {
                MC_type = MC_type|(1<<(cursor+0));
                int moffset=18570;
                int c_s_type=(int)(10000*Massa_gen-moffset + 2);
                MC_type = MC_type|(1<<(cursor+c_s_type));
        }
        cursor=16;
        if(Massa_gen>3) {
                MC_type = MC_type|(1<<(cursor+0));
                MC_type = MC_type|(1<<(cursor+2));
        }
        return  MC_type;
}

