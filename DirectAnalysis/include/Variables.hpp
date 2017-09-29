#ifndef VARIABLES_H
#define VARIABLES_H

#include "reweight.h"
#include "histUtils.h"
#include "TRandom3.h"

using namespace std;

int GetUnusedLayers(int hitbits){ return 7 - std::bitset<32>(hitbits & 0b1111111).count(); }

Reweighter ReweightInitializer(std::string galpropfilename="/storage/gpfs_ams/ams/users/fdimicco/Deutons/DirectAnalysis/include/CRDB_ProtonsAMS_R.galprop",float r_min=0.5, float r_max=100,float norm_at=1.05);

struct Variables{

	Reweighter reweighter;
	Reweighter reweighterHe;

	int 	   U_time;
        float 	   Latitude;
        float 	   Rcutoff;
	float 	   IGRFRcutoff;
        float      Livetime;
        int      PhysBPatt;
        float      R_pre;
        float      Beta_pre;
        int        CUTMASK;
        std::vector<float>  *       trtrack_edep=0;
        std::vector<float>  *       trtot_edep=0;
        std::vector<float>  *       Endep=0;
	float      BetaRICH_new;
	int        RICHmask_new;
	float      EdepECAL;
	int        NAnticluster;
	float	   NAnticluster_float;
	int        NTofClusters;
	int        NTofClustersusati;
	float	   NTofUsed=0;
	float      Rup;
	float      Rdown;
	float      R;
	float      Chisquare;
	float      Beta;
	float      BetaR;
	int        NTrackHits;                             
	int        Richtotused;
	float      Richtotused_float=0;
	float      RichPhEl;
	float      R_L1;
	int        hitbits;
	float	   qL1;
	float	   qInner;
	float	   qUtof;
	float	   qLtof;
	float 	   joinCutmask=0;
	float      diffR=0;
	float	   TOF_Up_Down=0;
	float	   Layernonusati=0;
	 
	//Discriminants
	float DistP=0;
	float DistD=0;
	float Likelihood=0;
	float BDTDiscr=0;

	//Summed E. deps.
	float EdepTOFU=0;
	float EdepTOFD=0;
	float EdepTrack=0;	
	float EdepL1=0;	
	
	//MC vars
	float	   Momento_gen;
	float	   Massa_gen;		
	float	   mcweight=0;


	Variables(){
		BDTreader();
		reweighter = ReweightInitializer();
		reweighterHe = ReweightInitializer("/storage/gpfs_ams/ams/users/fdimicco/Deutons/DirectAnalysis/include/CRDB_HeliumAMS_R.galprop",2,2000,2.05);
	}
	
	void ReadBranches(TTree * tree);
	void Update();
	void PrintCurrentState();
	void BDTreader();
	void Eval_Discriminants();
	bool IsFromNaF     (){ return (((int)joinCutmask>>11)==512&&BetaRICH_new>0);}
	bool IsFromAgl     (){ return (((int)joinCutmask>>11)==0&&BetaRICH_new>0);}
};

void Variables::BDTreader()
{
	TMVA::Tools::Instance();
	readerTOF = new TMVA::Reader( "V:Color:!Silent" );
	readerTOF->AddVariable("NAnticluster", &NAnticluster_float);
	readerTOF->AddVariable("Chisquare",&Chisquare);
	readerTOF->AddVariable("layernonusati := 7 - (hitbits&1) - (hitbits&2)*1/2 - (hitbits&4)*1/4 - (hitbits&8)*1/8 - (hitbits&16)*1/16 - (hitbits&32)*1/32 - (hitbits&64)*1/64",&Layernonusati);
	readerTOF->AddVariable("NTofUsed := NTofClusters - NTofClustersusati",&NTofUsed);
	readerTOF->AddVariable("diffR := TMath::Abs(Rup-Rdown)/R",&diffR);
	readerTOF->AddVariable("TOF_Up_Down := TMath::Abs(TOFEndep[2]+TOFEndep[3]-TOFEndep[0]-TOFEndep[1])", &TOF_Up_Down);
	readerTOF->BookMVA("BDTmethod", "/storage/gpfs_ams/ams/users/fdimicco/Deutons/DirectAnalysis/TMVA/weights/QualityTOF_BDT.weights.xml");
	readerTOF->BookMVA("TMVA::Types::kLikelihood", "/storage/gpfs_ams/ams/users/fdimicco/Deutons/DirectAnalysis/TMVA/weights/QualityTOF_Likelihood.weights.xml");

	TMVA::Tools::Instance();
	readerNaF = new TMVA::Reader( "V:Color:!Silent" );
	readerNaF->AddVariable("NAnticluster", &NAnticluster_float);
	readerNaF->AddVariable("Chisquare",&Chisquare);
	readerNaF->AddVariable("layernonusati := 7 - (hitbits&1) - (hitbits&2)*1/2 - (hitbits&4)*1/4 - (hitbits&8)*1/8 - (hitbits&16)*1/16 - (hitbits&32)*1/32 - (hitbits&64)*1/64",&Layernonusati);
	readerNaF->AddVariable("NTofUsed := NTofClusters - NTofClustersusati",&NTofUsed);
	readerNaF->AddVariable("diffR := TMath::Abs(Rup-Rdown)/R",&diffR);
	readerNaF->AddVariable("TOF_Up_Down := TMath::Abs(TOFEndep[2]+TOFEndep[3]-TOFEndep[0]-TOFEndep[1])", &TOF_Up_Down);
	readerNaF->AddVariable("Richtotused",&Richtotused_float);
	readerNaF->AddVariable("RichPhEl", &RichPhEl);
	readerNaF->BookMVA("BDTmethod", "/storage/gpfs_ams/ams/users/fdimicco/Deutons/DirectAnalysis/TMVA/weights/QualityNaF_BDT.weights.xml");
	readerNaF->BookMVA("TMVA::Types::kLikelihood", "/storage/gpfs_ams/ams/users/fdimicco/Deutons/DirectAnalysis/TMVA/weights/QualityNaF_Likelihood.weights.xml");

	TMVA::Tools::Instance();
	readerAgl = new TMVA::Reader( "V:Color:!Silent" );
	readerAgl->AddVariable("NAnticluster", &NAnticluster_float);
	readerAgl->AddVariable("Chisquare",&Chisquare);
	readerAgl->AddVariable("layernonusati := 7 - (hitbits&1) - (hitbits&2)*1/2 - (hitbits&4)*1/4 - (hitbits&8)*1/8 - (hitbits&16)*1/16 - (hitbits&32)*1/32 - (hitbits&64)*1/64",&Layernonusati);
	readerAgl->AddVariable("NTofUsed := NTofClusters - NTofClustersusati",&NTofUsed);
	readerAgl->AddVariable("diffR := TMath::Abs(Rup-Rdown)/R",&diffR);
	readerAgl->AddVariable("TOF_Up_Down := TMath::Abs(TOFEndep[2]+TOFEndep[3]-TOFEndep[0]-TOFEndep[1])", &TOF_Up_Down);
	readerAgl->AddVariable("Richtotused",&Richtotused_float);
	readerAgl->AddVariable("RichPhEl", &RichPhEl);
	readerAgl->BookMVA("BDTmethod", "/storage/gpfs_ams/ams/users/fdimicco/Deutons/DirectAnalysis/TMVA/weights/QualityNaF_BDT.weights.xml");
	readerAgl->BookMVA("TMVA::Types::kLikelihood", "/storage/gpfs_ams/ams/users/fdimicco/Deutons/DirectAnalysis/TMVA/weights/QualityNaF_Likelihood.weights.xml");

}


void Variables::ReadBranches(TTree * tree){

	 tree->SetBranchAddress("U_time"	 ,&U_time);
         tree->SetBranchAddress("Latitude"	 ,&Latitude);
         tree->SetBranchAddress("Rcutoff35"	 ,&IGRFRcutoff);
	 tree->SetBranchAddress("StoermerRcutoff",&Rcutoff);
         tree->SetBranchAddress("Livetime"	 ,&Livetime);
         tree->SetBranchAddress("PhysBPatt"	 ,&PhysBPatt);
         tree->SetBranchAddress("R"		 ,&R_pre);
         tree->SetBranchAddress("BetaRaw"	 ,&Beta_pre);
         tree->SetBranchAddress("CUTMASK"	 ,&CUTMASK);
         tree->SetBranchAddress("trtrack_edep"	 ,&trtrack_edep);
         tree->SetBranchAddress("trtot_edep"	 ,&trtot_edep);
         tree->SetBranchAddress("TOFEndep"	 ,&Endep);
	 tree->SetBranchAddress("BetaRICH"	 ,&BetaRICH_new);
	 tree->SetBranchAddress("RICHmask"	 ,&RICHmask_new);
	 tree->SetBranchAddress("EnergyECAL"	 ,&EdepECAL);
	 tree->SetBranchAddress("NAnticluster"	 ,&NAnticluster);
	 tree->SetBranchAddress("NTofClusters"	 ,&NTofClusters);
	 tree->SetBranchAddress("NTofClustersusati",&NTofClustersusati);
	 tree->SetBranchAddress("Rup"		 ,&Rup);
	 tree->SetBranchAddress("Rdown"		 ,&Rdown);
	 tree->SetBranchAddress("R"		 ,&R);
	 tree->SetBranchAddress("Chisquare"	 ,&Chisquare);
	 tree->SetBranchAddress("BetaHR"	 ,&Beta);
	 tree->SetBranchAddress("BetaOld"	 ,&BetaR);
	 tree->SetBranchAddress("NTrackHits"	 ,&NTrackHits );                  
	 tree->SetBranchAddress("Richtotused"	 ,&Richtotused );  
	 tree->SetBranchAddress("RichPhEl"	 ,&RichPhEl);  
	 tree->SetBranchAddress("R_L1"		 ,&R_L1);  
	 tree->SetBranchAddress("hitbits"	 ,&hitbits);  
	 tree->SetBranchAddress("qL1"		 ,&qL1);  
	 tree->SetBranchAddress("qInner"	 ,&qInner);  
	 tree->SetBranchAddress("qUtof"		 ,&qUtof);  
	 tree->SetBranchAddress("qLtof"		 ,&qLtof);  
	                                              
	 tree->SetBranchAddress("GenMomentum"	 ,&Momento_gen);  
	 tree->SetBranchAddress("GenMass"	 ,&Massa_gen);  		
}

void Variables::Update(){

	EdepTrack=0;
   	EdepTOFU=((*Endep)[0]+(*Endep)[1])/2;
   	EdepTOFD=((*Endep)[2]+(*Endep)[3])/2;
	for(int layer=1; layer<8; layer++) EdepTrack+=(*trtrack_edep)[layer];
	   EdepTrack=EdepTrack/7;
	DistP=100;
        DistD=100;
        joinCutmask=0;
	mcweight=0;
	
        RICHmask_new=(RICHmask_new&1023);
        joinCutmask=CUTMASK;
        joinCutmask=CUTMASK|(1<<10);
        joinCutmask=(int)joinCutmask|((RICHmask_new)<<11);	

	diffR=fabs(Rup-Rdown)/R;
	TOF_Up_Down = fabs(EdepTOFD - EdepTOFU);
	Layernonusati =	GetUnusedLayers(hitbits);

	NTofUsed = NTofClusters - NTofClustersusati;
        NAnticluster_float = NAnticluster;
	Richtotused_float= Richtotused;

	Eval_Discriminants();

	U_time-=1305200000; //time offset (in order to have small time stamp)	

	if(Momento_gen<1) mcweight=1; // in data Momento_gen=0
	else {
		if(Massa_gen<1) mcweight = reweighter.getWeight(fabs(Momento_gen));
		else if(Massa_gen>1&&Massa_gen<2) mcweight = 0.05*reweighter.getWeight(fabs(Momento_gen));
		else mcweight = reweighterHe.getWeight(fabs(Momento_gen));
	}	
}

void Variables::PrintCurrentState(){

	cout<<endl;
	cout<<endl;
	cout<<"***Current Values of Variables:***"<<endl;


	cout<<"U_time:                "<<U_time<<endl;
        cout<<"Latitude:              "<<Latitude<<endl;
        cout<<"Rcutoff:               "<<Rcutoff<<endl;
        cout<<"Livetime:              "<<Livetime<<endl;
        cout<<"PhysBPatt:             "<<PhysBPatt<<endl;
        cout<<"R_pre:                 "<<R_pre<<endl;
        cout<<"Beta_pre:              "<<Beta_pre<<endl;
        cout<<"CUTMASK:               "<<CUTMASK<<endl;
	cout<<"BetaRICH_new:          "<<BetaRICH_new<<endl;
	cout<<"RICHmask_new:          "<<RICHmask_new<<endl;
	cout<<"EdepECAL:              "<<EdepECAL<<endl;
	cout<<"NAnticluster:          "<<NAnticluster<<endl;
	cout<<"NTofClusters:          "<<NTofClusters<<endl;
	cout<<"NTofClustersusati:     "<<NTofClustersusati<<endl;
	cout<<"Rup:                   "<<Rup<<endl;
	cout<<"Rdown:                 "<<Rdown<<endl;
	cout<<"R:                     "<<R<<endl;
	cout<<"Chisquare:             "<<Chisquare<<endl;
	cout<<"Beta:                  "<<Beta<<endl;
	cout<<"BetaR:                 "<<BetaR<<endl;
	cout<<"NTrackHits:            "<<NTrackHits          <<endl;                   
	cout<<"Richtotused:           "<<Richtotused<<endl;
	cout<<"RichPhEl:              "<<RichPhEl<<endl;
	cout<<"R_L1:                  "<<R_L1<<endl;
	cout<<"hitbits:               "<<hitbits<<endl;
	cout<<"qL1:                   "<<qL1<<endl;
	cout<<"qInner:                "<<qInner<<endl;
	cout<<"qUtof:                 "<<qUtof<<endl;
	cout<<"qLtof:                 "<<qLtof<<endl;
	cout<<"Momento_gen:           "<<Momento_gen<<endl;
	cout<<"Massa_gen:	      "<<Massa_gen<<endl;		
	
	cout<<"EdepTOFU:	     "<<EdepTOFU<<endl;
	cout<<"EdepTOFD:             "<<EdepTOFD<<endl; 
	cout<<"EdepTrack:            "<<EdepTrack<<endl;
	
	cout<<"diffR:		      "<<diffR<<endl;
	cout<<"Layernonusati:         "<<Layernonusati<<endl;
	cout<<"TOF_Up_Down		"<<TOF_Up_Down<<endl;
	cout<<"NTofUsed			"<<NTofUsed<<endl;
	cout<<"DistP:		      "<<DistP<<endl;
	cout<<"DistD:                 "<<DistD<<endl; 
	cout<<"Likelihood:            "<<Likelihood<<endl;
	cout<<"BDTDiscr:            "<<BDTDiscr<<endl;
	
	cout<<"mcweight:            "<<mcweight<<endl;
	cout<<"joinCutmask:           "<<joinCutmask<<endl;




	cout<<"******"<<endl;
	cout<<endl;
        cout<<endl;

}

void Variables::Eval_Discriminants(){

	BDTDiscr = readerTOF->EvaluateMVA("BDTmethod");
	Likelihood = readerTOF->EvaluateMVA("TMVA::Types::kLikelihood");
	if(IsFromNaF()) {
		BDTDiscr = readerNaF->EvaluateMVA("BDTmethod");
		Likelihood = readerNaF->EvaluateMVA("TMVA::Types::kLikelihood");
	}
	if(IsFromAgl()) {
		BDTDiscr = readerAgl->EvaluateMVA("BDTmethod");
		Likelihood = readerAgl->EvaluateMVA("TMVA::Types::kLikelihood");
	}
	return;
}


float GetInverseRigidity (Variables * vars) {return 1/vars->R-1/vars->Momento_gen;}
float GetGenMomentum     (Variables * vars) {return vars->Momento_gen;}

float GetInverseEdepUToF (Variables * vars) {return 1/vars->EdepTOFU;}
float GetInverseEdepLToF (Variables * vars) {return 1/vars->EdepTOFD;}
float GetInverseEdepTrack(Variables * vars) {return 1/vars->EdepTrack;}
float GetBetaGen         (Variables * vars) {return pow((pow((vars->Momento_gen/vars->Massa_gen),2)/(pow((vars->Momento_gen/vars->Massa_gen),2)+1)),0.5);}
float GetInverseBetaTOF  (Variables * vars) {return 1/vars->Beta - 1/GetBetaGen(vars);}
float GetInverseBetaRICH (Variables * vars) {return 1/vars->BetaRICH_new - 1/GetBetaGen(vars);}
float GetBetaTOF         (Variables * vars) {return vars->Beta;}
float GetBetaRICH        (Variables * vars) {return vars->BetaRICH_new;}
float GetRecMassTOF	 (Variables * vars) {return (vars->R/vars->Beta)*pow((1-pow(vars->Beta,2)),0.5);}
float GetRecMassRICH     (Variables * vars) {return (vars->R/vars->BetaRICH_new)*pow((1-pow(vars->BetaRICH_new,2)),0.5);}

float GetRigidity (Variables * vars) {return vars->R;}


Reweighter ReweightInitializer(std::string galpropfilename, float r_min, float r_max, float norm_at){
        Histogram   mcFlux = makeLogUniform(500, r_min, r_max);
        Histogram dataFlux = loadGalpropFile(galpropfilename.c_str());
        dataFlux.multiply( mcFlux.at(norm_at) / dataFlux.getContent()[0] );
        Reweighter reweighter(mcFlux, dataFlux);
        return reweighter;
}

float GetUtofQ	(Variables * vars) {return vars->qUtof; }
float GetLtofQ	(Variables * vars) {return vars->qLtof; }
float GetInnerQ	(Variables * vars) {return vars->qInner;}
float GetL1Q (Variables * vars) {return vars->qL1;}
#endif
