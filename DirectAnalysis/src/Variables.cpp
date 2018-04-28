#include "Variables.hpp"
#include "histUtils.h"
#include "TRandom3.h"
#include <bitset>

int GetUnusedLayers(int hitbits){ return 7 - std::bitset<32>(hitbits & 0b1111111).count(); }
Reweighter ReweightInitializer(std::string galpropfilename="/afs/cern.ch/user/f/fdimicco/Work/Deutons/DirectAnalysis/include/CRDB_ProtonsAMS_R.galprop",float r_min=0.5, float r_max=100,float norm_at=1.05);


Variables::Variables(){
    BDTreader();
    reweighter = ReweightInitializer();
    reweighterHe = ReweightInitializer("/afs/cern.ch/user/f/fdimicco/Work/Deutons/DirectAnalysis/include/CRDB_HeliumAMS_R.galprop",2,2000,2.05);
}



void Variables::BDTreader()
{
	TMVA::Tools::Instance();
	readerTOF = new TMVA::Reader( "V:Color:!Silent" );

	TMVA::Tools::Instance();
	readerNaF = new TMVA::Reader( "V:Color:!Silent" );
	readerNaF->AddVariable("Richtotused",&Richtotused_float);	
    	readerNaF->AddVariable("RichPhEl",&RichPhEl);
    	readerNaF->AddVariable("RICHprob",&RICHprob);
    	readerNaF->AddVariable("RICHcollovertotal",&RICHcollovertotal);
     	readerNaF->AddVariable("RICHLipBetaConsistency",&RICHLipBetaConsistency);  
    	readerNaF->AddVariable("RICHTOFBetaConsistency",&RICHTOFBetaConsistency);  
    	readerNaF->AddVariable("RICHChargeConsistency",&RICHChargeConsistency);
	readerNaF->AddVariable("RICHPmts",&RICHPmts_float);
   	readerNaF->AddVariable("RICHgetExpected",&RICHgetExpected_float);
    	readerNaF->AddVariable("tot_hyp_p_uncorr",&tot_hyp_p_uncorr);
    	readerNaF->AddVariable("Bad_ClusteringRICH",&Bad_ClusteringRICH);
    	readerNaF->AddVariable("NSecondariesRICHrich",&NSecondariesRICHrich);
	readerNaF->BookMVA("BDTmethod", "/afs/cern.ch/user/f/fdimicco/Work/Deutons/DirectAnalysis/TMVA/Def/QualityNaF_BDT.weights.xml");

	TMVA::Tools::Instance();
	readerAgl = new TMVA::Reader( "V:Color:!Silent" );
    	readerAgl->AddVariable("Richtotused",&Richtotused_float);	
    	readerAgl->AddVariable("RichPhEl",&RichPhEl);
    	readerAgl->AddVariable("RICHprob",&RICHprob);
    	readerAgl->AddVariable("RICHcollovertotal",&RICHcollovertotal);
     	readerAgl->AddVariable("RICHLipBetaConsistency",&RICHLipBetaConsistency);  
    	readerAgl->AddVariable("RICHTOFBetaConsistency",&RICHTOFBetaConsistency);  
    	readerAgl->AddVariable("RICHChargeConsistency",&RICHChargeConsistency);
	readerAgl->AddVariable("RICHPmts",&RICHPmts_float);
   	readerAgl->AddVariable("RICHgetExpected",&RICHgetExpected_float);
    	readerAgl->AddVariable("tot_hyp_p_uncorr",&tot_hyp_p_uncorr);
    	readerAgl->AddVariable("Bad_ClusteringRICH",&Bad_ClusteringRICH);
    	readerAgl->AddVariable("NSecondariesRICHrich",&NSecondariesRICHrich);
	readerAgl->BookMVA("BDTmethod", "/afs/cern.ch/user/f/fdimicco/Work/Deutons/DirectAnalysis/TMVA/Def/QualityAgl_BDT.weights.xml");
}

void Variables::ResetVariables(){

    Run            = 0;
    Event          = 0;
    NEvent         = 0;
    U_time         = 0;
    NTracks        = 0;
    ThetaS         = 0;
    PhiS           = 0;
    Livetime       = 0;
    Latitude       = 0;
    PrescaleFactor = 0;

    // Cutoffs
    Rcutoff        = 0;
    IGRFRcutoff    = 0;

    // Bit fields
    JMembPatt      = 0;
    PhysBPatt      = 0;
    CUTMASK        = 0;
    RICHmask_new   = 0;

    // Counts
    NAnticluster_float = 0;
    NTofClusters       = 0;
    NTofClustersusati  = 0;
    NTofUsed           = 0;
    NTRDclusters       = 0;
    NAnticluster       = 0;
    NTRDSegments       = 0;
    NTrackHits         = 0;
    clustertrack       = 0;
    clustertottrack    = 0;


    // Track
    R_pre              = 0;
    Rup                = 0;
    Rdown              = 0;
    R                  = 0;
    R_L1               = 0;
    Chisquare          = 0;
    Chisquare_L1       = 0;
    Chisquare_y        = 0;
    Chisquare_L1_y     = 0;
    hitbits            = 0;
    FiducialVolume     = 0;	   
 
    // Tracker Charge
    qL1                = 0;
    qL1Status          = 0;
    qL2Status          = 0;
    qInner             = 0;
    qL2                = 0;	

    // TOF
    Beta               = 0;
    BetaR              = 0;
    Beta_pre           = 0;
    qUtof              = 0;
    qLtof              = 0;
    TOFchisq_s	=0;
    TOFchisq_t	=0;	
    NBadTOF     =0; 	

    // RICH 
    BetaRICH_new           = 0;
    Richtotused            = 0;
    Richtotused_float      = 0;
    RichPhEl               = 0;
    RICHprob               = 0;
    RICHPmts               = 0;
    RICHcollovertotal      = 0;
    RICHgetExpected        = 0;
    	
    RICHLipBetaConsistency =0;
    RICHTOFBetaConsistency =0;	
    RICHChargeConsistency  =0;	
    tot_hyp_p_uncorr	   =0;
    Bad_ClusteringRICH     =0;
    NSecondariesRICHrich   =0;	    

    //TRD
    TRDEdepovPath=0;
    TRDLikP=0;
    TRDLikD=0;
    TRDLike=0;
    EdepTRD=0;
	
    //MC vars
    Momento_gen            = 0;
    Massa_gen              = 0;
    mcweight               = 0;
    MCClusterGeantPids     = 0;
    Charge_gen             = 0;
    GenX  = 0; GenY  = 0; GenZ  = 0;
    GenPX = 0; GenPY = 0; GenPZ = 0;

    //Quality Discr.
    BDTDiscr = -1;
    Likelihood = -1;


}


void Variables::ReadBranches(TTree * tree){
	tree->SetBranchAddress("U_time"	 ,&U_time);
	tree->SetBranchAddress("NTracks"        ,&NTracks);
	tree->SetBranchAddress("Latitude"	 ,&Latitude);
	tree->SetBranchAddress("Rcutoff35"	 ,&IGRFRcutoff);
	tree->SetBranchAddress("StoermerRcutoff",&Rcutoff);
	tree->SetBranchAddress("Livetime"	 ,&Livetime);
	tree->SetBranchAddress("PhysBPatt"	 ,&PhysBPatt);
	tree->SetBranchAddress("R"		 ,&R);
	tree->SetBranchAddress("BetaR"	 ,&BetaR);
	tree->SetBranchAddress("CUTMASK"	 ,&CUTMASK);
	tree->SetBranchAddress("trtrack_edep"	 ,&trtrack_edep);
	tree->SetBranchAddress("trtot_edep"	 ,&trtot_edep);
	tree->SetBranchAddress("TOFEndep"	 ,&Endep);
	tree->SetBranchAddress("EdepTRD"	 ,&EdepTRD);
	tree->SetBranchAddress("BetaRICH_new"	 ,&BetaRICH_new);
	tree->SetBranchAddress("RICHmask"	 ,&RICHmask_new);
	tree->SetBranchAddress("EnergyECAL"	 ,&EdepECAL);
	tree->SetBranchAddress("NAnticluster"	 ,&NAnticluster);
	tree->SetBranchAddress("NTofClusters"	 ,&NTofClusters);
	tree->SetBranchAddress("NTofClustersusati",&NTofClustersusati);
	tree->SetBranchAddress("Rup"		 ,&Rup);
	tree->SetBranchAddress("Rdown"		 ,&Rdown);
	tree->SetBranchAddress("R"		 ,&R);
	tree->SetBranchAddress("Chisquare"	 ,&Chisquare);
	tree->SetBranchAddress("Beta"	 ,&Beta);
	tree->SetBranchAddress("NBadTOF"  ,&NBadTOF);
	tree->SetBranchAddress("NTrackHits"	 ,&NTrackHits );                  
	tree->SetBranchAddress("Richtotused"	 ,&Richtotused_float );  
	tree->SetBranchAddress("RichPhEl"	 ,&RichPhEl);  
	tree->SetBranchAddress("R_L1"		 ,&R_L1);  
	tree->SetBranchAddress("hitbits"	 ,&hitbits);  
	tree->SetBranchAddress("qL1"		 ,&qL1);  
	tree->SetBranchAddress("qL1Status"      ,&qL1Status);
	tree->SetBranchAddress("qL2Status"      ,&qL2Status);
	tree->SetBranchAddress("qL2"             ,&qL2);
	tree->SetBranchAddress("qInner"	 ,&qInner);  
	tree->SetBranchAddress("qUtof"		 ,&qUtof);  
	tree->SetBranchAddress("qLtof"		 ,&qLtof);  
	tree->SetBranchAddress("RICHprob",&RICHprob);
	tree->SetBranchAddress("RICHPmts",&RICHPmts_float);
	tree->SetBranchAddress("RICHcollovertotal",&RICHcollovertotal);
	tree->SetBranchAddress("RICHgetExpected",&RICHgetExpected_float);
	tree->SetBranchAddress("RICHLipBetaConsistency"		 ,&RICHLipBetaConsistency);  
	tree->SetBranchAddress("RICHTOFBetaConsistency"		 ,&RICHTOFBetaConsistency);  
	tree->SetBranchAddress("RICHChargeConsistency",&RICHChargeConsistency);
	tree->SetBranchAddress("tot_hyp_p_uncorr",&tot_hyp_p_uncorr);
	tree->SetBranchAddress("Bad_ClusteringRICH",&Bad_ClusteringRICH);
	tree->SetBranchAddress("NSecondariesRICHrich",&NSecondariesRICHrich);
	tree->SetBranchAddress("qL1InnerNoL2",&qL1InnerNoL2);
	tree->SetBranchAddress("TOFchisq_t",&TOFchisq_t);
	tree->SetBranchAddress("TOFchisq_s",&TOFchisq_s);

	tree->SetBranchAddress("TRDEdepovPath",&TRDEdepovPath);
	tree->SetBranchAddress("TRDLikP",&TRDLikP);
	tree->SetBranchAddress("TRDLikD",&TRDLikD);
	tree->SetBranchAddress("TRDLike",&TRDLike);
	tree->SetBranchAddress("FiducialVolume",&FiducialVolume);

	tree->SetBranchAddress("Momentum"	 ,&Momento_gen);  
	tree->SetBranchAddress("GenMass"	 ,&Massa_gen);  		
	tree->SetBranchAddress("MCClusterGeantPids",&MCClusterGeantPids); 	
	tree->SetBranchAddress("PrescaleFactor", &PrescaleFactor);
	
	// Derived Variables
	tree->SetBranchAddress("EdepTrack",&EdepTrack);
	tree->SetBranchAddress("EdepTOFU",&EdepTOFU);
	tree->SetBranchAddress("EdepTOFD",&EdepTOFD);

	tree->SetBranchAddress("joinCutmask",&joinCutmask);
	tree->SetBranchAddress("diffR",&diffR);
	tree->SetBranchAddress("TOF_Up_Down",&TOF_Up_Down);
	tree->SetBranchAddress("Layernonusati",&Layernonusati);
	tree->SetBranchAddress("NTofUsed",&NTofUsed);
	tree->SetBranchAddress("NAnticluster_float",&NAnticluster_float);
	tree->SetBranchAddress("RICHgetExpected_float",&RICHgetExpected_float);
        tree->SetBranchAddress("RICHPmts_float",&RICHPmts_float);

	tree->SetBranchAddress("BDTDiscr",&BDTDiscr);
	tree->SetBranchAddress("Likelihood",&Likelihood);
	tree->SetBranchAddress("U_time",&U_time);
	tree->SetBranchAddress("mcweight",&mcweight);

}

void Variables::RegisterTemplatesBranches(TTree * tree){
	tree->Branch("U_time"        ,&U_time);
	tree->Branch("StoermerRcutoff",&Rcutoff);
	tree->Branch("Livetime"        ,&Livetime);
	tree->Branch("R"               ,&R);
	tree->Branch("Beta"               ,&Beta);
	tree->Branch("BetaRICH_new"    ,&BetaRICH_new);
	tree->Branch("GenMass"         ,&Massa_gen);
        tree->Branch("U_time",&U_time);
	tree->Branch("Latitude",&Latitude);
        tree->Branch("mcweight",&mcweight);
	tree->Branch("joinCutmask",&joinCutmask);
	tree->Branch("qL1"		 ,&qL1);  
	tree->Branch("qL1Status"      ,&qL1Status);
	tree->Branch("qL2Status"      ,&qL2Status);
	tree->Branch("qL2"             ,&qL2);
	tree->Branch("qInner"	 	,&qInner);  
	tree->Branch("qUtof"		 ,&qUtof);  
	tree->Branch("qLtof"		 ,&qLtof);  
	tree->Branch("BDTDiscr",&BDTDiscr);
	tree->Branch("Likelihood",&Likelihood);
	tree->Branch("PrescaleFactor", &PrescaleFactor);
	tree->Branch("MCClusterGeantPids",&MCClusterGeantPids); 
}		

void Variables::RegisterBranches(TTree * tree){

	tree->Branch("U_time"	 ,&U_time);
	tree->Branch("NTracks"        ,&NTracks);
	tree->Branch("Latitude"	 ,&Latitude);
	tree->Branch("Rcutoff35"	 ,&IGRFRcutoff);
	tree->Branch("StoermerRcutoff",&Rcutoff);
	tree->Branch("Livetime"	 ,&Livetime);
	tree->Branch("PhysBPatt"	 ,&PhysBPatt);
	tree->Branch("R"		 ,&R);
	tree->Branch("BetaR"	 ,&BetaR);
	tree->Branch("CUTMASK"	 ,&CUTMASK);
	tree->Branch("trtrack_edep"	 ,&trtrack_edep);
	tree->Branch("trtot_edep"	 ,&trtot_edep);
	tree->Branch("TOFEndep"	 ,&Endep);
	tree->Branch("EdepTRD"	 ,&EdepTRD);
	tree->Branch("BetaRICH_new"	 ,&BetaRICH_new);
	tree->Branch("RICHmask"	 ,&RICHmask_new);
	tree->Branch("EnergyECAL"	 ,&EdepECAL);
	tree->Branch("NAnticluster"	 ,&NAnticluster);
	tree->Branch("NTofClusters"	 ,&NTofClusters);
	tree->Branch("NTofClustersusati",&NTofClustersusati);
	tree->Branch("Rup"		 ,&Rup);
	tree->Branch("Rdown"		 ,&Rdown);
	tree->Branch("R"		 ,&R);
	tree->Branch("Chisquare"	 ,&Chisquare);
	tree->Branch("Beta"	 ,&Beta);
	tree->Branch("NBadTOF"  ,&NBadTOF);
	tree->Branch("NTrackHits"	 ,&NTrackHits );                  
	tree->Branch("Richtotused"	 ,&Richtotused_float );  
	tree->Branch("RichPhEl"	 ,&RichPhEl);  
	tree->Branch("R_L1"		 ,&R_L1);  
	tree->Branch("hitbits"	 ,&hitbits);  
	tree->Branch("qL1"		 ,&qL1);  
	tree->Branch("qL1Status"      ,&qL1Status);
	tree->Branch("qL2Status"      ,&qL2Status);
	tree->Branch("qL2"             ,&qL2);
	tree->Branch("qL1InnerNoL2",&qL1InnerNoL2);
	tree->Branch("qInner"	 ,&qInner);  
	tree->Branch("qUtof"		 ,&qUtof);  
	tree->Branch("qLtof"		 ,&qLtof);  
	tree->Branch("RICHprob",&RICHprob);
	tree->Branch("RICHPmts",&RICHPmts_float);
	tree->Branch("RICHcollovertotal",&RICHcollovertotal);
	tree->Branch("RICHgetExpected",&RICHgetExpected_float);
	tree->Branch("RICHLipBetaConsistency"		 ,&RICHLipBetaConsistency);  
	tree->Branch("RICHTOFBetaConsistency"		 ,&RICHTOFBetaConsistency);  
	tree->Branch("RICHChargeConsistency",&RICHChargeConsistency);
	tree->Branch("tot_hyp_p_uncorr",&tot_hyp_p_uncorr);
	tree->Branch("Bad_ClusteringRICH",&Bad_ClusteringRICH);
	tree->Branch("NSecondariesRICHrich",&NSecondariesRICHrich);
	tree->Branch("PrescaleFactor", &PrescaleFactor);

	tree->Branch("TOFchisq_t",&TOFchisq_t);
	tree->Branch("TOFchisq_s",&TOFchisq_s);

	tree->Branch("TRDEdepovPath",&TRDEdepovPath);
	tree->Branch("TRDLikP",&TRDLikP);
	tree->Branch("TRDLikD",&TRDLikD);
	tree->Branch("TRDLike",&TRDLike);
	tree->Branch("FiducialVolume",&FiducialVolume);

	tree->Branch("Momentum"	 ,&Momento_gen);  
	tree->Branch("GenMass"	 ,&Massa_gen);  		
	tree->Branch("MCClusterGeantPids",&MCClusterGeantPids); 	

	// Derived Variables
	tree->Branch("EdepTrack",&EdepTrack);
	tree->Branch("EdepTOFU",&EdepTOFU);
	tree->Branch("EdepTOFD",&EdepTOFD);

	tree->Branch("joinCutmask",&joinCutmask);
	tree->Branch("diffR",&diffR);
	tree->Branch("TOF_Up_Down",&TOF_Up_Down);
	tree->Branch("Layernonusati",&Layernonusati);
	tree->Branch("NTofUsed",&NTofUsed);
	tree->Branch("NAnticluster_float",&NAnticluster_float);
	tree->Branch("RICHgetExpected_float",&RICHgetExpected_float);
        tree->Branch("RICHPmts_float",&RICHPmts_float);

	tree->Branch("BDTDiscr",&BDTDiscr);
	tree->Branch("Likelihood",&Likelihood);
	tree->Branch("U_time",&U_time);
	tree->Branch("mcweight",&mcweight);

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

    RICHmask_new=(RICHmask_new&2047);
    joinCutmask=CUTMASK;
    joinCutmask=CUTMASK|(1<<10);
    joinCutmask=(int)joinCutmask|((RICHmask_new)<<11);	

    diffR=fabs(Rup-Rdown)/R;
    TOF_Up_Down = fabs(EdepTOFD - EdepTOFU);
    Layernonusati =	GetUnusedLayers(hitbits);

    NTofUsed = NTofClusters - NTofClustersusati;
    NAnticluster_float = NAnticluster;
    Richtotused_float= Richtotused;
    RICHPmts_float = (float)RICHPmts;
    RICHgetExpected_float = (float) RICHgetExpected;	

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
    cout<<"EdepTRD:               "<<EdepTRD<<endl;
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

	Likelihood = 1;//readerTOF->EvaluateMVA("TMVA::Types::kLikelihood");

	if(IsFromNaF()) BDTDiscr = readerNaF->EvaluateMVA("BDTmethod");
	if(IsFromAgl()) BDTDiscr = readerAgl->EvaluateMVA("BDTmethod");
	
	return;
}


float GetInverseRigidity (Variables * vars) {return 1/vars->R-1/vars->Momento_gen;}
float GetGenMomentum     (Variables * vars) {return vars->Momento_gen;}

float GetInverseEdepUToF (Variables * vars) {return 1/vars->EdepTOFU;}
float GetInverseEdepLToF (Variables * vars) {return 1/vars->EdepTOFD;}
float GetInverseEdepTrack(Variables * vars) {return 1/vars->EdepTrack;}
float GetInverseEdepTRD  (Variables * vars) {return 1/vars->EdepTRD;}



float GetBetaGen         (Variables * vars) {return pow((pow((vars->Momento_gen/vars->Massa_gen),2)/(pow((vars->Momento_gen/vars->Massa_gen),2)+1)),0.5);}
float GetInverseBetaTOF  (Variables * vars) {return 1/vars->Beta - 1/GetBetaGen(vars);}
float GetInverseBetaRICH (Variables * vars) {return 1/vars->BetaRICH_new - 1/GetBetaGen(vars);}
float GetBetaTOF         (Variables * vars) {return vars->Beta;}
float GetBetaRICH        (Variables * vars) {return vars->BetaRICH_new;}
float GetRecMassTOF	 (Variables * vars) {return (vars->R/vars->Beta)*pow((1-pow(vars->Beta,2)),0.5);}
float GetRecMassRICH     (Variables * vars) {return (vars->R/vars->BetaRICH_new)*pow((1-pow(vars->BetaRICH_new,2)),0.5);}
float GetNegRecMassTOF	 (Variables * vars) {return (-vars->R/vars->Beta)*pow((1-pow(vars->Beta,2)),0.5);}
float GetNegRecMassRICH     (Variables * vars) {return (-vars->R/vars->BetaRICH_new)*pow((1-pow(vars->BetaRICH_new,2)),0.5);}


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
float GetL2Q (Variables * vars) {return vars->qL2;}
float GetL1InnerNoL2Q (Variables * vars) {return vars->qL1InnerNoL2;}
int   GetPIDatL1 (Variables * vars) {return static_cast<int> ((vars->MCClusterGeantPids)&255);}
int   GetPIDatL2 (Variables * vars) {return static_cast<int> ((vars->MCClusterGeantPids>>8)&255);}
int   GetPIDatL3 (Variables * vars) {return static_cast<int> ((vars->MCClusterGeantPids>>16)&255);}
float GetTRDePLikRatio    (Variables * vars) {return vars->TRDLikP/vars->TRDLike;}
float GetTRDDPLikRatio    (Variables * vars) {return vars->TRDLikD/vars->TRDLikP;}

float GetTOFSpatialChi (Variables * vars) {return vars->TOFchisq_s;}
float GetTOFTimeChi (Variables * vars) {return vars->TOFchisq_t;}


float GetTRDEdepovPath    (Variables * vars) {return vars->TRDEdepovPath;}

float GetLoweredBetaTOF  (Variables * vars) {	
	ResponseTOF->SetParameter(0,0.00347548);
	ResponseTOF->SetParameter(1,5.8474);
	return ResponseTOF->Eval(vars->Beta);				
}

float GetLoweredBetaNaF  (Variables * vars) {	
	ResponseNaF->SetParameter(0,-0.000859132);
	ResponseNaF->SetParameter(1,-30.5065);
	return ResponseNaF->Eval(vars->BetaRICH_new);	
}

float GetLoweredBetaAgl  (Variables * vars) {	
	ResponseAgl->SetParameter(0,4.28781e-05);
	ResponseAgl->SetParameter(1,67.8521);

	return ResponseAgl->Eval(vars->BetaRICH_new);				
}

float GetRICHBDT(Variables* vars) {	return vars->BDTDiscr;}
