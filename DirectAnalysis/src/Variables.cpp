#include "Variables.hpp"
#include "histUtils.h"
#include "TRandom3.h"
#include <bitset>
#include "GlobalPaths.h"


Double_t Variables::ChiXcut_X[] = {
			    0.5461929639082259,
                            0.6117777913465905,
                            0.7071657262317889,
                            0.782162253849283,
                            0.8542806062673918,
                            0.8705795905154753,
                            0.9041164064765199,
                            0.9330485467088381,
                            0.993719979859955,
                            1.058336569791599,
                            1.1705755610785775,
                            1.2624994283643205,
                            1.3702471441944273,
                            1.6141144911444367,
                            2.0250210803940845,
                            2.5565879564232303,
                            3.011592241466872,
                            3.6152599901105105,
                            4.593108484692199,
                            5.798796836339215,
                            7.555252339430445,
                            9.905947748010938,
                            13.65937777980889,
                            19.19436436007669,
                            23.187454671790675,
                            27.660527949179535,
                            31.374732184749277,
                            38.484163777366675,
                            60.41179653368272,
				120} ;

Double_t Variables::ChiXcut_Y[] = {15.712800278053805,
                            16.510284425867916,
                            17.695234763854792,
                            18.675191718987424,
                            19.678012360917403,
                            19.88313051150314,
                            19.518122513000353,
                            19.10754310056725,
                            18.263582972769445,
                            17.419622844971638,
                            15.95984038642349,
                            15.024612663033551,
                            14.36296938288989,
                            13.108086731397353,
                            11.647680434191704,
                            10.43805727730074,
                            9.479591563918955,
                            8.361423154217372,
                            7.584993561093857,
                            7.241819915603543,
                            6.9896955221752375,
                            6.7831425426772505,
                            6.507936118921474,
                            6.392245239888245,
                            6.3685081789704,
                            6.3448335019183055,
                            6.344209663260806,
                            6.366543087199277,
                            6.3865683081050015,
			6.4	};

Double_t Variables::ChiYcut_X[] = {  
	0.554988920241520,
	0.665321398928230,
	0.766094503809641,
	0.873268284254038,
	0.908899716158968,
	0.974650079710933,
	1.045112844951578,
	1.189765779872607,
	1.3019451667854716,
	1.3966580111639908,
	1.6236359122399298,
	1.9460478452088468,
	2.404375851115271,
	2.659348043437154,
	2.9710860026813823,
	3.2197752528619903,
	3.3855176480162017,
	3.6322634550753348,
	4.0169394677244465,
	4.625074192712959,
	5.43426215771492,
	5.9501497260948515,
	6.580995144088139,
	7.654857043964374,
	8.552185815608533,
	9.65086191304605,
	10.67383860718517,
	11.338151440439248,
	12.289778294130334,
	13.320434925356826,
	15.968890034503692,
	10.337539385333123,
	24.381694092636632,
	28.937574569441463,
	32.98652862434818,
	40.75800294332464,
	120
};

Double_t Variables::ChiYcut_Y[] = {
	14.975966562173461,  
	14.888192267502614,
	14.829676071055385,
	14.741901776384534,
	14.303030303030305,
	13.337513061650995,
	12.313479623824453,
	10.411703239289448,
	9.592476489028215, 
	9.153605015673982, 
	8.363636363636363, 
	8.012539184952981, 
	7.836990595611286, 
	7.866248693834902, 
	7.866248693834902, 
	7.544409613375131, 
	7.281086729362593, 
	7.017763845350053, 
	6.871473354231974, 
	6.725182863113898, 
	6.725182863113898, 
	6.725182863113898, 
	6.725182863113898, 
	6.725182863113898, 
	6.725182863113898, 
	6.637408568443052, 
	6.608150470219436, 
	6.491118077324974, 
	6.46185997910136,
	6.344827586206896, 
	6.286311389759667, 
	6.25705329153605,
	6.227795193312436, 
	6.227795193312436, 
	6.169278996865205, 
	6.081504702194359, 
	6};



int GetUnusedLayers(int hitbits){ return 7 - std::bitset<32>(hitbits & 0b1111111).count(); }
Reweighter ReweightInitializer(std::string galpropfilename=(workdir+"/include/CRDB_ProtonsAMS_R.galprop").c_str(),float r_min=0.5, float r_max=100,float norm_at=1.05);


Variables::Variables(){
    BDTreader();
    reweighter = ReweightInitializer();
    reweighterHe = ReweightInitializer((workdir+"include/CRDB_HeliumAMS_R.galprop").c_str(),2,2000,2.05);
   	
    Chi2Xcut = new TSpline3("Chi2Xcut", ChiXcut_X,ChiXcut_Y,30);
    Chi2Ycut = new TSpline3("Chi2Ycut", ChiYcut_X,ChiYcut_Y,37);


}


void Variables::BDTreader()
{
	TMVA::Tools::Instance();
	readerTOF = new TMVA::Reader( "V:Color:!Silent" );

	TMVA::Tools::Instance();
	readerNaF = new TMVA::Reader( "V:Color:!Silent" );
	readerNaF->AddVariable("Richtotused",&Richtotused);	
	readerNaF->AddVariable("RichPhEl",&RichPhEl);
    	readerNaF->AddVariable("RICHprob",&RICHprob);
	readerNaF->AddVariable("RICHcollovertotal",&RICHcollovertotal);
     	readerNaF->AddVariable("RICHLipBetaConsistency",&RICHLipBetaConsistency);  
	readerNaF->AddVariable("RICHTOFBetaConsistency",&RICHTOFBetaConsistency);  
    	readerNaF->AddVariable("RICHChargeConsistency",&RICHChargeConsistency);
	readerNaF->AddVariable("RICHPmts",&RICHPmts);
   	readerNaF->AddVariable("RICHgetExpected",&RICHgetExpected);
    	readerNaF->AddVariable("tot_hyp_p_uncorr",&tot_hyp_p_uncorr);
    	readerNaF->AddVariable("Bad_ClusteringRICH",&Bad_ClusteringRICH);
    	readerNaF->AddVariable("NSecondariesRICHrich",&NSecondariesRICHrich);
	//readerNaF->AddVariable("HitHValldir"  ,&HitHValldir);
    	//readerNaF->AddVariable("HitHVallrefl" ,&HitHVallrefl);
    	//readerNaF->AddVariable("HVBranchCheck:= (HitHValldir - HitHVoutdir) - (HitHVallrefl - HitHVoutrefl)"  ,&HVBranchCheck);
	readerNaF->AddVariable("HitHVoutdir"  ,&HitHVoutdir);
    	readerNaF->AddVariable("HitHVoutrefl" ,&HitHVoutrefl);
	readerNaF->AddSpectator("R", &R);
    	readerNaF->AddSpectator("BetaRICH_new", &BetaRICH_new);

	readerNaF->BookMVA("BDTmethod", (workdir+"/TMVA/New3/QualityNaF_BDT.weights.xml").c_str());


	TMVA::Tools::Instance();
	readerAgl = new TMVA::Reader( "V:Color:!Silent" );
    	readerAgl->AddVariable("Richtotused",&Richtotused);	
    	readerAgl->AddVariable("RichPhEl",&RichPhEl);
    	readerAgl->AddVariable("RICHprob",&RICHprob);
    	readerAgl->AddVariable("RICHcollovertotal",&RICHcollovertotal);
     	readerAgl->AddVariable("RICHLipBetaConsistency",&RICHLipBetaConsistency);  
    	readerAgl->AddVariable("RICHTOFBetaConsistency",&RICHTOFBetaConsistency);  
    	readerAgl->AddVariable("RICHChargeConsistency",&RICHChargeConsistency);
	readerAgl->AddVariable("RICHPmts",&RICHPmts);
   	readerAgl->AddVariable("RICHgetExpected",&RICHgetExpected);
    	readerAgl->AddVariable("tot_hyp_p_uncorr",&tot_hyp_p_uncorr);
    	readerAgl->AddVariable("Bad_ClusteringRICH",&Bad_ClusteringRICH);
    	readerAgl->AddVariable("NSecondariesRICHrich",&NSecondariesRICHrich);
	//readerAgl->AddVariable("HitHValldir"  ,&HitHValldir);
    	//readerAgl->AddVariable("HitHVallrefl" ,&HitHVallrefl);
    //	readerAgl->AddVariable("HVBranchCheck:= (HitHValldir - HitHVoutdir) - (HitHVallrefl - HitHVoutrefl)"  ,&HVBranchCheck);
	readerAgl->AddVariable("HitHVoutdir"  ,&HitHVoutdir);
    	readerAgl->AddVariable("HitHVoutrefl" ,&HitHVoutrefl);
	readerAgl->AddSpectator("R", &R);
    	readerAgl->AddSpectator("BetaRICH_new", &BetaRICH_new);

	readerAgl->BookMVA("BDTmethod", (workdir+"/TMVA/New3/QualityAgl_BDT.weights.xml").c_str());
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

    //RTI
    good_RTI=0;
    Rcutoff_RTI=0;
    Rcutoff_IGRFRTI =0;
    isinsaa=0;
    Livetime_RTI       = 0;
		
    //Tracking Efficiency
    theta_track=0;
    phi_track=0;
    entrypointcoo[0]=0;	
    entrypointcoo[1]=0;	
    entrypointcoo[2]=0;	
    beta_SA=0;
    betapatt_SA=0;	
    beta_ncl_SA=0;
    beta_chiT_SA=0;
    Trd_chi_SA=0;	
    qUtof_SA=0;
    qLtof_SA=0;			
    qTrd_SA=0;	
    EdepECAL=0;
    trd_int_inside_tracker=0;

    //L1 pick-up Efficiency
    exthit_closest_q	=0;     
    exthit_closest_status	=0;
    hitdistfromint  =0;


    // Bit fields
    JMembPatt      = 0;
    PhysBPatt      = 0;
    CUTMASK        = 0;
    RICHmask_new   = 0;
    P_standard_sel = 0;

    // Counts
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
    patty	       =0;
    FiducialVolume     = 0;	   
    R_sec	       = 0;

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

    HitHValldir	    =0;
    HitHVallrefl    =0;
    HitHVoutdir     =0;
    HitHVoutrefl    =0;
    HVBranchCheck   =0;

    //TRD
    TRDEdepovPath=0;
    TRDLikP=0;
    TRDLikD=0;
    TRDLike=0;
    EdepTRD=0;
	
    //MC vars
    Momento_gen            = 0;
    Momento_gen_UTOF       = 0;
    Momento_gen_LTOF       = 0;
    Momento_gen_RICH       = 0;
    Massa_gen              = 0;
    mcweight               = 0;
    MCClusterGeantPids     = 0;
    Charge_gen             = 0;
    GenX  = 0; GenY  = 0; GenZ  = 0;
    GenPX = 0; GenPY = 0; GenPZ = 0;

    //Quality Discr.
    BDTDiscr = -1;
    Likelihood = -1;

    //checks on variables
    beta_ncl =  0;
    chisqcn  = -1; 
    chisqtn  = -1; 
    nTrTracks= -1;
    sumclsn  = -1; 




}

void Variables::ReadBranches(TTree * tree){
	tree->SetBranchAddress("U_time"	 ,&U_time);
	tree->SetBranchAddress("NTracks"        ,&NTracks);
	tree->SetBranchAddress("Latitude"	 ,&Latitude);
	tree->SetBranchAddress("Rcutoff35"	 ,&IGRFRcutoff);
	tree->SetBranchAddress("StoermerRcutoff",&Rcutoff);
	tree->SetBranchAddress("Livetime"	 ,&Livetime);
	tree->SetBranchAddress("PhysBPatt"	 ,&PhysBPatt);
	tree->SetBranchAddress("Livetime_RTI"        ,&Livetime_RTI);
	tree->SetBranchAddress("Rcutoff_RTI",&Rcutoff_RTI);
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
	tree->SetBranchAddress("Chisquare_y"       ,&Chisquare_y);
	tree->SetBranchAddress("Beta"	 ,&Beta);
	tree->SetBranchAddress("NBadTOF"  ,&NBadTOF);
	tree->SetBranchAddress("NTrackHits"	 ,&NTrackHits );                  
	tree->SetBranchAddress("Richtotused"	 ,&Richtotused );  
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
	tree->SetBranchAddress("RICHPmts",&RICHPmts);
	tree->SetBranchAddress("RICHcollovertotal",&RICHcollovertotal);
	tree->SetBranchAddress("RICHgetExpected",&RICHgetExpected);
	tree->SetBranchAddress("RICHLipBetaConsistency"		 ,&RICHLipBetaConsistency);  
	tree->SetBranchAddress("RICHTOFBetaConsistency"		 ,&RICHTOFBetaConsistency);  
	tree->SetBranchAddress("RICHChargeConsistency",&RICHChargeConsistency);
	tree->SetBranchAddress("tot_hyp_p_uncorr",&tot_hyp_p_uncorr);
	tree->SetBranchAddress("HitHValldir" ,&HitHValldir); 
	tree->SetBranchAddress("HitHVallrefl",&HitHVallrefl);
	tree->SetBranchAddress("HitHVoutdir",&HitHVoutdir);
	tree->SetBranchAddress("HitHVoutrefl",&HitHVoutrefl);
	tree->SetBranchAddress("HVBranchCheck",&HVBranchCheck);

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
	tree->SetBranchAddress("NAnticluster",&NAnticluster);
	tree->SetBranchAddress("RICHgetExpected",&RICHgetExpected);
        tree->SetBranchAddress("RICHPmts",&RICHPmts);

	tree->SetBranchAddress("BDTDiscr",&BDTDiscr);
	tree->SetBranchAddress("Likelihood",&Likelihood);
	tree->SetBranchAddress("U_time",&U_time);
	tree->SetBranchAddress("mcweight",&mcweight);

}

void Variables::RegisterTemplatesBranches(TTree * tree){
	tree->Branch("U_time"        ,&U_time);
	tree->Branch("Livetime_RTI"        ,&Livetime_RTI);
	tree->Branch("Rcutoff_RTI",&Rcutoff_RTI);
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
	tree->Branch("Livetime_RTI"        ,&Livetime_RTI);
	tree->Branch("Rcutoff_RTI",&Rcutoff_RTI);
	
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
	tree->Branch("Chisquare_y"       ,&Chisquare_y);
	tree->Branch("Beta"	 ,&Beta);
	tree->Branch("NBadTOF"  ,&NBadTOF);
	tree->Branch("NTrackHits"	 ,&NTrackHits );                  
	tree->Branch("Richtotused"	 ,&Richtotused );  
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
	tree->Branch("RICHPmts",&RICHPmts);
	tree->Branch("RICHcollovertotal",&RICHcollovertotal);
	tree->Branch("RICHgetExpected",&RICHgetExpected);
	tree->Branch("RICHLipBetaConsistency"		 ,&RICHLipBetaConsistency);  
	tree->Branch("RICHTOFBetaConsistency"		 ,&RICHTOFBetaConsistency);  
	tree->Branch("RICHChargeConsistency",&RICHChargeConsistency);
	tree->Branch("tot_hyp_p_uncorr",&tot_hyp_p_uncorr);
	tree->Branch("HitHValldir" ,&HitHValldir); 
	tree->Branch("HitHVallrefl",&HitHVallrefl);
	tree->Branch("HitHVoutdir",&HitHVoutdir);
	tree->Branch("HitHVoutrefl",&HitHVoutrefl);
	tree->Branch("HVBranchCheck",&HVBranchCheck);
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
	tree->Branch("NAnticluster",&NAnticluster);
	tree->Branch("RICHgetExpected",&RICHgetExpected);
        tree->Branch("RICHPmts",&RICHPmts);

	tree->Branch("BDTDiscr",&BDTDiscr);
	tree->Branch("Likelihood",&Likelihood);
	tree->Branch("U_time",&U_time);
	tree->Branch("mcweight",&mcweight);

}


void Variables::Update(){

    EdepTrack=0;
	
    //EdepTOFU=((*Endep)[0]+(*Endep)[1])/2;
    //EdepTOFD=((*Endep)[2]+(*Endep)[3])/2;
    //for(int layer=1; layer<8; layer++) EdepTrack+=(*trtrack_edep)[layer];
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
    HVBranchCheck = (HitHValldir - HitHVoutdir) - (HitHVallrefl - HitHVoutrefl);	
    Eval_Discriminants();

    U_time-=1305200000; //time offset (in order to have small time stamp)	

    if(Momento_gen==0) mcweight=1; // in data Momento_gen=0
    else {
			if(isinf(reweighter.getWeight(fabs(Momento_gen)))) mcweight = 0;
			else {
				mcweight = reweighter.getWeight(fabs(Momento_gen));
				if(Massa_gen>1&&Massa_gen<2) mcweight *= 0.05;
				if(Massa_gen>2&&Massa_gen<4) mcweight = reweighterHe.getWeight(fabs(Momento_gen));
			}
	}	

}

void Variables::PrintCurrentState(){

    cout<<endl;
    cout<<endl;
    cout<<"***Current Values of Variables:***"<<endl;


	cout<<"Massa_gen		"<<Massa_gen<<endl;
	cout<<"Momento_gen		"<<Momento_gen<<endl;
	cout<<"Charge_gen		"<<Charge_gen<<endl;
	cout<<"mcweight			"<<mcweight<<endl;
	cout<<"PrescaleFactor		"<<PrescaleFactor	    <<endl;
	cout<<"P_standard_sel           " <<P_standard_sel            <<endl;
	cout<<"Rcutoff                  " <<Rcutoff               <<endl;
        cout<<"Rcutoff_RTI              " <<Rcutoff_RTI           <<endl;
	cout<<"Rcutoff_IGRFRTI          " <<Rcutoff_IGRFRTI       <<endl;
	cout<<"Livetime_RTI             " <<Livetime_RTI          <<endl;
	cout<<"good_RTI                 " <<good_RTI              <<endl;
	cout<<"isinsaa                  " <<isinsaa               <<endl;
	cout<<"NTRDSegments             " <<NTRDSegments          <<endl;
	cout<<"NTRDclusters             " <<NTRDclusters          <<endl;
	cout<<"NTrackHits               " <<NTrackHits            <<endl;
	cout<<"NTracks                  " <<NTracks               <<endl;
	cout<<"NAnticluster             " <<NAnticluster          <<endl;
	cout<<"PhysBPatt                " <<PhysBPatt             <<endl;
	cout<<"JMembPatt                " <<JMembPatt             <<endl;
	cout<<"CUTMASK                  " <<CUTMASK               <<endl;
	cout<<"beta_SA                  " <<beta_SA               <<endl;
	cout<<"qUtof_SA                 " <<qUtof_SA              <<endl;
	cout<<"qLtof_SA                 " <<qLtof_SA              <<endl;
	cout<<"EdepECAL                 " <<EdepECAL              <<endl;
	cout<<"exthit_closest_q         " <<exthit_closest_q      <<endl;
	cout<<"hitdistfromint           " <<hitdistfromint        <<endl;
	cout<<"R                        " <<R                     <<endl;
	cout<<"Chisquare                " <<Chisquare             <<endl;
	cout<<"Chisquare_y              " <<Chisquare_y           <<endl;
	cout<<"hitbits                  " <<hitbits               <<endl;
	cout<<"FiducialVolume           " <<FiducialVolume        <<endl;
	cout<<"qL1                      " <<qL1                   <<endl;
	cout<<"qL1Status                " <<qL1Status             <<endl;
	cout<<"qL2                      " <<qL2                   <<endl;
	cout<<"qL2Status                " <<qL2Status             <<endl;
	cout<<"qInner                   " <<qInner                <<endl;
	cout<<"qUtof                    " <<qUtof                 <<endl;
	cout<<"qLtof                    " <<qLtof                 <<endl;
	cout<<"Beta                     " <<Beta                  <<endl;
	cout<<"TOFchisq_s               " <<TOFchisq_s            <<endl;
	cout<<"TOFchisq_t               " <<TOFchisq_t            <<endl;
	cout<<"BetaRICH_new             " <<BetaRICH_new          <<endl;
	cout<<"RICHmask_new             " <<RICHmask_new          <<endl;
	cout<<"Richtotused              " <<Richtotused           <<endl;
	cout<<"RichPhEl                 " <<RichPhEl              <<endl;
	cout<<"RICHprob                 " <<RICHprob              <<endl;
	cout<<"RICHgetExpected          " <<RICHgetExpected      	<<endl;  
	cout<<"RICHTOFBetaConsistency   " <<RICHTOFBetaConsistency<<endl;
	cout<<"BDTDiscr                 " <<BDTDiscr              <<endl;
	      
              
             
    cout<<"******"<<endl;
    cout<<endl;
    cout<<endl;


}

void Variables::Eval_Discriminants(){

	Likelihood = 1;//readerTOF->EvaluateMVA("TMVA::Types::kLikelihood");

	if(IsFromNaF()&&BDTDiscr==-1) BDTDiscr = readerNaF->EvaluateMVA("BDTmethod");
	if(IsFromAgl()&&BDTDiscr==-1) BDTDiscr = readerAgl->EvaluateMVA("BDTmethod");
	
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
float GetBetaMeas        (Variables * vars) { if(!(vars->IsFromNaF()) && !(vars->IsFromAgl())) return GetBetaTOF(vars); else return GetBetaRICH(vars); }
float GetBetaTOF         (Variables * vars) {return vars->Beta;}
float GetBetaRICH        (Variables * vars) {return vars->BetaRICH_new;}
float GetRecMassTOF	 (Variables * vars) {return (vars->R/vars->Beta)*pow((1-pow(vars->Beta,2)),0.5);}
float GetRecMassRICH     (Variables * vars) {return (vars->R/vars->BetaRICH_new)*pow((1-pow(vars->BetaRICH_new,2)),0.5);}
float GetNegRecMassTOF	 (Variables * vars) {return (-vars->R/vars->Beta)*pow((1-pow(vars->Beta,2)),0.5);}
float GetNegRecMassRICH     (Variables * vars) {return (-vars->R/vars->BetaRICH_new)*pow((1-pow(vars->BetaRICH_new,2)),0.5);}


float GetRigidity (Variables * vars) {return vars->R;}
float GetRigiditySecondTrack (Variables * vars) {return vars->R_sec;}



Reweighter ReweightInitializer(std::string galpropfilename, float r_min, float r_max, float norm_at){
	Histogram   mcFlux = makeLogUniform(500, r_min, r_max);
        Histogram dataFlux = loadGalpropFile(galpropfilename.c_str());
        dataFlux.multiply( mcFlux.at(norm_at) / dataFlux.getContent()[0] );
        Reweighter reweighter(mcFlux, dataFlux);
        return reweighter;
}

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
	ResponseTOF->SetParameter(2,0);
	return ResponseTOF->Eval(vars->Beta);				
}

float GetLoweredBetaNaF  (Variables * vars) {	
	ResponseNaF->SetParameter(0,-0.000859132);
	ResponseNaF->SetParameter(1,-30.5065);
	ResponseNaF->SetParameter(2,0);
	return ResponseNaF->Eval(vars->BetaRICH_new);	
}

float GetLoweredBetaAgl  (Variables * vars) {	
	ResponseAgl->SetParameter(0,4.28781e-05);
	ResponseAgl->SetParameter(1,67.8521);
	ResponseAgl->SetParameter(2,0);
	
	return ResponseAgl->Eval(vars->BetaRICH_new);				
}

float GetRICHBDT(Variables* vars) {	return vars->BDTDiscr;}
float GetRupdown(Variables* vars) {  return fabs(vars->Rup - vars->Rdown)/vars->R;}
float GetChisquareX(Variables * vars) {return vars->Chisquare;}
float GetChisquareY(Variables * vars) {return vars->Chisquare_y;}

float GetEdepECAL(Variables * vars) {return vars->EdepECAL;}

float GetMomentumProxy(Variables *vars) {
	float momproxy = -1;
	if(vars->beta_SA>0 && vars->beta_SA<0.8) momproxy = (1/0.983) * vars->beta_SA/pow((1-pow(vars->beta_SA,2)),0.5);
	if(vars->beta_SA>0 &&vars->EdepECAL>0.45) momproxy = pow((vars->EdepECAL/0.3),(float)(1/1.1));

	return momproxy;	
}

float GetBetaSlowTOF     (Variables * vars) { return GetBetaGen(vars) - GetBetaTOF(vars);}

float GetBetaSlowRICH    (Variables * vars) { return GetBetaGen(vars) - GetBetaRICH(vars);}

float GetRigSlow          (Variables * vars) { return  vars->Momento_gen - GetRigidity(vars); } 


float GetBetaGen_SlowDownTOF(Variables * vars) {

	float beta = GetBetaGen(vars);
	ResponseTOF->SetParameter(0,0.00347548);
	ResponseTOF->SetParameter(1,5.8474);
	ResponseTOF->SetParameter(2,0);
	return beta-(0.9/vars->Massa_gen)*ResponseTOF->Eval(beta);				
}

float GetBetaGen_SlowDownNaF(Variables * vars) {

	float beta = GetBetaGen(vars);
	ResponseNaF->SetParameter(0,-0.000859132);
	ResponseNaF->SetParameter(1,-30.5065);
	ResponseNaF->SetParameter(2,0);
	return beta-(0.9/vars->Massa_gen)*ResponseNaF->Eval(beta);				
}

float GetBetaGen_SlowDownAgl(Variables * vars) {

	float beta = GetBetaGen(vars);
	ResponseAgl->SetParameter(0,4.28781e-05);
	ResponseAgl->SetParameter(1,67.8521);
	ResponseAgl->SetParameter(2,0);
	return beta-(0.9/vars->Massa_gen)*ResponseNaF->Eval(beta);				
}



// Tests of the Variables

float GetNToFClusters(Variables * vars) { return vars->beta_ncl;}
float GetUtofQ	(Variables * vars) {return vars->qUtof; }
float GetLtofQ	(Variables * vars) {return vars->qLtof; }
float GetTofChisqcn (Variables * vars) {return vars->chisqcn; }
float GetTofChisqtn  (Variables * vars) {return vars->chisqtn; }
float GetNTracks (Variables * vars) {return vars->nTrTracks; }
float GetTofOnTime (Variables * vars) {return vars->sumclsn; } 

