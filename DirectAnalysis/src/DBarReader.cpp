// g++ -c DBarReader.cpp -std=c++14 -I/afs/cern.ch/user/k/kostams/data_reduction/include/ -I`root-config --incdir`


#include "TTree.h"
#include "Variables.hpp"
#include "DBarReader.h"

void DBarReader::Init(){

    rtiInfo 	  = new RTIInfo();
    ntpSHeader    = new NtpSHeader(); 
    ntpHeader     = new NtpHeader(); 
    ntpMCHeader   = new NtpMCHeader();          
    ntpTrd        = new NtpTrd();
    ntpTof        = new NtpTof();
    ntpTracker    = new NtpTracker();
    ntpRich       = new NtpRich();
    ntpEcal       = new NtpEcal();
    ntpAnti       = new NtpAnti();
    ntpStandAlone = new NtpStandAlone();
}


UInt_t DBarReader::getPackedLayers_1to4() {
    UInt_t id1 = ntpMCHeader->hit_pid[0];
    UInt_t id2 = ntpMCHeader->hit_pid[1];
    UInt_t id3 = ntpMCHeader->hit_pid[2];
    UInt_t id4 = ntpMCHeader->hit_pid[3];

    id1 = id1 < 127 ? id1 : 0;
    id2 = id2 < 127 ? id2 : 0;
    id3 = id3 < 127 ? id3 : 0;
    id4 = id4 < 127 ? id4 : 0;

    return id1 + (id2 << 8) + (id3 << 16) + (id4 << 24);
}

int countBits(int n) {
    int counter = 0;
    while(n) {
        counter += n % 2;
        n >>= 1;
    }
    return counter;
}

bool DBarReader::minTOF(){

    return (ntpTof->z_nhit>=3)&&(ntpTof->z_like>0)&&((ntpTof->clsn[0]+ntpTof->clsn[2])<=4);	

}

bool DBarReader::goldenTOF(){

    if (!((ntpTof->trk_ncl==4)&&(ntpTof->beta_patt==0)&&(ntpTof->beta_ncl==4))) return false; 
    if ( (ntpTof->flag) != 0   ) return false;
    if (  ntpTof->chisqcn > 10 || ntpTof->chisqcn < 0     ) return false;
    if (  ntpTof->chisqtn > 10 || ntpTof->chisqtn < 0     ) return false;
    return true;	
}

int DBarReader::RICHmaskConverter(){
    int richMASK;
    if(ntpRich->selection<0) richMASK = -ntpRich->selection;
    else richMASK = 0;
    if(ntpRich->is_naf) richMASK=richMASK|1024;;
    return richMASK;
}




void DBarReader::FillVariables(int NEvent, Variables * vars){

    Tree->GetEntry(NEvent);

    vars->ResetVariables();

    ////////////////////// EVENT INFORMATION ///////////////////////////////////////
    vars->Run               = ntpSHeader->run;
    vars->Event             = ntpSHeader->event;
    vars->NEvent            = NEvent;
    vars->U_time            = ntpHeader->utc_time; 
    vars->Livetime          = ntpHeader->livetime;
    vars->Latitude          = ntpHeader->thetam;
    vars->ThetaS            = ntpHeader->thetas;
    vars->PhiS              = ntpHeader->phis;
    vars->PrescaleFactor    = ntpHeader->pres_weight;

    // Stroemer cutoff is in the tracker data
    vars->Rcutoff = ntpTracker->stoermer_cutoff[0];
    vars->Rcutoff_RTI  =   rtiInfo->cf[0][2][1]; //stoermer
    // vars->Rcutoff_RTI  =   rtiInfo->cf[1][2][1]; //IGRF
    
    vars->Livetime_RTI =   rtiInfo->lf;    
    vars->good_RTI   = rtiInfo->good;
    vars->isinsaa = rtiInfo->isinsaa;

    vars->NAnticluster      = ntpHeader->nanti;
    vars->NTofClusters      = ntpHeader->ntofh; // NOTE: That sees to be an error
    vars->NTofClustersusati = ntpTof->beta_ncl;
    vars->NTrackHits        = ntpHeader->ntrrechit;
    vars->NTracks           = ntpHeader->ntrtrack;
    vars->NTRDclusters      = ntpHeader->ntrdcluster;
    vars->NTRDSegments      = ntpHeader->ntrdsegment[1];

    ////////////////////// MONTE CARLO INFO ////////////////////////////////////////
    if(isMC) {
	    vars->Charge_gen   = ntpMCHeader->charge;
	    vars->Massa_gen     = ntpMCHeader->mass;
	    vars->Momento_gen = ntpMCHeader->momentum[0];
	    vars->Momento_gen_UTOF = ntpMCHeader->momentum[6];
	    vars->Momento_gen_LTOF = ntpMCHeader->momentum[16];
	    vars->Momento_gen_RICH = ntpMCHeader->momentum[18];

	    vars->GenX = ntpMCHeader->coo[0][0];  vars->GenPX = ntpMCHeader->dir[0];;
	    vars->GenY = ntpMCHeader->coo[0][1];  vars->GenPY = ntpMCHeader->dir[1];;
	    vars->GenZ = ntpMCHeader->coo[0][2];  vars->GenPZ = ntpMCHeader->dir[2];;
	    vars->MCClusterGeantPids = getPackedLayers_1to4();
    }
    if(!isMC) {
   		if(((ntpHeader->trigpatt & 0x2) != 0) && ((ntpHeader->sublvl1&0x3E) ==0)) vars->PrescaleFactor*=100;
    }	    

    /////////////////////////////////// UNBIAS ////////////////////////////////////////
    vars->PhysBPatt = ntpHeader->sublvl1;
    vars->JMembPatt = ntpHeader->trigpatt;

    bool goodChi2 =  (ntpTracker->chisqn[1][0] < vars->Chi2Xcut->Eval(abs(ntpTracker->rig[1])) &&
			ntpTracker->chisqn[1][1] < vars->Chi2Ycut->Eval(abs(ntpTracker->rig[1])));	
    
    /*bool goodChi2 =  (ntpTracker->chisqn[1][0] < 10 &&
			ntpTracker->chisqn[1][1] < 10);	
   */


    /////////////////////////////// PRESELECTION CUTMASK //////////////////////////////////	
   
    if( ((ntpHeader->trigpatt & 0x2) != 0) && ((ntpHeader->sublvl1&0x3E) !=0) )      vars->CUTMASK |= 1 << 0;
    if( minTOF()                          )  vars->CUTMASK |= 1 << 1;
    if( (ntpHeader->trigpatt & 0x2) != 0)  vars->CUTMASK |= 1 << 2;
    if( ntpTracker->rig[1] != 0.0          )  vars->CUTMASK |= 1 << 3;
    if( goodChi2                          )  vars->CUTMASK |= 1 << 4;  
    if( goldenTOF()                       )  vars->CUTMASK |= 1 << 5;  
                                                                // 6
    if( ntpHeader->nparticle == 1  && vars->NTracks == 1 )  vars->CUTMASK |= 1 << 7;
    if( ntpTracker->rig[4] != 0.0          )  vars->CUTMASK |= 1 << 8;

    //////////////////////////////  Tracking Efficiency /////////////////////////////////////

    vars->theta_track     =ntpStandAlone->theta;;
    vars->phi_track       =ntpStandAlone->phi;
    vars->entrypointcoo[0]=ntpStandAlone->coo[0];	
    vars->entrypointcoo[1]=ntpStandAlone->coo[1];	
    vars->entrypointcoo[2]=ntpStandAlone->coo[2];	
    vars->beta_SA	  =ntpStandAlone->beta;
    vars->betapatt_SA     =ntpStandAlone->beta_patt;	
    vars->qUtof_SA	  =(ntpStandAlone->beta_q_lay[0]+ntpStandAlone->beta_q_lay[1])/2;
    vars->qLtof_SA	  =(ntpStandAlone->beta_q_lay[2]+ntpStandAlone->beta_q_lay[3])/2;
    vars->qTrd_SA	  =ntpStandAlone->trd_q;	
    vars->EdepECAL	  =ntpEcal->energyE[0];


    //////////////////////////////  L1 PICK-UP Efficiency /////////////////////////////////////
    vars->exthit_int[0]		=ntpStandAlone->exthit_int[0][0];	          
    vars->exthit_int[1]		=ntpStandAlone->exthit_int[0][1];	          
    vars->exthit_int[2]		=ntpStandAlone->exthit_int[0][2];	          
    vars->exthit_closest_coo[0]	=ntpStandAlone->exthit_closest_coo[0][0];
    vars->exthit_closest_coo[1]	=ntpStandAlone->exthit_closest_coo[0][1];
    vars->exthit_closest_coo[2]	=ntpStandAlone->exthit_closest_coo[0][2];
    vars->exthit_closest_q		=ntpStandAlone->exthit_closest_q[0]	 ;     
    vars->exthit_closest_status	=ntpStandAlone->exthit_closest_status[0];
    vars->exthit_largest_coo[0]	=ntpStandAlone->exthit_largest_coo[0][0];
    vars->exthit_largest_coo[1]	=ntpStandAlone->exthit_largest_coo[0][1];
    vars->exthit_largest_coo[2]	=ntpStandAlone->exthit_largest_coo[0][2];
    vars-> exthit_largest_q		=ntpStandAlone->exthit_largest_q[0]	 ;     
    vars-> exthit_largest_status	=ntpStandAlone->exthit_largest_status[0];


    /////////////////////////////// TRACKER ////////////////////////////////////
    
    vars->R     = ntpTracker->rig[1]; // 1 -- Inner tracker
    vars->Rup   = ntpTracker->rig[2]; // 2 -- Upper inner tracker
    vars->Rdown = ntpTracker->rig[3]; // 3 -- Lower inner tracker
    vars->R_L1  = ntpTracker->rig[4]; // 4 -- Inner + L1
    vars->R_noMS= ntpTracker->rig[9]; // 9 -- Inner tracker NoMS


    vars->Chisquare         = ntpTracker->chisqn[1][0]; // 1 = Inner      , 0 = X side
    vars->Chisquare_L1      = ntpTracker->chisqn[4][0]; // 4 = L1 + Inner , 0 = X side
    vars->Chisquare_y       = ntpTracker->chisqn[1][1]; // 1 = Inner      , 1 = Y side
    vars->Chisquare_L1_y    = ntpTracker->chisqn[4][1]; // 4 = L1 + Inner , 1 = Y side
    vars->hitbits           = ntpTracker->pattxy; 
    vars->FiducialVolume    = ntpTracker->GetPatternInsideTracker();   

    vars->qL1               = ntpTracker->q_lay[1][0];
    vars->qL1Status         = ntpTracker->q_lay_status[1][0];
    vars->qL2               = ntpTracker->q_lay[1][1];
    vars->qL2Status         = ntpTracker->q_lay_status[1][1];
    vars->qInner            = ntpTracker->q_inn;
    vars->clustertottrack   = ntpHeader->ntrrechit;
    vars->clustertrack      = countBits(vars->hitbits);
    vars->qL1InnerNoL2	    = (ntpTracker->q_lay[1][0]+ntpTracker->q_lay[1][0]+ntpTracker->q_lay[1][0]+ntpTracker->q_lay[1][0]+ntpTracker->q_lay[1][0]+ntpTracker->q_lay[1][0]+ntpTracker->q_lay[1][0])/7;	

    vars->trtrack_edep = new std::vector<float>;
    vars->trtot_edep   = new std::vector<float>;
    for(int il=1;il<=9;il++) {
        vars->trtrack_edep->push_back(ntpTracker->edep_lay[1][il-1][0]);
        vars->trtot_edep->push_back(ntpTracker->edep_lay[1][il-1][0]);
    }
    
    /////////////////////////////// TOF ////////////////////////////////////
    
    // TODO: proper averaging inclding errors?
    vars->qUtof   = ( ntpTof->q_lay[0] + ntpTof->q_lay[1] ) / 2.0;
    vars->qLtof   = ( ntpTof->q_lay[2] + ntpTof->q_lay[3] ) / 2.0;

    vars->BetaR   = ntpTof->evgeni_beta;
    vars->Beta    = ntpTof->beta;
    
    vars->Endep = new std::vector<float>;
    for(int il=0;il<4;il++) {
       //cout<<vars->Endep<<endl; 
       vars->Endep->push_back(ntpTof->edep[il][0]);
    }
    vars->TOFchisq_s = ntpTof->chisqcn;
    vars->TOFchisq_t = ntpTof->chisqtn;
    vars->NBadTOF    = 0;
    for(int il=0;il<4;il++) { 
	if(ntpTof->flagp[il]!=0) vars->NBadTOF++;
    }	
    /////////////////////////////// TRD ////////////////////////////////////
    vars->TRDLike= ntpTrd->trdk_like_e[0];
    vars->TRDLikP= ntpTrd->trdk_like_p[0];
    vars->TRDLikD= ntpTrd->trdk_like_d[0];	    
    float Edep=0;
    float path=0;
    for(int il=0;il<20;il++) {		
	Edep+=ntpTrd->trdk_ampl[il];
	path+=ntpTrd->trdk_path[il];
    }
    vars->TRDEdepovPath = Edep/path;
    vars->EdepTRD = Edep;	

    /////////////////////////////// RICH ////////////////////////////////////
    vars->BetaRICH_new      = ntpRich->beta_corrected;
    vars->RICHmask_new      = RICHmaskConverter();
    vars->Richtotused       = ntpHeader->nrichhit-ntpRich->nhit;
    if(ntpRich->np_uncorr>0) vars->RichPhEl = ntpRich->np_exp_uncorr/ntpRich->np_uncorr; else vars->RichPhEl = 0;
    vars->RICHprob          = ntpRich->prob;
    vars->RICHPmts          = ntpRich->npmt;
    if(ntpRich->tot_p_uncorr>0) vars->RICHcollovertotal = ntpRich->np_uncorr/ntpRich->tot_p_uncorr; else vars->RICHcollovertotal=0;
    vars->RICHgetExpected   = ntpRich->np_exp_uncorr;

    vars->RICHLipBetaConsistency = fabs(ntpRich->lip_beta-ntpRich->beta_corrected);	
    vars->RICHTOFBetaConsistency = fabs(ntpRich->beta_corrected - ntpTof->beta)/ntpRich->beta_corrected;
    vars->RICHChargeConsistency  = ntpRich->q_consistency;
    vars->tot_hyp_p_uncorr	 = ntpRich->tot_hyp_p_uncorr[1]; 	
    vars->Bad_ClusteringRICH=0;
	 for (int is=0; is<10; is++) if ((ntpRich->clus_mean[is]-ntpRich->beta)>0.01) vars->Bad_ClusteringRICH++;
    vars->NSecondariesRICHrich=0;
	 for (int is=0; is< 5; is++) if ((ntpRich->pmt_np_uncorr[is]>5)&&(ntpRich->pmt_dist[is]>3.5)) vars->NSecondariesRICHrich++;

   vars->HitHValldir   =ntpRich->tot_hyp_hit_uncorr[0][0];
   vars->HitHVallrefl  =ntpRich->tot_hyp_hit_uncorr[0][1];
   vars->HitHVoutdir   =ntpRich->tot_hyp_hit_uncorr[1][0];
   vars->HitHVoutrefl  =ntpRich->tot_hyp_hit_uncorr[1][1];   

  


 
}
DBarReader::DBarReader(TTree * tree, bool _isMC, TTree * tree_RTI) {
    Init();
    Tree = tree;
    Tree_RTI = tree_RTI	;   

    if(Tree){
    Tree->SetBranchAddress( "SHeader" , &ntpSHeader     );
    Tree->SetBranchAddress( "Header"  , &ntpHeader     );
    Tree->SetBranchAddress( "Trd"     , &ntpTrd        );
    Tree->SetBranchAddress( "Tof"     , &ntpTof        );
    Tree->SetBranchAddress( "Tracker" , &ntpTracker    );
    Tree->SetBranchAddress( "Rich"    , &ntpRich       );
    Tree->SetBranchAddress( "Ecal"   , &ntpEcal       );
//  Tree->SetBranchAddress( "Anti"   , &ntpAnti       );
    Tree->SetBranchAddress( "SA"     , &ntpStandAlone );
    if(Tree_RTI){
	Tree_RTI->SetBranchAddress( "RTIInfo" , &rtiInfo  );		
    	Tree_RTI->BuildIndex("SHeader.utime");
	Tree->AddFriend(Tree_RTI);
    }	

    isMC = _isMC;
    if (isMC) Tree->SetBranchAddress("MCHeader",&ntpMCHeader);
    cout<<"************ TOT ENTRIES ***************"<<endl;
    cout<<Tree->GetEntries()<<endl;
    }	
}

DBarReader::DBarReader(TTree * tree, bool _isMC) {
    Init();
    Tree = tree;
    if(Tree){
	    Tree->SetBranchAddress( "RTIInfo" , &rtiInfo  );	
	    Tree->SetBranchAddress( "SHeader" , &ntpSHeader     );
	    Tree->SetBranchAddress( "Header"  , &ntpHeader     );
	    Tree->SetBranchAddress( "Trd"     , &ntpTrd        );
	    Tree->SetBranchAddress( "Tof"     , &ntpTof        );
	    Tree->SetBranchAddress( "Tracker" , &ntpTracker    );
	    Tree->SetBranchAddress( "Rich"    , &ntpRich       );
	    Tree->SetBranchAddress( "Ecal"   , &ntpEcal       );
	    //  Tree->SetBranchAddress( "Anti"   , &ntpAnti       );
	    Tree->SetBranchAddress( "SA"     , &ntpStandAlone );

	    isMC = _isMC;
	    if (isMC) Tree->SetBranchAddress("MCHeader",&ntpMCHeader);
	    cout<<"************ TOT ENTRIES ***************"<<endl;
	    cout<<Tree->GetEntries()<<endl;
    }	
}

