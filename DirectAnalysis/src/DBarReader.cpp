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

    ntpCompact    = new NtpCompact();	
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


float GetQTOF(float QL1, float QL2){
	if(QL1>0&&QL2>0) return (QL1+QL2)/2;
	else if(QL1>0) return QL1;
	else return QL2;
}

bool DBarReader::minTOF_Cpct(){
   if (!(ntpCompact->tof_beta_patt==0)) return false; 
   else return true;	
}

bool DBarReader::minTOF(){

   if (!(ntpTof->beta_patt==0)) return false; 
// else return (ntpTof->z_nhit>=3)&&(ntpTof->z_like>0)&&((ntpTof->clsn[0]+ntpTof->clsn[2])<=4);	
   else return true;  	
}

bool DBarReader::goldenTOF(){

//   if ( (ntpTof->flag) != 0   ) return false;
    if (  ntpTof->chisqcn > 5 || ntpTof->chisqcn < 0     ) return false;
    if (  ntpTof->chisqtn > 5 || ntpTof->chisqtn < 0     ) return false;
    return true;	
}



bool DBarReader::goldenTOF_Cpct(){

  // if ( (ntpTof->flag) != 0   ) return false; ?? Is it already in compact selection?
    if (  ntpCompact->tof_chisqcn > 5 || ntpCompact->tof_chisqcn < 0     ) return false;
    if (  ntpCompact->tof_chisqtn > 5 || ntpCompact->tof_chisqtn < 0     ) return false;
    return true;	
}


int DBarReader::RICHmaskConverter(){
    int richMASK;
    if(ntpRich->selection<0) richMASK = -ntpRich->selection;
    else richMASK = 0;
    if(ntpRich->is_naf) richMASK=richMASK|1024;;
    return richMASK;
}

int DBarReader::RICHmaskConverter_Cpt(){
    int richMASK;
    if(ntpCompact->rich_select<0) richMASK = -ntpCompact->rich_select;
    else richMASK = 0;
    if(ntpCompact->rich_select==1) richMASK=richMASK|1024;;
    return richMASK;
}

Long64_t DBarReader::ProtonCandidateSelection() {
  if ( (!ntpHeader)||(!ntpTof)||(!ntpTracker)||(!ntpTrd)||(!ntpRich) ) return -1;
  int n = 0;
  double rms = 0;
  double q_utof = (ntpTof) ? ntpTof->GetQ(n,rms,0x3) : 0;
  double q_ltof = (ntpTof) ? ntpTof->GetQ(n,rms,0xc) : 0;
  bool selection[22] = {
    ntpHeader->IsFTCP0(),
    ntpHeader->IsChargedPhysTrigger(),
    ntpTof->beta_ncl==4,
    ntpTof->beta>0.3,
    q_utof>0,
    (q_ltof>0.5)&&(q_ltof<2),
    ntpTof->chisqcn<100.,
    (ntpTof->chisqtn>0)&&(ntpTof->chisqtn<10.),
    ntpHeader->ntrtrack==1,
    ((ntpTracker->patty&0x2)!=0)&&((ntpTracker->patty&0xc)!=0)&&((ntpTracker->patty&0x30)!=0)&&((ntpTracker->patty&0xc0)!=0),
    (ntpTracker->q_inn[0]>0.7)&&(ntpTracker->q_inn[0]<1.5),
    (ntpTracker->pattxy&0x2)!=0,
    (ntpTracker->chisqn[7][0]>0)&&(ntpTracker->chisqn[7][0]<10),
    (ntpTracker->chisqn[7][1]>0)&&(ntpTracker->chisqn[7][1]<10),
    (ntpTof->clsn[0]+ntpTof->clsn[2])<=4,
    (ntpTrd->trdk_is_align_ok)&&(ntpTrd->trdk_is_calib_ok)&&(ntpTrd->trdk_like_nhit[0]>15)&&(ntpTrd->trdk_like_valid[0]>0),
    ntpTrd->chisq<10,
    ntpRich->correction_status==0,
    ntpRich->selection>0,
    ntpRich->nhit>2,
    fabs(ntpTof->beta-ntpRich->beta)<0.1*ntpRich->beta,
    fabs(ntpTracker->rig[7][1])>0.8
  };
  Long64_t pattern = 0;
  for (int i=0; i<22; i++) if (!selection[i]) pattern |= (0x1LL<<i);
  return pattern;
}



void DBarReader::FillVariables(int NEvent, Variables * vars){

    int w = Tree->GetEntry(NEvent);
    if(!(w>0)) { return;}
	
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

    vars->P_standard_sel    = ProtonCandidateSelection(); 	


    // Stroemer cutoff is in the tracker data
    vars->Rcutoff = ntpTracker->stoermer_cutoff[0];
    vars->Rcutoff_RTI  =   rtiInfo->cf[0][2][1]; //stoermer
    vars->Rcutoff_IGRFRTI = rtiInfo->cf[1][2][1]; //IGRF 
    // vars->Rcutoff_RTI  =   rtiInfo->cf[1][2][1]; //IGRF
    vars->Livetime_RTI =   rtiInfo->lf;    
    vars->good_RTI   = rtiInfo->good;
    vars->isinsaa = rtiInfo->isinsaa;


   cout<<"CUTOFF: "<<setprecision(5)<<vars->Rcutoff_IGRFRTI<<" "<<rtiInfo->cf[1][2][1]<<endl;
 

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
	    vars->Momento_gen_cpct = ntpMCHeader->momentum[0];
	

	    vars->GenX = ntpMCHeader->coo[0][0];  vars->GenPX = ntpMCHeader->dir[0];;
	    vars->GenY = ntpMCHeader->coo[0][1];  vars->GenPY = ntpMCHeader->dir[1];;
	    vars->GenZ = ntpMCHeader->coo[0][2];  vars->GenPZ = ntpMCHeader->dir[2];;
	    vars->MCClusterGeantPids = getPackedLayers_1to4();
    }
    if(!isMC) {
   //		if(((ntpHeader->trigpatt & 0x2) != 0) && ((ntpHeader->sublvl1&0x3E) ==0)) vars->PrescaleFactor*=100;
    }	    

    /////////////////////////////////// UNBIAS ////////////////////////////////////////
    vars->PhysBPatt = ntpHeader->sublvl1;
    vars->JMembPatt = ntpHeader->trigpatt;

/*    bool pct
			ntpTracker->chisqn[1][1] < vars->Chi2Ycut->Eval(abs(ntpTracker->rig[1][1])));	
    */
    bool goodChi2 =  (ntpTracker->chisqn[6][0] < 10 &&
			ntpTracker->chisqn[6][1] < 10);	
   


    /////////////////////////////// PRESELECTION CUTMASK //////////////////////////////////	
   
    if( ((ntpHeader->trigpatt & 0x2) != 0) && ((ntpHeader->sublvl1&0x3E) !=0) )      vars->CUTMASK |= 1 << 0;
    if( minTOF()                          )  vars->CUTMASK |= 1 << 1;
    if( (ntpHeader->trigpatt & 0x2) != 0)  vars->CUTMASK |= 1 << 2;
    if( ntpTracker->rig[1][1] != 0.0          )  vars->CUTMASK |= 1 << 3;
    if( goodChi2                          )  vars->CUTMASK |= 1 << 4;  
    if( goldenTOF()                       )  vars->CUTMASK |= 1 << 5;  
                                                                // 6
    if( vars->NTracks == 1 )  vars->CUTMASK |= 1 << 7;
    if( ntpTracker->rig[4][1] != 0.0          )  vars->CUTMASK |= 1 << 8;

    //////////////////////////////  Tracking Efficiency /////////////////////////////////////

    vars->trd_int_inside_tracker = ntpTrd->GetPatternInsideTracker();
    vars->theta_track     =ntpStandAlone->theta;;
    vars->phi_track       =ntpStandAlone->phi;
    vars->entrypointcoo[0]=ntpStandAlone->coo[0];	
    vars->entrypointcoo[1]=ntpStandAlone->coo[1];	
    vars->entrypointcoo[2]=ntpStandAlone->coo[2];	
    vars->beta_SA	  =ntpStandAlone->beta;
    vars->betapatt_SA     =ntpStandAlone->beta_patt;	
    vars->beta_ncl_SA	  =ntpStandAlone->beta_ncl;
    vars->beta_chiT_SA    =ntpStandAlone->beta_chisqtn;
    vars->Trd_chi_SA      =ntpStandAlone->trd_chisq;     	   
    vars->qUtof_SA	  =GetQTOF(ntpStandAlone->beta_q_lay[0],ntpStandAlone->beta_q_lay[1]);
    vars->qLtof_SA	  =GetQTOF(ntpStandAlone->beta_q_lay[2],ntpStandAlone->beta_q_lay[3]);
    vars->qTrd_SA	  =ntpStandAlone->trd_q;	
    vars->EdepECAL	  =ntpEcal->energyD[0];


    //////////////////////////////  L1 PICK-UP Efficiency /////////////////////////////////////

    vars->exthit_closest_q		=ntpStandAlone->exthit_closest_q[0][0]	 ;     
    vars->exthit_closest_status	=ntpStandAlone->exthit_closest_status[0];
    vars->hitdistfromint = pow(pow(ntpStandAlone->exthit_closest_coo[0][0] - ntpStandAlone->exthit_int[0][0],2)+ pow(ntpStandAlone->exthit_closest_coo[0][1] - ntpStandAlone->exthit_int[0][1],2),0.5); 

    /////////////////////////////// TRACKER ////////////////////////////////////
    
    vars->R     = ntpTracker->rig[1][1]; // 1 -- Inner tracker Kalman
    vars->Rup   = ntpTracker->rig[2][1]; // 2 -- Upper inner tracker
    vars->Rdown = ntpTracker->rig[3][1]; // 3 -- Lower inner tracker
    vars->R_L1  = ntpTracker->rig[4][1]; // 4 -- Inner + L1
    vars->R_noMS= ntpTracker->rig[7][1]; // 9 -- Inner tracker NoMS
    vars->RInner= ntpTracker->rig[6][1];  // 7 -- Inner tracker
    vars->R_sec = ntpTracker-> sec_inn_rig; // rig secondary track;

      vars->Chisquare         = ntpTracker->chisqn[1][0]; // 1 = Inner      , 0 = X side
    vars->Chisquare_L1      = ntpTracker->chisqn[4][0]; // 4 = L1 + Inner , 0 = X side
    vars->Chisquare_y       = ntpTracker->chisqn[1][1]; // 1 = Inner      , 1 = Y side
    vars->Chisquare_L1_y    = ntpTracker->chisqn[4][1]; // 4 = L1 + Inner , 1 = Y side
    vars->Chisquare_Inner   = ntpTracker->chisqn[6][0]; // 1 = Inner      , 0 = X side
    vars->Chisquare_Inner_y = ntpTracker->chisqn[6][1]; // 4 = Inner      , 1 = Y side
    	
    vars->hitbits           = ntpTracker->pattxy; 
     vars->patty           = ntpTracker->patty; 
    vars->FiducialVolume    = ntpTracker->GetPatternInsideTracker();   

    vars->qL1               = ntpTracker->q_lay[1][0];
    vars->qL1Status         = ntpTracker->q_clu_status[1][0];
    vars->qL1Status_SA      = ((ntpStandAlone->exthit_closest_status[0])&0x10013D);
    vars->qL2               = ntpTracker->q_lay[1][1];
    vars->qL2Status         = ntpTracker->q_clu_status[1][1];
    vars->qInner            = ntpTracker->q_inn[0];
    vars->clustertottrack   = ntpHeader->ntrrechit;
    vars->clustertrack      = countBits(vars->hitbits);
    int layerswithhit=0;
    vars->qL1InnerNoL2=0;
    for(int i=0;i<7;i++) {
	    if(i!=1){
		    if(ntpTracker->q_lay[1][i]>0) layerswithhit++;
		    vars->qL1InnerNoL2 += 	ntpTracker->q_lay[1][i];	
	    }	
    }
    vars->qL1InnerNoL2/= layerswithhit;

    vars->trtrack_edep = new std::vector<float>;
    vars->trtot_edep   = new std::vector<float>;
    for(int il=1;il<=9;il++) {
        vars->trtrack_edep->push_back(ntpTracker->edep_lay[1][il-1][0]);
        vars->trtot_edep->push_back(ntpTracker->edep_lay[1][il-1][0]);
    }
    
    /////////////////////////////// TOF ////////////////////////////////////
    
    int n = 0;
    double rms = 0;
  
    vars->qUtof   =GetQTOF(ntpTof->q_lay[0],ntpTof->q_lay[1]);
    vars->qLtof   =GetQTOF(ntpTof->q_lay[2],ntpTof->q_lay[3]);


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
    vars->Richtothits       = ntpHeader->nrichhit;
    vars->Richtotused       = ntpHeader->nrichhit-ntpRich->nhit;
    if(ntpRich->np_uncorr>0) vars->RichPhEl = ntpRich->np_exp_uncorr/ntpRich->np_uncorr; else vars->RichPhEl = 0;
    vars->RICHprob          = ntpRich->prob;
    vars->RICHPmts          = ntpRich->npmt;
    if(ntpRich->tot_p_uncorr>0) vars->RICHcollovertotal = ntpRich->np_uncorr/ntpRich->tot_p_uncorr; else vars->RICHcollovertotal=0;
    vars->RICHgetExpected   = ntpRich->np_exp_uncorr;
	vars->RichPhEl_tot = ntpRich->np_uncorr;
	vars->RichPhEl_ring = ntpRich->tot_p_uncorr;
	
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

 

   //////////////////////// CHECKS on VARIABLES ///////////////////////////
    vars-> beta_ncl =  ntpTof->beta_ncl;
    vars-> chisqcn  =  ntpTof->chisqcn; 
    vars-> chisqtn  =  ntpTof->chisqtn; 
    vars-> nTrTracks=  ntpHeader->ntrtrack;
    vars-> sumclsn  =  ntpTof->clsn[0] + ntpTof->clsn[2] ; 

}



void DBarReader::FillCompact(int NEvent, Variables * vars, float massgen, float momentumscaling){
	int w = Tree_Cpct->GetEntry(NEvent);
	if(!(w>0)) { cout<<"I/O error event: "<<NEvent<<endl; return;}
	vars->ResetVariables();

	///////////////////// MONTE CARLO //////////////////////////
	if(isMC) {
	    vars->mcweight =1;

	    vars->Massa_gen     = massgen;

	    if(vars->Massa_gen<2) vars->Charge_gen=1;
	    else vars->Charge_gen=2;	

	    vars->Momento_gen 	   = momentumscaling*ntpCompact->mc_momentum;
	    vars->Momento_gen_cpct = momentumscaling*ntpCompact->mc_momentum;
	    vars->Momento_gen_UTOF = momentumscaling*ntpMCHeader->momentum[6];
	    vars->Momento_gen_LTOF = momentumscaling*ntpMCHeader->momentum[16];
	    vars->Momento_gen_RICH = momentumscaling*ntpMCHeader->momentum[18];


	    vars->GenX = ntpMCHeader->coo[0][0];  vars->GenPX = ntpMCHeader->dir[0];;
	    vars->GenY = ntpMCHeader->coo[0][1];  vars->GenPY = ntpMCHeader->dir[1];;
	    vars->GenZ = ntpMCHeader->coo[0][2];  vars->GenPZ = ntpMCHeader->dir[2];;
	    vars->MCClusterGeantPids = getPackedLayers_1to4();
	}

	//cout<<NEvent<<" "<<vars->Momento_gen<<endl;
        vars->Event             = ntpSHeader->event;
  
	////////////////////// EVENT INFORMATION ///////////////////////////////////////
	vars->PrescaleFactor    = 1; 
	vars->P_standard_sel    = 0; 	

		
	// Stroemer cutoff is in the tracker data
	vars->Rcutoff = ntpCompact->trk_stoermer[1];
	vars->Rcutoff_RTI  =   rtiInfo->cf[0][2][1]; //stoermer
	vars->Rcutoff_IGRFRTI = rtiInfo->cf[1][2][1]; //IGRF 
	// vars->Rcutoff_RTI  =   rtiInfo->cf[1][2][1]; //IGRF

   
	vars->Livetime_RTI =   rtiInfo->lf;    
	vars->good_RTI   = rtiInfo->good;
	vars->isinsaa = rtiInfo->isinsaa;
	
	long long  int status = ntpCompact->status; 

	vars->NTRDSegments = (status%1000000000)/100000000; 
	
	vars->NTRDclusters = (status%100000000)/1000000; 
	
	vars->NTrackHits  = (status%100000)/10000; 
	
	vars->NTracks  = (status%10000)/1000; 
	
	vars->NAnticluster   = (status%100)/10; 
	
	//cout<<vars->NAnticluster<<" "<<vars->NTracks<<endl;	
	/////////////////////////////////// UNBIAS ////////////////////////////////////////
	vars->PhysBPatt = ntpCompact->sublvl1;
	vars->JMembPatt = ntpCompact->trigpatt;

//	bool goodChi2 =    (ntpCompact->trk_chisqn[0] < vars->Chi2Xcut->Eval(abs(ntpCompact->rig[1])) && ntpCompact->trk_chisqn[1] < vars->Chi2Ycut->Eval(abs(ntpCompact->rig[1])));	
	  
	bool goodChi2 =  (ntpCompact->trk_chisqn[3][0]< 10 &&
			ntpCompact->trk_chisqn[3][1]< 10);	


	/////////////////////////////// PRESELECT:wION CUTMASK //////////////////////////////////	

	if( ((ntpCompact->trigpatt & 0x2) != 0) && ((ntpCompact->sublvl1&0x3E) !=0) )      vars->CUTMASK |= 1 << 0;
	if( minTOF_Cpct()                          )  vars->CUTMASK |= 1 << 1; //minTOF already in compact selection
	if( (ntpCompact->trigpatt & 0x2) != 0)  vars->CUTMASK |= 1 << 2;
	if( ntpCompact->trk_kal_rig[1] != 0.0          )  vars->CUTMASK |= 1 << 3;
	if( goodChi2                          )  vars->CUTMASK |= 1 << 4;	
	if( goldenTOF_Cpct()                       )  vars->CUTMASK |= 1 << 5;  
	// 6
	if( vars->NTracks == 1 )  vars->CUTMASK |= 1 << 7;
	if( false  )  vars->CUTMASK |= 1 << 8;

	//////////////////////////////  Tracking Efficiency /////////////////////////////////////

	vars->beta_SA	  =ntpCompact->sa_tof_beta;
	vars->qUtof_SA	  =GetQTOF(ntpCompact->sa_tof_q_lay[0],ntpCompact->sa_tof_q_lay[1]);
	vars->qLtof_SA	  =GetQTOF(ntpCompact->sa_tof_q_lay[2],ntpCompact->sa_tof_q_lay[3]);
	vars->EdepECAL	  =ntpCompact->sa_ecal_edepd;
	vars->beta_chiT_SA    =ntpCompact->sa_tof_chisqtn;
        vars->Trd_chi_SA      =ntpCompact->sa_trd_chi2;     	   
        vars->beta_ncl_SA	  = ntpCompact->sa_tof_beta_ncl;  //already in selection
	vars->trd_int_inside_tracker = ntpCompact->sa_trd_fiducial; 
	vars->qTrd_SA         =ntpCompact->sa_trd_q;
        vars->Trd_chi_SA      =ntpCompact->sa_trd_chi2;		
	//////////////////////////////  L1 PICK-UP Efficiency /////////////////////////////////////
	vars->exthit_closest_q		=ntpCompact ->sa_exthit_ql1 ;     
	vars->hitdistfromint = pow(pow(ntpCompact->sa_exthit_dl1[0],2)+ pow(ntpCompact->sa_exthit_dl1[1],2),0.5);


	/////////////////////////////// TRACKER ////////////////////////////////////

	vars->R     = momentumscaling*ntpCompact->trk_kal_rig[1]; // 1 -- Inner tracker (Kalman)
	vars->R_L1  = momentumscaling*ntpCompact->trk_rig[1]; // 4 -- Inner + L1
        vars->R_noMS= momentumscaling*ntpCompact->trk_rig[4]; // 9 -- Inner tracker NoMS
        vars->RInner= momentumscaling*ntpCompact->trk_rig[3]; // 6 -- Inner tracker 

	vars->Chisquare         = ntpCompact->trk_kal_chisqn[0]; // 1 = Inner      , 0 = X side
	vars->Chisquare_L1      = ntpCompact->trk_chisqn[1][0]; // 4 = L1 + Inner , 0 = X side
	vars->Chisquare_y       = ntpCompact->trk_kal_chisqn[1]; // 1 = Inner      , 1 = Y side
	vars->Chisquare_L1_y    = ntpCompact->trk_chisqn[1][1]; // 4 = L1 + Inner , 1 = Y side
	vars->Chisquare_Inner   = ntpCompact->trk_chisqn[3][0];
    	vars->Chisquare_Inner_y = ntpCompact->trk_chisqn[3][1];

	//if(vars->R!=0) cout<<vars->Momento_gen<<" "<<vars->R<<" "<<vars->Momento_gen_cpct<<endl;

	vars->hitbits           = ntpCompact->trk_pattxy; 
	vars->patty		= ntpCompact->trk_patty;
	vars->FiducialVolume    = ntpCompact->trk_fiducial[1];   

	vars->qL1               = ntpCompact->trk_q_lay[0];
	if(vars->qL1>0)         vars->qL1Status         = 0; else   vars->qL1Status         = 1; 
	vars->qL1Status_SA     = ((ntpCompact->sa_exthit_status_l1)&0x10013D); 
	vars->qL2               = ntpCompact->trk_q_lay[1];
	if(vars->qL2>0)         vars->qL2Status         = 0; else   vars->qL2Status         = 1; 
	vars->qInner            = ntpCompact->trk_q_inn;

	int layerswithhit=0;
	vars->qL1InnerNoL2=0;
	for(int i=0;i<7;i++) {
		if(i!=1){
			if(ntpCompact->trk_q_lay[i]>0) layerswithhit++;
			vars->qL1InnerNoL2 += 	ntpCompact->trk_q_lay[i];	
		}	
	}
	vars->qL1InnerNoL2/= layerswithhit;


	/////////////////////////////// TOF ////////////////////////////////////

	vars->qUtof	  =GetQTOF(ntpCompact->tof_q_lay[0],ntpCompact->tof_q_lay[1]);
	vars->qLtof	  =GetQTOF(ntpCompact->tof_q_lay[2],ntpCompact->tof_q_lay[3]);
	vars->Beta    = ntpCompact-> tof_beta;

	vars->TOFchisq_s = ntpCompact->tof_chisqcn;
	vars->TOFchisq_t = ntpCompact->tof_chisqtn;

	/////////////////////////////// RICH ////////////////////////////////////
	vars->BetaRICH_new      = ntpCompact->rich_beta;
	vars->RICHmask_new      = RICHmaskConverter_Cpt();

	int richstatus = ntpCompact->rich_status;
	int nrich_hits, rich_pmts, rich_usedhits=0;	

	rich_usedhits= ((int)richstatus/10000);
	richstatus -= ((int)richstatus/10000)*10000;

	rich_pmts=((int)richstatus/100);
	richstatus -= ((int)richstatus/100)*100;	 

	nrich_hits = richstatus;	
        vars->Richtothits       = nrich_hits;
	vars->Richtotused       = nrich_hits-rich_usedhits;
	if(ntpCompact->rich_np>0) vars->RichPhEl = ntpCompact->rich_np_exp/ntpCompact->rich_np; else vars->RichPhEl = 0; 
	vars->RICHprob          = ntpCompact->rich_prob;
	vars->RICHPmts          = rich_pmts;
	if(ntpCompact->rich_np>0) vars->RICHcollovertotal = ntpCompact->rich_np/ntpCompact->rich_tot_np; else vars->RICHcollovertotal=0; 
	vars->RICHgetExpected   = ntpCompact->rich_np_exp;
	vars->RichPhEl_tot = ntpCompact->rich_tot_np;
	vars->RichPhEl_ring = ntpCompact->rich_np;
	

	vars->RICHTOFBetaConsistency = fabs(ntpCompact->rich_beta - ntpCompact->tof_beta)/ntpCompact->rich_beta;

	vars->BDTDiscr = ntpCompact->rich_bdt;

	//////////////////////// CHECKS on VARIABLES ///////////////////////////
	vars-> chisqcn  =  vars->TOFchisq_s;
	vars-> chisqtn  =  vars->TOFchisq_t;
	vars-> nTrTracks=  vars->NTracks; 

}



DBarReader::DBarReader(TTree * tree, bool _isMC, TTree * tree_RTI, TTree * tree_cpct) {
	Init();
	Tree = tree;
	Tree_RTI = tree_RTI	;   
	Tree_Cpct = tree_cpct; 

	cout<<Tree<<" "<<Tree_RTI<<" "<<Tree_Cpct<<endl;

	if(Tree){
		Tree->SetAutoFlush( 0 );
		Tree->SetBranchAddress( "SHeader" , &ntpSHeader     );
		Tree->SetBranchAddress( "Header"  , &ntpHeader     );
		Tree->SetBranchAddress( "Trd"     , &ntpTrd        );
		Tree->SetBranchAddress( "Tof"     , &ntpTof        );
		Tree->SetBranchAddress( "Tracker" , &ntpTracker    );
		Tree->SetBranchAddress( "Rich"    , &ntpRich       );
		Tree->SetBranchAddress( "Ecal"   , &ntpEcal       );
		//  Tree->SetBranchAddress( "Anti"   , &ntpAnti       );
		Tree->SetBranchAddress( "SA"     , &ntpStandAlone );
		if (_isMC) Tree->SetBranchAddress("MCHeader",&ntpMCHeader);

	}
	if(Tree_Cpct){
		Tree->SetAutoFlush( 0 );
		Tree_Cpct->SetAutoFlush( 0 );
		Tree_Cpct->SetBranchAddress( "SHeader" , &ntpSHeader    );
		Tree_Cpct->SetBranchAddress( "Compact" , &ntpCompact  );
		//Tree->BuildIndex("SHeader.run","SHeader.event");
		//Tree_Cpct->BuildIndex("SHeader.run","SHeader.event");
		//if(_isMC) if(Tree) Tree_Cpct->AddFriend(Tree);				
	}	
	if(!_isMC){
		if(Tree_RTI){
			Tree_RTI->SetAutoFlush( 0 );
			Tree_RTI->SetBranchAddress( "RTIInfo" , &rtiInfo  );		
			Tree_RTI->BuildIndex("SHeader.utime");
			//Tree_Cpct->BuildIndex("SHeader.utime");
			if(Tree)  {
			//	Tree->BuildIndex("SHeader.utime");
				Tree->AddFriend(Tree_RTI);
			}
			Tree_Cpct->AddFriend(Tree_RTI);
		}	
	}

	isMC = _isMC;
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

