#include <vector>
#include <string>
#include <iostream>

#ifndef _PGTRACK_
#define _PGTRACK_
#include "TrTrack.h"
#endif

#include "selezioni.h"
#include "tkdcards.h"
#include "amschain.h"
#include <Tofrec02_ihep.h>
#include "Variables.hpp"


int getRICHmask(AMSEventR * ev){
    int cutmask=0;
    if( ! ev->pRichRing(0) )  return 0;

    RichRingR & ring = *(ev->pRichRing(0));
    ParticleR * particle =  ev->pParticle(0);

    float cut_prob=0.01;                        //  Kolmogorov test probability
    float cut_pmt=3;                            //  number of pmts
    float cut_collovertotal=0.4;                //  ring photoelctrons / total photoelectrons in the event
    float cut_chargeconsistency=10;              //  hit/PMT charge consistency test
    float cut_betaconsistency[2]={0.01,0.005};  //  beta_lip vs beta_ciemat consistency ([0]=NaF, [1]=aerogel)
    float cut_expphe[2]={1,2};                  //  expected number of photoelectrons   ([0]=NaF, [1]=aerogel)
    float cut_aerogelexternalborder=3500.;      //  aerogel external border (r**2)
    float cut_aerogel_nafborder[2]={17.,19.};   //  aerogel/NaF border                  ([0]=NaF, [1]=aerogel)

    int nbadtiles=5;                           
    int kbadtile[nbadtiles];
    kbadtile[0]=3;
    kbadtile[1]=7;
    kbadtile[2]=87;
    kbadtile[3]=100;
    kbadtile[4]=108;   //  tiles with bad beta recosntruction

    if(!ring.IsGood() || !ring.IsClean())           cutmask |= 1 << 0;
    if(ring.getProb() < cut_prob)                   cutmask |= 1 << 1;
    if(ring.getPMTs() < cut_pmt)                    cutmask |= 1 << 2;
    if(ring.getPhotoElectrons()/RichHitR::getCollectedPhotoElectrons() < cut_collovertotal) 
                                                    cutmask |= 1 << 3;
    if(ring.getPMTChargeConsistency() > cut_chargeconsistency)  
                                                    cutmask |= 1 << 4;

    float x=ring.getTrackEmissionPoint()[0];
    float y=ring.getTrackEmissionPoint()[1];

    if(ring.IsNaF()) {
        if(ring.getExpectedPhotoelectrons() < cut_expphe[0])    cutmask |= 1 << 5;
        if(ring.getBetaConsistency() > cut_betaconsistency[0])  cutmask |= 1 << 6;
        if(max(abs(x),abs(y)) > cut_aerogel_nafborder[0])       cutmask |= 1 << 7;
    } else {
        if(ring.getExpectedPhotoelectrons() < cut_expphe[1])    cutmask |= 1 << 5;
        if(ring.getBetaConsistency()>cut_betaconsistency[1])    cutmask |= 1 << 6;
        if(x*x+y*y > cut_aerogelexternalborder)                 cutmask |= 1 << 7;
        if(max(abs(x),abs(y)) < cut_aerogel_nafborder[1])       cutmask |= 1 << 7;
        for(int kbad=0;kbad<nbadtiles;kbad++) {
            if(ring.getTileIndex()==kbadtile[kbad])             cutmask |= 1 << 8;
        }
    }
    if(ring.IsNaF())                                            cutmask |= 1 << 9;

    // Old legacy selections
    int totalHits = ev->NRichHit();
    int hotspots = 0;
    int used = ev->pRichRing(0)->Used;
    for(int i = 0; i < ev->NRichHit(); i++) {
        RichHitR* Hit= ev->pRichHit(i);
        if(Hit->IsCrossed()) hotspots++;
    }
    float phelectronFrac = ring.getExpectedPhotoelectrons()/ring.getPhotoElectrons();
    bool pbetaoff = false;
    if(particle) {
        if(particle->pBeta()) {
            float pbeta = particle->pBeta()->Beta;
            float rbeta = ring.getBeta();
            pbetaoff = (pbeta - rbeta) / pbeta;
        }
    }

    if(ev->NRichRing() > 1)                         cutmask |= 1 << 10; 
    if(ring.getHits() < 5)                          cutmask |= 1 << 11;
    if(phelectronFrac < 0.4 || phelectronFrac > 2)  cutmask |= 1 << 12;
    if(ring.getProb() < 0.2)                        cutmask |= 1 << 13;
    if(totalHits - used - hotspots > 5)             cutmask |= 1 << 14;
    if(pbetaoff)                                    cutmask |= 1 << 15;

    return cutmask;
}



int main(int argc, char * argv[])
{
    if(argc < 2) return -1;


	AMSChain* ch= new AMSChain;
    ch->Add(argv[1]);
	int entries = ch->GetEntries();
	std::cout << "Processing " << entries << " events." << std::endl;

	TFile * File = new TFile(argv[2], "RECREATE");
	TTree * measure_stuff= new TTree("parametri_geo","parametri_geo");
    Variables * vars = new Variables;
    vars->RegisterBranches(measure_stuff); 

    /// --- Begin https://twiki.cern.ch/twiki/bin/view/AMS/PHeFluxStandardSelection
    AMSSetupR::RTI::UseLatest(6); 
	TkDBc::UseFinal();

    TRMCFFKEY_DEF::ReadFromFile = 0;
    TRFITFFKEY_DEF::ReadFromFile = 0;
    TRFITFFKEY.magtemp = 0;
    // --- End
    
    // --- Begin RTI stuff
    AMSSetupR setup;
    AMSSetupR::RTI::UseLatest(6);
    AMSSetupR::RTI rti;
    // --- End
    
    bool richflag = false;

	for(int ii=0;ii<entries;ii++)
    { 
        vars->ResetVariables();

        AMSEventR * ev = ch->GetEvent();
        ev->SetDefaultMCTuningParameters();


        if( ii % 10000 == 0){
            std::cout << "Processed " << ii << " out of " << entries << "\n";
            // Output your testing here
        }		

        if(ev->fHeader.Zenith() > 40) continue;

        ////////////////////// EVENT INFORMATION ///////////////////////////////////////
        vars->Run    = ev->Run();
        vars->Event  = ev->Event();
        vars->NEvent = ii;
        vars->U_time = ev->UTime();

        ////////////////////// MONTE CARLO INFO ////////////////////////////////////////
        MCEventgR * mc = ev->GetPrimaryMC();
        if(mc){
            TofMCPar::MCtuneDT = -87.0;
            TofMCPar::MCtuneST =  10.0;	

            vars->GenCharge   = mc->Charge;
            vars->GenMass     = mc->Mass;
            vars->GenMomentum = mc->Momentum;
            vars->GenX = mc->Coo[0];  vars->GenPX = mc->Dir[0];;
            vars->GenY = mc->Coo[1];  vars->GenPY = mc->Dir[1];;
            vars->GenZ = mc->Coo[2];  vars->GenPZ = mc->Dir[2];;
        }    

        ////////////////////// GEOG. VARIABLES ////////////////////////////////////////
        vars->Zenith   = ev->fHeader.Zenith();
        vars->Livetime = ev->LiveTime();
        vars->Latitude = ev->fHeader.ThetaM;
        vars->ThetaS   = ev->fHeader.ThetaS;
        vars->PhiS     = ev->fHeader.PhiS;
        if (AMSEventR::GetRTI(rti,ev->UTime())==0) {
            vars->Rcutoff35 = rti.cfi[2][1];
            vars->Rcutoff40 = rti.cfi[3][1];
        }
        if(ev->pParticle(0)){
            double dirTheta_ams= ev->pParticle(0)->Theta;
            double dirPhi_ams= ev->pParticle(0)->Phi;
            AMSDir amsdir(dirTheta_ams, dirPhi_ams);
            double cutoff =0;
            int momOrR = 1; //==0 GeV/c  ; ==1 GV
            int positive = 1; //==0 particle.mon == 1 (+) , == 2 (-)
            int ss = ev->pParticle(0)->GetStoermerCutoff(cutoff,momOrR,positive ,amsdir);
            vars->StoermerRcutoff=cutoff;
        }
        //////////////////////////////////////////////////////////////////////////////////		


        /////////////////////////////////// UNBIAS ////////////////////////////////////////
        Level1R* trig = ev->pLevel1(0);
        if(!trig) continue;

        vars->PhysBPatt = trig->PhysBPatt;
        vars->JMembPatt = trig->JMembPatt;

        if(ev->pRichRing(0) && richflag) {
            ev->pRichRing(0)->switchDynCalibration(); 
            richflag = true;
        }
        /////////////////////////////////////////////////////////////////////////////////////


        //////////////////////////////////// CUTMASK ////////////////////////////////////////
        bool MinTOF[2];
        minimumbiasTOF( ev, MinTOF);

        if( minimumbiasTRIGG(ev) )       vars->CUTMASK |= 1 << 0;
        if( MinTOF[0] )                  vars->CUTMASK |= 1 << 1;
        if( minimumbiasTRD(ev))          vars->CUTMASK |= 1 << 2;
        if(minimumbiasTRACKER(ev,3))     vars->CUTMASK |= 1 << 3;
        if(goldenTRACKER(ev,0,3))        vars->CUTMASK |= 1 << 4;
        if(goldenTOF_new(ev))            vars->CUTMASK |= 1 << 5;
        if(goldenTRD(ev,0,3))            vars->CUTMASK |= 1 << 6;
        if(OneParticle(ev))              vars->CUTMASK |= 1 << 7;
        if(minimumbiasTRACKER(ev,5))     vars->CUTMASK |= 1 << 8;
        if(goldenTRACKER(ev,0,5))        vars->CUTMASK |= 1 << 9;
        if(ev->IsTrackPickingUpNoise())  vars->CUTMASK |= 1 << 10;



        ////////////////////////////////////////////////////////////////////////////////////////


        //////////////////////  RIGIDITY //////////////////////////////////////////////////////
        TrTrackR* Tr = ev->pTrTrack(0);
        if(Tr){
            int fitID1=Tr->iTrTrackPar(1,1,1);
            int fitID2=Tr->iTrTrackPar(1,2,1);
            int fitID3=Tr->iTrTrackPar(1,3,1);
            int fitID5=Tr->iTrTrackPar(1,5,1);

            if(Tr->ParExists(fitID1)) vars->Rup   = Tr->GetRigidity(fitID1); 
            if(Tr->ParExists(fitID2)) vars->Rdown = Tr->GetRigidity(fitID2); 
            if(Tr->ParExists(fitID3)) vars->R     = Tr->GetRigidity(fitID3); 
            if(Tr->ParExists(fitID5)) vars->R_L1  = Tr->GetRigidity(fitID5); 

            if(Tr->ParExists(fitID3)) vars->Chisquare   = Tr->GetChisq(fitID3); else vars->Chisquare    = 1e7;
            if(Tr->ParExists(fitID5)) vars->Chisquare_L1= Tr->GetChisq(fitID5); else vars->Chisquare_L1 = 1e7;
        }
        //////////////////////////////////////////////////////////////////////////////////////////


        ///////////////////////////  E. DEP. /////////////////////////////////////////////////////

        /////// TRACK EDEP on track
        if(minimumbiasTRACKER(ev,3)) {
            vars->NTrackHits = Tr->NTrRecHit();
            for(int i=0; i < Tr->NTrRecHit(); i++){
                TrRecHitR * hit = ev->pTrTrack(0)->pTrRecHit(i);
                int ilay = hit->GetLayerJ()-1;
                TrClusterR* cluster = hit->GetYCluster();
                vars->trtrack_edep[ilay]  = cluster->GetEdep();
                vars->clustertrack++;
            }
        }

        /////// TRACK EDEP tot
        for (int i=0; i<ev->NTrCluster(); i++) {
            TrClusterR* cluster = ev->pTrCluster(i);
            int ilay = cluster->GetLayerJ() - 1;
            if(cluster->GetSide()==1) {
                vars->trtot_edep[ilay] += cluster->GetEdep(); 
                vars->clustertottrack++;
            }
        }
        /////// TOF EDEP
        for(int j=0; j<ev->NTofCluster(); j++)
        {
            int layer = ev->pTofCluster(j)->Layer - 1;
            vars->TOFEndepR[layer] += ev->pTofCluster(j)->Edep;
        }

        for(int j=0; j<ev->NTofClusterH(); j++)
        {
            int layer = ev->pTofClusterH(j)->Layer;
            vars->TOFEndep[layer] += ev->pTofClusterH(j)->GetEdep();
        }

        ////// TRD EDEP
        TrdTrackR * trdtrack = ev->pTrdTrack(0);
        if(trdtrack){
            for(int j = 0; j < trdtrack->NTrdSegment(); j++) {
                for(int i = 0; i < trdtrack->pTrdSegment(j)->NTrdCluster(); i++) {
                    vars->EdepTRD += trdtrack->pTrdSegment(j)->pTrdCluster(i)->EDep;
                    vars->NTRDclusters++;
                }
            }
        }


        /////// ECAL EDEP
        int NhitECAL = ev->NEcalHit();
        vars->EnergyECAL=-100;
        if(ev->NEcalShower() == 1) {
            EcalShowerR* show = ev->pEcalShower(0); 
            vars->EnergyECAL = show->GetCorrectedEnergy(2, 2);
        }
        ////////////////////////////////////////////////////////////////////////////////////////////


        ///////////////////////////////// BETA TOF ////////////////////////////////////////////////
        if(ev->pBeta(0))  vars->BetaOld  = ev->pBeta(0)->Beta;

        float bestBeta = 0;
        BetaHR * betaH = ev->pBetaH(0);
        if(betaH) { 
            vars->BetaRaw = betaH->GetBeta();

            TofRecH::BuildOpt = 0;
            TofRecH::ReBuild();
            vars->BetaHR  = betaH->GetBeta();
            bestBeta = vars->BetaHR;

            if(mc){
                betaH->DoMCtune();
                vars->BetaHRS = betaH->GetBeta();
                bestBeta = vars->BetaHRS;
            }
        }
        /////////////////////////////////////////////////////////////////////////////////////////////


        ////////////////////////////// CHARGE /////////////////////////////////////////////////////////
        for(int i = 0; i < ev->nTrRecHit(); i++) {
            TrRecHitR * hit = ev->pTrRecHit(i);
            if(hit->GetLayerJ() != 1) continue;

            float q = hit->GetQ(1);
            if(q > vars->qL1Max) {
                vars->qL1Max = q;
                vars->qL1MaxStatus = hit->GetQStatus();
            }
        }


        if(minimumbiasTRACKER(ev,3)) {
            vars->qL1       = Tr->GetLayerJQ(1, vars->BetaRaw);
            vars->qL1Status = Tr->GetLayerJQStatus(1);
            vars->qInner= Tr->GetInnerQ(vars->BetaRaw);
        }
        int   utoflay, ltoflay=0;
        float utofrms, ltofrms=0;
        if(ev->pBetaH(0)) vars->qUtof = ev->pBetaH(0)->GetQ(utoflay,utofrms,2,TofClusterHR::DefaultQOptIonW,1100,0,vars->R);
        if(ev->pBetaH(0)) vars->qLtof = ev->pBetaH(0)->GetQ(ltoflay,ltofrms,2,TofClusterHR::DefaultQOptIonW,11  ,0,vars->R);
        /////////////////////////////////////////////////////////////////////////////////////////////////

        vars->RICHmask = getRICHmask(ev); 

        //////////////////////////////// BETA RICH /////////////////////////////////////////////////
        if(ev->NRichRing() > 0){
            vars->BetaRICH=ev->pRichRing(0)->getBeta();
            int totali = ev->NRichHit();
            int hotspots=0;
            int usate=ev->pRichRing(0)->Used;
            for(int i=0;i<ev->NRichHit();i++)
            {
                RichHitR* Hit= ev->pRichHit(i);
                if(Hit->IsCrossed()) hotspots++;
            }
            vars->Richtotused=totali-usate-hotspots;
            vars->RichPhEl=ev->pRichRing(0)->getExpectedPhotoelectrons()/ev->pRichRing(0)->getPhotoElectrons();
        }
        ///////////////////////////////////////////////////////////////////////////////////////////////	


        /////////////////////////////// LIKELIHOOD VARIABLES //////////////////////////////////////////
        // for next production: add the info about tracker hits being XY or only Y 
        if(Tr) vars->hitbits = Tr->GetHitBits();

        vars->NAnticluster = ev->NAntiCluster();
        vars->NTRDSegments = ev->NTrdSegment();
        vars->NTofClusters = ev->NTofCluster();
        if(ev->pBetaH(0))
            vars->NTofClustersusati = ev->pBetaH(0)->NTofClusterH();
        ///////////////////////////////////////////////////////////////////////////////////////


        measure_stuff->Fill();
    }

	File->Write();
	File->Close();

	return 0;
}

