struct Variables{

	int 	   U_time;
        float 	   Latitude;
        float 	   Rcutoff;
        float      Livetime;
        int        PhysBPatt;
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
	int        NTofClusters;
	int        NTofClustersusati;
	float      Rup;
	float      Rdown;
	float      R;
	float      Chisquare;
	float      Beta;
	float      BetaR;
	int        NTrackHits;                             
	int        Richtotused;
	float      RichPhEl;
	float      R_L1;
	int        hitbits;
	float	   qL1;
	float	   qInner;
	float	   qUtof;
	float	   qLtof;


	//Discriminants
	float DistP=0;
	float DistD=0;
	float Likelihood=0;
	int   joinCutmask=0;

	float EdepTOFU=0;
	float EdepTOFD=0;
	float EdepTrack=0;	
	
	//MC vars
	float	   Momento_gen;
	float	   Massa_gen;		
	float	   mcweight=0;

	void ReadBranches(TTree * tree);
	void Update();
	void PrintCurrentState();
	void FillwithAnalysisVariables(TNtuple *grandezzesepd,bool isMC);
	void FillwithRawVariables(TNtuple *trig,bool isMC);
	void FillwithChargeInfos(TNtuple *Q,bool isMC);
	
};



void Variables::ReadBranches(TTree * tree){

	 tree->SetBranchAddress("U_time"	,&U_time);
         tree->SetBranchAddress("Latitude"	,&Latitude);
         tree->SetBranchAddress("Rcutoff35"	,&Rcutoff);
         tree->SetBranchAddress("Livetime"	,&Livetime);
         tree->SetBranchAddress("PhysBPatt"	,&PhysBPatt);
         tree->SetBranchAddress("R"		,&R_pre);
         tree->SetBranchAddress("BetaRaw"	,&Beta_pre);
         tree->SetBranchAddress("CUTMASK"	,&CUTMASK);
         tree->SetBranchAddress("trtrack_edep"	,&trtrack_edep);
         tree->SetBranchAddress("trtot_edep"	,&trtot_edep);
         tree->SetBranchAddress("TOFEndep"	,&Endep);
	 tree->SetBranchAddress("BetaRICH"	,&BetaRICH_new);
	 tree->SetBranchAddress("RICHmask"	,&RICHmask_new);
	 tree->SetBranchAddress("EnergyECAL"	,&EdepECAL);
	 tree->SetBranchAddress("NAnticluster"	,&NAnticluster);
	 tree->SetBranchAddress("NTofClusters"	,&NTofClusters);
	 tree->SetBranchAddress("NTofClustersusati",&NTofClustersusati);
	 tree->SetBranchAddress("Rup"		,&Rup);
	 tree->SetBranchAddress("Rdown"		,&Rdown);
	 tree->SetBranchAddress("R"		,&R);
	 tree->SetBranchAddress("Chisquare"	,&Chisquare);
	 tree->SetBranchAddress("BetaHR"	,&Beta);
	 tree->SetBranchAddress("BetaOld"	,&BetaR);
	 tree->SetBranchAddress("NTrackHits"	,&NTrackHits );                  
	 tree->SetBranchAddress("Richtotused"	,&Richtotused );  
	 tree->SetBranchAddress("RichPhEl"	,&RichPhEl);  
	 tree->SetBranchAddress("R_L1"		,&R_L1);  
	 tree->SetBranchAddress("hitbits"	,&hitbits);  
	 tree->SetBranchAddress("qL1"		,&qL1);  
	 tree->SetBranchAddress("qInner"	,&qInner);  
	 tree->SetBranchAddress("qUtof"		,&qUtof);  
	 tree->SetBranchAddress("qLtof"		,&qLtof);  
                                                      
	                                              
	 tree->SetBranchAddress("GenMomentum"	,&Momento_gen);  
	 tree->SetBranchAddress("GenMass"	,&Massa_gen);  		


}

void Variables::Update(){

	EdepTrack=0;
   	EdepTOFU=((*Endep)[0]+(*Endep)[1])/2;
   	EdepTOFD=((*Endep)[2]+(*Endep)[3])/2;
	for(int layer=1; layer<8; layer++) EdepTrack+=(*trtrack_edep)[layer];
	   EdepTrack=EdepTrack/7;
	DistP=100;
        DistD=100;
        Likelihood=0;
        joinCutmask=0;
	mcweight=0;
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
	
	cout<<"DistP:		      "<<DistP<<endl;
	cout<<"DistD:                 "<<DistD<<endl; 
	cout<<"Likelihood:            "<<Likelihood<<endl;
	cout<<"mcweight:            "<<mcweight<<endl;
	cout<<"joinCutmask:           "<<joinCutmask<<endl;


	cout<<"******"<<endl;
	cout<<endl;
        cout<<endl;

}



void Variables::FillwithAnalysisVariables(TNtuple *grandezzesepd,bool isMC){
	if(isMC) grandezzesepd->Fill(R,Beta,(*trtrack_edep)[0],Massa_gen,joinCutmask,PhysBPatt,EdepTOFU,EdepTrack,EdepTOFD,Momento_gen,BetaRICH_new,Likelihood,mcweight,DistD,DistP);
	else	  grandezzesepd->Fill(R,Beta,(*trtrack_edep)[0],joinCutmask,Latitude,PhysBPatt,EdepTOFU,EdepTrack,EdepTOFD,Rcutoff,BetaRICH_new,Likelihood,DistD,DistP);
	return;
}

void Variables::FillwithRawVariables(TNtuple *trig,bool isMC){
	if(isMC) trig->Fill(Massa_gen,Momento_gen,R_L1,R,Beta_pre,joinCutmask,(*trtrack_edep)[0],EdepTOFU,EdepTOFD,EdepTrack,BetaRICH_new,EdepECAL,PhysBPatt,mcweight);
	else     trig->Fill(U_time,Latitude,Rcutoff,R_L1,R,Beta_pre,joinCutmask,(*trtrack_edep)[0],EdepTOFU,EdepTOFD,EdepTrack,BetaRICH_new,EdepECAL,PhysBPatt,Livetime);
	return;
}
               
void Variables::FillwithChargeInfos(TNtuple *Q,bool isMC){
	if(isMC) Q->Fill(R,Beta,Massa_gen,joinCutmask,BetaRICH_new,DistD,DistP,Likelihood,R,qL1,qInner,qUtof,qLtof,Momento_gen );
	else     Q->Fill(R,Beta,joinCutmask,BetaRICH_new,DistD,DistP,Likelihood,R,qL1,qInner,qUtof,qLtof);
	return;
}


