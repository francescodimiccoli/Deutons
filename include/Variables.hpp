struct Variables{

	int 	   U_time;
        float 	   Latitude;
        float 	   Rcutoff;
	float 	   IGRFRcutoff;
        float      Livetime;
        float      PhysBPatt;
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
	float 	   joinCutmask=0;

	//Discriminants
	float DistP=0;
	float DistD=0;
	float Likelihood=0;

	//Summed E. deps.
	float EdepTOFU=0;
	float EdepTOFD=0;
	float EdepTrack=0;	
	float EdepL1=0;	
	
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
	
	void ReadAnalysisBranches(TNtuple * tree);
	void AnalysisVariablseReset();
	void PrintCurrentAnalysisState();
};



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
        Likelihood=0;
        joinCutmask=0;
	mcweight=0;
	
        RICHmask_new=(RICHmask_new&1023);
        joinCutmask=CUTMASK;
        joinCutmask=CUTMASK|(1<<10);
        joinCutmask=(int)joinCutmask|((RICHmask_new)<<11);	


	U_time-=1305200000; //time offset (in order to have small time stamp)	
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
	else	  grandezzesepd->Fill(R,Beta,(*trtrack_edep)[0],joinCutmask,Latitude,PhysBPatt,EdepTOFU,EdepTrack,EdepTOFD,IGRFRcutoff,BetaRICH_new,Likelihood,DistD,DistP,Rcutoff);
	return;
}

void Variables::FillwithRawVariables(TNtuple *trig,bool isMC){
	if(isMC) trig->Fill(Massa_gen,Momento_gen,R_L1,R,Beta_pre,joinCutmask,(*trtrack_edep)[0],EdepTOFU,EdepTOFD,EdepTrack,BetaRICH_new,U_time-Timebeg,PhysBPatt,mcweight);
	else     trig->Fill(U_time,Latitude,Rcutoff,IGRFRcutoff,R,Beta_pre,joinCutmask,(*trtrack_edep)[0],EdepTOFU,EdepTOFD,EdepTrack,BetaRICH_new,U_time-Timebeg,PhysBPatt,Livetime);
	return;
}
               
void Variables::FillwithChargeInfos(TNtuple *Q,bool isMC){
	if(isMC) Q->Fill(R,Beta,qL1,Massa_gen,joinCutmask,PhysBPatt,qUtof,qInner,qLtof,Momento_gen,BetaRICH_new,Likelihood,mcweight,DistD,DistP );
	else     Q->Fill(R,Beta,qL1,joinCutmask,Latitude,PhysBPatt,qUtof,qInner,qLtof,IGRFRcutoff,BetaRICH_new,Likelihood,DistD,DistP,Rcutoff);
	return;
}


void Variables::ReadAnalysisBranches(TNtuple * tree){

         tree->SetBranchAddress("R"		 ,&R);
         tree->SetBranchAddress("Beta"	         ,&Beta);
         tree->SetBranchAddress("EdepL1"	 ,&EdepL1);
	 tree->SetBranchAddress("Massa_gen"	 ,&Massa_gen);
	 tree->SetBranchAddress("Cutmask"	 ,&joinCutmask);
         tree->SetBranchAddress("PhysBPatt"	 ,&PhysBPatt);
         tree->SetBranchAddress("EdepTOF"	 ,&EdepTOFU);
	 tree->SetBranchAddress("EdepTOFU"        ,&EdepTOFU);
	 tree->SetBranchAddress("EdepTrack"	 ,&EdepTrack);	
		  	
 	 tree->SetBranchAddress("EdepTOFD"	 ,&EdepTOFD);
         tree->SetBranchAddress("Momentogen"	 ,&Momento_gen);
         tree->SetBranchAddress("BetaRICH_new"	 ,&BetaRICH_new);
	 tree->SetBranchAddress("LDiscriminant"	 ,&Likelihood);
	 tree->SetBranchAddress("mcweight"	 ,&mcweight);
         tree->SetBranchAddress("Dist5D"	 ,&DistD);
         tree->SetBranchAddress("Dist5D_P"	 ,&DistP);
	 tree->SetBranchAddress("Latitude"	 ,&Latitude);
 	 tree->SetBranchAddress("Rcutoff"        ,&Rcutoff);
 	 tree->SetBranchAddress("IGRFRcutoff"    ,&IGRFRcutoff);

	 tree->SetBranchAddress("qUtof"	 ,&qUtof);
	 tree->SetBranchAddress("qLtof"  ,&qLtof);
	 tree->SetBranchAddress("qInner" ,&qInner);
	 tree->SetBranchAddress("qL1" 	 ,&qL1);
	   	
	 return;
}

void Variables::AnalysisVariablseReset(){

	R	    =0;           
	R_L1        =0;
        Beta        =0; 
        EdepL1      =0; 
        Massa_gen   =0; 
        joinCutmask =0; 
        PhysBPatt   =0; 
        EdepTOFU    =0; 
        EdepTrack   =0; 
	            
        EdepTOFD    =0; 
        Momento_gen =0; 
        BetaRICH_new=0; 
        Likelihood  =0; 
	mcweight    =1; 
        DistD       =0; 
        DistP       =0; 
        Latitude    =0; 
        Rcutoff     =0; 
        IGRFRcutoff =0; 
	qL1         =0;
	qUtof       =0;
	qLtof	    =0;
	qInner	    =0;
	return;
}

void Variables::PrintCurrentAnalysisState(){
	cout<<endl;
	cout<<endl;
	cout<<"***Current Values of Analysis Variables:***"<<endl;

	cout<< "R: "	        <<R           	<<endl; 
	cout<< "Beta: "	 	<<Beta        	<<endl;         
	cout<< "EdepL1: "       <<EdepL1      	<<endl; 
	cout<< "Massa_gen: "    <<Massa_gen   	<<endl; 
	cout<< "Cutmask: "      <<(int)joinCutmask 	<<endl; 
	cout<< "PhysBPatt: "    <<(int)PhysBPatt  	<<endl; 
	cout<< "EdepTOF: "      <<EdepTOFU    	<<endl; 
	cout<< "EdepTrack: "    <<EdepTrack   	<<endl; 

	cout<< "EdepTOFD: "     <<EdepTOFD    	<<endl; 
	cout<< "Momentogen: "   <<Momento_gen 	<<endl; 
	cout<< "BetaRICH_new: " <<BetaRICH_new	<<endl; 
	cout<< "LDiscriminant: "<<Likelihood  	<<endl; 
	cout<< "mcweight: "     <<mcweight    	<<endl; 
	cout<< "Dist5D: "       <<DistD       	<<endl; 
	cout<< "Dist5D_P: "     <<DistP       	<<endl; 
	cout<< "Latitude: "     <<Latitude    	<<endl; 
	cout<< "Rcutoff: "      <<Rcutoff     	<<endl; 
	cout<< "IGRFRcutoff: "  <<IGRFRcutoff 	<<endl; 
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

float GetSmearedBetaTOF  (Variables * vars) {   
	if(vars->Massa_gen>0){
		float time = 1.2/(vars->Beta*3e-4);
		time = time + ToFsmearShift + Rand->Gaus(0,ToFsmearSigma);
		return 1.2/(time*3e-4);
	}
	else return GetBetaTOF(vars);
}

float GetSmearedRecMassTOF(Variables * vars) { 
	if(vars->Massa_gen>0){		
		float smearbeta=GetSmearedBetaTOF(vars);   
		return (vars->R/smearbeta)*pow((1-pow(smearbeta,2)),0.5);
	}
	else return GetRecMassTOF(vars);
}




