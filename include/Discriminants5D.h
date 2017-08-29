int GetUnusedLayers(int hitbits){ return 7 - std::bitset<32>(hitbits & 0b1111111).count(); }


void Likelihood(Variables * vars){

	float var[9]= {1,1,1,1,1,1,1,1,1};
	var[0]=vars->NAnticluster;
	var[1]=(vars->NTofClusters-vars->NTofClustersusati);
	if(vars->R!=0) var[2]=(fabs((vars->Rup-vars->Rdown)/vars->R));
        else var[2]=1;
	var[3]=GetUnusedLayers(vars->hitbits);
	var[4]=1; //former "fuoriX", not present in the new DSTs ->to be added
	var[5]=vars->Chisquare;
	
	if(!ApplyCuts("IsOnlyFromToF",vars)){
		var[6]=vars->Richtotused;
		var[7]=vars->RichPhEl;
	}
	double Ltrue=1;
        double Lfalse=1;
        if(ApplyCuts("IsOnlyFromToF",vars))
	for(int m=0; m<6; m++){
                Lfalse=Lfalse*Bkgnd[m]->Eval(var[m]);
                Ltrue=Ltrue*Signal[m]->Eval(var[m]);
		}
	else if(ApplyCuts("IsFromNaF",vars))
		for(int m=0; m<8; m++){
                Lfalse=Lfalse*BkgndNaF[m]->Eval(var[m]);
                Ltrue=Ltrue*SignalNaF[m]->Eval(var[m]);
        	}
	else if(ApplyCuts("IsFromAgl",vars))
		for(int m=0; m<8; m++){
			Lfalse=Lfalse*BkgndAgl[m]->Eval(var[m]);
			Ltrue=Ltrue*SignalAgl[m]->Eval(var[m]);
		}

	if(Ltrue<0) Ltrue=0;
	if(Lfalse<0) Lfalse=0;

	
	vars->Likelihood=Ltrue/(Ltrue+Lfalse);
	
	return;
}

float Remove_SlowDown(float meas, float gen, float corr){
	float MEAS = 1/meas;
	float GEN  = 1/gen;
	GEN = MEAS/(corr+1);
	return 1/GEN;

}

void Eval_DistanceFromP(Variables * vars){
	float RGDT=0.35;
        float BT=0.01;
        float distR=0;
        float distB=0;
        float distETOFU=0;
        float distETOFD=0;
        float distETrack=0;
        float distETRD=0;
        float DistTOFU=0;
        float DistTOFD=0;
        float DistTrack=0;
        float DistTRD=0;
        float Dist=0;
	int DR1=0;
        int DR2=0;
	float Dist1=1000000;
        float Dist2=1000000;
	float step = 0.1;
	
	// before slow down rig and vel.
	float RG;   
	float BETA,BETANAF,BETAAGL;
	float Rmin=0;
	//cout<<endl;
	for(int z=0; z<1e6; z++) {
		BT=protons->Eval(RGDT);
		
		//slow down correction
		RG      = Remove_SlowDown(vars->R,    RGDT, Rigmean_p->Eval(RGDT)); 
		BETA    = Remove_SlowDown(vars->Beta, BT,   betamean_p->Eval(BT));
		if(ApplyCuts("IsFromNaF",vars)) BETANAF = Remove_SlowDown(vars->BetaRICH_new, BT, betaNaFmean_p->Eval(BT));
		if(ApplyCuts("IsFromAgl",vars)) BETAAGL = Remove_SlowDown(vars->BetaRICH_new, BT, betaAglmean_p->Eval(BT));
	
		distR=(RGDT-RG)/(pow(RGDT,2)*Rig_p->Eval(RGDT));
		distB=(BT-BETA)/(pow(BT,2)*beta_p->Eval(BT));
		if(ApplyCuts("IsFromNaF",vars)) distB=(BT-BETANAF)/(pow(BT,2)*betaNaF_p->Eval(BT));
		if(ApplyCuts("IsFromAgl",vars)) distB=(BT-BETAAGL)/(pow(BT,2)*betaAgl_p->Eval(BT));
		distETOFU=(1/EdepTOFbeta_p->Eval(BT)-vars->EdepTOFU)/(pow(1/EdepTOFbeta_p->Eval(BT),2)*etofu_p->Eval(BT));
		distETrack=(1/EdepTrackbeta_p->Eval(BT)-vars->EdepTrack)/(pow(1/EdepTrackbeta_p->Eval(BT),2)*etrack_p->Eval(BT));
		distETOFD=(1/EdepTOFbeta_p->Eval(BT)-vars->EdepTOFD)/(pow(1/EdepTOFbeta_p->Eval(BT),2)*etofd_p->Eval(BT));

		Dist=pow(pow(distR,2)+pow(distB,2)+pow(distETrack,2)+pow(distETOFU,2)+pow(distETOFD,2),0.5);
		if(Dist<Dist1) {
			DR1=0;
			Dist1=Dist;	
			Rmin=RGDT;
		}else DR1++;

		//break conditions
		if(ApplyCuts("IsOnlyFromToF",vars)){
			if(DR1>30) break;
		}
		else{
			if(fabs(Dist-Dist2)<0.001) DR2++;
			if(DR2>4||DR1>300) break;
		}
		RGDT+=step;
		Dist2=Dist;
		if(z>1e5) std::cout<<"cazzo"<<std::endl;
	}
	vars->DistP=Dist1;
	return;

}




void Eval_DistanceFromD(Variables * vars){
	float RGDT=0.35;
        float BT=0.01;
        float distR=0;
        float distB=0;
        float distETOFU=0;
        float distETOFD=0;
        float distETrack=0;
        float distETRD=0;
        float DistTOFU=0;
        float DistTOFD=0;
        float DistTrack=0;
        float DistTRD=0;
        float Dist=0;
	int DR1=0;
        int DR2=0;
	float Dist1=1000000;
        float Dist2=1000000;
	float step = 0.1;
	
	// before slow down rig and vel.
	float RG;   
	float BETA,BETANAF,BETAAGL;
	float Rmin=0;
	//cout<<endl;
	for(int z=0; z<1e6; z++) {
		BT=deutons->Eval(RGDT);
		
		//slow down correction
		RG      = Remove_SlowDown(vars->R,    RGDT, Rigmean_d->Eval(RGDT)); 
		BETA    = Remove_SlowDown(vars->Beta, BT,   betamean_d->Eval(BT));
		if(ApplyCuts("IsFromNaF",vars)) BETANAF = Remove_SlowDown(vars->BetaRICH_new, BT, betaNaFmean_d->Eval(BT));
		if(ApplyCuts("IsFromAgl",vars)) BETAAGL = Remove_SlowDown(vars->BetaRICH_new, BT, betaAglmean_d->Eval(BT));
	
		distR=(RGDT-RG)/(pow(RGDT,2)*Rig_d->Eval(RGDT));
		distB=(BT-BETA)/(pow(BT,2)*beta_d->Eval(BT));
		if(ApplyCuts("IsFromNaF",vars)) distB=(BT-BETANAF)/(pow(BT,2)*betaNaF_d->Eval(BT));
		if(ApplyCuts("IsFromAgl",vars)) distB=(BT-BETAAGL)/(pow(BT,2)*betaAgl_d->Eval(BT));
		distETOFU=(1/EdepTOFbeta_d->Eval(BT)-vars->EdepTOFU)/(pow(1/EdepTOFbeta_d->Eval(BT),2)*etofu_d->Eval(BT));
		distETrack=(1/EdepTrackbeta_d->Eval(BT)-vars->EdepTrack)/(pow(1/EdepTrackbeta_d->Eval(BT),2)*etrack_d->Eval(BT));
		distETOFD=(1/EdepTOFbeta_d->Eval(BT)-vars->EdepTOFD)/(pow(1/EdepTOFbeta_d->Eval(BT),2)*etofd_d->Eval(BT));

		Dist=pow(pow(distR,2)+pow(distB,2)+pow(distETrack,2)+pow(distETOFU,2)+pow(distETOFD,2),0.5);
		if(Dist<Dist1) {
			DR1=0;
			Dist1=Dist;	
			Rmin=RGDT;
		}else DR1++;

		//break conditions
		if(ApplyCuts("IsOnlyFromToF",vars)){
			if(DR1>30) break;
		}
		else{
			if(fabs(Dist-Dist2)<0.001) DR2++;
			if(DR2>4||DR1>300) break;
		}
		RGDT+=step;
		Dist2=Dist;
		if(z>1e5) std::cout<<"cazzo"<<std::endl;
	}
	vars->DistD=Dist1;
	return;

}



