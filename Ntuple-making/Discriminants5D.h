int GetUnusedLayers(int hitbits){ return 7 - std::bitset<32>(hitbits & 0b1111111).count(); }

bool IsNaF(Variables * vars){ if(((vars->joinCutmask)>>11)==512) return true; else return false;}

bool IsAgl(Variables * vars){ if(((vars->joinCutmask)>>11)==0 && vars->BetaRICH_new>0) return true; else return false;}

bool IsPreselected(Variables * vars){if((vars->joinCutmask&187)==187) return true; else return false;  }

void Likelihood(Variables * vars){

	float var[9]= {1,1,1,1,1,1,1,1,1};
	var[0]=vars->NAnticluster;
	var[1]=(vars->NTofClusters-vars->NTofClustersusati);
	if(vars->R!=0) var[2]=(fabs(vars->Rup-vars->Rdown)/vars->R);
        else var[2]=1;
	var[3]=GetUnusedLayers(vars->hitbits);
	var[4]=1; //former "fuoriX", not present in the new DSTs ->to be added
	var[5]=vars->Chisquare;
	
	if(IsNaF(vars)||IsAgl(vars)){
		var[6]=vars->Richtotused;
		var[7]=vars->RichPhEl;
	}
	double Ltrue=1;
        double Lfalse=1;
        if(!IsNaF(vars)&&!IsAgl(vars))
	for(int m=0; m<6; m++){
                Lfalse=Lfalse*Bkgnd[m]->Eval(var[m]);
                Ltrue=Ltrue*Signal[m]->Eval(var[m]);
	}
	else if(IsNaF(vars))
		for(int m=0; m<8; m++){
                Lfalse=Lfalse*BkgndNaF[m]->Eval(var[m]);
                Ltrue=Ltrue*SignalNaF[m]->Eval(var[m]);
        }
	else if(IsAgl(vars))
                for(int m=0; m<8; m++){
                Lfalse=Lfalse*BkgndAgl[m]->Eval(var[m]);
                Ltrue=Ltrue*SignalAgl[m]->Eval(var[m]);
        }
       

	vars->Likelihood=Ltrue/(Ltrue+Lfalse);
	
	return;
}


void Eval_Distance(Variables * vars,bool forprotons=true){
	float RGDT=0.01;
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
	float step = 0.075;
	float RG = vars->R;
	for(int z=0; z<1e6; z++) {
		if(forprotons) BT=protons->Eval(RGDT); else BT=deutons->Eval(RGDT);
		distR=(RGDT-RG)/(pow(RGDT,2)*Rig->Eval(RGDT));
		distB=(BT-vars->Beta)/(pow(BT,2)*beta->Eval(BT));
		if(IsNaF(vars)) distB=(BT-vars->BetaRICH_new)/(pow(BT,2)*betaNaF->Eval(BT));
		if(IsAgl(vars)) distB=(BT-vars->BetaRICH_new)/(pow(BT,2)*betaAgl->Eval(BT));
		distETOFU=(EdepTOFbeta->Eval(BT)-vars->EdepTOFU)/(pow(EdepTOFbeta->Eval(BT),2)*etofu->Eval(BT));
		distETrack=(EdepTrackbeta->Eval(BT)-vars->EdepTrack)/(pow(EdepTrackbeta->Eval(BT),2)*etrack->Eval(BT));
		distETOFD=(EdepTOFbeta->Eval(BT)-vars->EdepTOFD)/(pow(EdepTOFbeta->Eval(BT),2)*etofd->Eval(BT));

		Dist=pow(pow(distR,2)+pow(distB,2)+pow(distETrack,2)+pow(distETOFU,2)+pow(distETOFD,2),0.5);

		if(Dist<Dist1) {
			DR1=0;
			Dist1=Dist;	
		}else DR1++;

		//break conditions
		if(!IsNaF(vars)&&!IsAgl(vars)){
			if(DR1>25) break;
		}
		else{
			if(fabs(Dist-Dist2)<0.001) DR2++;
			if(DR2>4||DR1>300) break;
		}
		RGDT+=step;
		Dist2=Dist;
		if(z>1e5) std::cout<<"cazzo"<<std::endl;
	}
	if(forprotons) vars->DistP=Dist1;
	else 	       vars->DistD=Dist1;   			
	return;

}

