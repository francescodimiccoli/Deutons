int GetUnusedLayers(int hitbits){ return 7 - std::bitset<32>(hitbits & 0b1111111).count(); }


void Likelihood(Variables * vars){

	float var[9]= {0,0,0,0,0,0,0,0,0};
	var[0]=vars->NAnticluster;
	var[1]=(vars->NTofClusters-vars->NTofClustersusati);
	if(vars->R!=0) var[2]=(fabs(vars->Rup-vars->Rdown)/vars->R);
        else var[2]=1;
	var[3]=GetUnusedLayers(vars->hitbits);
	var[4]=1; //former "fuoriX", not present in the new DSTs ->to be added
	var[5]=vars->Chisquare;
	double Ltrue=1;
        double Lfalse=1;
	for(int m=0; m<6; m++){
		Lfalse=Lfalse*Bkgnd[m]->Eval(var[m]);
                Ltrue=Ltrue*Signal[m]->Eval(var[m]);	
	}	
	vars->Likelihood=Ltrue/(Ltrue+Lfalse);
	return;


}
