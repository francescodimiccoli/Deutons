#include "TemplateFITbetasmear.h"

TSpline3 * GetSplineFromHisto(TH1F * Graph, Binning bins){
	
    
        double X[Graph->GetNbinsX()];
        double Y[Graph->GetNbinsX()];
   	cout<<"********* EXPOSURE SPLINE **********"<<endl; 
        for(int i=0; i<Graph->GetNbinsX(); i++){        X[i]=bins.RigBinCent(i); Y[i]=Graph->GetBinContent(i+1); }
        TSpline3 *Exposure = new TSpline3("Exposure",X,Y,Graph->GetNbinsX());
        return Exposure;
}

void TFit::RegularizeTemplateError(){
	float integral = Templ_P->Integral();
	float mincontent=9999999;
	for(int i=0;i<Templ_P->GetNbinsX();i++)
		if(Templ_P->GetBinContent(i+1)>0&&Templ_P->GetBinContent(i+1)<mincontent){
			mincontent = Templ_P->GetBinContent(i+1);
	}
	Templ_P->Scale(1/mincontent);
	if(Templ_P->Integral(Templ_P->FindBin(3),Templ_P->FindBin(7))>0)
	{
		float content = Templ_P->Integral(Templ_P->FindBin(3),Templ_P->FindBin(7))/(Templ_P->FindBin(7)-Templ_P->FindBin(3));
		for(int i=0;i<Templ_P->GetNbinsX();i++) {
			if(Templ_P->GetBinContent(i+1)==0) 
				Templ_P->SetBinError(i+1,pow(content,0.5));
		}
	}
	Templ_P->Scale(integral/Templ_P->Integral());
}

void TemplateFIT::Eval_TransferFunction(){
	cout<<"Transf:"<<endl;	
	for(int bin=0;bin<fits.size();bin++){
		cout<<"Transf: "<< fits[bin][0][5]->DataPrim<<endl;
		TH1F * transferfunction = (TH1F *) fits[bin][0][5]->DataPrim->Clone();
		transferfunction->Sumw2();
		transferfunction->Divide(fits[bin][0][5]->DataAmbient);
		TransferFunction.push_back(transferfunction);
	}
	return;
}

void TemplateFIT::Fill(TTree * treeMC,TTree * treeDT, Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars) ){
	
	cout<<basename.c_str()<<" Filling ... (Data)"<< endl;
	vars->ReadBranches(treeDT);

	for(int i=0;i<treeDT->GetEntries()/FRAC;i++){
		UpdateProgressBar(i, treeDT->GetEntries()/FRAC);
		treeDT->GetEvent(i);
		vars->Update();
		FillEventByEventData(vars,var,discr_var);
	}
	
	cout<<basename.c_str()<<" Filling ... (MC Protons)"<< endl;
	vars->ReadBranches(treeMC);

	for(int i=0;i<treeMC->GetEntries();i++){
		if(i%(int)FRAC!=0) continue;
		UpdateProgressBar(i, treeMC->GetEntries());
		treeMC->GetEvent(i);
		vars->Update();

		if(BadEvSim) BadEvSim->LoadEvent(vars);
		FillEventByEventMC(vars,var,discr_var);
	}

	return;
}


void TemplateFIT::FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars)){
	int kbin;
	kbin = 	bins.GetBin(discr_var(vars));
	if(!ApplyCuts("IsData",vars)) return;
							
	if(ApplyCuts((cut+"&RigSafetyCut").c_str(),vars)&&kbin>=0){
		//vars->PrintCurrentState();
		for(int i=0;i<systpar.steps;i++)
			for(int j=0;j<systpar.steps;j++){
				if(fits[kbin][i][j]->Data){
				if(ApplyCuts(cut,vars)) fits[kbin][i][j]->DataAmbient->Fill(var(vars),vars->PrescaleFactor);		
				if(ApplyCuts((cut+"&"+cutoff).c_str(),vars)) fits[kbin][i][j]->Data->Fill(var(vars),vars->PrescaleFactor);		
				if(ApplyCuts(cutprimary,vars)) 		     fits[kbin][i][j]->DataPrim->Fill(var(vars),vars->PrescaleFactor);
					}
				}
	}
	return;	
}


void TemplateFIT::SimpleExtractPrimaries(){
        for(int bin=0;bin<bins.size();bin++)
                for(int sigma=0;sigma<systpar.steps;sigma++)
                        for(int shift=0;shift<systpar.steps;shift++){
                               cout<<"basename: "<<fits[bin][sigma][shift]->DataPrim<<" "<<fits[bin][sigma][shift]->Data<<" "<<fits[bin][sigma][shift]->DataAmbient<<endl;
				fits[bin][sigma][shift]->Templ_DPrim=(TH1F*) fits[bin][sigma][shift]->Templ_D->Clone();
                                fits[bin][sigma][shift]->Templ_PPrim=(TH1F*) fits[bin][sigma][shift]->Templ_P->Clone();
				fits[bin][sigma][shift]->Templ_HePrim=(TH1F*) fits[bin][sigma][shift]->Templ_He->Clone();

                                fits[bin][sigma][shift]->Templ_PPrim->Multiply(fits[bin][sigma][shift]->DataPrim);
                                fits[bin][sigma][shift]->Templ_DPrim->Multiply(fits[bin][sigma][shift]->DataPrim);
                                fits[bin][sigma][shift]->Templ_HePrim->Multiply(fits[bin][sigma][shift]->DataPrim);

                                fits[bin][sigma][shift]->Templ_PPrim->Divide(fits[bin][sigma][shift]->DataAmbient);
                                fits[bin][sigma][shift]->Templ_DPrim->Divide(fits[bin][sigma][shift]->DataAmbient);
		                fits[bin][sigma][shift]->Templ_HePrim->Divide(fits[bin][sigma][shift]->DataAmbient);
				
                }
} 
 
float Systpar::GetShiftStart(){
	if(mode==0) return (-shift/5000)*100;   //delta T.o.F. (%)
	if(mode==1) return (-shift*1e-4)*0.75*100;   //delta (1/beta) (%)
	if(mode==2) return (-shift/(10e4*acos(1/(1.2*0.98))))*100; //delta angle (%)			///(1-1/(1.2*cos( ((acos(1/(1.2)) *10e4)+shift)/10e4)))*100; 
	if(mode==3) return -shift*1e-5*100;     //delta (1/beta) (%)
}

float Systpar::GetShiftEnd(){
	if(mode==0) return (shift/5000)*100;
	if(mode==1) return (shift*1e-4)*100;
	if(mode==2) return (shift/(10e4*acos(1/(1.2*0.98))))*100;
	if(mode==3) return shift*1e-5*100;
}

float Systpar::GetSigmaStart(){
	if(mode==0) return  0;  //sigma T.o.F.(%)
	if(mode==1) return  ((0.0284/2 - 0.0284/2*sigma/100)/0.0318)*100; //delta (sigma(1/beta)) (%)
	if(mode==2) return  0; //delta (sigma angle) %
	if(mode==3) return  ((1.15694e-03*(1-sigma/100.) - 1.36634e-03)/1.36634e-03)*100; //delta (sigma(1/beta)) (%)
}

float Systpar::GetSigmaEnd(){
	if(mode==0) return ((2*sigma/5000)/0.0318)*100;
	if(mode==1) return ((0.0284/2 + 0.0284/2*sigma/100)/0.0318)*100;
	if(mode==2) return /*sigma angolo da prop. errori di 0.3% di errore su beta*/ (2*sigma/(100000*1/sqrt(1-pow(1/(1.2*0.98),2))*((1/1.2)/(0.98*0.98))*0.003))*100  ;
	if(mode==3) return ((1.15694e-03*(1+sigma/100.)-1.36634e-03)/1.36634e-03)*100;
}



float TemplateFIT::SmearBetaRICH(float Beta, float stepsigma, float stepshift){
	float angle;
	angle= acos(1/(1.2*Beta))*10e4;
	float shiftstart=-systpar.shift;
	angle = angle + (shiftstart+(2*systpar.shift/(float)systpar.steps)*stepshift) + Rand->Gaus(0,(float)((2*systpar.sigma/systpar.steps)*stepsigma));
	return 1/(1.2*cos(angle/10e4));
}

float TemplateFIT::SmearBetaRICH_v2(float Beta, float Beta_gen, float stepsigma, float stepshift){
	if(stepsigma==0 && stepshift==5) return Beta;
	float delta = (1/Beta - 1/Beta_gen)/1.36634e-03;
	float sigma = 1.15694e-03*(1-systpar.sigma/100.);//1.15694e-03*(1-systpar.sigma/100.);
	sigma+=(float)((2*fabs(1.15694e-03-sigma)/(float)systpar.steps)*stepsigma);
	float shift = 0;//- systpar.shift*1e-5 + (2*systpar.shift*1e-5/(float)systpar.steps)*stepshift;		
	return 1/(1/Beta_gen +sigma*delta - shift );
}

float TemplateFIT::SmearRRICH_R(float R, float R_gen, float stepsigma, float stepshift){
	if(stepsigma==0 && stepshift==5) return R;
	float delta = (1/R - 1/R_gen)/9.16394e-03;
	float sigma = 9.16394e-03*(1-systpar.shift/100.);//1.15694e-03*(1-systpar.sigma/100.);
	sigma+=(float)((2*fabs(9.16394e-03-sigma)/(float)systpar.steps)*stepshift);
	float shift = 0;//- systpar.shift*1e-5 + (2*systpar.shift*1e-5/(float)systpar.steps)*stepshift;		
	return 1/(1/R_gen +sigma*delta - shift );
}


float TemplateFIT::SmearBeta(float Beta,float M_gen, float stepsigma, float stepshift,float R){

	float time = 1.2/(Beta*3e-4);
	float shiftstart=-systpar.shift/M_gen;

	float sigma = ((2*systpar.sigma/systpar.steps)*stepsigma);
	float shift = (shiftstart+(2*systpar.shift/(float)systpar.steps)*stepshift/M_gen -7);
	float smeartime;
	float norm = Rand->Rndm();
	if (norm<=2024.44/3165.1211) 
		smeartime = shift + Rand->Gaus(0,sigma);
	else if (norm<=((2024.44+1136.55)/3165.1211)) 
		smeartime = shift*1.001 + Rand->Gaus(0,1.72834*sigma);
	else 
		smeartime = shift*1.02747 +Rand->Gaus(0,3.84716*sigma);

	time = time + smeartime;
	return 1.2/(time*3e-4);
}

float TemplateFIT::SmearBeta_v2(float Beta, float Beta_gen,float M_gen, TF1* MCmodel,  float stepsigma, float stepshift){

	if(stepsigma==0 && stepshift==5) return Beta;



	float beta_inv = 1/Beta;
	float sigmasmear = 0.0284/2;
	float shiftstart=-systpar.shift/M_gen;
	float sigmastart = sigmasmear - sigmasmear*systpar.sigma/100;
	float sigma = sigmastart +  ((2*systpar.sigma/(systpar.steps))*stepsigma)/100*sigmasmear;
	float shift = (- systpar.shift*1e-4 + (2*systpar.shift*1e-4/(float)systpar.steps)*stepshift)/M_gen;	
	float norm = Rand->Rndm();
	
	float MCsigma = (0.0318/0.0335)*MCmodel->Eval(1/Beta_gen);

	if (norm<=0.6667) {
		beta_inv = beta_inv + shift + Rand->Gaus(0,sigma);
		}
	else if (norm<=(0.6667+0.332)) {
		float sigma2 = sqrt(pow(1.25*sigma,2)+(1.25*1.25-1)*pow(MCsigma,2));
		beta_inv = beta_inv + shift + Rand->Gaus(0,sigma2);
		}
	else  {
		float sigma3 = sqrt(pow(2.29*sigma,2)+(2.29*2.29-1)*pow(MCsigma,2));
		beta_inv = beta_inv + shift + Rand->Gaus(0.02747,sigma3);
		}
//	cout<<"Smearing: "<<stepshift<<" "<<stepsigma<<" MC sigma: "<<MCsigma<< " :"<<shift<<" "<<sigma<<": "<<Beta<<" "<<1/beta_inv<<endl;
	return 1/(beta_inv);
		
}



void TemplateFIT::FillEventByEventMC(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars)){


	std::string cutP=cut+"&IsProtonMC";
	std::string cutD=cut+"&IsPureDMC";
	std::string cutHe=cut+"&IsProtonMC";
	
	if((ApplyCuts(cutP,vars)||ApplyCuts(cutD,vars)||ApplyCuts(cutHe,vars))){
	for(int i=0;i<systpar.steps;i++)
			for(int j=0;j<systpar.steps;j++){
		
				
				float betasmear=0;
				float Rsmear = vars->R;

				if(ApplyCuts("IsOnlyFromToF",vars)) 	
					if(systpar.mode==0) betasmear = SmearBeta(vars->Beta,vars->Massa_gen,(float)i,(float)j,vars->R);    
					if(systpar.mode==1) betasmear = SmearBeta_v2(vars->Beta,GetBetaGen(vars),vars->Massa_gen,MCmodel,(float)i,(float)j);    
				else if(ApplyCuts("IsFromNaF",vars))	 
					if(systpar.mode==2) betasmear = SmearBetaRICH(vars->BetaRICH_new,(float)i,(float)j);
                                        if(systpar.mode==3) { 
						betasmear = SmearBetaRICH_v2(vars->BetaRICH_new,GetBetaGen(vars),(float)i,(float)j); 
						Rsmear = SmearRRICH_R(vars->R,GetGenMomentum(vars),(float)i,(float)j);
						}
				else if(ApplyCuts("IsFromAgl",vars))  	
					if(systpar.mode==2) betasmear = SmearBetaRICH(vars->BetaRICH_new,(float)i,(float)j);
                                        if(systpar.mode==3) {
						betasmear = SmearBetaRICH_v2(vars->BetaRICH_new,GetBetaGen(vars),(float)i,(float)j); 
						Rsmear = SmearRRICH_R(vars->R,GetGenMomentum(vars),(float)i,(float)j);
						}
					
				float mctotalweight = 
						   vars->mcweight 
						   * vars->GetCutoffCleaningWeight(GetRFromBeta(bins.getParticle().getMass(),betasmear),vars->Momento_gen,BetacutoffCut) 
						   * vars->GetTimeDepWeight(vars->R)
						   ;
 
				int kbin;
	
				if(!IsFitNoise){
				   if(bins.IsUsingBetaEdges()){	 
					kbin = bins.GetBin(betasmear);
					float mass=0;
					mass = Rsmear/betasmear * pow((1-pow(betasmear,2)),0.5);
					if(ApplyCuts((cutP+"&RigSafetyCut").c_str(),vars)&&kbin>=0) {  fits[kbin][i][j]->Templ_P->Fill(mass,mctotalweight);}
					if(ApplyCuts((cutD+"&RigSafetyCut_D").c_str(),vars)&&kbin>=0)  {  fits[kbin][i][j]->Templ_D->Fill(mass,mctotalweight);}
					if(ApplyCuts(cutHe,vars)&&kbin>=0) fits[kbin][i][j]->Templ_He->Fill((2.793/0.938)*mass,mctotalweight);
					}
				    else {
				   	kbin = bins.GetBin(discr_var(vars));
					if(ApplyCuts((cutP+"&RigSafetyCut").c_str(),vars)&&kbin>=0)   fits[kbin][i][j]->Templ_P->Fill(betasmear,mctotalweight);
                                        if(ApplyCuts((cutD).c_str(),vars)&&kbin>=0)   fits[kbin][i][j]->Templ_D->Fill(betasmear,mctotalweight);
                                        if(ApplyCuts(cutHe,vars)&&kbin>=0) fits[kbin][i][j]->Templ_He->Fill(betasmear,mctotalweight);

				   }	

				}
				else{
					if(bins.IsUsingBetaEdges()) {
						kbin = bins.GetBin(betasmear);
						float mass = Rsmear/betasmear * pow((1-pow(betasmear,2)),0.5);		

						if(ApplyCuts((cutP+"&RigSafetyCut").c_str(),vars)&&kbin>=0)  fits[kbin][i][j]->Templ_P->Fill(mass,mctotalweight);		
						if(ApplyCuts((cutD+"&RigSafetyCut_D").c_str(),vars)&&kbin>=0) fits[kbin][i][j]->Templ_D->Fill(mass,vars->mcweight);
						if(ApplyCuts(cutHe,vars)&&kbin>=0) fits[kbin][i][j]->Templ_He->Fill((2.793/0.938)*mass,vars->mcweight);

						float betabad = betasmear;
						if(BadEvSim) {betabad=BadEvSim->SimulateBadEvents(betasmear); 
							kbin = bins.GetBin(betabad);
						}
						float mass_bad = Rsmear/betabad * pow((1-pow(betabad,2)),0.5);
						if(ApplyCuts((cutP+"&RigSafetyCut").c_str(),vars)&&kbin>=0)  fits[kbin][i][j]->Templ_Noise->Fill(mass_bad,mctotalweight);
					}					
					else {
						cout<<"ECCO"<<endl;
						if(BadEvSim) betasmear=BadEvSim->SimulateBadEvents(betasmear);
						kbin = bins.GetBin(discr_var(vars));
						if(ApplyCuts((cutP+"&RigSafetyCut").c_str(),vars)&&kbin>=0)  fits[kbin][i][j]->Templ_P->Fill(betasmear,mctotalweight);
	                                        if(ApplyCuts((cutD).c_str(),vars)&&kbin>=0)  fits[kbin][i][j]->Templ_D->Fill(betasmear,vars->mcweight);
        	                                if(ApplyCuts(cutHe,vars)&&kbin>=0) fits[kbin][i][j]->Templ_He->Fill(betasmear,vars->mcweight);
					}
				}
			}
	}
	return;	
}



void TemplateFIT::Save(){

	for(int bin=0;bin<bins.size();bin++){ 
		finalhistos.Add(fits[bin][0][5]->Data);
		finalhistos.Add(fits[bin][0][5]->DataPrim);
		finalhistos.Add(fits[bin][0][5]->DataAmbient);


		finalhistos.writeObjsInFolder((basename + "/Bin "+ to_string(bin)+"/Data").c_str(),false);
	}

	for(int bin=0;bin<bins.size();bin++){ 
		for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){

				finalhistos.Add(fits[bin][i][j]->Templ_P);
		}
		finalhistos.writeObjsInFolder((basename + "/Bin "+ to_string(bin)+"/TemplateP").c_str(),false);
	}

	for(int bin=0;bin<bins.size();bin++){ 
		for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){

				finalhistos.Add(fits[bin][i][j]->Templ_D);
		}
		finalhistos.writeObjsInFolder((basename + "/Bin "+ to_string(bin)+"/TemplateD").c_str(),false);
	}

	for(int bin=0;bin<bins.size();bin++){ 
		for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){

				finalhistos.Add(fits[bin][i][j]->Templ_He);
		}
		finalhistos.writeObjsInFolder((basename + "/Bin "+ to_string(bin)+"/TemplateHe").c_str(),false);
	}

	for(int bin=0;bin<bins.size();bin++){ 
		for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){

				finalhistos.Add(fits[bin][i][j]->Templ_Noise);
		}
		finalhistos.writeObjsInFolder((basename + "/Bin "+ to_string(bin)+"/TemplateNoise").c_str(),false);
	}

	return;
}

void TemplateFIT::SaveFitResults(FileSaver finalhistos){

	
	//save best
	for(int bin=0;bin<bins.size();bin++){
		TH1F * OriginalP=(TH1F*)fits[bin][0][5]->Templ_P->Clone();

		TH1F* BestP;
		if(!IsLocalFit) {
				BestP=(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_P->Clone();
				BestP->SetName(("Best #chi^{2}: " + to_string(BestChiSquare->i) + "_" + to_string(BestChiSquare->j)).c_str());
		}	
		else {
			BestChi * bestlocal = new BestChi();
			bestlocal->FindMinimum(TFitChisquare[bin]);
			BestP=(TH1F*)fits[bin][bestlocal->i][bestlocal->j]->Templ_P->Clone();	
			BestP->SetName(("Best #chi^{2}: " + to_string(bestlocal->i) + "_" + to_string(bestlocal->j)).c_str());
	
		}
		TH1F * BestPPrim =(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_PPrim->Clone();
		OriginalP->SetName("Original Proton MC ");
		BestPPrim->SetName("Best #chi^{2} Overcutoff Proton MC ");

		finalhistos.Add(OriginalP);
		finalhistos.Add(BestP);
		finalhistos.Add(BestPPrim);
		finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesP/Bin"+to_string(bin)).c_str());

		TH1F * OriginalD=(TH1F*)fits[bin][0][5]->Templ_D->Clone();
		TH1F * BestD;
		if(!IsLocalFit){
				 BestD=(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_D->Clone();
				BestD->SetName(("Best #chi^{2}: " + to_string(BestChiSquare->i) + "_" + to_string(BestChiSquare->j)).c_str());
			}
		else {
			BestChi * bestlocal = new BestChi();
			bestlocal->FindMinimum(TFitChisquare[bin]);
			BestD=(TH1F*)fits[bin][bestlocal->i][bestlocal->j]->Templ_D->Clone();	
			BestD->SetName(("Best #chi^{2}: " + to_string(bestlocal->i) + "_" + to_string(bestlocal->j)).c_str());
		}
		TH1F * BestDPrim =(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_DPrim->Clone();
		OriginalD->SetName("Original Deuton MC ");
		BestDPrim->SetName("Best #chi^{2} Overcutoff Deuton MC ");
		finalhistos.Add(OriginalD);
		finalhistos.Add(BestD);
		finalhistos.Add(BestDPrim);
		finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesD/Bin"+to_string(bin)).c_str());

		TH1F * OriginalHe=(TH1F*)fits[bin][0][5]->Templ_He->Clone();
		TH1F * BestHe=(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_He->Clone();
		TH1F * BestHePrim =(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_HePrim->Clone();
		OriginalHe->SetName("Original Tritium MC ");
		BestHe->SetName("Best #chi^{2} Tritium MC ");
		BestHePrim->SetName("Best #chi^{2} Overcutoff Triton MC ");
		finalhistos.Add(OriginalHe);
		finalhistos.Add(BestHe);
		finalhistos.Add(BestHePrim);
		finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesHe/Bin"+to_string(bin)).c_str());

		TH1F * OriginalNoise=(TH1F*)fits[bin][0][5]->Templ_Noise->Clone();
		TH1F * BestNoise=(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_Noise->Clone();
		OriginalNoise->SetName("Original Noise Template MC ");
		BestNoise->SetName("Best #chi^{2} Noise Template MC ");
		finalhistos.Add(OriginalNoise);
		finalhistos.Add(BestNoise);
		finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesNoise/Bin"+to_string(bin)).c_str());

	}


	//save data
	for(int bin=0;bin<bins.size();bin++){ 
		finalhistos.Add(fits[bin][0][5]->Data);
		finalhistos.Add(fits[bin][0][5]->DataPrim);
		finalhistos.Add(fits[bin][0][5]->DataAmbient);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/Data/Bin"+to_string(bin)).c_str());	
	}
	
	//save all
	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(!(i==0&&i==3)) if(fits[bin][i][j]->Templ_P) finalhistos.Add(fits[bin][i][j]->Templ_P);
			}	
	finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesP/Bin"+to_string(bin)).c_str());
	}

	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(!(i==0&&i==3))if(fits[bin][i][j]->Templ_D) finalhistos.Add(fits[bin][i][j]->Templ_D);
				}
	finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesD/Bin"+to_string(bin)).c_str());
	}
	
	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(!(i==0&&i==3))if(fits[bin][i][j]->Templ_He) finalhistos.Add(fits[bin][i][j]->Templ_He);
				}
	finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesHe/Bin"+to_string(bin)).c_str());
	}

	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(!(i==0&&i==3)) if(fits[bin][i][j]->Templ_Noise) finalhistos.Add(fits[bin][i][j]->Templ_Noise);
				}
	finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesNoise/Bin"+to_string(bin)).c_str());
	}

	//save fits
	for(int bin=0;bin<bins.size();bin++){
		if(BestChiSquare&&fits[bin][BestChiSquare->i][BestChiSquare->j]->Tfit){
			TH1F * FIT;
			if(fits[bin][BestChiSquare->i][BestChiSquare->j]->Tfit_outcome!=-1){ 
				FIT=(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Tfit-> GetPlot();
				FIT->SetName(("Fraction Fit bin" + to_string(bin)).c_str());
				finalhistos.Add(FIT);
				finalhistos.writeObjsInFolder((basename+"/Fit Results/FractionFits/Bin"+to_string(bin)).c_str());
			}
		}
	}

	for(int i=0;i<bins.size();i++)
		if(DCountsSpread[i]) 
		{		
			finalhistos.Add(DCountsSpread[i]);
		}
	finalhistos.writeObjsInFolder((basename+"/Fit Results/Spreads/DCounts").c_str());
	
	for(int i=0;i<bins.size();i++)
		if(TFitChisquare[i]) { 
			TFitChisquare[i]->Scale(fits[i][BestChiSquare->i][BestChiSquare->j]->ndf);	
			finalhistos.Add(TFitChisquare[i]);
		}
		finalhistos.writeObjsInFolder((basename+"/Fit Results/Spreads/ChiSquare").c_str());	

	for(int i=0;i<bins.size();i++)
		if(WeightedDCounts[i]) finalhistos.Add(WeightedDCounts[i]);
		finalhistos.writeObjsInFolder((basename+"/Fit Results/Spreads/Weighted D counts").c_str());	
	

	if(HeContModel) finalhistos.Add(HeContModel);
	if(HeContError) finalhistos.Add(HeContError);
	if(MeasuredHeContRatio) finalhistos.Add(MeasuredHeContRatio);
	if(MCHeContRatio) finalhistos.Add(MCHeContRatio);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/HeliumContamination").c_str());	

	if(Global_ChiSquare) {
		Global_ChiSquare->Scale(fits[0][BestChiSquare->i][BestChiSquare->j]->ndf);	
		finalhistos.Add(Global_ChiSquare);
		TH1D * BestSigma = (TH1D*) Global_ChiSquare->ProjectionX("Best Sigma",BestChiSquare->j+1,BestChiSquare->j+1);
		finalhistos.Add(BestSigma);
		TH1D * BestShift = (TH1D*) Global_ChiSquare->ProjectionY("Best Shift",BestChiSquare->i+1,BestChiSquare->i+1);
		finalhistos.Add(BestShift);
		BestChiSquares->Scale(fits[0][BestChiSquare->i][BestChiSquare->j]->ndf);
		OriginalChiSquares->Scale(fits[0][BestChiSquare->i][BestChiSquare->j]->ndf);
		finalhistos.Add(BestChiSquares   ); 
		finalhistos.Add(OriginalChiSquares);
	}

	finalhistos.Add(StatErrorP);
	finalhistos.Add(StatErrorD);
	finalhistos.Add(StatErrorT);
	finalhistos.Add(SystError);
	finalhistos.Add(ProtonCounts);
	finalhistos.Add(DeuteronCounts);
	finalhistos.Add(TritiumCounts);
	finalhistos.Add(ProtonCountsPrim);
	finalhistos.Add(DeuteronCountsPrim);

	finalhistos.writeObjsInFolder((basename+"/Fit Results/").c_str());

}


void TemplateFIT::SumUpMassDistrib(FileSaver finalhistos){

	TH1F * SummedMass = (TH1F*)fits[0][0][5]->Data ->Clone(); 
	for(int bin=0; bin<bins.size()-2;bin++){
		SummedMass->Add(fits[bin][0][5]->Data);
	}
	finalhistos.Add(SummedMass);
	finalhistos.writeObjsInFolder((basename + "/SummedData").c_str());
	return;
}

void TemplateFIT::SetFitRangeByQuantiles(float quant_min,float quant_max){
	for(int bin=0;bin<bins.size();bin++) 
		for(int i=0;i<systpar.steps;i++) 
			for(int j=0;j<systpar.steps;j++){
				Double_t q[2]= {quant_min,quant_max};
				Double_t x[2];
				fits[bin][i][j]->Data->GetQuantiles(2,x,q);
				fits[bin][i][j]->fitrangemin = x[0];
				fits[bin][i][j]->fitrangemax = x[1];
			}
	return;
}


void Pre_Scale(TH1F * Data, TH1F * PHisto,TH1F * DHisto, TH1F * HeHisto){
	if(DHisto>0&&PHisto>0) DHisto->Scale(0.02*PHisto->Integral()/DHisto->Integral());
	if(HeHisto>0&&PHisto>0) HeHisto->Scale(0.005*PHisto->Integral()/HeHisto->Integral());
	return;
}

void CleanTemplates(TH1F * Histo){
	for(int i=0;i<Histo->GetNbinsX();i++)
		if(Histo->GetBinContent(i+1)>0&&Histo->GetBinContent(i)==0&&Histo->GetBinContent(i+1)==0&&Histo->GetBinContent(i+2)==0){
			Histo->SetBinContent(i+1,0);
			Histo->SetBinError(i+1,0);
		}
}

float CalculateAmountOfHighMassComponent(TH1F * Data, TH1F* First, TH1F * Third, float massref, float base_constrain){
		TH1F * third = (TH1F*) Third->Clone();
		int bin_norm = Data->FindBin(massref); 
		float datapeak = (Data->GetBinContent(bin_norm-2)+Data->GetBinContent(bin_norm-1)+Data->GetBinContent(bin_norm)+Data->GetBinContent(bin_norm+1)+Data->GetBinContent(bin_norm+2))/5;
		third->Scale(datapeak/Third->GetBinContent(bin_norm));
		third->Sumw2();
		cout<<"Component: "<<third->Integral()/Data->Integral()<<endl;
		if(third->Integral()/Data->Integral()>0){
		if(First->GetBinContent(bin_norm)==0||First->GetBinContent(bin_norm+1)==0||First->GetBinContent(bin_norm-1)==0)
			return third->Integral()/Data->Integral();
		}
		else 
		return base_constrain;
}

float GetChiSquare(TH1 * result, TH1 * data, float min, float max){

	TH1 * Result = (TH1*) result->Clone();
	TH1 * Data   = (TH1*) data->Clone();

	int binmin = Data->FindBin(min);
	int binmax = Data->FindBin(4.5);
	float chi = 0;
	float err = 0;
	int ndf =binmax-binmin;
	for(int i=binmin; i<binmax; i++){
		if(Data->GetBinContent(i+1)>=1&&Result->GetBinContent(i+1)>=0.4){
			err = pow(pow(Data->GetBinError(i+1),2) + pow(Result->GetBinError(i+1),2),0.5);
			chi += pow((Data->GetBinContent(i+1) - Result->GetBinContent(i+1)),2)/pow(err,2);
		}
	}
	ndf = ndf-3;
	chi /= ndf;
	cout<<"chi: "<<chi<<" ndf: "<<ndf<<endl;
	return chi;

}


void Do_TemplateFIT(TFit * Fit,float fitrangemin,float fitrangemax,float constrain_min[], float constrain_max[], bool isfitnoise, bool highmasstailconstrain, bool IsFitPrim){
	if(isfitnoise) cout<<"********** FIT NOISE MODE **********"<<endl;
	if(isfitnoise) Fit->Templ_Noise->Scale(0.01*Fit ->  Templ_P->Integral()/Fit ->  Templ_Noise->Integral());
	TObjArray *Tpl;
	Tpl = new TObjArray(3);


	if(!IsFitPrim) {
		if(Fit ->  Templ_P)  if(Fit ->  Templ_P ->GetEntries()>0) Tpl -> Add( Fit ->  Templ_P );
		if(Fit ->  Templ_D)  if(Fit ->  Templ_D ->GetEntries()>0) Tpl -> Add( Fit ->  Templ_D );
        	if(Fit ->  Templ_He) if(Fit ->  Templ_He->GetEntries()>0) Tpl -> Add( Fit ->  Templ_He);
	}
	else {
		if(Fit ->  Templ_PPrim)  if(Fit ->  Templ_PPrim ->GetEntries()>0) Tpl -> Add( Fit ->  Templ_PPrim );
		if(Fit ->  Templ_DPrim)  if(Fit ->  Templ_DPrim ->GetEntries()>0) Tpl -> Add( Fit ->  Templ_DPrim );
        	if(Fit ->  Templ_HePrim) if(Fit ->  Templ_HePrim->GetEntries()>0) Tpl -> Add( Fit ->  Templ_HePrim);
	}


	if(isfitnoise) if(Fit ->  Templ_Noise) if(Fit ->  Templ_Noise->GetEntries()>50) Tpl -> Add( Fit ->  Templ_Noise);

	float min=fitrangemin;
	float max=fitrangemax;
	cout<<	Fit -> Data<<" "<<Fit -> Templ_P<<" "<<Fit -> Templ_D<<" "<<endl;

	bool fitcondition = (Fit -> Data->GetEntries()>100)&&(Fit -> Templ_P->GetEntries()>100);
	cout<<"fit"<<endl;
	cout<<	Fit -> Data->GetEntries()<<" "<<Fit -> Templ_P->GetEntries()<<" "<<Fit -> Templ_D->GetEntries()<<endl;

	if(fitcondition){	
		cout<<"Conditions for fit OK!"<<endl;	
		
		if(!IsFitPrim) Fit -> Tfit = new TFractionFitter(Fit -> Data, Tpl ,"q");
		else   	Fit -> Tfit = new TFractionFitter(Fit -> DataPrim, Tpl ,"q");
		int n = Fit -> Data->Integral() / Fit -> Templ_P->Integral();
		Fit -> Templ_P-> Scale(n);
		Fit -> Templ_D-> Scale(n);
		Fit -> Templ_He-> Scale(n);

		cout<<"fitting..."<<endl;			
		cout<<"Constrainsts: "<<endl;
		cout<<constrain_min[0]<<" "<<constrain_max[0]<<endl;	
		cout<<constrain_min[1]<<" "<<constrain_max[1]<<endl;	
		cout<<constrain_min[2]<<" "<<constrain_max[2]<<endl;	
		Fit -> Tfit -> SetRangeX(Fit -> Data -> FindBin(min), Fit -> Data -> FindBin(max));
		
		Fit -> Tfit -> Constrain(1, constrain_min[0] ,constrain_max[0]);
                Fit -> Tfit -> Constrain(2, constrain_min[1] ,constrain_max[1]);
	 	Fit -> Tfit -> Constrain(3, constrain_min[2] ,constrain_max[2]);
		if(highmasstailconstrain) 
			Fit -> Tfit -> Constrain(3, 0.1*CalculateAmountOfHighMassComponent(Fit -> Data,Fit ->  Templ_P,Fit ->  Templ_Noise,3.3,constrain_min[2]),1.2*CalculateAmountOfHighMassComponent(Fit -> Data,Fit ->  Templ_P,Fit ->  Templ_Noise,3.1,constrain_max[2]));	 
	        if(isfitnoise&&highmasstailconstrain) 
			Fit -> Tfit -> Constrain(4, 0.1*CalculateAmountOfHighMassComponent(Fit -> Data,Fit ->  Templ_P,Fit ->  Templ_Noise,4.5,0.001),1.2*CalculateAmountOfHighMassComponent(Fit -> Data,Fit ->  Templ_P,Fit ->  Templ_Noise,4.5,0.001));	 


		//fitting
		if(Fit -> Tfit ) 
		{
			Fit -> Tfit_outcome = -1;
			try { Fit -> Tfit_outcome = Fit -> Tfit -> Fit(); } catch(const std::invalid_argument& e) 
			{ cout << "Failed TFractionFItting. Giving up.\n"; }
		}
		for(int fit_attempt=0; fit_attempt<20; fit_attempt++) {
			cout<<"fit attempt: "<<fit_attempt<<endl;
			if(Fit -> Tfit_outcome == 0) break;
			else {
				cout<<fit_attempt<<endl;
				Fit -> Tfit -> SetRangeX(Fit -> Data -> FindBin((1+0.01*fit_attempt)*min), Fit -> Data -> FindBin((1+0.01*fit_attempt)*max));

				Fit -> Tfit_outcome = -1;
				try { Fit -> Tfit_outcome = Fit -> Tfit -> Fit(); } catch(const std::invalid_argument& e) 
				{ cout << "Failed TFractionFItting. Giving up.\n"; }
			}
		}
		//

		if(Fit -> Tfit_outcome==0){
			TH1F * Result = (TH1F *) Fit-> Tfit -> GetPlot();
			float itot= Result->Integral();
			double w1,e1 = 0;
			double w2,e2 = 0;
			double w3,e3 = 0;
			double w4,e4 = 0;

			Fit -> Tfit ->GetResult(0,w1,e1);
			Fit -> Tfit ->GetResult(1,w2,e2);
			if(Fit-> Templ_He) Fit -> Tfit ->GetResult(2,w3,e3);
			if(isfitnoise) Fit -> Tfit ->GetResult(3,w4,e4);

			float i1=1;
			float i2=1;
			float i3=1;
			float i4=1;

			if(!IsFitPrim){
				 i1 = Fit-> Templ_P  ->Integral(Fit->Templ_P -> FindBin(min), Fit->Templ_P -> FindBin(max));
				 i2 = Fit-> Templ_D  ->Integral(Fit->Templ_D -> FindBin(min), Fit->Templ_D -> FindBin(max));
				if(Fit-> Templ_He) i3 = Fit-> Templ_He ->Integral(Fit->Templ_He -> FindBin(min), Fit->Templ_He -> FindBin(max));
				if(isfitnoise) i4 = Fit-> Templ_Noise ->Integral(Fit->Templ_Noise -> FindBin(min), Fit->Templ_Noise -> FindBin(max));
			}

			else{
				i1 = Fit-> Templ_PPrim  ->Integral(Fit->Templ_P -> FindBin(min), Fit->Templ_P -> FindBin(max));
				i2 = Fit-> Templ_DPrim  ->Integral(Fit->Templ_D -> FindBin(min), Fit->Templ_D -> FindBin(max));
				if(Fit-> Templ_HePrim) i3 = Fit-> Templ_HePrim ->Integral(Fit->Templ_HePrim -> FindBin(min), Fit->Templ_HePrim -> FindBin(max));
				if(isfitnoise) i4 = Fit-> Templ_Noise ->Integral(Fit->Templ_Noise -> FindBin(min), Fit->Templ_Noise -> FindBin(max));
			}

	

			Fit ->ContribP= w1;
			Fit ->ContribD= w2;
			Fit ->ContribHe=w3;
			Fit ->ContribNoise=w4;

			Fit ->errP= e1;
			Fit ->errD= e2;
			Fit ->errHe=e3;
			Fit ->errNoise=e4;

			Fit ->wheightP= w1*itot/i1;
			Fit ->wheightD= w2*itot/i2;
			Fit->wheightHe= w3*itot/i3;
			Fit->wheightNoise= w4*itot/i4;

			float Cov00 = Fit -> Tfit->GetFitter()->GetCovarianceMatrixElement(0,0);
			float Cov11 = Fit -> Tfit->GetFitter()->GetCovarianceMatrixElement(1,1);
			float Cov22 = Fit -> Tfit->GetFitter()->GetCovarianceMatrixElement(2,2);
			float Cov01 = Fit -> Tfit->GetFitter()->GetCovarianceMatrixElement(0,1);
			float Cov02 = Fit -> Tfit->GetFitter()->GetCovarianceMatrixElement(0,2);
			float Cov12 = Fit -> Tfit->GetFitter()->GetCovarianceMatrixElement(1,2);
		
			cout<<endl;
			cout<<"Cov Matrix:"<<endl;
			cout<<Cov00<<" "<<Cov01<<" "<<Cov02<<endl;
			cout<<Cov01<<" "<<Cov11<<" "<<Cov12<<endl;
			cout<<Cov02<<" "<<Cov12<<" "<<Cov22<<endl;

			float df1dw1 = (w2+w3)/pow(w1+w2+w3,2);
                        float df1dw2 = (-w1)/pow(w1+w2+w3,2);
                        float df1dw3 = (-w1)/pow(w1+w2+w3,2);

                        float df2dw1 = (-w2)/pow(w1+w2+w3,2);
                        float df2dw2 = (w1+w3)/pow(w1+w2+w3,2);
                        float df2dw3 = (-w2)/pow(w1+w2+w3,2);

                        float df3dw1 = (-w3)/pow(w1+w2+w3,2);
                        float df3dw2 = (-w3)/pow(w1+w2+w3,2);
                        float df3dw3 = (w1+w2)/pow(w1+w2+w3,2);
			
			cout<<"Derivatives:"<<endl;
			cout<<df1dw1<<" "<<df1dw2<<" "<<df1dw3<<endl;
			cout<<df2dw1<<" "<<df2dw2<<" "<<df2dw3<<endl;
			cout<<df3dw1<<" "<<df3dw2<<" "<<df3dw3<<endl;

			cout<<endl;

	

			if(!IsFitPrim){
				Fit ->  Templ_P  -> Scale(Fit ->wheightP);
				Fit ->  Templ_D  -> Scale(Fit ->wheightD);
				Fit ->  Templ_He  -> Scale(Fit ->wheightHe);
			}
			else{
				Fit ->  Templ_PPrim  -> Scale(Fit ->wheightP);
				Fit ->  Templ_DPrim  -> Scale(Fit ->wheightD);
				Fit ->  Templ_HePrim  -> Scale(Fit ->wheightHe);
			}

			if(isfitnoise) Fit ->  Templ_Noise -> Scale(Fit ->wheightNoise);

			TH1F * Sum = (TH1F *)Fit ->  Templ_P->Clone();
			Sum -> Add(Fit ->  Templ_D);
			Sum -> Add(Fit ->  Templ_He);
			if(isfitnoise) Sum -> Add(Fit ->  Templ_Noise);			
			for(int n=0;n<Sum->GetNbinsX();n++) if(Sum->GetBinContent(n+1)==0) Sum->SetBinError(n+1,1);
	
			if(( df1dw1*df1dw1*Cov00 + df1dw2*df1dw2*Cov11 + df1dw3*df1dw3*Cov22 +2*df1dw1*df1dw2*Cov01 +2*df1dw1*df1dw3*Cov02+ 2*df1dw2*df1dw3*Cov12)>0){
				Fit -> StatErrP =  pow(( df1dw1*df1dw1*Cov00 + df1dw2*df1dw2*Cov11 + df1dw3*df1dw3*Cov22 +2*df1dw1*df1dw2*Cov01 +2*df1dw1*df1dw3*Cov02+ 2*df1dw2*df1dw3*Cov12)/2,0.5);
				Fit -> StatErrD =  pow(( df2dw1*df2dw1*Cov00 + df2dw2*df2dw2*Cov11 + df2dw3*df2dw3*Cov22 +2*df2dw1*df2dw2*Cov01 +2*df2dw1*df2dw3*Cov02+ 2*df2dw2*df2dw3*Cov12)/2,0.5);
				Fit -> StatErrT =  pow(( df3dw1*df3dw1*Cov00 + df3dw2*df3dw2*Cov11 + df3dw3*df3dw3*Cov22 +2*df3dw1*df3dw2*Cov01 +2*df3dw1*df3dw3*Cov02+ 2*df3dw2*df3dw3*Cov12)/2,0.5);
			}
			else{

				Fit -> StatErrP = e1;
				Fit -> StatErrD = e2;			
				Fit -> StatErrT = e3;
			}
			cout<<"fract: "<<w1<<" "<<w2<<" "<<w3<<endl;
			cout<<endl;
			cout<<"err: "  <<e1<<" "<<e2<<" "<<e3<<endl;
			cout<<"er c:"  <<Fit -> StatErrP<<" "<<Fit -> StatErrD<<" "<<Fit -> StatErrT<<endl;	
			
			//Fit -> ChiSquare = Fit -> Tfit -> GetChisquare()/(float) (Fit ->  Tfit -> GetNDF());
			Fit -> ChiSquare = GetChiSquare(Sum,Fit->Data,min,max);
			Fit->ndf = Fit ->  Data->FindBin(4.5)-Fit ->Data->FindBin(0.85) -3;
			Fit -> DCounts = Fit ->  Templ_D -> Integral();
			Fit -> PCounts = Fit ->  Templ_P -> Integral();
			Fit -> TCounts = Fit ->  Templ_He -> Integral();
			cout<<"HOFINITO: "<<Fit ->  Templ_DPrim<<" "<<Fit ->  Templ_PPrim<<" ndf: "<<Fit->ndf<<endl;
			Fit -> DCountsPrim = Fit ->  Templ_DPrim -> Integral();
			Fit -> PCountsPrim = Fit ->  Templ_PPrim -> Integral();
		

	}
	
		else{
			cout<<"Fit not converged: returning original templates"<<endl; 
			Fit -> PCounts = Fit -> Data -> Integral();
			Fit -> DCounts = Fit -> Data -> Integral();
			Fit -> PCountsPrim = Fit -> Data -> Integral();
			Fit -> DCountsPrim = Fit -> Data -> Integral();
		
			Fit ->wheightP= 1;
			Fit ->wheightD= 1;
			Fit ->wheightHe= 1;
			Fit ->StatErrP= 0.5;
			Fit ->StatErrD= 0.5;
			Fit ->StatErrT= 0.5;
	

		}
	}
	else{
		cout<<"Fit conditions not OK: returning original templates"<<endl; 
		Fit -> PCounts = Fit -> Data -> Integral();
		Fit -> DCounts = Fit -> Data -> Integral();
		Fit -> PCountsPrim = Fit -> Data -> Integral();
		Fit -> DCountsPrim = Fit -> Data -> Integral();
		
		Fit ->wheightP= 1;
		Fit ->wheightD= 1;
		Fit ->wheightHe= 1;
		Fit ->StatErrP= 0.5;
		Fit ->StatErrD= 0.5;
		Fit ->StatErrT= 0.5;


	} 
	return;
}

float EvalFitProbability(float chi){
	TF1 * CumulativeChiSquare=new TF1("Chi","exp(-x)*x^-0.5",0.001,1000);
	if(chi>0.01&&chi<100)
		return CumulativeChiSquare->Integral(chi,100);
	else return 0;
	//if(chi>4) return 0;
	//else return 1;
}

void  TemplateFIT::RebinAll(int f){
	for(int bin=0;bin<bins.size();bin++){
	 for(int sigma=0;sigma<systpar.steps;sigma++)
		 for(int shift=0;shift<systpar.steps;shift++){
			cout<<bin<<" "<<sigma<<" "<<shift<<" "<< fits[bin][sigma][shift] -> DataPrim->GetNbinsX()<<endl;
			 fits[bin][sigma][shift] ->  Templ_P->Rebin(f);
			 fits[bin][sigma][shift] ->  Templ_D->Rebin(f);
			 fits[bin][sigma][shift] ->  Templ_He->Rebin(f);
			 fits[bin][sigma][shift] -> Data->Rebin(f);
			 fits[bin][sigma][shift] -> DataPrim->Rebin(f);			
			 fits[bin][sigma][shift] -> DataAmbient->Rebin(f);			
		
		 }
	}
};





void TemplateFIT::ExtractCounts(FileSaver finalhistos,int force_shift){

	SimpleExtractPrimaries();
				
	for(int bin=0;bin<bins.size();bin++){

		for(int sigma=0;sigma<systpar.steps;sigma++){
			for(int shift=0;shift<systpar.steps;shift++){

			// MAIN Template Fit 
				if(!fitDisabled) {
					cout<<endl;
					cout<<"Bin: "<<bin<<": "<<sigma<<" "<<shift<<endl;
					fits[bin][sigma][shift]->RegularizeTemplateError();	
					if((systpar.mode==3)) fits[bin][sigma][shift]->Templ_D->Smooth();
					Do_TemplateFIT(fits[bin][sigma][shift],fits[bin][sigma][shift]->fitrangemin,fits[bin][sigma][shift]->fitrangemax,constrainmin,constrainmax,IsFitNoise,highmassconstrain);


				}
			}
		}
		// Histograms for systematic error evaluation
		TH2F * dcountsspread;
		TH2F * tfitchisquare;
	
		std::string xaxisname;
		std::string yaxisname;
		
		if(systpar.mode==0) { xaxisname = "#Delta #sigma_{T} (%)";
			      yaxisname = "#Delta (T) (%)";}
		if(systpar.mode==1) { xaxisname = "#Delta #sigma_{1/#beta} (%)";
			      yaxisname = "#Delta (1/#beta) (%)";}
		if(systpar.mode==2) { xaxisname = "#Delta #sigma_{#theta_{C}} (%)";
			      yaxisname = "#Delta (#theta_{C}) (%)";}
		if(systpar.mode==3) { xaxisname = "#Delta #sigma_{1/#beta} (%)";
			      yaxisname = "#Delta (1/#beta) (%)";}
		

		dcountsspread = new TH2F(("DCountsSpread Bin " +to_string(bin)).c_str(),("DCountsSpread Bin " +to_string(bin)+";"+xaxisname+";"+yaxisname).c_str(),systpar.steps,systpar.GetSigmaStart(),systpar.GetSigmaEnd(),systpar.steps,systpar.GetShiftStart(),systpar.GetShiftEnd());

		tfitchisquare = new TH2F(("ChiSquare Bin " +to_string(bin)).c_str(),("ChiSquare Bin " +to_string(bin)+";"+xaxisname+";"+yaxisname).c_str(),systpar.steps,systpar.GetSigmaStart(),systpar.GetSigmaEnd(),systpar.steps,systpar.GetShiftStart(),systpar.GetShiftEnd());


		for(int sigma=0;sigma<systpar.steps;sigma++){
			for(int shift=0;shift<systpar.steps;shift++){
				dcountsspread->SetBinContent(sigma+1,shift+1,fits[bin][sigma][shift]->DCounts);
			
				if(fits[bin][sigma][shift]->ChiSquare>0)
					 tfitchisquare->SetBinContent(sigma+1,shift+1,fits[bin][sigma][shift]->ChiSquare);	
				else  tfitchisquare->SetBinContent(sigma+1,shift+1,5000);
				tfitchisquare->SetBinError(sigma+1,shift+1,0.2);	
				}
		}
	
		
	

		DCountsSpread.push_back(dcountsspread);
		TFitChisquare.push_back(tfitchisquare);
	
	}
	Global_ChiSquare = (TH2F*) TFitChisquare[0]->Clone("Global ChiSquare");
	for(int i=0;i<TFitChisquare.size();i++) Global_ChiSquare->Add(TFitChisquare[i]);
	BestChiSquare = new BestChi();
	BestChiSquare->FindMinimum(Global_ChiSquare);

	for(int bin=0;bin<bins.size();bin++){
	BestChi * bestlocal = new BestChi();
	bestlocal->FindMinimum(TFitChisquare[bin]);
	TH1F * weighteddcounts  = new TH1F(("WeightedCountsBin " +to_string(bin)).c_str(),("WeightedCountsBin " +to_string(bin)+";Counts;").c_str(),100,0.5*fits[bin][bestlocal->i][bestlocal->j]->DCounts,2*fits[bin][bestlocal->i][bestlocal->j]->DCounts);
		for(int sigma=0;sigma<systpar.steps;sigma++)
			for(int shift=0;shift<systpar.steps;shift++)
				if(fits[bin][sigma][shift]->DCounts&&fits[bin][0][5]->DCounts)
					weighteddcounts -> Fill(fits[bin][sigma][shift]->DCounts,EvalFitProbability(fits[bin][sigma][shift]->ChiSquare));

		WeightedDCounts.push_back(weighteddcounts);
			
	}


	EvalFinalParameters();
	EvalFinalErrors();
	CalculateFinalPDCounts();	
	return;
}


void TemplateFIT::EvalFinalParameters(){
	for(int bin=0;bin<bins.size();bin++){

		BestChiSquares     ->SetBinContent(bin+1,TFitChisquare[bin]->GetBinContent(BestChiSquare->i+1,BestChiSquare->j+1));
                if(TFitChisquare[bin]->GetBinContent(1,6)<1000) OriginalChiSquares ->SetBinContent(bin+1,TFitChisquare[bin]->GetBinContent(1,6));
		BestChiSquares     ->SetBinError(bin+1,0.25);	
                OriginalChiSquares ->SetBinError(bin+1,0.25);
	}
	return;
}

void TemplateFIT::EvalFinalErrors(){

	for(int bin=0;bin<bins.size();bin++){

		if(fits[bin][BestChiSquare->i][BestChiSquare->j]->StatErrP>0 
	        && fits[bin][BestChiSquare->i][BestChiSquare->j]->StatErrD>0
		&& fits[bin][BestChiSquare->i][BestChiSquare->j]->DCounts>0){
			StatErrorP -> SetBinContent(bin+1,0.5*fits[bin][BestChiSquare->i][BestChiSquare->j]->StatErrP);
			StatErrorD -> SetBinContent(bin+1,0.5*fits[bin][BestChiSquare->i][BestChiSquare->j]->StatErrD);
			StatErrorT -> SetBinContent(bin+1,0.5*fits[bin][BestChiSquare->i][BestChiSquare->j]->StatErrT);
			SystError -> SetBinContent(bin+1,WeightedDCounts[bin]->GetStdDev()/fits[bin][BestChiSquare->i][BestChiSquare->j]->DCounts);
			StatErrorP -> SetBinError(bin+1,0);
			StatErrorD -> SetBinError(bin+1,0);
			StatErrorT -> SetBinError(bin+1,0);
			SystError -> SetBinError(bin+1,0);
		}

	}
	return;

}


void TemplateFIT::CalculateFinalPDCounts(){

	HeContError = new TH1F ("HeContError","HeContError",bins.size(),0,bins.size());	
	MeasuredHeContRatio = new TH1F("Measured He over P","Measured He over P",bins.size(),0,bins.size());

	for(int bin=0;bin<bins.size();bin++){

		float staterrP= StatErrorP -> GetBinContent(bin+1) * fits[bin][BestChiSquare->i][BestChiSquare->j]->PCounts;
		float staterrD= StatErrorD -> GetBinContent(bin+1) * fits[bin][BestChiSquare->i][BestChiSquare->j]->PCounts;
		float staterrT= StatErrorT -> GetBinContent(bin+1) * fits[bin][BestChiSquare->i][BestChiSquare->j]->PCounts;
	
		float conterrP= 0.1 * fits[bin][BestChiSquare->i][BestChiSquare->j]->TCounts;
		float conterrD= 0.1 * fits[bin][BestChiSquare->i][BestChiSquare->j]->TCounts;
		float conterrT= 0.1 * fits[bin][BestChiSquare->i][BestChiSquare->j]->TCounts;

		float systerrP= SystError -> GetBinContent(bin+1) * fits[bin][BestChiSquare->i][BestChiSquare->j]->PCounts;
		float systerrD= SystError -> GetBinContent(bin+1) * fits[bin][BestChiSquare->i][BestChiSquare->j]->DCounts;
		float systerrT= SystError -> GetBinContent(bin+1) * fits[bin][BestChiSquare->i][BestChiSquare->j]->TCounts;

		float CountsP 	= fits[bin][BestChiSquare->i][BestChiSquare->j]->PCounts - 8.5 * fits[bin][BestChiSquare->i][BestChiSquare->j]->TCounts;
	
		float CountsD;	
		if(IsLocalFit)  CountsD      = WeightedDCounts[bin]->GetMean() - 1.5 * fits[bin][BestChiSquare->i][BestChiSquare->j]->TCounts;
		else CountsD      = fits[bin][BestChiSquare->i][BestChiSquare->j]->DCounts - 1.5 * fits[bin][BestChiSquare->i][BestChiSquare->j]->TCounts;

		float CountsT      = fits[bin][BestChiSquare->i][BestChiSquare->j]->TCounts;
		float CountsP_prim = fits[bin][BestChiSquare->i][BestChiSquare->j]->PCountsPrim - 8.5 * fits[bin][BestChiSquare->i][BestChiSquare->j]->TCounts;
		float CountsD_prim = fits[bin][BestChiSquare->i][BestChiSquare->j]->DCountsPrim - 1.5 * fits[bin][BestChiSquare->i][BestChiSquare->j]->TCounts;

		HeContError -> SetBinContent(bin+1,conterrD);

		ProtonCounts->SetBinContent(bin+1,CountsP);
		ProtonCounts->SetBinError(bin+1,pow(pow(staterrP,2)+pow(systerrP,2)+pow(conterrP,2),0.5));

		DeuteronCounts->SetBinContent(bin+1,CountsD);
		DeuteronCounts->SetBinError(bin+1,pow(pow(staterrD,2)+pow(systerrD,2)+pow(conterrD,2),0.5));		

		TritiumCounts->SetBinContent(bin+1,CountsT);
		TritiumCounts->SetBinError(bin+1,pow(pow(staterrT,2)+pow(systerrT,2)+pow(conterrT,2),0.5));		


		if(fits[bin][BestChiSquare->i][BestChiSquare->j]->PCountsPrim>0){
			ProtonCountsPrim->SetBinContent(bin+1,fits[bin][BestChiSquare->i][BestChiSquare->j]->PCountsPrim);
			ProtonCountsPrim->SetBinError(bin+1,pow(pow(staterrP,2)+pow(systerrP,2)+pow(conterrP,2),0.5) * ProtonCountsPrim->GetBinContent(bin+1)/ProtonCounts->GetBinContent(bin+1) );
		}
		if(fits[bin][BestChiSquare->i][BestChiSquare->j]->DCountsPrim>0){
			DeuteronCountsPrim->SetBinContent(bin+1,fits[bin][BestChiSquare->i][BestChiSquare->j]->DCountsPrim);
			DeuteronCountsPrim->SetBinError(bin+1,pow(pow(staterrD,2)+pow(systerrD,2)+pow(conterrD,2),0.5)* DeuteronCountsPrim->GetBinContent(bin+1)/DeuteronCounts->GetBinContent(bin+1) );		
		}
	}	

}

TH1F * TemplateFIT::Eval_MCHeContRatio(std::string name){
	TH1F * HeCountRatio = new TH1F(name.c_str(),name.c_str(),bins.size(),0,bins.size());
	for(int bin=0;bin<bins.size();bin++){
                TH1F * HeliumCounts = (TH1F*) fits[bin][0][5]->Templ_He->Clone();
                TH1F * ProtonCounts = (TH1F*) fits[bin][0][5]->Templ_P->Clone();	
		if(ProtonCounts->Integral()>0){
			HeCountRatio->SetBinContent(bin+1,HeliumCounts->Integral()/ProtonCounts->Integral());
			HeCountRatio->SetBinError(bin+1,pow(HeliumCounts->Integral(),0.5)/ProtonCounts->Integral());
		}
	}
	return HeCountRatio;
}

void TemplateFIT::Eval_ContError(){

	return;
}
