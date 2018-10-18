#include "TemplateFITbetasmear.h"

void TemplateFIT::Eval_TransferFunction(){
	for(int bin=0;bin<fits.size();bin++){
		TH1F * transferfunction = (TH1F *) fits[bin][0][5]->DataPrim->Clone();
		transferfunction->Sumw2();
		transferfunction->Divide(fits[bin][0][5]->Data);
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
	if((int)vars->U_time!=ActualTime) {
		UpdateZoneLivetime(vars->Livetime,vars->Rcutoff,ExposureTime,bins);
		ActualTime=vars->U_time;
	}
	if(ApplyCuts(cut,vars)&&kbin>0){
		//cout<<"Bin; "<<kbin<<endl;
		//
		//vars->PrintCurrentState();
		for(int i=0;i<systpar.steps;i++)
			for(int j=0;j<systpar.steps;j++){
				fits[kbin][i][j]->Data->Fill(var(vars),vars->PrescaleFactor);		
				if(ApplyCuts(cutprimary,vars)) fits[kbin][i][j]->DataPrim->Fill(var(vars),vars->PrescaleFactor);
				}
	}
	return;	
}

float TemplateFIT::SmearBetaRICH(float Beta, float stepsigma, float stepshift){
	float angle;
	angle= acos(1/(1.2*Beta))*10e4;
	float shiftstart=-systpar.shift;
	angle = angle + (shiftstart+(2*systpar.shift/(float)systpar.steps)*stepshift) + Rand->Gaus(0,(float)((2*systpar.sigma/systpar.steps)*stepsigma));
	return 1/(1.2*cos(angle/10e4));
}


float TemplateFIT::SmearBeta(float Beta, float stepsigma, float stepshift,float R){

	float time = 1.2/(Beta*3e-4);
	float shiftstart=-systpar.shift;

	float tailcontrolfactor=1;	//migration tail fixing
/*	if(R>2.3&&R<2.5) tailcontrolfactor *=1.1;
	if(R>2.5&&R<2.7) tailcontrolfactor *=1.1;
	if(R>2.7&&R<3)   tailcontrolfactor *=1.1;
	if(R>3)          tailcontrolfactor *=1.1;
*/
	if(R>2.5) tailcontrolfactor = 135./90.;
	
	float smeartime = (shiftstart+(2*systpar.shift/(float)systpar.steps)*stepshift) + Rand->Gaus(0,(float) tailcontrolfactor*((2*systpar.sigma/systpar.steps)*stepsigma));
	time = time + smeartime;
	return 1.2/(time*3e-4);

}



void TemplateFIT::FillEventByEventMC(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars)){

	std::string cutP=cut+"&IsProtonMC";
	std::string cutD=cut+"&IsPureDMC";
	std::string cutHe=cut+"&IsPureTMC";
	//std::string cutD=cut+"&IsDeutonMC";
	//std::string cutHe=cut+"&IsHeliumMC";


	cutHe.erase(cutHe.find("IsMinimumBias&IsLooseCharge1&"),14);
	cutHe = "IsMinimumBias&" + cutHe;  //releasing cut for more stat. in Tritium templates

	if((ApplyCuts(cutP,vars)||ApplyCuts(cutD,vars)||ApplyCuts(cutHe,vars))){
		for(int i=0;i<systpar.steps;i++)
			for(int j=0;j<systpar.steps;j++){
				float betasmear;
				float mctotalweight = vars->mcweight; // Latweighter->GetWeight(fabs(vars->R));
				if(!ApplyCuts("IsOnlyFromToF",vars)) betasmear = SmearBetaRICH(vars->BetaRICH_new,(float)i,(float)j); 
				else 	{ 
						betasmear = SmearBeta(vars->Beta,(float)i,(float)j,vars->R);    
						if(betasmear>0.7776) mctotalweight = vars->mcweight / Latweighter->GetWeight(fabs(vars->R));
				}
				int kbin;
	
				if(!IsFitNoise){
				   if(bins.IsUsingBetaEdges()){	 
					kbin = bins.GetBin(betasmear);

					float mass = vars->R/betasmear * pow((1-pow(betasmear,2)),0.5);

					if(ApplyCuts(cutP,vars)&&kbin>0)  fits[kbin][i][j]->Templ_P->Fill(mass,mctotalweight);		
					if(ApplyCuts(cutD,vars)&&kbin>0)  fits[kbin][i][j]->Templ_D->Fill(mass,vars->mcweight);
					if(ApplyCuts(cutHe,vars)&&kbin>0) fits[kbin][i][j]->Templ_He->Fill(mass,vars->mcweight);
					}
				    else {
				   	kbin = bins.GetBin(discr_var(vars));
					if(ApplyCuts(cutP,vars)&&kbin>0)  fits[kbin][i][j]->Templ_P->Fill(betasmear,mctotalweight);
                                        if(ApplyCuts(cutD,vars)&&kbin>0)  fits[kbin][i][j]->Templ_D->Fill(betasmear,vars->mcweight);
                                        if(ApplyCuts(cutHe,vars)&&kbin>0) fits[kbin][i][j]->Templ_He->Fill(betasmear,vars->mcweight);

				   }	

				}
				else{
					if(bins.IsUsingBetaEdges()) {
						kbin = bins.GetBin(betasmear);
						float mass = vars->R/betasmear * pow((1-pow(betasmear,2)),0.5);		

						if(ApplyCuts(cutP,vars)&&kbin>0)  fits[kbin][i][j]->Templ_P->Fill(mass,mctotalweight);		
						if(ApplyCuts(cutD,vars)&&kbin>0)  fits[kbin][i][j]->Templ_D->Fill(mass,vars->mcweight);
						if(ApplyCuts(cutHe,vars)&&kbin>0) fits[kbin][i][j]->Templ_He->Fill(mass,vars->mcweight);

						float betabad = betasmear;
						if(BadEvSim) {betabad=BadEvSim->SimulateBadEvents(betasmear); 
							kbin = bins.GetBin(betabad);
						}
						float mass_bad = vars->R/betabad * pow((1-pow(betabad,2)),0.5);
						if(ApplyCuts(cutP,vars)&&kbin>0)  fits[kbin][i][j]->Templ_Noise->Fill(mass_bad,mctotalweight);
					}					
					else {
						if(BadEvSim) betasmear=BadEvSim->SimulateBadEvents(betasmear);
						kbin = bins.GetBin(discr_var(vars));
						if(ApplyCuts(cutP,vars)&&kbin>0)  fits[kbin][i][j]->Templ_P->Fill(betasmear,mctotalweight);
	                                        if(ApplyCuts(cutD,vars)&&kbin>0)  fits[kbin][i][j]->Templ_D->Fill(betasmear,vars->mcweight);
        	                                if(ApplyCuts(cutHe,vars)&&kbin>0) fits[kbin][i][j]->Templ_He->Fill(betasmear,vars->mcweight);
					}
				}
			}
	}
	return;	
}



void TemplateFIT::Save(){

	finalhistos.Add(ExposureTime);
	finalhistos.writeObjsInFolder((basename + "/ExposureTime").c_str(),false);	

	for(int bin=0;bin<bins.size();bin++){ 
		finalhistos.Add(fits[bin][0][5]->Data);
		finalhistos.Add(fits[bin][0][5]->DataPrim);

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


	for(int bin=0;bin<bins.size();bin++){
	TH1F * OriginalP=(TH1F*)fits[bin][0][5]->Templ_P->Clone();
	TH1F * BestP=(TH1F*)fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_P->Clone();
	OriginalP->SetName("Original Proton MC ");
	BestP->SetName("Best #chi^{2} Mod. Proton MC ");
	finalhistos.Add(OriginalP);
	finalhistos.Add(BestP);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesP/Bin"+to_string(bin)).c_str());
	
	TH1F * OriginalD=(TH1F*)fits[bin][0][5]->Templ_D->Clone();
	TH1F * BestD=(TH1F*)fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_D->Clone();
	OriginalD->SetName("Original Deuton MC ");
	BestD->SetName("Best #chi^{2} Mod. Deuton MC ");
	finalhistos.Add(OriginalD);
	finalhistos.Add(BestD);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesD/Bin"+to_string(bin)).c_str());

	TH1F * OriginalHe=(TH1F*)fits[bin][0][5]->Templ_He->Clone();
	TH1F * BestHe=(TH1F*)fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_He->Clone();
	OriginalHe->SetName("Original Tritium MC ");
	BestHe->SetName("Best #chi^{2} Tritium MC ");
	finalhistos.Add(OriginalHe);
	finalhistos.Add(BestHe);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesHe/Bin"+to_string(bin)).c_str());

	TH1F * OriginalNoise=(TH1F*)fits[bin][0][5]->Templ_Noise->Clone();
	TH1F * BestNoise=(TH1F*)fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_Noise->Clone();
	OriginalNoise->SetName("Original Noise Template MC ");
	BestNoise->SetName("Best #chi^{2} Noise Template MC ");
	finalhistos.Add(OriginalNoise);
	finalhistos.Add(BestNoise);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesNoise/Bin"+to_string(bin)).c_str());
	
	}

	for(int bin=0;bin<bins.size();bin++){ 
		finalhistos.Add(fits[bin][0][5]->Data);
		finalhistos.Add(fits[bin][0][5]->DataPrim);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/Data/Bin"+to_string(bin)).c_str());	
	}
	
	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(!(i==0&&i==5))finalhistos.Add(fits[bin][i][j]->Templ_P);
			}	
	finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesP/Bin"+to_string(bin)).c_str());
	}

	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(!(i==0&&i==5))finalhistos.Add(fits[bin][i][j]->Templ_D);
				}
	finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesD/Bin"+to_string(bin)).c_str());
	}
	
	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(!(i==0&&i==5))finalhistos.Add(fits[bin][i][j]->Templ_He);
				}
	finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesHe/Bin"+to_string(bin)).c_str());
	}

	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(!(i==0&&i==5))finalhistos.Add(fits[bin][i][j]->Templ_Noise);
				}
	finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesNoise/Bin"+to_string(bin)).c_str());
	}
	for(int bin=0;bin<bins.size();bin++){
		if(BestChiSquare[bin]&&fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Tfit){
			TH1F * FIT;
			if(fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Tfit_outcome!=-1){ 
				FIT=(TH1F*)fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Tfit-> GetPlot();
				FIT->SetName(("Fraction Fit bin" + to_string(bin)).c_str());
				finalhistos.Add(FIT);
				finalhistos.writeObjsInFolder((basename+"/Fit Results/FractionFits/Bin"+to_string(bin)).c_str());
			}
		}
	}


	for(int i=0;i<bins.size();i++) 
                finalhistos.Add(TransferFunction[i]);
        finalhistos.writeObjsInFolder((basename + "/Fit Results/TrasnferFunctions/").c_str());

	for(int i=0;i<bins.size();i++){
		finalhistos.Add(DCountsSpread[i]);
		finalhistos.writeObjsInFolder((basename+"/Fit Results/Spreads/DCounts").c_str());
	}
	for(int i=0;i<bins.size();i++){
		finalhistos.Add(DErrSpread[i]);
		finalhistos.writeObjsInFolder((basename+"/Fit Results/Spreads/DErrs").c_str());
	}
	
	for(int i=0;i<bins.size();i++){
		finalhistos.Add(TFitChisquare[i]);
		finalhistos.writeObjsInFolder((basename+"/Fit Results/Spreads/ChiSquare").c_str());	
	}
	for(int i=0;i<bins.size();i++){
		finalhistos.Add(WeightedDCounts[i]);
		finalhistos.writeObjsInFolder((basename+"/Fit Results/Spreads/Weighted D counts").c_str());	
	}

	if(HeContModel) finalhistos.Add(HeContModel);
	if(HeContError) finalhistos.Add(HeContError);
	if(MeasuredHeContRatio) finalhistos.Add(MeasuredHeContRatio);
	if(MCHeContRatio) finalhistos.Add(MCHeContRatio);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/HeliumContamination").c_str());	

	finalhistos.Add(BestFitSigma);
	finalhistos.Add(BestFitShift);
	finalhistos.Add(StatErrorP);
	finalhistos.Add(StatErrorD);
	finalhistos.Add(StatErrorT);
	finalhistos.Add(SystError);
	finalhistos.Add(ProtonCounts);
	finalhistos.Add(DeuteronCounts);
	finalhistos.Add(TritiumCounts);
	finalhistos.Add(ProtonCountsPrim);
	finalhistos.Add(DeuteronCountsPrim);

	finalhistos.Add(BestChiSquares   ); 
	finalhistos.Add(OriginalChiSquares);

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
				cout<<"bin: "<<bin<<": "<<fits[bin][i][j]->fitrangemin<<" "<<fits[bin][i][j]->fitrangemax<<endl;
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

float GetChiSquare(TH1 * Result, TH1 * Data, float min, float max){

	int binmin = Data->FindBin(min);
	int binmax = Data->FindBin(max);
	float chi = 0;
	float err = 0;
	for(int i=binmin; i<binmax; i++){
		if(Data->GetBinContent(i+1)>0&&Result->GetBinContent(i+1)>0){
			err = pow(pow(Data->GetBinError(i+1),2) + pow(Result->GetBinError(i+1),2),0.5);
			chi += pow((Data->GetBinContent(i+1) - Result->GetBinContent(i+1)),2)/pow(err,2);
		}
	}
	chi /= ((binmax-binmin)-3);
	cout<<"chi: "<<chi<<endl;
	return chi;

}


void Do_TemplateFIT(TFit * Fit,float fitrangemin,float fitrangemax,float constrain_min[], float constrain_max[], bool isfitnoise, bool highmasstailconstrain ){
Pre_Scale(Fit->Data, Fit ->  Templ_P,Fit ->  Templ_D,Fit ->  Templ_He);
	//CleanTemplates(Fit ->  Templ_P);
	if(isfitnoise) cout<<"********** FIT NOISE MODE **********"<<endl;
	if(isfitnoise) Fit->Templ_Noise->Scale(0.01*Fit ->  Templ_P->Integral()/Fit ->  Templ_Noise->Integral());
	TObjArray *Tpl;
	Tpl = new TObjArray(3);
	if(Fit ->  Templ_P)  if(Fit ->  Templ_P ->GetEntries()>0) Tpl -> Add( Fit ->  Templ_P );
	if(Fit ->  Templ_D)  if(Fit ->  Templ_D ->GetEntries()>0) Tpl -> Add( Fit ->  Templ_D );
        if(Fit ->  Templ_He) if(Fit ->  Templ_He->GetEntries()>0) Tpl -> Add( Fit ->  Templ_He);
	if(isfitnoise) if(Fit ->  Templ_Noise) if(Fit ->  Templ_Noise->GetEntries()>50) Tpl -> Add( Fit ->  Templ_Noise);

	float min=fitrangemin;
	float max=fitrangemax;

	bool fitcondition = (Fit -> Data->Integral()>50)&&(Fit -> Templ_P->Integral()>50) &&(Fit -> Templ_D->Integral()>50);
	cout<<	Fit -> Data->Integral()<<" "<<Fit -> Templ_P->Integral()<<" "<<Fit -> Templ_D->Integral()<<endl;
	
if(fitcondition){	
		Fit -> Tfit = new TFractionFitter(Fit -> Data, Tpl ,"q");
		Fit -> Tfit -> SetRangeX(Fit -> Data -> FindBin(min), Fit -> Data -> FindBin(max));
		
		Fit -> Tfit -> Constrain(1, constrain_min[0] ,constrain_max[0]);
                Fit -> Tfit -> Constrain(2, constrain_min[1] ,constrain_max[1]);
	 	Fit -> Tfit -> Constrain(3, constrain_min[2] ,constrain_max[2]);
		if(highmasstailconstrain) 
			Fit -> Tfit -> Constrain(3, 0.1*CalculateAmountOfHighMassComponent(Fit -> Data,Fit ->  Templ_P,Fit ->  Templ_Noise,3.3,constrain_min[2]),1.2*CalculateAmountOfHighMassComponent(Fit -> Data,Fit ->  Templ_P,Fit ->  Templ_Noise,3.1,constrain_max[2]));	 
	        if(isfitnoise&&highmasstailconstrain) 
			Fit -> Tfit -> Constrain(4, 0.1*CalculateAmountOfHighMassComponent(Fit -> Data,Fit ->  Templ_P,Fit ->  Templ_Noise,4.5,0.001),1.2*CalculateAmountOfHighMassComponent(Fit -> Data,Fit ->  Templ_P,Fit ->  Templ_Noise,4.5,0.001));	 
			
		if(Fit -> Tfit ) Fit -> Tfit_outcome = Fit -> Tfit -> Fit();

		for(int fit_attempt=0; fit_attempt<20; fit_attempt++) {
			cout<<"fit attempt: "<<fit_attempt<<endl;
			if(Fit -> Tfit_outcome == 0) break;
			else {
				cout<<fit_attempt<<endl;
				Fit -> Tfit -> SetRangeX(Fit -> Data -> FindBin((1+0.01*fit_attempt)*min), Fit -> Data -> FindBin((1+0.01*fit_attempt)*max));
				Fit -> Tfit_outcome = Fit -> Tfit -> Fit();
			}
		}

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
		
	

			float i1 = Fit-> Templ_P  ->Integral(Fit->Templ_P -> FindBin(min), Fit->Templ_P -> FindBin(max));
			float i2 = Fit-> Templ_D  ->Integral(Fit->Templ_D -> FindBin(min), Fit->Templ_D -> FindBin(max));
			float i3=1;
			if(Fit-> Templ_He) i3 = Fit-> Templ_He ->Integral(Fit->Templ_He -> FindBin(min), Fit->Templ_He -> FindBin(max));
			float i4=1;
			if(isfitnoise) i4 = Fit-> Templ_Noise ->Integral(Fit->Templ_Noise -> FindBin(min), Fit->Templ_Noise -> FindBin(max));


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

			cout<<"fract: "<<w1<<" "<<w2<<" "<<w3<<" "<<w4<<endl;
			cout<<"err: "  <<e1<<" "<<e2<<" "<<e3<<" "<<e4<<endl;


			Fit ->  Templ_P  -> Scale(Fit ->wheightP);
			Fit ->  Templ_D  -> Scale(Fit ->wheightD);
			if(Fit-> Templ_He) Fit ->  Templ_He -> Scale(Fit ->wheightHe);
			if(isfitnoise) Fit ->  Templ_Noise -> Scale(Fit ->wheightNoise);
			float Cov01=0;

			Cov01= Fit -> Tfit->GetFitter()->GetCovarianceMatrixElement(0,1);

			TH1F * Sum = (TH1F *)Fit ->  Templ_P->Clone();
			Sum -> Add(Fit ->  Templ_D);
			Sum -> Add(Fit ->  Templ_He);
			if(isfitnoise) Sum -> Add(Fit ->  Templ_Noise);			

			Fit -> StatErrP = e1;
			Fit -> StatErrD = e2;
			Fit -> StatErrT = e3;
			
	//		Fit -> ChiSquare = Fit -> Tfit -> GetChisquare()/(float) (Fit ->  Tfit -> GetNDF());
			Fit -> ChiSquare = GetChiSquare(Sum,Fit->Data,min,max);
			Fit -> DCounts = Fit ->  Templ_D -> Integral();
			Fit -> PCounts = Fit ->  Templ_P -> Integral();
			Fit -> TCounts = Fit ->  Templ_He -> Integral();


		}
		else{
			Fit -> PCounts = Fit -> Data -> Integral();
			Fit -> DCounts = Fit -> Data -> Integral();
			Fit ->wheightP= 1;
			Fit ->wheightD= 1;
		}
	}
	return;
}

float EvalFitProbability(float chi){
	TF1 * CumulativeChiSquare=new TF1("Chi","exp(-x)*x^-0.5",0.001,10);
	if(chi>0)
		return CumulativeChiSquare->Integral(chi,100);
	else return 0;
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
		 }
	}
};

void TemplateFIT::ExtractCounts(FileSaver finalhistos){
	Eval_TransferFunction();
	for(int bin=0;bin<bins.size();bin++){

		// Template Fits
		for(int sigma=0;sigma<systpar.steps;sigma++){
			for(int shift=0;shift<systpar.steps;shift++){

				if(!fitDisabled) {
					cout<<"Bin: "<<bin<<": "<<sigma<<" "<<shift<<endl;

					Do_TemplateFIT(fits[bin][sigma][shift],fits[bin][sigma][shift]->fitrangemin,fits[bin][sigma][shift]->fitrangemax,
					constrainmin,constrainmax,IsFitNoise,highmassconstrain);
				}
	
				fits[bin][sigma][shift]->Templ_DPrim=(TH1F*) fits[bin][sigma][shift]->Templ_D->Clone();
				fits[bin][sigma][shift]->Templ_DPrim->Multiply(TransferFunction[bin]);
				if(fits[bin][sigma][shift]->wheightD==0) fits[bin][sigma][shift]->DCountsPrim=0;	
				else{
					if(!fitDisabled) fits[bin][sigma][shift]->DCountsPrim = fits[bin][sigma][shift]->Templ_DPrim->Integral();
					else { 
						fits[bin][sigma][shift]->DCountsPrim = fits[bin][sigma][shift]->DataPrim->Integral();
						fits[bin][sigma][shift]->DCounts     = fits[bin][sigma][shift]->DataPrim->Integral();
					}
				}
				fits[bin][sigma][shift]->Templ_PPrim=(TH1F*) fits[bin][sigma][shift]->Templ_P->Clone();
				fits[bin][sigma][shift]->Templ_PPrim->Multiply(TransferFunction[bin]);
				if(fits[bin][sigma][shift]->wheightD==0) fits[bin][sigma][shift]->PCountsPrim=0;
				else{
					if(!fitDisabled) fits[bin][sigma][shift]->PCountsPrim = fits[bin][sigma][shift]->Templ_PPrim->Integral();
					else { 
						fits[bin][sigma][shift]->PCountsPrim = fits[bin][sigma][shift]->DataPrim->Integral();
						fits[bin][sigma][shift]->PCounts     = fits[bin][sigma][shift]->DataPrim->Integral();
						fits[bin][sigma][shift]->TCounts     = fits[bin][sigma][shift]->DataPrim->Integral();
					}
				}
			}
		}
		// Histograms for systematic error evaluation
		TH2F * dcountsspread = new TH2F(("DCountsSpread Bin " +to_string(bin)).c_str(),("DCountsSpread Bin " +to_string(bin)).c_str(),systpar.steps,0,2*systpar.sigma,systpar.steps,-systpar.shift,systpar.shift);
		TH2F * derrspread = new TH2F(("DerrSpread Bin " +to_string(bin)).c_str(),("DerrSpread Bin " +to_string(bin)).c_str(),systpar.steps,0,2*systpar.sigma,systpar.steps,-systpar.shift,systpar.shift);
		TH2F * tfitchisquare = new TH2F(("ChiSquare Bin " +to_string(bin)).c_str(),("ChiSquare Bin " +to_string(bin)).c_str(),systpar.steps,0,2*systpar.sigma,systpar.steps,-systpar.shift,systpar.shift);
		TH1F * weighteddcounts  = new TH1F(("Weighted Counts Bin " +to_string(bin)).c_str(),("Weighted Counts Bin " +to_string(bin)).c_str(),35,0.4*fits[bin][0][4]->DCounts,1.7*fits[bin][0][4]->DCounts);
		BestChi * MinimumChi = new BestChi();



		for(int sigma=0;sigma<systpar.steps;sigma++){
			for(int shift=0;shift<systpar.steps;shift++){
				dcountsspread->SetBinContent(sigma+1,shift+1,fits[bin][sigma][shift]->DCounts);
				derrspread->SetBinContent(sigma+1,shift+1,fits[bin][sigma][shift]->StatErrD);///(fits[bin][sigma][shift]->DCounts/fits[bin][sigma][shift]->PCounts));
				if(fits[bin][sigma][shift]->ChiSquare>0)
					tfitchisquare->SetBinContent(sigma+1,shift+1,fits[bin][sigma][shift]->ChiSquare);
				else  tfitchisquare->SetBinContent(sigma+1,shift+1,500);
				if(fits[bin][sigma][shift]->DCounts&&fits[bin][0][4]->DCounts){
					if(	fits[bin][sigma][shift]->DCounts>0.1*fits[bin][0][4]->DCounts &&
							fits[bin][sigma][shift]->DCounts<2*fits[bin][0][4]->DCounts  ){
							weighteddcounts -> Fill(fits[bin][sigma][shift]->DCounts,EvalFitProbability(fits[bin][sigma][shift]->ChiSquare));
						}
				}
				cout<<sigma<<" "<<shift<<endl;
			}
		}
		
		MinimumChi->FindMinimum(tfitchisquare,derrspread);	


		BestChiSquare.push_back(MinimumChi);			
		DCountsSpread.push_back(dcountsspread);
		DErrSpread.push_back(derrspread);
		TFitChisquare.push_back(tfitchisquare);
		WeightedDCounts.push_back(weighteddcounts);


	}

	
	EvalFinalParameters();
	EvalFinalErrors();
	CalculateFinalPDCounts();	
	return;
}

void TemplateFIT::EvalFinalParameters(){
	for(int bin=0;bin<bins.size();bin++){
		float BinSTD=0;
		for(int i=0;i<systpar.steps;i++){
			for(int j=0;j<systpar.steps;j++){
				BinSTD+=pow(TFitChisquare[bin]->GetBinContent(i+1,j+1)-BestChiSquare[bin]->chimin,2);
			}
		}
		BinSTD=pow(BinSTD,0.5)/(systpar.steps*systpar.steps);

		BestFitSigma->SetBinContent(bin+1,(float)((2*systpar.sigma/systpar.steps)*BestChiSquare[bin]->i));
		BestFitShift->SetBinContent(bin+1,-systpar.shift+(2*systpar.shift/(float)systpar.steps)*BestChiSquare[bin]->j);

		
	
		if(BinSTD>0){
		BestFitSigma->SetBinError(bin+1,pow(pow(4/BinSTD,2)+pow(1.5*(2*systpar.sigma/(float)systpar.steps),2),0.5));
		BestFitShift->SetBinError(bin+1,pow(pow(4/BinSTD,2)+pow(1.5*(2*systpar.shift/(float)systpar.steps),2),0.5));
		}

		if(BestChiSquare[bin]->chimin<1000) 	       BestChiSquares     ->SetBinContent(bin+1,BestChiSquare[bin]->chimin);
                if(TFitChisquare[bin]->GetBinContent(1,6)<1000)OriginalChiSquares ->SetBinContent(bin+1,TFitChisquare[bin]->GetBinContent(1,6));
		BestChiSquares     ->SetBinError(bin+1,0.25);	
                OriginalChiSquares ->SetBinError(bin+1,0.25);
	}
	return;
}

void TemplateFIT::EvalFinalErrors(){

	for(int bin=0;bin<bins.size();bin++){

		if(fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts>0){
			StatErrorP -> SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->StatErrP);
			StatErrorD -> SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->StatErrD);
			StatErrorT -> SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->StatErrT);
			SystError -> SetBinContent(bin+1,WeightedDCounts[bin]->GetStdDev()/fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts*fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->ChiSquare);
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
	Eval_ContError();
	StatErrorP -> Smooth(1);
	StatErrorD -> Smooth(1);
	StatErrorT -> Smooth(1);

	SystError -> Smooth(1);

	for(int bin=0;bin<bins.size();bin++){

	float staterrP= StatErrorP -> GetBinContent(bin+1) * fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCounts;
		float staterrD= StatErrorD -> GetBinContent(bin+1) * fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCounts;
		float staterrT= StatErrorT -> GetBinContent(bin+1) * fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCounts;
	
		float conterrP= 0.3 * fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->TCounts;
		float conterrD= 0.3 * fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->TCounts;
		float conterrT= 0.3 * fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->TCounts;

		float systerrP= SystError -> GetBinContent(bin+1) * fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCounts;
		float systerrD= SystError -> GetBinContent(bin+1) * fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts;
		float systerrT= SystError -> GetBinContent(bin+1) * fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->TCounts;

		float CountsP      = fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCounts - 8.5 * fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->TCounts;
		float CountsD      = fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts - 1.5 * fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->TCounts;
		float CountsT      = fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->TCounts;
		float CountsP_prim = fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCountsPrim - 8.5 * fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->TCounts;
		float CountsD_prim = fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCountsPrim - 1.5 * fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->TCounts;

		HeContError -> SetBinContent(bin+1,conterrD);

		ProtonCounts->SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCounts);
		ProtonCounts->SetBinError(bin+1,pow(pow(staterrP,2)+pow(systerrP,2)+pow(conterrP,2),0.5));

		DeuteronCounts->SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts);
		DeuteronCounts->SetBinError(bin+1,pow(pow(staterrD,2)+pow(systerrD,2)+pow(conterrD,2),0.5));		

		TritiumCounts->SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->TCounts);
		TritiumCounts->SetBinError(bin+1,pow(pow(staterrT,2)+pow(systerrT,2)+pow(conterrT,2),0.5));		


		if(fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCountsPrim>0){
			ProtonCountsPrim->SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCountsPrim);
			ProtonCountsPrim->SetBinError(bin+1,pow(pow(staterrP,2)+pow(systerrP,2)+pow(conterrP,2),0.5) * ProtonCountsPrim->GetBinContent(bin+1)/ProtonCounts->GetBinContent(bin+1) );
		}
		if(fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCountsPrim>0){
			DeuteronCountsPrim->SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCountsPrim);
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
