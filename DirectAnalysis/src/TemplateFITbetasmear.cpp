#include "TemplateFITbetasmear.h"


TCanvas * DrawTemplateFit(int bin, TH1F * Data, TH1F* BestP,TH1F * BestD,TH1F *BestHe,TH1F * Residuals){

	TCanvas *c1 = new TCanvas(("Bin: "+to_string(bin)).c_str());
	c1->SetCanvasSize(800,800);
	c1->Divide(1,2,0,0.1);
	c1->cd(1);
	gPad->SetLogy();
	gPad->SetTickx();
	gPad->SetTicky();

	TH1F * sum = (TH1F*) BestP->Clone("sum");
        sum->Add(BestD);
        if(BestHe) if(BestHe->Integral()>0) sum->Add(BestHe);
        sum->SetLineColor(1);
        sum->SetLineWidth(2);


        Data->SetMarkerStyle(8);
        //Data->GetXaxis()->SetRangeUser(0.5,6.5);
        Data->Draw("same");
        sum->Draw("histsame");
        BestP->SetLineColor(2);
        BestD->SetLineColor(4);
        if(BestHe) if(BestHe->Integral()>0) BestHe->SetLineColor(3);
        BestP->SetLineWidth(2);
        BestD->SetLineWidth(2);
        if(BestHe) if(BestHe->Integral()>0) BestHe->SetLineWidth(2);


        BestP->Draw("histsame");
        BestD->Draw("histsame");
        if(BestHe) if(BestHe->Integral()>0) BestHe->Draw("histsame");

	c1->cd(2);
        gPad->SetTickx();
        gPad->SetTicky();
        gPad->SetGridy();


	Residuals->SetMarkerStyle(8);
        Residuals->SetMarkerColor(1);
        //Residuals->GetXaxis()->SetRangeUser(0.5,6.5);
        Residuals->Draw("P");


	return c1;

}



void AdjustErrorsforFullStat(TH1F* histo,int ntimebins){
	for(int i=0;i<histo->GetNbinsX();i++) histo->SetBinError(i+1,ntimebins*histo->GetBinError(i+1));
	return;
};




float EqualizeTemplate(TH1F * data, TH1F * templ){

	TH1F * DT = new TH1F("DT","DT",data->GetNbinsX(),data->GetBinLowEdge(1),data->GetBinLowEdge(data->GetNbinsX()+1));
	TH1F * MC = new TH1F("MC","MC",data->GetNbinsX(),data->GetBinLowEdge(1),data->GetBinLowEdge(data->GetNbinsX()+1));

	for(int i=0;i<data->GetNbinsX();i++){
		DT->SetBinContent(i+1,data->GetBinContent(i+1));
		DT->SetBinError(i+1,data->GetBinError(i+1));
		MC->SetBinContent(i+1,templ->GetBinContent(i+1));
		MC->SetBinError(i+1,templ->GetBinError(i+1));
	}
/*
	TF1 * fitdata = new TF1("fitdata","gaus",0,3);
	TF1 * fittemp = new TF1("fittemp","gaus",0,3);
*/
	for(int i=0;i<data->GetNbinsX();i++)
	if(data->GetBinLowEdge(i+1)<0.7||data->GetBinLowEdge(i+1)>1.08){
		DT->SetBinContent(i+1,0);
		MC->SetBinContent(i+1,0);
	}	
	//DT->Fit("fitdata","","",0.75,1.07);
	//MC->Fit("fittemp","","",0.75,1.07);


	float meanDT=/*fitdata->GetParameter(1);*/DT->GetMean();
	float meanTP=/*fittemp->GetParameter(1);*/MC->GetMean();
	cout<<"Shift: "<<meanTP<<" "<<meanDT<<endl;
	float shift =  meanTP - meanDT;

	return shift;
}





TH1F * TemplateFIT::BuildTemplateFromExternal(TH1F* core, TH1F *reference, TH1F* external){
	float External = external->Integral();
	float Core = core->Integral();
	float Reference = reference->Integral();
	
	TH1F * templ = (TH1F*) core->Clone();
	templ->SetName(reference->GetName());

	TH1F* a = (TH1F*) external->Clone("a");
	TH1F* b = (TH1F*) reference->Clone("b");
	a->Scale(1/External);
	b->Scale(1/External);
	a->Add(b,-1);

	templ->Scale(1/Core);
	templ->Add(a);		
	
	
	templ->Scale(Core);
	//for(int i=0;i<templ->GetNbinsX();i++) if(templ->GetBinContent(i+1)<0) templ->SetBinContent(i+1,0);

	return templ;
}


TSpline3 * GetSplineFromHisto(TH1F * Graph, Binning bins){
	
    
        double X[Graph->GetNbinsX()];
        double Y[Graph->GetNbinsX()];
   	cout<<"********* EXPOSURE SPLINE **********"<<endl; 
        for(int i=0; i<Graph->GetNbinsX(); i++){        X[i]=bins.RigBinCent(i); Y[i]=Graph->GetBinContent(i+1); }
        TSpline3 *Exposure = new TSpline3("Exposure",X,Y,Graph->GetNbinsX());
        return Exposure;
}

void TFit::AddSysttoTempl(){

	for(int i=0;i<Templ_P->GetNbinsX();i++){
		Templ_P->SetBinError(i+1,2*Templ_P->GetBinError(i+1));
		Templ_D->SetBinError(i+1,2*Templ_D->GetBinError(i+1));
		Templ_He->SetBinError(i+1,2*Templ_He->GetBinError(i+1));
	}	

}

void TFit::AdjoustProtonTail(float factor){

	for(int i=0;i<Templ_P->GetNbinsX();i++){
		if(Templ_P->GetBinCenter(i+1)>1.4){
		Templ_P->SetBinContent(i+1,factor*Templ_P->GetBinContent(i+1));
		Templ_P->SetBinError(i+1,factor*Templ_P->GetBinError(i+1));
		}
	}	

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
			if(Templ_P->GetBinContent(i+1)==0) Templ_P->SetBinError(i+1,pow(content,0.5));
			if(Templ_P->GetBinCenter(i+1)>4) Templ_P->SetBinError(i+1,0.5*Templ_P->GetBinError(i+1));
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
	//cout<<vars->Event<<endl;
	//vars->PrintCurrentState();
	if(ApplyCuts((cut+"&RigSafetyCut").c_str(),vars)/*&&kbin>=0*/){
		if(!(kbin<0)){
		//if(vars->BetaRICH_new>0) cout<<vars->Event<<" "<<vars->BetaRICH_new<<endl;
		for(int i=0;i<systpar.steps;i++)
			for(int j=0;j<systpar.steps;j++){
				if(fits[kbin][i][j]->Data){
				if(ApplyCuts(cut,vars)) 	fits[kbin][i][j]->DataAmbient->Fill(var(vars),vars->PrescaleFactor);		
				if(ApplyCuts((cut+"&"+cutoff).c_str(),vars)) fits[kbin][i][j]->Data->Fill(var(vars),vars->PrescaleFactor);		
				if(ApplyCuts(cutprimary,vars)) 		     fits[kbin][i][j]->DataPrim->Fill(var(vars),vars->PrescaleFactor);
					}
				}
		}
	}
	return;	
}


void TemplateFIT::SimpleExtractPrimaries(){
        for(int bin=0;bin<bins.size();bin++)
                for(int sigma=0;sigma<systpar.steps;sigma++)
                        for(int shift=0;shift<systpar.steps;shift++){
				if(fits[bin][sigma][shift]->Templ_D){
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
	if(mode==3) return  ((1.28694e-03*(1-sigma/100.) - 1.36634e-03)/1.36634e-03)*100; //delta (sigma(1/beta)) (%)
}

float Systpar::GetSigmaEnd(){
	if(mode==0) return ((2*sigma/5000)/0.0318)*100;
	if(mode==1) return ((0.0284/2 + 0.0284/2*sigma/100)/0.0318)*100;
	if(mode==2) return /*sigma angolo da prop. errori di 0.3% di errore su beta*/ (2*sigma/(100000*1/sqrt(1-pow(1/(1.2*0.98),2))*((1/1.2)/(0.98*0.98))*0.003))*100  ;
	if(mode==3) return ((1.28694e-03*(1+sigma/100.)-1.36634e-03)/1.36634e-03)*100;
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

	float sigma = 1.28694e-03*(1-systpar.sigma/100.); //1.15694e-03*(1-systpar.sigma/100.);
	sigma+=(float)((2*fabs(1.28694e-03-sigma)/(float)systpar.steps)*stepsigma);

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

	float masstocharge=0;
	if(M_gen<2) masstocharge=M_gen/1;
	else masstocharge = M_gen/2;

	float beta_inv = 1/Beta;
	float sigmasmear = 0.0284/2;
	float shiftstart=-systpar.shift/masstocharge;
	float sigmastart = sigmasmear - sigmasmear*systpar.sigma/100;
	float sigma = sigmastart +  ((2*systpar.sigma/(systpar.steps))*stepsigma)/100*sigmasmear;
	float shift = (- systpar.shift*1e-4 + (2*systpar.shift*1e-4/(float)systpar.steps)*stepshift)/masstocharge;	
	float norm = Rand->Rndm();
	
	float MCsigma;
	if(M_gen<2) MCsigma = /*(0.0318/0.0335)*/MCmodel->Eval(1/Beta_gen);
	else MCsigma = (0.0155/0.0169)*MCmodel->Eval(1/Beta_gen);
	
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
	//IsPhysTrig(vars) && IsDownGoing(vars) && IsGoodTrack(vars) && IsGoodChi2(vars) && IsCharge1Track(vars) && IsGoodKalman(vars)
	//if(vars->Event==717) cout<<"ECCO! "<<vars->R<<" "<<vars->Beta<<" "<<vars->Chisquare<<" "<<vars->Chisquare_y<<" "<< ApplyCuts("IsPositive&IsPhysTrig&IsGoodKalman",vars)<<endl;
	if((ApplyCuts(cutP,vars)||ApplyCuts(cutD,vars)||ApplyCuts(cutHe,vars))){

	int imin,imax,jmin,jmax;
	imin=0; jmin=0; imax=systpar.steps; jmax=systpar.steps;
	if(IsExtern) {imin=0; jmin=5; imax=1; jmax=6;}
	
	for(int i=imin;i<imax;i++)
			for(int j=jmin;j<jmax;j++){
		
	
				float betasmear=0;
				float Rsmear = vars->R;

				if(!isrich) 	
					if(systpar.mode==0) betasmear = SmearBeta(vars->Beta,vars->Massa_gen,(float)i,(float)j,vars->R);    
					if(systpar.mode==1) betasmear = SmearBeta_v2(vars->Beta,GetBetaGen_cpct(vars),vars->Massa_gen,MCmodel,(float)i,(float)j);    
				else if(ApplyCuts("IsFromNaF",vars))	 
					if(systpar.mode==2) betasmear = SmearBetaRICH(vars->BetaRICH_new,(float)i,(float)j);
                                        if(systpar.mode==3) { 
						betasmear = SmearBetaRICH_v2(vars->BetaRICH_new,GetBetaGen_cpct(vars),(float)i,(float)j); 
						Rsmear = SmearRRICH_R(vars->R,GetGenMomentum(vars),(float)i,(float)j);
						}
				else if(ApplyCuts("IsFromAgl",vars))  	
					if(systpar.mode==2) betasmear = SmearBetaRICH(vars->BetaRICH_new,(float)i,(float)j);
                                        if(systpar.mode==3) {
						betasmear = SmearBetaRICH_v2(vars->BetaRICH_new,GetBetaGen_cpct(vars),(float)i,(float)j); 
						Rsmear = SmearRRICH_R(vars->R,GetGenMomentum(vars),(float)i,(float)j);
						}
					
				float mctotalweight = 
						   vars->mcweight 
						   * vars->GetCutoffCleaningWeight(GetRFromBeta(bins.getParticle().getMass(),betasmear),GetGenMomentum(vars),BetacutoffCut) 
						   * vars->GetTimeDepWeight(vars->R)
						   ;
 
				int kbin;
	
				if(!IsFitNoise){
				   if(bins.IsUsingBetaEdges()){	 
					kbin = bins.GetBin(betasmear);
					float scaledR;
					float mass;
					float rigcut;
					if(vars->Massa_gen<2) rigcut=0.55;
					else rigcut=1.05; 
					scaledR=Rsmear*template1scalefactor1;
					mass = scaledR/betasmear * pow((1-pow(betasmear,2)),0.5);
					if(ApplyCuts(cutP,vars)&&kbin>=0&&scaledR>=rigcut) { 
				//		if(i==0&&j==5) cout<<vars->R<<" "<<vars->RInner<<" "<<betasmear<<" "<<vars->Event<<endl;

					 fits[kbin][i][j]->Templ_P->Fill(mass,mctotalweight);}
					scaledR=Rsmear*template1scalefactor2;
					mass = scaledR/betasmear * pow((1-pow(betasmear,2)),0.5);
					if(ApplyCuts(cutD,vars)&&kbin>=0&&scaledR>=rigcut)  {fits[kbin][i][j]->Templ_D->Fill(mass,mctotalweight);}
					scaledR=Rsmear*template1scalefactor3;
					mass = scaledR/betasmear * pow((1-pow(betasmear,2)),0.5);
					if(ApplyCuts(cutHe,vars)&&kbin>=0&&scaledR>=rigcut) fits[kbin][i][j]->Templ_He->Fill(mass,mctotalweight);
					}
				    else {
				   	kbin = bins.GetBin(discr_var(vars));
					if(ApplyCuts(cutP,vars)&&kbin>=0)   fits[kbin][i][j]->Templ_P->Fill(betasmear,mctotalweight);
                                        if(ApplyCuts((cutD).c_str(),vars)&&kbin>=0)   fits[kbin][i][j]->Templ_D->Fill(betasmear,mctotalweight);
                                        if(ApplyCuts(cutHe,vars)&&kbin>=0) fits[kbin][i][j]->Templ_He->Fill(betasmear,mctotalweight);

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

	return;
}



void TemplateFIT::SaveFitResults(FileSaver finalhistos){

	cout<<"Saving best"<<endl;	
	//save best
	for(int bin=0;bin<bins.size();bin++){
		TH1F * OriginalP=(TH1F*)fits[bin][0][5]->Templ_P->Clone();
		TH1F* BestP;
		if(!IsLocalFit&&!IsLocalConstrainedFit) {
				BestP=(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_P->Clone();
				BestP->SetName(("Best #chi^{2}: " + to_string(BestChiSquare->i) + "_" + to_string(BestChiSquare->j)).c_str());
			}	
		else {
			BestChi * bestlocal = new BestChi();
			if(IsLocalConstrainedFit)  bestlocal->FindMinimum(TFitChisquare[bin],BestChiSquare->i);
			else bestlocal->FindMinimum(TFitChisquare[bin]);
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
		if(!IsLocalFit&!IsLocalConstrainedFit){
				 BestD=(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_D->Clone();
				BestD->SetName(("Best #chi^{2}: " + to_string(BestChiSquare->i) + "_" + to_string(BestChiSquare->j)).c_str());
			}
		else {
			BestChi * bestlocal = new BestChi();
			if(IsLocalConstrainedFit)  bestlocal->FindMinimum(TFitChisquare[bin],BestChiSquare->i);
			else bestlocal->FindMinimum(TFitChisquare[bin]);
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
		
		TH1F * BestHe=0x0;
		TH1F * OriginalHe=(TH1F*)fits[bin][0][5]->Templ_He->Clone();
		BestHe=(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_He->Clone();
		TH1F * BestHePrim =(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_HePrim->Clone();
		OriginalHe->SetName("Original Tritium MC ");
		BestHe->SetName("Best #chi^{2} Tritium MC ");
		BestHePrim->SetName("Best #chi^{2} Overcutoff Triton MC ");
		finalhistos.Add(OriginalHe);
		finalhistos.Add(BestHe);
		finalhistos.Add(BestHePrim);
		finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesHe/Bin"+to_string(bin)).c_str());

		TH1F * Residuals = (TH1F*) fits[bin][0][5]->Data->Clone(("Residuals " + to_string(bin)).c_str());
		TH1F * model = (TH1F*) BestP->Clone("model");
		model->Add(BestD);
		if(BestHe) if(BestHe->Integral()>0) model->Add(BestHe);
	//	for(int i=0;i<model->GetNbinsX();i++) model->SetBinError(i+1,2*model->GetBinError(i+1));
		Residuals->Add(model,-1);
		for(int i=0;i<Residuals->GetNbinsX();i++) {
                if(Residuals->GetBinError(i+1)>0 && fits[bin][0][5]->Data->GetBinContent(i+1)>0 && model->GetBinContent(i+1)>0 ){
                        Residuals->SetBinContent(i+1,Residuals->GetBinContent(i+1)/Residuals->GetBinError(i+1));
                        Residuals->SetBinError(i+1,1);
                }
                else {
                        Residuals->SetBinContent(i+1,0);
                        Residuals->SetBinError(i+1,0);
                        }
        	}
		
		TGraphErrors * ResidualGraph = new TGraphErrors();
		
		std::string binlow = to_string(bins.RigTOIBins()[bin]);
		std::string binhigh = to_string(bins.RigTOIBins()[bin+1]);
		binlow.erase(4,15);
		binhigh.erase(4,15);
		std::string title = (binlow + " < R < " + binhigh + " GV").c_str();
		ResidualGraph->SetName(("Residuals: " + title).c_str());
		ResidualGraph->SetTitle(("Residuals: " + title).c_str());
	
		int n=0;
		for(int i=0;i<Residuals->GetNbinsX();i++) if(Residuals->GetBinContent(i+1)!=0) { 
			ResidualGraph->SetPoint(n,Residuals->GetBinCenter(i+1),Residuals->GetBinContent(i+1));
			ResidualGraph->SetPointError(n,0,Residuals->GetBinError(i+1));
			n++;
		}
		finalhistos.Add(ResidualGraph);
		finalhistos.writeObjsInFolder((basename+"/Fit Results/Residuals").c_str());
	
		finalhistos.Add(DrawTemplateFit(bin,fits[bin][0][5]->Data,BestP,BestD,BestHe,Residuals));
		finalhistos.writeObjsInFolder((basename+"/Fit Results/DrawFits").c_str());

	}


	




	cout<<"Saving data"<<endl;	
	
	//save data
	for(int bin=0;bin<bins.size();bin++){ 
		finalhistos.Add(fits[bin][0][5]->Data);
		finalhistos.Add(fits[bin][0][5]->DataPrim);
		finalhistos.Add(fits[bin][0][5]->DataAmbient);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/Data/Bin"+to_string(bin)).c_str());	
	}
	////////////////////////////////////////////////////////	
/*	cout<<"Saving all"<<endl;	
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
*/	/////////////////////////////////////////////////////////
	cout<<"Saving fits"<<endl;	
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
		BestSigma   = ConvertBinnedHisto(BestSigma,"Best Sigma",bins,false);
		finalhistos.Add(BestSigma);
		TH1D * BestShift = (TH1D*) Global_ChiSquare->ProjectionY("Best Shift",BestChiSquare->i+1,BestChiSquare->i+1);
		BestShift   = ConvertBinnedHisto(BestShift,"Best Shift",bins,false);
		finalhistos.Add(BestShift);
		BestChiSquares->Scale(fits[0][BestChiSquare->i][BestChiSquare->j]->ndf);
		OriginalChiSquares->Scale(fits[0][BestChiSquare->i][BestChiSquare->j]->ndf);

		BestChiSquares   = ConvertBinnedHisto((TH1F*)BestChiSquares,"Best ChiSquare",bins,false); 
		OriginalChiSquares   = ConvertBinnedHisto((TH1F*)OriginalChiSquares,"Original ChiSquare",bins,false); 
	
		finalhistos.Add(BestChiSquares   ); 
		finalhistos.Add(OriginalChiSquares);
	}

	ProtonCounts   = ConvertBinnedHisto(ProtonCounts,"Proton Counts",bins,true); 
	DeuteronCounts = ConvertBinnedHisto(DeuteronCounts,"Deuteron Counts",bins,true); 
	TritiumCounts = ConvertBinnedHisto(TritiumCounts,"Tritium Counts",bins,true); 


	TH1F * RatioDP = (TH1F*) DeuteronCounts->Clone("RatioPD");
	RatioDP->Divide(ProtonCounts);

	finalhistos.Add(StatErrorP);
	finalhistos.Add(StatErrorD);
	finalhistos.Add(StatErrorT);
	finalhistos.Add(SystError);
	finalhistos.Add(ProtonCounts);
	finalhistos.Add(DeuteronCounts);
	finalhistos.Add(TritiumCounts);
	finalhistos.Add(RatioDP);
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

float GetChiSquare(TH1 * result, TH1 * data, float min, float max,bool prioritizetail){

	TH1 * Result = (TH1*) result->Clone();
	TH1 * Data   = (TH1*) data->Clone();

	int binmin = Data->FindBin(min);
	int binmax = Data->FindBin(max);
	float chi = 0;
	float err = 0;
	int ndf =0;
	if(prioritizetail) cout<<"PRIORITIZING TAIL!"<<endl;
	//core
	for(int i=binmin; i<data->FindBin(4); i++){
			if(Result->GetBinError(i+1)<=0||Data->GetBinError(i+1)<=0) continue;
			float tailmult=1;
			if(Data->GetBinCenter(i+1)>3&&Data->GetBinCenter(i+1)<5&&prioritizetail) tailmult=4;
			if(Data->GetBinCenter(i+1)>5&&Data->GetBinCenter(i+1)<10&&prioritizetail) tailmult=0.5;
			err = pow(pow(Data->GetBinError(i+1),2) + pow(Result->GetBinError(i+1),2),0.5);
			chi += tailmult * pow((Data->GetBinContent(i+1) - Result->GetBinContent(i+1)),2)/pow(err,2);
			ndf++;
	}

	float totalcontent_DT=0;
	float totalcontent_MC=0;
	float totalerror=0;
	for(int i=Data->FindBin(4); i<binmax; i++){
		totalcontent_DT+=Data->GetBinContent(i+1);		
		totalcontent_MC+=Result->GetBinContent(i+1);		
		totalerror+=pow(Data->GetBinError(i+1),2)+pow(Result->GetBinError(i+1),2);		
	}
	chi+=pow(totalcontent_DT-totalcontent_MC,2)/totalerror;
	ndf++;
	ndf = ndf-3;
	chi /= ndf;
	cout<<"chi: "<<chi<<" ndf: "<<ndf<<endl;
	return chi;

}


void Do_TemplateFIT(TFit * Fit,float fitrangemin,float fitrangemax,float constrain_min[], float constrain_max[], bool isfitnoise, bool prioritizetail, bool IsFitPrim){
	if(isfitnoise) cout<<"********** FIT NOISE MODE **********"<<endl;
	if(isfitnoise) Fit->Templ_Noise->Scale(0.01*Fit ->  Templ_P->Integral()/Fit ->  Templ_Noise->Integral());
	TObjArray *Tpl;
	Tpl = new TObjArray(3);


	if(!IsFitPrim) {
		if(Fit ->  Templ_P)  if(Fit ->  Templ_P ->Integral()>0) Tpl -> Add( Fit ->  Templ_P );
		if(Fit ->  Templ_D)  if(Fit ->  Templ_D ->Integral()>0) Tpl -> Add( Fit ->  Templ_D );
        	if(Fit ->  Templ_He) if(Fit ->  Templ_He->Integral()>0) Tpl -> Add( Fit ->  Templ_He);
	}
	else {
		if(Fit ->  Templ_PPrim)  if(Fit ->  Templ_PPrim ->Integral()>0) Tpl -> Add( Fit ->  Templ_PPrim );
		if(Fit ->  Templ_DPrim)  if(Fit ->  Templ_DPrim ->Integral()>0) Tpl -> Add( Fit ->  Templ_DPrim );
        	if(Fit ->  Templ_HePrim) if(Fit ->  Templ_HePrim->Integral()>0) Tpl -> Add( Fit ->  Templ_HePrim);
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
			if(Fit-> Templ_He->Integral()>0) Sum -> Add(Fit ->  Templ_He);
			if(isfitnoise) Sum -> Add(Fit ->  Templ_Noise);			
			for(int n=0;n<Sum->GetNbinsX();n++) if(Sum->GetBinContent(n+1)==0) Sum->SetBinError(n+1,1);
	
			if(( df1dw1*df1dw1*Cov00 + df1dw2*df1dw2*Cov11 + df1dw3*df1dw3*Cov22 +2*df1dw1*df1dw2*Cov01 +2*df1dw1*df1dw3*Cov02+ 2*df1dw2*df1dw3*Cov12)>0){
				Fit -> StatErrP =  pow(( df1dw1*df1dw1*Cov00 + df1dw2*df1dw2*Cov11 + df1dw3*df1dw3*Cov22 +2*df1dw1*df1dw2*Cov01 +2*df1dw1*df1dw3*Cov02+ 2*df1dw2*df1dw3*Cov12)/2,0.5);//w1;
				Fit -> StatErrD =  pow(( df2dw1*df2dw1*Cov00 + df2dw2*df2dw2*Cov11 + df2dw3*df2dw3*Cov22 +2*df2dw1*df2dw2*Cov01 +2*df2dw1*df2dw3*Cov02+ 2*df2dw2*df2dw3*Cov12)/2,0.5);//w2;
				Fit -> StatErrT =  pow(( df3dw1*df3dw1*Cov00 + df3dw2*df3dw2*Cov11 + df3dw3*df3dw3*Cov22 +2*df3dw1*df3dw2*Cov01 +2*df3dw1*df3dw3*Cov02+ 2*df3dw2*df3dw3*Cov12)/2,0.5);//w3;
			}
			else{

				Fit -> StatErrP = 0.5*e1/w1;
				Fit -> StatErrD = 0.5*e2/w2;			
				Fit -> StatErrT = 0.5*e3/w3;
			}
			cout<<"fract: "<<w1<<" "<<w2<<" "<<w3<<endl;
			cout<<endl;
			cout<<"err: "  <<e1<<" "<<e2<<" "<<e3<<endl;
			cout<<"er c:"  <<Fit -> StatErrP<<" "<<Fit -> StatErrD<<" "<<Fit -> StatErrT<<endl;	
			
			//Fit -> ChiSquare = Fit -> Tfit -> GetChisquare()/(float) (Fit ->  Tfit -> GetNDF());
			Fit -> ChiSquare = GetChiSquare(Sum,Fit->Data,min,max,prioritizetail);
			Fit->ndf = Fit ->  Data->FindBin(4.5)-Fit ->Data->FindBin(0.85) -3;
			Fit -> DCounts = Fit ->  Templ_D -> Integral();
			Fit -> PCounts = Fit ->  Templ_P -> Integral();
			if(Fit-> Templ_He->Integral()>0) Fit -> TCounts = Fit ->  Templ_He -> Integral();
			else Fit -> TCounts = 0;
			//Fit -> DCountsPrim = Fit ->  Templ_DPrim -> Integral();
			//Fit -> PCountsPrim = Fit ->  Templ_PPrim -> Integral();
			cout<<"HOFINITO: "<<Fit ->  Templ_DPrim<<" "<<Fit ->  Templ_PPrim<<" ndf: "<<Fit->ndf<<endl;
		

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

float EvalFitProbability(float chi,float chimin){

/*	TF1 * CumulativeChiSquare=new TF1("Chi","exp(-x)*x^-0.5",0.001,10000);
	if(chi>0.01&&chi<10000)
		return CumulativeChiSquare->Integral(chi,10000);
	else return 0;
*/
	
	return 1/fabs(chimin-chi+0.5);

}

void  TemplateFIT::RebinAll(int f){
	for(int bin=0;bin<bins.size();bin++){
	 for(int sigma=0;sigma<systpar.steps;sigma++)
		 for(int shift=0;shift<systpar.steps;shift++){
			for(int i=0;i<fits[bin][sigma][shift] -> Data->GetNbinsX();i++){
			if(fits[bin][sigma][shift] -> Data->GetBinContent(i+1)<10)
			 fits[bin][sigma][shift] -> Data->SetBinContent(i+1,0);
			 fits[bin][sigma][shift] -> Data->SetBinError(i+1,0);
			}
		 }
	}
};


float FindConstraintD(TH1F * dt, TH1F * mcp, TH1F * mcd){
	float corrp = mcp->Integral()/mcp->Integral(mcp->FindBin(0),mcp->FindBin(0.938));  
	float corrd = mcd->Integral()/mcd->Integral(mcd->FindBin(1.875),mcd->FindBin(4));  


	float counts_p=dt->Integral(dt->FindBin(0),dt->FindBin(0.938));
	float counts_d=dt->Integral(dt->FindBin(1.875),dt->FindBin(4));

	return corrd*counts_d/(corrp*counts_p);
}

void TemplateFIT::ExtractCounts(FileSaver finalhistos,int force_shift){

	SimpleExtractPrimaries();

	for(int bin=0;bin<bins.size();bin++){
		
		TH1F * DD = (TH1F*) fits[bin][0][5]->Data->Clone("tempdata");
		TH1F * TT = (TH1F*) fits[bin][0][5]->Templ_P->Clone("temptemp");
	
		float eq = 1;
		if(bin<=14&&(systpar.mode==3)) eq= EqualizeTemplate(DD, TT);			
		TH1F * DD2 = (TH1F*) fits[bin][0][6]->Data->Clone("tempdata2");


			
		for(int sigma=0;sigma<systpar.steps;sigma++){
			for(int shift=0;shift<systpar.steps;shift++){
		 	float constrain_D = FindConstraintD (DD2,fits[bin][0][5]->Templ_P,fits[bin][0][5]->Templ_D);
	
			cout<<"************ CONSTRAIN FIT: "<<constrain_D<<endl;
	
	
			// MAIN Template Fit 
				if(!fitDisabled) {
					cout<<endl;
					cout<<"Bin: "<<bin<<": "<<sigma<<" "<<shift<<endl;
					fits[bin][sigma][shift]->RegularizeTemplateError();	
					if(!isrich) fits[bin][sigma][shift]->AddSysttoTempl();
					if((systpar.mode==2)&&bin>4&&bin<9) fits[bin][sigma][shift]->AdjoustProtonTail(0.3);
					bool prioritizetail=false;	
					/*if(bin<=4&&(systpar.mode==3)) fits[bin][sigma][shift]->Templ_P = SimpleShiftHisto(fits[bin][sigma][shift]->Templ_P,-2.7*eq);
					if(bin<=4&&(systpar.mode==3)) fits[bin][sigma][shift]->Templ_D = SimpleShiftHisto(fits[bin][sigma][shift]->Templ_D,+2.7*eq);
					if(bin<=4&&(systpar.mode==3)) fits[bin][sigma][shift]->Templ_He = SimpleShiftHisto(fits[bin][sigma][shift]->Templ_He,-eq);*/
					cout<<"Eq: "<<eq<<endl; 
					//if(bin>=6&&bin<8&&(systpar.mode==3))  fits[bin][sigma][shift]->Templ_D = SimpleShiftHisto(fits[bin][sigma][shift]->Templ_D,-0.06);
					//if(bin>=3&&bin<8&&(systpar.mode==3))  fits[bin][sigma][shift]->Templ_He = SimpleShiftHisto(fits[bin][sigma][shift]->Templ_He,+0.06);
					//if((systpar.mode==2)) fits[bin][sigma][shift]->Templ_D = SimpleShiftHisto(fits[bin][sigma][shift]->Templ_D,-0.06);
					//if((systpar.mode==3)&&bin>4) fits[bin][sigma][shift]->Templ_D = SimpleShiftHisto(fits[bin][sigma][shift]->Templ_D,-0.03);
					/// NOISE AGL			
					//if(bin>=6&&bin<=8&&(systpar.mode==3)) AddConstantNoise(fits[bin][sigma][shift]->Data,fits[bin][sigma][shift]->Templ_P);	
				
			
					//Empirical constrain AGL
			/*		float normtime=1;
					BestChi * besty1 = new BestChi();
					BestChi * besty2 = new BestChi();
					BestChi * besty3 = new BestChi();
					if(bin>5&&(systpar.mode==3)){ 
						besty1->FindMinimum(TFitChisquare[5]);
						besty2->FindMinimum(TFitChisquare[4]);
						besty3->FindMinimum(TFitChisquare[3]);
						normtime = fits[5][besty1->i][besty1->j]->DCounts/(fits[5][besty1->i][besty1->j]->PCounts)
							 + fits[4][besty2->i][besty2->j]->DCounts/(fits[4][besty2->i][besty2->j]->PCounts)
							 + fits[3][besty3->i][besty3->j]->DCounts/(fits[3][besty3->i][besty3->j]->PCounts);
						normtime/=3;
						cout<<"*********************** NORMTIME: "<<normtime<<endl;
						normtime/=0.008124;

					}
					if(bin==6&&(systpar.mode==3)) constrainmin[1]*=normtime*1.12;
					if(bin==7&&(systpar.mode==3)) constrainmin[1]*=normtime*1.10;
					if(bin==8&&(systpar.mode==3)) constrainmin[1]*=normtime*1.05;*/
				/* 	if(bin>5&&bin<=8&&(systpar.mode==3)&&IsLocalConstrainedFit){
						BestChi * besty2 = new BestChi();
						besty2->FindMinimum(TFitChisquare[5]);
						BestChi * besty1 = new BestChi();
						besty1->FindMinimum(TFitChisquare[2]);
						float y1=log(fits[1][besty1->i][besty1->j]->DCounts/fits[1][besty1->i][besty1->j]->PCounts);
						
						float deriv = fabs((y2-y1)/(x2-x1));
						float mind = exp(y3-deriv*fabs(x4-x3));
						constrainmin[1]=mind; 	
						cout<<" CONSTRAINMIN "<< fabs(x4-x3) <<" "<<deriv <<" "<<mind  <<endl;		
					}*/

					float constrain_deut=constrainmin[1];
			
					
 			/*		if(systpar.mode==2&&bin<7){
					constrainmin[1]=0.95*constrain_D;
					constrainmax[1]=1.3*constrain_D;}

					if(systpar.mode==2&&bin>=7&&bin<9){
					constrainmin[1]=0.9*constrain_D;
					constrainmax[1]=1.3*constrain_D;}
				       

					if(systpar.mode==3&&bin>3&&bin<7){
					constrainmin[1]=0.95*constrain_D;
					constrainmax[1]=1.3*constrain_D;}
			*/
				Do_TemplateFIT(fits[bin][sigma][shift],fits[bin][sigma][shift]->fitrangemin,fits[bin][sigma][shift]->fitrangemax,constrainmin,constrainmax,IsFitNoise,prioritizetail);
		    		constrainmin[1]=constrain_deut; 
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
	Global_ChiSquare = (TH2F*) TFitChisquare[TFitChisquare.size()-1]->Clone("Global ChiSquare");
	for(int i=TFitChisquare.size();i<TFitChisquare.size()-1;i--) Global_ChiSquare->Add(TFitChisquare[i]);
	BestChiSquare = new BestChi();
	BestChiSquare->FindMinimum(Global_ChiSquare);

	for(int bin=0;bin<bins.size();bin++){
	BestChi * bestlocal = new BestChi();
	if(IsLocalConstrainedFit)  bestlocal->FindMinimum(TFitChisquare[bin],BestChiSquare->i);
	else bestlocal->FindMinimum(TFitChisquare[bin]);

	BestChi * besterror = new BestChi();
	besterror->FindMinimum(TFitChisquare[bin]);

	TH1F * weighteddcounts  = new TH1F(("WeightedCountsBin " +to_string(bin)).c_str(),("WeightedCountsBin " +to_string(bin)+";Counts;").c_str(),100,0.5*fits[bin][bestlocal->i][bestlocal->j]->DCounts,2*fits[bin][bestlocal->i][bestlocal->j]->DCounts);
		for(int sigma=0;sigma<systpar.steps;sigma++)
			for(int shift=0;shift<systpar.steps;shift++)
				if(fits[bin][sigma][shift]->DCounts&&fits[bin][0][5]->DCounts){
					float totcounts=fits[bin][sigma][shift]->DCounts+fits[bin][sigma][shift]->Templ_P->Integral(fits[bin][sigma][shift]->Templ_P->FindBin(2),fits[bin][sigma][shift]->Templ_P->FindBin(3));
					weighteddcounts -> Fill(totcounts,EvalFitProbability(fits[bin][sigma][shift]->ChiSquare,TFitChisquare[bin]->GetBinContent(besterror->i+1,besterror->j+1) ));
					
					}
		WeightedDCounts.push_back(weighteddcounts);
			
	}


	EvalFinalParameters();
	EvalFinalErrors();
	CalculateFinalPDCounts();	
	return;
}


void TemplateFIT::EvalFinalParameters(){
	for(int bin=0;bin<bins.size();bin++){

		if(!IsLocalFit&&!IsLocalConstrainedFit) BestChiSquares     ->SetBinContent(bin+1,TFitChisquare[bin]->GetBinContent(BestChiSquare->i+1,BestChiSquare->j+1));
		else{
			BestChi * bestlocal = new BestChi();
			if(IsLocalConstrainedFit)  bestlocal->FindMinimum(TFitChisquare[bin],BestChiSquare->i);
			else bestlocal->FindMinimum(TFitChisquare[bin]);
			BestChiSquares     ->SetBinContent(bin+1,TFitChisquare[bin]->GetBinContent(bestlocal->i+1,bestlocal->j+1));
		}
		if(TFitChisquare[bin]->GetBinContent(1,6)<1000) OriginalChiSquares ->SetBinContent(bin+1,TFitChisquare[bin]->GetBinContent(1,6));
		BestChiSquares     ->SetBinError(bin+1,0.25);	
                OriginalChiSquares ->SetBinError(bin+1,0.25);
	}
	return;
}

void TemplateFIT::EvalFinalErrors(){

	for(int bin=0;bin<bins.size();bin++){
		BestChi * bestlocal = new BestChi();
	if(IsLocalConstrainedFit)  bestlocal->FindMinimum(TFitChisquare[bin],BestChiSquare->i);
	else bestlocal->FindMinimum(TFitChisquare[bin]);
				
		if(fits[bin][bestlocal->i][bestlocal->j]->StatErrP>0 
	        && fits[bin][bestlocal->i][bestlocal->j]->StatErrD>0
		&& fits[bin][bestlocal->i][bestlocal->j]->DCounts>0){
			StatErrorP -> SetBinContent(bin+1,1.2*fits[bin][bestlocal->i][bestlocal->j]->StatErrP);
			StatErrorD -> SetBinContent(bin+1,max(1.2*fits[bin][bestlocal->i][bestlocal->j]->StatErrD,0.002));
			StatErrorT -> SetBinContent(bin+1,1.2*fits[bin][bestlocal->i][bestlocal->j]->StatErrT);
			if(IsLocalConstrainedFit) SystError -> SetBinContent(bin+1,std::max(WeightedDCounts[bin]->GetStdDev()/fits[bin][bestlocal->i][bestlocal->j]->DCounts,0.045));
			else SystError -> SetBinContent(bin+1,0.5*WeightedDCounts[bin]->GetStdDev()/fits[bin][bestlocal->i][bestlocal->j]->DCounts);
			StatErrorP -> SetBinError(bin+1,0);
			StatErrorD -> SetBinError(bin+1,0);
			StatErrorT -> SetBinError(bin+1,0);
			SystError -> SetBinError(bin+1,0);
		}

	}

	StatErrorP->Smooth();
	StatErrorD->Smooth();
	StatErrorT->Smooth();


	return;

}


void TemplateFIT::CalculateFinalPDCounts(){

	HeContError = new TH1F ("HeContError","HeContError",bins.size(),0,bins.size());	
	MeasuredHeContRatio = new TH1F("Measured He over P","Measured He over P",bins.size(),0,bins.size());

	for(int bin=0;bin<bins.size();bin++){

		BestChi * bestlocal = new BestChi();
		if(IsLocalConstrainedFit)  bestlocal->FindMinimum(TFitChisquare[bin],BestChiSquare->i);
		else bestlocal->FindMinimum(TFitChisquare[bin]);

		float staterrP= StatErrorP -> GetBinContent(bin+1) * fits[bin][bestlocal->i][bestlocal->j]->PCounts;
		float staterrD= StatErrorD -> GetBinContent(bin+1) * fits[bin][bestlocal->i][bestlocal->j]->DCounts;
		float staterrT= StatErrorT -> GetBinContent(bin+1) * fits[bin][bestlocal->i][bestlocal->j]->PCounts;
	
		float conterrP= 0.1 * fits[bin][bestlocal->i][bestlocal->j]->TCounts;
		float conterrD= 0.1 * fits[bin][bestlocal->i][bestlocal->j]->TCounts;
		float conterrT= 0.1 * fits[bin][bestlocal->i][bestlocal->j]->TCounts;

		float systerrP= SystError -> GetBinContent(bin+1) * fits[bin][bestlocal->i][bestlocal->j]->PCounts;
		float systerrD= SystError -> GetBinContent(bin+1) * fits[bin][bestlocal->i][bestlocal->j]->DCounts;
		float systerrT= SystError -> GetBinContent(bin+1) * fits[bin][bestlocal->i][bestlocal->j]->TCounts;

		float CountsP; 	
		if(IsLocalFit||IsLocalConstrainedFit)  CountsP      = fits[bin][bestlocal->i][bestlocal->j]->PCounts - 8.5 * fits[bin][bestlocal->i][bestlocal->j]->TCounts;
		else CountsP      = fits[bin][BestChiSquare->i][BestChiSquare->j]->PCounts- 8.5 * fits[bin][bestlocal->i][bestlocal->j]->TCounts;


		float CountsD;	
		if(IsLocalFit||IsLocalConstrainedFit)  CountsD      = fits[bin][bestlocal->i][bestlocal->j]->DCounts - 1.5 * fits[bin][bestlocal->i][bestlocal->j]->TCounts;
		else CountsD      = fits[bin][BestChiSquare->i][BestChiSquare->j]->DCounts - 1.5 * fits[bin][bestlocal->i][bestlocal->j]->TCounts;

		float CountsT      = fits[bin][bestlocal->i][bestlocal->j]->TCounts;
		float CountsP_prim = fits[bin][bestlocal->i][bestlocal->j]->PCountsPrim - 8.5 * fits[bin][bestlocal->i][bestlocal->j]->TCounts;
		float CountsD_prim = fits[bin][bestlocal->i][bestlocal->j]->DCountsPrim - 1.5 * fits[bin][bestlocal->i][bestlocal->j]->TCounts;

		HeContError -> SetBinContent(bin+1,conterrD);

		ProtonCounts->SetBinContent(bin+1,CountsP);
		ProtonCounts->SetBinError(bin+1,pow(pow(staterrP,2)+pow(systerrP,2)+pow(conterrP,2),0.5));

		DeuteronCounts->SetBinContent(bin+1,CountsD);
		DeuteronCounts->SetBinError(bin+1,pow(pow(staterrD,2)+pow(systerrD,2)+pow(conterrD,2),0.5));		

		TritiumCounts->SetBinContent(bin+1,CountsT);
		TritiumCounts->SetBinError(bin+1,pow(pow(staterrT,2)+pow(systerrT,2)+pow(conterrT,2),0.5));		


		if(fits[bin][bestlocal->i][bestlocal->j]->PCountsPrim>0){
			ProtonCountsPrim->SetBinContent(bin+1,fits[bin][bestlocal->i][bestlocal->j]->PCountsPrim);
			ProtonCountsPrim->SetBinError(bin+1,pow(pow(staterrP,2)+pow(systerrP,2)+pow(conterrP,2),0.5) * ProtonCountsPrim->GetBinContent(bin+1)/ProtonCounts->GetBinContent(bin+1) );
		}
		if(fits[bin][bestlocal->i][bestlocal->j]->DCountsPrim>0){
			DeuteronCountsPrim->SetBinContent(bin+1,fits[bin][bestlocal->i][bestlocal->j]->DCountsPrim);
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
