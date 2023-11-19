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

	TH1F * residuals;
	if(Residuals->GetEntries()>0) residuals=(TH1F*)Residuals->Clone();
	else	{
		residuals=(TH1F*)Data->Clone("residuals");
		TH1F* model=(TH1F*) BestP->Clone("model");
		model->Add(BestD);
		model->Add(BestHe);
		residuals->Add(model,-1);
		for(int i=0;i<model->GetNbinsX();i++){
			if(residuals->GetBinError(i+1)>0&&Data->GetBinContent(i+1)>0 && model->GetBinContent(i+1)>0){
			residuals->SetBinContent(i+1,residuals->GetBinContent(i+1)/residuals->GetBinError(i+1));
			residuals->SetBinError(i+1,1);
			}
			else{
			residuals->SetBinContent(i+1,0);
			residuals->SetBinError(i+1,1);
			}
		}
	}

	residuals->SetMarkerStyle(8);
        residuals->SetMarkerColor(1);
        //residuals->GetXaxis()->SetRangeUser(0.5,6.5);
        residuals->Draw("P");


	return c1;

}

float EvalFitProbability(float chi,float chimax){

/*	TF1 * CumulativeChiSquare=new TF1("Chi","exp(-x)*x^-0.5",0.001,10000);
	if(chi>0.01&&chi<10000)
		return CumulativeChiSquare->Integral(chi,10000);
	else return 0;
*/
	
	return fabs(chi-chimax)/1000;
}

double EvalFitProbability(float chi,int x, int y, float centroid_x, float centroid_y, float chimin){

	if(x==1) return 0;	
	float centroidweight=1;
//	return 1/(chi * centroidweight );
        chi/=chimin;
	double f=pow(2.71,-fabs(chi-1)/2);

	return f;	

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

TFit* TFit::Clone(){
	TFit * fit = new TFit();

	if(Templ_P 	) fit->  Templ_P 		=(TH1F*)Templ_P 	->Clone();
	if(Templ_D 	) fit->  Templ_D 		=(TH1F*)Templ_D 	->Clone();
	if(Templ_DPrim	) fit->  Templ_DPrim		=(TH1F*)Templ_DPrim	->Clone();
	if(Templ_PPrim	) fit->  Templ_PPrim		=(TH1F*)Templ_PPrim	->Clone();
	if(Templ_He	) fit->  Templ_He		=(TH1F*)Templ_He	->Clone();
	if(Templ_HePrim	) fit->  Templ_HePrim		=(TH1F*)Templ_HePrim	->Clone();
	if(Templ_Noise	) fit->  Templ_Noise		=(TH1F*)Templ_Noise	->Clone(); 	
	if(Data 	) fit->  Data 			=(TH1F*)Data 		->Clone();
	if(DataAmbient 	) fit->  DataAmbient 		=(TH1F*)DataAmbient 	->Clone();
	if(DataPrim	) fit->  DataPrim		=(TH1F*)DataPrim	->Clone();	

	fit->fitrangemin = fitrangemin;
	fit->fitrangemax = fitrangemax;
	fit->DCounts = DCounts;
	fit->PCounts = PCounts;
	fit->TCounts = TCounts;
	fit -> PCountsPrim =PCountsPrim;
	fit -> DCountsPrim =DCountsPrim;

	return fit;
}



	
void TFit::RegularizeTemplateError(){
/*	for(int i=0;i<Templ_P->GetNbinsX();i++) {
		Templ_P->SetBinError(i+1,pow(Templ_P->GetBinContent(i+1),0.5));
	}
	for(int i=0;i<Templ_D->GetNbinsX();i++) {
		Templ_D->SetBinError(i+1,pow(Templ_D->GetBinContent(i+1),0.5));
	}*/
	for(int i=0;i<Templ_He->GetNbinsX();i++) {
		Templ_He->SetBinError(i+1,pow(Templ_He->GetBinContent(i+1),0.5));
	}
	float integral = Templ_P->Integral();
	float mincontent=9999999;
	for(int i=0;i<Templ_P->GetNbinsX();i++)
		if(Templ_P->GetBinContent(i+1)>0&&Templ_P->GetBinContent(i+1)<mincontent){
			mincontent = Templ_P->GetBinContent(i+1);
	}
	Templ_P->Scale(1/mincontent);
	float maximum= Templ_P->GetBinCenter(Templ_P->GetMaximumBin());
	if(Templ_P->Integral(Templ_P->FindBin(3/maximum),Templ_P->FindBin(7/maximum))>0)
	{
		float content = Templ_P->Integral(Templ_P->FindBin(3/maximum),Templ_P->FindBin(7/maximum))/(Templ_P->FindBin(7/maximum)-Templ_P->FindBin(3/maximum));
		for(int i=0;i<Templ_P->GetNbinsX();i++) {
			if(Templ_P->GetBinContent(i+1)==0||Templ_P->GetBinError(i+1)==0) Templ_P->SetBinError(i+1,pow(content,0.5));
			Templ_P->SetBinError(i+1,(2/maximum)*Templ_P->GetBinError(i+1));
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

//void AdjoustTail(float factor){

//}

void TemplateFIT::FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars)){
	int kbin;
	if(!ApplyCuts("IsData",vars)) return;
	//cout<<vars->Event<<endl;
	//vars->PrintCurrentState();
	if(ApplyCuts((cut+"").c_str(),vars)/*&&kbin>=0*/){
		if(bins.IsUsingBetaEdges()){
			kbin = 	bins.GetBin(discr_var(vars));

			if(!(kbin<0)){
				//if(vars->BetaRICH_new>0) cout<<vars->Event<<" "<<vars->BetaRICH_new<<endl;
				if(fits[kbin][0][5]->Data){
					if(ApplyCuts(cut,vars)) 	fits[kbin][0][5]->DataAmbient->Fill(var(vars),vars->PrescaleFactor);		
					if(ApplyCuts((cut+"&"+cutoff).c_str(),vars)) fits[kbin][0][5]->Data->Fill(var(vars),vars->PrescaleFactor);		
					if(ApplyCuts(cutprimary,vars)) 		     fits[kbin][0][5]->DataPrim->Fill(var(vars),vars->PrescaleFactor);
				}
			}


		}
		else{
			kbin = 	bins.GetBin(discr_var(vars));
			if(!(kbin<0)){
				if(fits[kbin][0][5]->Data){
					if(ApplyCuts(cut,vars)) 	fits[kbin][0][5]->DataAmbient->Fill(var(vars),vars->PrescaleFactor);		
					if(ApplyCuts((cut+"&"+cutoff).c_str(),vars)) fits[kbin][0][5]->Data->Fill(var(vars),vars->PrescaleFactor);		
					if(ApplyCuts(cutprimary,vars)) 		     fits[kbin][0][5]->DataPrim->Fill(var(vars),vars->PrescaleFactor);
				}
			}



		}
		if(vars->R>50) {
			if(!isrich)	BetaResolutionDT->Fill(1/vars->Beta);
			else BetaResolutionDT->Fill(1/vars->BetaRICH_new);
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
	return ShiftStart;
}

float Systpar::GetShiftEnd(){
	return ShiftEnd;
}

float Systpar::GetSigmaStart(){
	return SigmaStart;
}

float Systpar::GetSigmaEnd(){
	return SigmaEnd;
}


float TemplateFIT::SmearBetaRICH_v2(int kbin, float Beta, float Beta_gen, float HighRsigmaMC, float HighRsigmaDT, float relativeshift, float stepsigma, float stepshift){
	float sigmaMC=HighRsigmaMC;
	float sigmaDT=HighRsigmaDT;


	if(kbin<=0) kbin =0;

	float delta = (1/Beta - 1/Beta_gen);


	float shift; 
	float sigma; 


	float shiftstart = -relativeshift - 0.0004*sigmaDT/0.00114; //-relativeshift - 5*fabs(relativeshift);
	float sigmastart = sigmaDT - fabs(sigmaDT-0.4*sigmaDT);


	if(stepsigma==-1 && stepshift==-1){
		shift= -relativeshift;
		sigma=sigmaDT;

	} 

	else{
		shift = shiftstart+2*(fabs(relativeshift-shiftstart)/systpar.steps)*stepshift; 
		sigma = sigmastart+2*(fabs(sigmaDT-sigmastart)/systpar.steps)*stepsigma;
		shift= -relativeshift;

	}

	return 1/(1/Beta_gen +(sigma/sigmaMC)*(delta + shift) );
}




void TemplateFIT::BuildRegularizedTemplates(int bin,float centroid_x,float centroid_y, TH2F* Chi){
/*	if(regularizationfile>0){
		Reg_sigma = (TF1*)regularizationfile->Get(("Parameters/"+basename+"/Sigma/Model_"+to_string(bin)).c_str());	
		Reg_shift = (TF1*)regularizationfile->Get(("Parameters/"+basename+"/Shift/Model_"+to_string(bin)).c_str());	
		cout<<"********* Regularization functions: "<<("Parameters/"+basename+"/Sigma/Model_"+to_string(bin)).c_str()<<" "<<Reg_sigma<<" "<<Reg_shift<<"**********"<<endl;

		if(Reg_sigma>0){
			cout<< "****************************+ REGULARIZATION PARAMETERS *****************************" <<endl;
			cout<<"Sigma: "<<Reg_sigma->Eval(time)<<endl;
			cout<<"Shift: "<<Reg_shift->Eval(time)<<endl;
			cout<< "*************************************************************************************" <<endl;

			int sigma1= 	Reg_sigma->Eval(time);
			int sigma2= 	Reg_sigma->Eval(time)+1;
			int shift1= 	Reg_shift->Eval(time);
			int shift2= 	Reg_shift->Eval(time)+1;
		}
	}
*/
	TH1F * TemplP = (TH1F*) fits[bin][0][5]->Templ_P->Clone("fitbest_P");
	TH1F * TemplD = (TH1F*) fits[bin][0][5]->Templ_D->Clone("fitbest_D"); 
	TemplP->Reset();
	float w=0;
	float wl=0;
	float x,y;
	float chimin=0;
	for(int sigma=0;sigma<systpar.steps;sigma++)
		for(int shift=0;shift<systpar.steps;shift++)
		if(fits[bin][sigma][shift]->DCounts<15*fits_best[bin]->DCounts && fits[bin][sigma][shift]->DCounts>0.006*fits_best[bin]->DCounts)
		{
			x=sigma+1;y=shift+1;
			wl=EvalFitProbability(Chi->GetBinContent(sigma+1,shift+1),x,y,centroid_x,centroid_y,chimin);
			w+=EvalFitProbability(Chi->GetBinContent(sigma+1,shift+1),x,y,centroid_x,centroid_y,chimin);
			TemplP->Add(fits[bin][sigma][shift]->Templ_P,wl);
		}
	TemplP->Scale(1/w);
	
	for(int i=0;i<fits[bin][0][5]->Templ_P->GetNbinsX();i++) {
		float se=1.3;
		if(fits[bin][0][5]->Templ_P->GetBinCenter(i+1)<1) se=2;
		TemplP->SetBinError(i+1,se*fits[bin][0][5]->Templ_P->GetBinError(i+1));
	}

	TemplD->Reset();
	w=0;
	for(int sigma=0;sigma<systpar.steps;sigma++)
		for(int shift=0;shift<systpar.steps;shift++)
		if(fits[bin][sigma][shift]->DCounts<15*fits_best[bin]->DCounts && fits[bin][sigma][shift]->DCounts>0.006*fits_best[bin]->DCounts)
		{
			x=sigma+1;y=shift+1;
			wl=EvalFitProbability(Chi->GetBinContent(sigma+1,shift+1),x,y,centroid_x,centroid_y,chimin);
			w+=EvalFitProbability(Chi->GetBinContent(sigma+1,shift+1),x,y,centroid_x,centroid_y,chimin);
			TemplD->Add(fits[bin][sigma][shift]->Templ_D,wl);
		}
	TemplD->Scale(1/w);
	for(int i=0;i<fits[bin][0][5]->Templ_D->GetNbinsX();i++) {
		int k = fits_best[bin]->Templ_P->FindBin(fits_best[bin]->Templ_D->GetBinCenter(i+1)/2.);
		//TemplD->SetBinError(i+1,(fits_best[bin]->Templ_P->GetBinError(k)/fits_best[bin]->Templ_P->GetBinContent(k))*fits_best[bin]->Templ_D->GetBinContent(i+1));
	}

	fits_best[bin]->Templ_P = (TH1F*) TemplP->Clone();
	fits_best[bin]->Templ_D = (TH1F*) TemplD->Clone();

	fits_best[bin] -> DCounts = fits_best[bin] ->  Templ_D -> Integral();
	fits_best[bin] -> PCounts = fits_best[bin] ->  Templ_P -> Integral();
	fits_best[bin] -> TCounts = fits_best[bin] ->  Templ_He -> Integral();
	
	return;
}



float TemplateFIT::SmearBetaTOF_v2(float massagen, float R, float Beta, float Beta_gen, float HighRsigmaMC, float HighRsigmaDT, float relativeshift, float stepsigma, float stepshift){
	float sigmaMC=HighRsigmaMC;
	float sigmaDT=HighRsigmaDT;

	
	if(massagen<1){
		slowdownmodel->SetParameter(0, 1.47911e+00);
		slowdownmodel->SetParameter(1,-3.67686e+00);
		slowdownmodel->SetParameter(2, 5.59880e+00);
		slowdownmodel->SetParameter(3,-3.05480e+00);
		slowdownmodel->SetParameter(4, 6.55016e-01);
	}	

	else{
		slowdownmodel->SetParameter(0, 1.00972e+00);
		slowdownmodel->SetParameter(1,-2.04976e+00);
		slowdownmodel->SetParameter(2, 3.45564e+00);
		slowdownmodel->SetParameter(3,-1.76667e+00);
		slowdownmodel->SetParameter(4, 3.52331e-01);
	}	


	Sigmamodel->SetParameter(0,1.23);
	Sigmamodel->SetParameter(1,-0.115);
	Sigmamodel->SetParameter(2,-0.084);



	float shift; 
	float sigma; 


	float shiftstart = -relativeshift - 4*fabs(relativeshift);
	float sigmastart = sigmaDT - fabs(sigmaDT-0.8*sigmaDT);


	if(stepsigma==-1 && stepshift==-1){
		shift= -relativeshift;
		sigma=sigmaDT;

	} 

	else{
		sigma = sigmastart+2*(fabs(sigmaDT-sigmastart)/systpar.steps)*stepsigma;//sigmastart + 2*fabs(sigmaDT-sigmaMC)/(systpar.steps)*(stepsigma+0.5);
		shift= -relativeshift;

	}

	richcore->SetRange(0.5,1.2);

	sigma = max((double)-0.75*(1/Beta_gen)+1.75,0.)*sigma;
	//sigma = Sigmamodel->Eval(1/Beta_gen) * sigma;	

	//core
	richcore->SetNpx(1e3);
	
	richcore->SetParameter(0,1.987e4);
	richcore->SetParameter(1,1+shift);
	richcore->SetParameter(2,sigma);
	richcore->SetParameter(3,0);
	richcore->SetParameter(4,0);
	richcore->SetParameter(5,1.064e4);
	richcore->SetParameter(6,1.284);
	richcore->SetParameter(7,4.433);
	richcore->SetParameter(8,3.543);
	richcore->SetParameter(9,0.000625);
	richcore->SetParameter(10,0.025);

	float betasmear = 1/(slowdownmodel->Eval(1/Beta_gen) + (richcore->GetRandom()-1)    );
//	if(R/betasmear * pow((1-pow(betasmear,2)),0.5)<0.7 && betasmear<0.85) cout<<R/betasmear * pow((1-pow(betasmear,2)),0.5)<<" "<<R/Beta * pow((1-pow(Beta,2)),0.5)<<" "<<betasmear<<" "<<Beta<<" "<<Beta_gen<<endl; 

	//if(R/betasmear * pow((1-pow(betasmear,2)),0.5)<0.7 && betasmear<0.85) cout<<fabs(betasmear-Beta)/Beta<<endl;
	if(fabs(betasmear-Beta)/Beta>0.10 && betasmear>Beta) return Beta;
	return betasmear;
}


float TemplateFIT::SmearBetaNaF_v2(int kbin, float Beta, float Beta_gen, float HighRsigmaMC, float HighRsigmaDT, float relativeshift, float stepsigma, float stepshift){
	float sigmaMC=HighRsigmaMC;
	float sigmaDT=HighRsigmaDT;


	Sigmamodel->SetParameter(0,6.18);
	Sigmamodel->SetParameter(1,-10.07);
	Sigmamodel->SetParameter(2,4.90);


	if(kbin<=0) kbin =0;

	float norm = Rand->Rndm();

	float delta = (1/Beta - 1/Beta_gen);

	/*	float sigma = (sigmaDT)*(1-systpar.sigma/100.); 
		sigma+=(float)((2*fabs(sigmaDT-sigma)/(float)systpar.steps)*stepsigma);

		float shift = -relativeshift - systpar.shift*1e-5 + (2*systpar.shift*1e-5/(float)systpar.steps)*stepshift;		
		*/

	float shift; 
	float sigma; 


	float shiftstart = -relativeshift - 0.0004*sigmaDT/0.00114; //-relativeshift - 5*fabs(relativeshift);
	float sigmastart = sigmaDT - fabs(sigmaDT-0.8*sigmaDT);


	if(stepsigma==-1 && stepshift==-1){
		shift= -relativeshift;
		sigma=sigmaDT;

	} 

	else{
		//	shift = shiftstart+2*(fabs(relativeshift-shiftstart)/systpar.steps)*stepshift;  //shiftstart + 10*fabs(relativeshift)/(systpar.steps)*stepshift; 
		sigma = sigmastart+2*(fabs(sigmaDT-sigmastart)/systpar.steps)*stepsigma;//sigmastart + 2*fabs(sigmaDT-sigmaMC)/(systpar.steps)*(stepsigma+0.5);
		shift= -relativeshift;

	}

	richcore->SetRange(0.97,1.2);
	
	sigma = Sigmamodel->Eval(1/Beta_gen)*sigma;//max((double)-0.6*(1/Beta_gen)+1.6,0.)*sigma;
	
	//core
	richcore->SetNpx(1e3);
	
	richcore->SetParameter(0,1798);
	richcore->SetParameter(1,1+shift);
	richcore->SetParameter(2,sigma);
	richcore->SetParameter(3,2.276e4);
	richcore->SetParameter(4,-31.3);
	richcore->SetParameter(5,1160);
	richcore->SetParameter(6,1.491);
	richcore->SetParameter(7,67.08);
	richcore->SetParameter(8,2.541);
	richcore->SetParameter(9,0.0002029);
	richcore->SetParameter(10,0.001101);

	return 1/(1/Beta_gen + (richcore->GetRandom()-1)    );


	//return 1/(1/Beta_gen +(sigma/sigmaMC)*(delta + shift) );



}


float TemplateFIT::SmearBetaAgl_v2(int kbin, float Beta, float Beta_gen, float HighRsigmaMC, float HighRsigmaDT, float relativeshift, float stepsigma, float stepshift){
	float sigmaMC=HighRsigmaMC;
	float sigmaDT=HighRsigmaDT;

	Sigmamodel->SetParameter(0,-75.42);
	Sigmamodel->SetParameter(1,153.42);
	Sigmamodel->SetParameter(2,-77);



	float delta = (1/Beta - 1/Beta_gen);

	/*	float sigma = (sigmaDT)*(1-systpar.sigma/100.); 
		sigma+=(float)((2*fabs(sigmaDT-sigma)/(float)systpar.steps)*stepsigma);

		float shift = -relativeshift - systpar.shift*1e-5 + (2*systpar.shift*1e-5/(float)systpar.steps)*stepshift;		
		*/

	float shift; 
	float sigma; 


	float shiftstart = -relativeshift - 0.0004*sigmaDT/0.00114; //-relativeshift - 5*fabs(relativeshift);
	float sigmastart = sigmaDT - fabs(sigmaDT-0.8*sigmaDT);


	if(stepsigma==-1 && stepshift==-1){
		shift= -relativeshift;
		sigma=sigmaDT;

	} 

	else{
		//	shift = shiftstart+2*(fabs(relativeshift-shiftstart)/systpar.steps)*stepshift;  //shiftstart + 10*fabs(relativeshift)/(systpar.steps)*stepshift; 
		sigma = sigmastart+2*(fabs(sigmaDT-sigmastart)/systpar.steps)*stepsigma;//sigmastart + 2*fabs(sigmaDT-sigmaMC)/(systpar.steps)*(stepsigma+0.5);
		shift= -relativeshift;

	}
	sigma = Sigmamodel->Eval(1/Beta_gen)*sigma;
	
	richcore->SetRange(0.97,1.06);
	//core
	richcore->SetParameter(0,38852.35);
	richcore->SetParameter(1,1+shift);
	richcore->SetParameter(2,sigma);
	richcore->SetParameter(3,251600.7);
	richcore->SetParameter(4,-168.157);
	richcore->SetParameter(5,17962.76);
	richcore->SetParameter(6,1.599928);
	richcore->SetParameter(7,676.7463);
	richcore->SetParameter(8,3.005362);
	richcore->SetParameter(9,-0.0001519769);
	richcore->SetParameter(10,-0.0001153038);


	return 1/(1/Beta_gen + (richcore->GetRandom()-1)    );

	//return 1/(1/Beta_gen +(sigma/sigmaMC)*(delta + shift) );


}



float TemplateFIT::SmearBeta(float Beta, float Beta_gen,float  HighRsigmaMC, float HighRsigmaDT, float relativeshift , TF1* MCmodel,  float stepsigma, float stepshift){


	float beta_inv = 1/Beta;
	float sigmasmear = pow(pow(HighRsigmaDT,2)-pow(HighRsigmaMC,2),0.5);

	float shift;
        float sigma;

	float shiftstart = -relativeshift - 4*fabs(relativeshift);
	float sigmastart = 0;


	if(stepsigma==-1 && stepshift==-1){
		shift = -relativeshift;
		sigma = sigmasmear;
	}

	else{
		shift = -relativeshift; //shiftstart + 10*fabs(relativeshift)/(systpar.steps)*stepshift; 
		sigma = sigmastart + 1.2*fabs(sigmasmear-sigmastart)/(systpar.steps)*stepsigma;
	}
	
	float norm = Rand->Rndm();
	
	float MCsigma = MCmodel->Eval(1/Beta_gen);


	beta_inv = beta_inv + shift + Rand->Gaus(0,sigma);

	return 1/(beta_inv);
		
}





float TemplateFIT::SmearBeta_v2(float Beta, float Beta_gen,float  HighRsigmaMC, float HighRsigmaDT, float relativeshift , TF1* MCmodel,  float stepsigma, float stepshift){


	float beta_inv = 1/Beta;
	float sigmasmear = pow(pow(HighRsigmaDT,2)-pow(HighRsigmaMC,2),0.5);

/*	float shiftstart=-systpar.shift;
	float sigmastart = sigmasmear - sigmasmear*systpar.sigma/100;
	float sigma = sigmastart +  ((2*systpar.sigma/(systpar.steps))*stepsigma)/100*sigmasmear;
	float shift = -relativeshift - (systpar.shift*1e-5 + (2*systpar.shift*1e-5/(float)systpar.steps)*stepshift);	
*/

	float shift;
        float sigma;

	float shiftstart = -relativeshift - 4*fabs(relativeshift);
	float sigmastart = 0;


	if(stepsigma==-1 && stepshift==-1){
		shift = -relativeshift;
		sigma = sigmasmear;
	}

	else{
		shift = -relativeshift; //shiftstart + 10*fabs(relativeshift)/(systpar.steps)*stepshift; 
		sigma = sigmastart + 1.2*fabs(sigmasmear-sigmastart)/(systpar.steps)*stepsigma;
	}
	
	float norm = Rand->Rndm();
	
	float MCsigma = MCmodel->Eval(1/Beta_gen);


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
	return 1/(beta_inv);

		
}



void TemplateFIT::FillEventByEventMC(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars)){
	if(vars->R>100) return;
	if((ApplyCuts(cutP,vars)||ApplyCuts(cutD,vars)||ApplyCuts(cutHe,vars))){
		int imin,imax,jmin,jmax;
		imin=0; jmin=0; imax=systpar.steps; jmax=systpar.steps;
		if(IsExtern) {imin=0; jmin=5; imax=1; jmax=6;}
		int kbin=0;

		for(int i=0;i<systpar.steps;i++)
			for(int j=0;j<systpar.steps;j++)
			{

				if(!((i==0&&j==5)||(j==5))) continue;
				float betasmear=0;
				float Rsmear = vars->R;

				if(!isrich){ 	
					if(systpar.mode==1) {
						float ss,si1,si2;
						if(useMCtuning) {	
							if(vars->Massa_gen<2)	{
								ss=-0.0024; si1=0.0329; si2=0.0428; 
								//	betasmear = SmearBeta_v2(vars->Beta,GetBetaGen_cpct(vars),si1,si2,ss,MCmodel,i,j);   
								betasmear = SmearBetaTOF_v2(vars->Massa_gen,vars->R,vars->Beta,GetBetaGen_cpct(vars),si1,si2,ss,i,j); 
							}

							else  {
								ss=-0.0001; si1=1.66514e-02; si2=1.91982e-02; 
								betasmear = SmearBeta(vars->Beta,GetBetaGen_cpct(vars),si1,si2,ss,MCmodel,i,j);   

							}
						}
					}

				}
				else {
					if(ApplyCuts("IsFromNaF",vars))	 {
						if(systpar.mode==3) { 
							float ss,si1,si2;
							if(useMCtuning) {	
								if(vars->Massa_gen<2)	{
									ss=-0.00014; si1=0.00334; si2=0.002786; 
									betasmear = SmearBetaNaF_v2(kbin,vars->BetaRICH_new,GetBetaGen_cpct(vars),si1,si2,ss,i,j); 
								}
								else  {
									ss=0.00005; si1=0.00213; si2=0.00233; 
									betasmear = SmearBetaRICH_v2(kbin,vars->BetaRICH_new,GetBetaGen_cpct(vars),si1,si2,ss,i,j); 
								}
							}

							Rsmear = vars->R;
						}
					}
					if(ApplyCuts("IsFromAgl",vars))  {
						if(systpar.mode==3) {
							float ss,si1,si2;
							if(useMCtuning) {	
								if(vars->Massa_gen<2)	{
									ss=0.00006; si1=0.0012369; si2=0.0009792274; 
									betasmear = SmearBetaAgl_v2(kbin,vars->BetaRICH_new,GetBetaGen_cpct(vars),si1,si2,ss,i,j); 
								}
								else  {
									ss=-0.0001; si1=0.000733; si2=0.000760; 
									betasmear = SmearBetaRICH_v2(kbin,vars->BetaRICH_new,GetBetaGen_cpct(vars),si1,si2,ss,i,j); 

								}
							}
							Rsmear = vars->R;
						}
					}
				}

				float mctotalweight =1;
				mctotalweight = vars->mcweight;

				if(useMCreweighting){
					mctotalweight = vars->mcweight;
					//* vars->GetCutoffCleaningWeight(GetRFromBeta(bins.getParticle().getMass(),betasmear),GetGenMomentum(vars),BetacutoffCut) 
					if(!isrich) mctotalweight*=vars->GetTimeDepWeight(vars->R);
				}
				if(!IsFitNoise){
					if(bins.IsUsingBetaEdges()){	 
						kbin = bins.GetBin(betasmear);
						float scaledR;
						float mass;
						float rigcut=0;
						scaledR=Rsmear*template1scalefactor1;
						mass = scaledR/betasmear * pow((1-pow(betasmear,2)),0.5);
				//			if(ApplyCuts(cutD,vars))	cout<<"************ Accepted: "<<mass<<" "<<vars->Beta<<endl;

						if(ApplyCuts(cutP,vars)&&kbin>=0&&scaledR>=rigcut) {fits[kbin][i][5]->Templ_P->Fill(mass,mctotalweight);}
						scaledR=Rsmear*template1scalefactor2;
						mass = scaledR/betasmear * pow((1-pow(betasmear,2)),0.5);
						if(ApplyCuts(cutD,vars)&&kbin>=0&&scaledR>=rigcut)  {fits[kbin][i][5]->Templ_D->Fill(mass,mctotalweight);}
						scaledR=Rsmear*template1scalefactor3;
						mass = scaledR/betasmear * pow((1-pow(betasmear,2)),0.5);
						if(ApplyCuts(cutHe,vars)&&kbin>=0&&scaledR>=rigcut) fits[kbin][i][5]->Templ_He->Fill(mass,mctotalweight);
					}
					else {
						kbin = bins.GetBin(discr_var(vars));

						float scaledR;
						float mass;
						float rigcut=0;
						scaledR=Rsmear*template1scalefactor1;
						mass = scaledR/betasmear * pow((1-pow(betasmear,2)),0.5);
						if(ApplyCuts(cutP,vars)&&kbin>=0)   fits[kbin][i][5]->Templ_P->Fill(mass,mctotalweight);
						scaledR=Rsmear*template1scalefactor2;
						mass = scaledR/betasmear * pow((1-pow(betasmear,2)),0.5);
						if(ApplyCuts((cutD).c_str(),vars)&&kbin>=0)   fits[kbin][i][5]->Templ_D->Fill(mass,mctotalweight);
						scaledR=Rsmear*template1scalefactor3;
						mass = scaledR/betasmear * pow((1-pow(betasmear,2)),0.5);
						if(ApplyCuts(cutHe,vars)&&kbin>=0) fits[kbin][i][5]->Templ_He->Fill(mass,mctotalweight);

					}	

				}
			}
	}
	if(ApplyCuts(cutP,vars))
		if(vars->R>50) {
			if(!isrich)	{ 
				BetaResolutionMC->Fill(1/vars->Beta);
				float ss,si1,si2;
				if(vars->Massa_gen<2)	{ss=-0.0024; si1=0.0329; si2=0.0428; 
					BetaResolutionMC_tuned->Fill(1/SmearBetaTOF_v2(0.981,vars->R,vars->Beta,GetBetaGen_cpct(vars),si1,si2,ss,-1,-1)  );
				}
				else  {ss=-0.0001; si1=1.66514e-02; si2=1.91982e-02; 
					BetaResolutionMC_tuned->Fill(1/ SmearBeta(vars->Beta,GetBetaGen_cpct(vars),si1,si2,ss,MCmodel,-1,-1)   );
				}
			}
			else {
				BetaResolutionMC->Fill(1/vars->BetaRICH_new);
				if(ApplyCuts("IsFromNaF",vars))  {
					float ss,si1,si2;
					if(vars->Massa_gen<2)	{ss=-0.00014; si1=0.00334; si2=0.002786; 
						BetaResolutionMC_tuned->Fill(1/SmearBetaNaF_v2(10,vars->BetaRICH_new,GetBetaGen_cpct(vars),si1,si2,ss,-1,-1));
					}
					else  {ss=0.00005; si1=0.00213; si2=0.00233; 
						BetaResolutionMC_tuned->Fill(1/SmearBetaRICH_v2(10,vars->BetaRICH_new,GetBetaGen_cpct(vars),si1,si2,ss,-1,-1));
					}
				}

				if(ApplyCuts("IsFromAgl",vars))  {
					float ss,si1,si2;
					if(vars->Massa_gen<2)	{ss=0.00006; si1=0.0012369; si2=0.0009792274; 
						BetaResolutionMC_tuned->Fill(1/SmearBetaAgl_v2(10,vars->BetaRICH_new,GetBetaGen_cpct(vars),si1,si2,ss,-1,-1));
					}
					else  {ss=-0.0001; si1=0.000733; si2=0.000760; 
						BetaResolutionMC_tuned->Fill(1/SmearBetaRICH_v2(10,vars->BetaRICH_new,GetBetaGen_cpct(vars),si1,si2,ss,-1,-1));
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

	for(int i=0;i<bins.size();i++){
		if(BetaResolutionDT) finalhistos.Add(BetaResolutionDT);
		if(BetaResolutionMC) finalhistos.Add(BetaResolutionMC);
		if(BetaResolutionMC_tuned) finalhistos.Add(BetaResolutionMC_tuned);
	
	finalhistos.writeObjsInFolder((basename + "/BetaRes").c_str(),false);
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

	// Best è il migliore regolarizzato
	// Best2 + il migliore locale
	
	cout<<"Saving best"<<endl;	
	//save best
	for(int bin=0;bin<bins.size();bin++){
		TH1F * OriginalP=(TH1F*)fits[bin][0][5]->Templ_P->Clone();
		TH1F* BestP;
		TH1F* BestP2;
		if(!IsLocalFit&&!IsLocalConstrainedFit) {
				BestP=(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_P->Clone();
				BestP->SetName("Best #chi^{2}: "/* + to_string(BestChiSquare->i) + "_" + to_string(BestChiSquare->j)).c_str()*/);
			}	
		else {
			BestChi * bestlocal = new BestChi();
			if(IsLocalConstrainedFit)  bestlocal->FindMinimum(TFitChisquare[bin],BestChiSquare->i);
			else bestlocal->FindMinimum(TFitChisquare[bin]);
			BestP2=(TH1F*)fits_best[bin]->Templ_P->Clone();	
			BestP=(TH1F*)fits[bin][bestlocal->i][bestlocal->j]->Templ_P->Clone();	
			BestP->SetName(("BestP: " + to_string(bestlocal->i) + "_" + to_string(bestlocal->j)).c_str());
			BestP2->SetName("BestP: ");

		}
		TH1F * BestPPrim =(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_PPrim->Clone();
		OriginalP->SetName("Original Proton MC ");
		BestPPrim->SetName("Best #chi^{2} Overcutoff Proton MC ");
		finalhistos.Add(OriginalP);
		finalhistos.Add(BestP);
		finalhistos.Add(BestP2);
		finalhistos.Add(BestPPrim);
		finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesP/Bin"+to_string(bin)).c_str());
		TH1F * OriginalD=(TH1F*)fits[bin][0][5]->Templ_D->Clone();
		TH1F * BestD;
		TH1F * BestD2;
		if(!IsLocalFit&!IsLocalConstrainedFit){
				 BestD=(TH1F*)fits_best[bin]->Templ_D->Clone();
				 BestD2=(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_D->Clone();
			//	BestD->SetName("Best #chi^{2}: "/* + to_string(BestChiSquare->i) + "_" + to_string(BestChiSquare->j)).c_str()*/);
				BestD2->SetName("BestD: "/* + to_string(BestChiSquare->i) + "_" + to_string(BestChiSquare->j)).c_str()*/);
			}
		else {
			BestChi * bestlocal = new BestChi();
			if(IsLocalConstrainedFit)  bestlocal->FindMinimum(TFitChisquare[bin],BestChiSquare->i);
			else bestlocal->FindMinimum(TFitChisquare[bin]);
			BestD=(TH1F*)fits_best[bin]->Templ_D->Clone();	
			BestD2=(TH1F*)fits[bin][bestlocal->i][bestlocal->j]->Templ_D->Clone();	
			BestD->SetName(("BestD " + to_string(bestlocal->i) + "_" + to_string(bestlocal->j)).c_str());
			BestD2->SetName("BestD"/* + to_string(BestChiSquare->i) + "_" + to_string(BestChiSquare->j)).c_str()*/);
		}
		TH1F * BestDPrim =(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_DPrim->Clone();
		OriginalD->SetName("Original Deuton MC ");
		BestDPrim->SetName("Best #chi^{2} Overcutoff Deuton MC ");
		finalhistos.Add(OriginalD);
		finalhistos.Add(BestD);
		finalhistos.Add(BestD2);
		finalhistos.Add(BestDPrim);
		finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesD/Bin"+to_string(bin)).c_str());
		
		TH1F * BestHe=0x0;
		TH1F * OriginalHe=(TH1F*)fits[bin][0][5]->Templ_He->Clone();
		BestHe=(TH1F*)fits_best[bin]->Templ_He->Clone();
		TH1F * BestHePrim =(TH1F*)fits[bin][BestChiSquare->i][BestChiSquare->j]->Templ_HePrim->Clone();
		OriginalHe->SetName("Original Tritium MC ");
		BestHe->SetName("Best Tritium MC ");
		BestHePrim->SetName("Best #chi^{2} Overcutoff Triton MC ");
		finalhistos.Add(OriginalHe);
		finalhistos.Add(BestHe);
		finalhistos.Add(BestHePrim);
		finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesHe/Bin"+to_string(bin)).c_str());

		if(IsFitNoise){
			TH1F * BestNoise=0x0;
			BestNoise=(TH1F*)fits_best[bin]->Templ_Noise->Clone();
			BestNoise->SetName("Best Noise ");
			finalhistos.Add(BestNoise);
			finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesNoise/Bin"+to_string(bin)).c_str());
		}

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
	
	//	finalhistos.Add(DrawTemplateFit(bin,fits[bin][0][5]->Data,fits[bin][0][5]->Templ_P,fits[bin][0][5]->Templ_D,BestHe,Residuals));
		finalhistos.Add(DrawTemplateFit(bin,fits[bin][0][5]->Data,BestP,BestD,BestHe,Residuals));
		finalhistos.writeObjsInFolder((basename+"/Fit Results/DrawFits").c_str());

	}

	if(adjousttail) {
		TH1F* dummy=new TH1F();
		for(int i=0;i<tailFT.size();i++)
			finalhistos.Add(DrawTemplateFit(i,tailFT[i]->Data,tailFT[i]->Templ_P,tailFT[i]->Templ_D,tailFT[i]->Templ_He,dummy));
			for(int i=0;i<chitail.size();i++) 
			finalhistos.Add(chitail[i]);
			finalhistos.writeObjsInFolder((basename+"/Fit Results/TailFT").c_str());

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
	cout<<"Saving all"<<endl;	

	//save all
/*	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(fits[bin][i][j]->Templ_P) {
					fits[bin][i][j]->Templ_P->SetName(("P_"+to_string(bin)+"_"+to_string(i)+"_"+to_string(j)).c_str());
					finalhistos.Add(fits[bin][i][j]->Templ_P);
				}
			}
	finalhistos.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesP/Bin"+to_string(bin)).c_str());
	}

	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(fits[bin][i][j]->Templ_D) {
					fits[bin][i][j]->Templ_D->SetName(("D_"+to_string(bin)+"_"+to_string(i)+"_"+to_string(j)).c_str());
					finalhistos.Add(fits[bin][i][j]->Templ_D);
					}
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
				if(!(i==0&&i==3))if(fits[bin][i][j]->Residual) finalhistos.Add(fits[bin][i][j]->Residual);
				}
	finalhistos.writeObjsInFolder((basename+"/Fit Results/Residual/Bin"+to_string(bin)).c_str());
	}
*/	
	/////////////////////////////////////////////////////////
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
	
	if(BetaResolutionDT) finalhistos.Add(BetaResolutionDT);
	if(BetaResolutionMC) finalhistos.Add(BetaResolutionMC);
	if(BetaResolutionMC_tuned) finalhistos.Add(BetaResolutionMC_tuned);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/BetaResolutions").c_str());


	for(int i=0;i<bins.size();i++)
		if(TFitChisquare[i]) { 
			TFitChisquare[i]->Scale(fits[i][BestChiSquare->i][BestChiSquare->j]->ndf);	
			finalhistos.Add(TFitChisquare[i]);
		}
		finalhistos.writeObjsInFolder((basename+"/Fit Results/Spreads/ChiSquare").c_str());	

	for(int i=0;i<bins.size();i++)
		if(WeightedDCounts_unb[i]) finalhistos.Add(WeightedDCounts_unb[i]);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/WeightedDCounts_Unb").c_str());	
	for(int i=0;i<bins.size();i++)
		if(WeightedPCounts_unb[i]) finalhistos.Add(WeightedPCounts_unb[i]);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/WeightedPCounts_Unb").c_str());	
	
	for(int i=0;i<bins.size();i++)
		if(WeightedDCounts[i]) finalhistos.Add(WeightedDCounts[i]);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/WeightedDCounts").c_str());	
	for(int i=0;i<bins.size();i++)
		if(WeightedPCounts[i]) finalhistos.Add(WeightedPCounts[i]);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/WeightedPCounts").c_str());	
	for(int i=0;i<bins.size();i++)
		if(WeightedShift[i]) finalhistos.Add(WeightedShift[i]);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/WeightedShift").c_str());	
	for(int i=0;i<bins.size();i++)
		if(WeightedSigma[i]) finalhistos.Add(WeightedSigma[i]);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/WeightedSigma").c_str());	
		
		
	

	if(HeContModel) finalhistos.Add(HeContModel);
	if(HeContError) finalhistos.Add(HeContError);
	if(MeasuredHeContRatio) finalhistos.Add(MeasuredHeContRatio);
	if(MCHeContRatio) finalhistos.Add(MCHeContRatio);
	finalhistos.writeObjsInFolder((basename+"/Fit Results/HeliumContamination").c_str());	
	
	if(BetaResolutionDT && BetaResolutionMC && BetaResolutionMC_tuned) {
		TCanvas * BetaRes = new TCanvas("BetaRes HR");
		BetaResolutionDT->SetMarkerStyle(8);
		BetaResolutionDT->SetMarkerColor(1);
		BetaResolutionDT->SetLineColor(1);
		BetaResolutionDT->Scale(1/BetaResolutionDT->Integral());
		BetaResolutionMC->SetLineColor(2);
		BetaResolutionMC->Scale(1/BetaResolutionMC->Integral());

		BetaResolutionMC_tuned->SetLineColor(2);
		BetaResolutionMC_tuned->Scale(1/BetaResolutionMC_tuned->Integral());

		BetaResolutionDT->Draw("P");
		BetaResolutionMC->Draw("histsame");
		BetaResolutionMC_tuned->Draw("histsame");
		finalhistos.Add(BetaRes);			
		finalhistos.writeObjsInFolder((basename+"/Fit Results").c_str());	
	}


	if(Global_ChiSquare) {
		Global_ChiSquare->Scale(fits[0][BestChiSquare->i][BestChiSquare->j]->ndf);	
		finalhistos.Add(Global_ChiSquare);
		BestSigma   = ConvertBinnedHisto((TH1F*)BestSigma,"Best Sigma",bins,false); 
		BestShift   = ConvertBinnedHisto((TH1F*)BestShift,"Best Shift",bins,false); 
		finalhistos.Add(new TGraphErrors(BestSigma));
		finalhistos.Add(new TGraphErrors(BestShift));
		finalhistos.Add(BestSigma);
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
	ProtonCounts_unb   = ConvertBinnedHisto(ProtonCounts_unb,"Proton Counts_unb",bins,true); 
	DeuteronCounts_unb = ConvertBinnedHisto(DeuteronCounts_unb,"Deuteron Counts_unb",bins,true); 
	ProtonTail = ConvertBinnedHisto(ProtonTail,"Proton Tail",bins,true); 
	TritiumCounts = ConvertBinnedHisto(TritiumCounts,"Tritium Counts",bins,true); 

	TH1F * ContaminationD = (TH1F*) TritiumCounts->Clone("ContaminationD");
	ContaminationD->Divide(DeuteronCounts);
	ContaminationD->Scale(1.5);

	TH1F * RatioDP = (TH1F*) DeuteronCounts->Clone("RatioPD");
	RatioDP->Divide(ProtonCounts);

	TH1F * RatioDP_unb = (TH1F*) DeuteronCounts_unb->Clone("RatioPD_unb");
	RatioDP_unb->Divide(ProtonCounts_unb);



	finalhistos.Add(StatErrorP);
	finalhistos.Add(StatErrorD);
	finalhistos.Add(StatErrorT);
	finalhistos.Add(SystError);
	finalhistos.Add(SystErrorP);
	finalhistos.Add(ProtonCounts);
	finalhistos.Add(DeuteronCounts);
	finalhistos.Add(SystError_unb);
	finalhistos.Add(SystErrorP_unb);
	finalhistos.Add(ProtonCounts_unb);
	finalhistos.Add(DeuteronCounts_unb);
	finalhistos.Add(TritiumCounts);
	finalhistos.Add(ProtonTail);
	finalhistos.Add(RatioDP);
	finalhistos.Add(RatioDP_unb);
	finalhistos.Add(ContaminationD);
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
	float tailmult=1;
	//core
	for(int i=binmin; i<binmax/*data->FindBin(4)*/; i++){
			if(Result->GetBinError(i+1)<=0||Data->GetBinError(i+1)<=0) continue;
			err = pow(pow(Data->GetBinError(i+1),2) + pow(Result->GetBinError(i+1),2),0.5);
			chi += pow((Data->GetBinContent(i+1) - Result->GetBinContent(i+1)),2)/pow(err,2);
			ndf++;
	}
	/*
	float totalcontent_DT=0;
	float totalcontent_MC=0;
	float totalerror=0;
	for(int i=Data->FindBin(4); i<binmax; i++){
		totalerror+=pow(Data->GetBinError(i+1),2)+pow(Result->GetBinError(i+1),2);		
		if(Result->GetBinError(i+1)<=0||Data->GetBinError(i+1)<=0) continue;
		totalcontent_DT+=Data->GetBinContent(i+1);		
		totalcontent_MC+=Result->GetBinContent(i+1);		
	}
	chi+=pow(totalcontent_DT-totalcontent_MC,2)/totalerror;
	ndf++;*/
	if(prioritizetail){
		chi=0;
		ndf=0;
		for(int i=binmin; i<binmax; i++){
			if(Result->GetBinCenter(i+1)>4) tailmult=4;
			if(Result->GetBinError(i+1)<=0||Data->GetBinError(i+1)<=0) continue;
			err = pow(pow(Data->GetBinError(i+1),2) + pow(Result->GetBinError(i+1),2),0.5);
			chi += tailmult*pow((Data->GetBinContent(i+1) - Result->GetBinContent(i+1)),2)/pow(err,2);
			ndf++;
		}

	}	

	ndf = ndf-3;
//	chi /= ndf;
	cout<<"chi: "<<chi<<" ndf: "<<ndf<<endl;
	return chi;

}


void Do_TemplateFIT(TFit * Fit,float fitrangemin,float fitrangemax,float constrain_min[], float constrain_max[], bool isfitnoise, bool prioritizetail, bool IsFitPrim){
	cout<<	Fit -> Data->GetEntries()<<" "<<Fit -> Templ_P->GetEntries()<<" "<<Fit -> Templ_D->GetEntries()<<endl;

	Fit->Templ_Noise = (TH1F*)Fit -> Templ_P->Clone();
	for(int i=0;i<Fit->Templ_Noise->GetNbinsX();i++){

		if(Fit->Templ_Noise->GetBinCenter(i+1)>1) Fit->Templ_Noise->SetBinContent(i+1, 100*exp(-0.6*Fit->Templ_Noise->GetBinCenter(i+1)));
		else Fit->Templ_Noise->SetBinContent(i+1, 100*exp(-0.6*1));
		Fit->Templ_Noise->SetBinError(i+1,0.1*Fit->Templ_Noise->GetBinContent(i+1));
	}

	TObjArray *Tpl;
	Tpl = new TObjArray(3);

		

	if(!IsFitPrim) {

		if(Fit ->  Templ_P) 	Fit -> Templ_P-> Scale(Fit -> Templ_P->GetEntries()/Fit -> Templ_P->Integral());
		if(Fit ->  Templ_D) 	Fit -> Templ_D-> Scale(Fit -> Templ_D->GetEntries()/Fit -> Templ_D->Integral());
		if(Fit ->  Templ_He)	Fit -> Templ_He-> Scale(Fit -> Templ_He->GetEntries()/Fit -> Templ_He->Integral());
	
		if(Fit ->  Templ_P)  if(Fit ->  Templ_P ->Integral()>0) Tpl -> Add( Fit ->  Templ_P );
		if(Fit ->  Templ_D)  if(Fit ->  Templ_D ->Integral()>0) Tpl -> Add( Fit ->  Templ_D );
		if(Fit ->  Templ_He) if(Fit ->  Templ_He->Integral()>0) Tpl -> Add( Fit ->  Templ_He);
	}
	else {
		if(Fit ->  Templ_PPrim)  if(Fit ->  Templ_PPrim ->Integral()>0) Tpl -> Add( Fit ->  Templ_PPrim );
		if(Fit ->  Templ_DPrim)  if(Fit ->  Templ_DPrim ->Integral()>0) Tpl -> Add( Fit ->  Templ_DPrim );
		if(Fit ->  Templ_HePrim) if(Fit ->  Templ_HePrim->Integral()>0) Tpl -> Add( Fit ->  Templ_HePrim);
	}


	if(isfitnoise) if(Fit ->  Templ_Noise) Tpl -> Add( Fit ->  Templ_Noise);

	float min=fitrangemin;
	float max=fitrangemax;
	cout<<	Fit -> Data<<" "<<Fit -> Templ_P<<" "<<Fit -> Templ_D<<" "<<Fit -> Templ_He<<endl;

	bool fitcondition = (Fit -> Data->GetEntries()>100)&&(Fit -> Templ_P->GetEntries()>100);
	cout<<"fit"<<endl;
	cout<<	Fit -> Data->GetEntries()<<" "<<Fit -> Templ_P->GetEntries()<<" "<<Fit -> Templ_D->GetEntries()<<endl;

	if(fitcondition){	
		cout<<"Conditions for fit OK!"<<endl;	
		
		if(!IsFitPrim) Fit -> Tfit = new TFractionFitter(Fit -> Data, Tpl ,"q");
		else   	Fit -> Tfit = new TFractionFitter(Fit -> DataPrim, Tpl ,"q");
		cout<<"fitting..."<<endl;			
		cout<<"Constraints: "<<endl;
		cout<<constrain_min[0]<<" "<<constrain_max[0]<<endl;	
		cout<<constrain_min[1]<<" "<<constrain_max[1]<<endl;	
		cout<<constrain_min[2]<<" "<<constrain_max[2]<<endl;	
		Fit -> Tfit -> SetRangeX(Fit -> Data -> FindBin(min), Fit -> Data -> FindBin(max));
		
		Fit -> Tfit -> Constrain(1, constrain_min[0] ,constrain_max[0]);
                Fit -> Tfit -> Constrain(2, constrain_min[1] ,constrain_max[1]);
	 	Fit -> Tfit -> Constrain(3, constrain_min[2] ,constrain_max[2]);
	 	if(isfitnoise) Fit -> Tfit -> Constrain(4, 0.000005 ,0.1);



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
				Fit -> Tfit -> SetRangeX(Fit -> Data -> FindBin((1-0.03*fit_attempt)*min), Fit -> Data -> FindBin((1+0.03*fit_attempt)*max));

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
			
			Fit->Residual = (TH1F*)Fit->Data->Clone(Sum->GetName());
			Fit->Residual ->Add(Sum,-1);
			for(int l=0;l<Fit->Residual->GetNbinsX();l++){
				if(Fit->Data->GetBinError(l+1)>0){
					Fit->Residual->SetBinContent(l+1,pow(Fit->Residual->GetBinContent(l+1)/Fit->Data->GetBinError(l+1),2));		
					Fit->Residual->SetBinError(l+1,1);
				}
				else{
					Fit->Residual->SetBinContent(l+1,0);
					Fit->Residual->SetBinContent(l+1,1);
				}
			}
	
			if(( df1dw1*df1dw1*Cov00 + df1dw2*df1dw2*Cov11 + df1dw3*df1dw3*Cov22 +2*df1dw1*df1dw2*Cov01 +2*df1dw1*df1dw3*Cov02+ 2*df1dw2*df1dw3*Cov12)>0){
				Fit -> StatErrP =  pow(( df1dw1*df1dw1*Cov00 + df1dw2*df1dw2*Cov11 + df1dw3*df1dw3*Cov22 +2*df1dw1*df1dw2*Cov01 +2*df1dw1*df1dw3*Cov02+ 2*df1dw2*df1dw3*Cov12)/2,0.5)/w1;
				Fit -> StatErrD =  pow(( df2dw1*df2dw1*Cov00 + df2dw2*df2dw2*Cov11 + df2dw3*df2dw3*Cov22 +2*df2dw1*df2dw2*Cov01 +2*df2dw1*df2dw3*Cov02+ 2*df2dw2*df2dw3*Cov12)/2,0.5)/w2;
				Fit -> StatErrT =  pow(( df3dw1*df3dw1*Cov00 + df3dw2*df3dw2*Cov11 + df3dw3*df3dw3*Cov22 +2*df3dw1*df3dw2*Cov01 +2*df3dw1*df3dw3*Cov02+ 2*df3dw2*df3dw3*Cov12)/2,0.5)/w3;
			}
			else{

				Fit -> StatErrP = e1/w1;
				Fit -> StatErrD = e2/w2;
				Fit -> StatErrT = e3/w3;

			}
			cout<<"fract: "<<w1<<" "<<w2<<" "<<w3<<endl;
			cout<<endl;
			if(isfitnoise) cout<<"frac noise: "<<w4<<endl;
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



void SharpenDistrib(TH1F* Distrib, float factor){

	float sigma= Distrib->GetStdDev();
	float mean = Distrib->GetBinCenter(Distrib->GetMaximumBin());

	float sigma_r = sigma*pow(pow(1-factor,2)/(1-pow(1-factor,2)),0.5);

TF1 * G = new TF1("G","gaus",0,10);	
	G->SetParameter(0,1);
	G->SetParameter(1,mean);
	G->SetParameter(2,sigma_r);

	for(int i=0;i<Distrib->GetNbinsX();i++){
		float Mid = mean + 4*sigma;
		if(Distrib->GetBinCenter(i+1)<Mid)
			Distrib->SetBinContent(i+1,Distrib->GetBinContent(i+1)*G->Eval(Distrib->GetBinCenter(i+1)));
		else
			Distrib->SetBinContent(i+1,Distrib->GetBinContent(i+1)*G->Eval(Mid));
	}
	return;

}




void TemplateFIT::ExtractCounts(FileSaver finalhistos,int force_shift){

	SimpleExtractPrimaries();
	//magnitude fine tuning
	float mag=0;
	int sigmabest[4]={0,0,0,0};
	float chibest=9e7; 
	float magbest[4]={0,0,0,0};

	if(sharpen) {
		for(int bin=0;bin<bins.size();bin++)	
			 for(int sigma=0;sigma<systpar.steps;sigma++){
			SharpenDistrib(fits[bin][sigma][5]->Templ_P,sharpenfactor);
			//if(bin<(bins.size()-2))	SharpenDistrib(fits[bin][sigma][5]->Templ_P,sharpenfactor);
			SharpenDistrib(fits[bin][sigma][5]->Templ_D,sharpenfactor);

		}
	}		
	
	if(adjoustpeak){
	    for(int bin=0;bin<bins.size();bin++){	
		 for(int sigma=0;sigma<systpar.steps;sigma++){
		if(iso==0) fits[bin][sigma][5]->Templ_P =SimpleShiftHisto(fits[bin][sigma][5]->Templ_P,magpeak);
		if(iso==1) fits[bin][sigma][5]->Templ_D =SimpleShiftHisto(fits[bin][sigma][5]->Templ_D,magpeak);
				}
		}
	}

	if(adjousttail){
		mag=Mag;
		for(int bin=0;bin<4;bin++) {
			if((Midbin-2+bin)>=bins.size()) continue;
			for(int sigma=0;sigma<systpar.steps;sigma++) tailFT.push_back(fits[Midbin-2+bin][sigma][5]->Clone());
			chibest=9e7;
			int sigmaindex=0;
			for(int k=(systpar.steps+20)*bin;k<(systpar.steps+20)*(bin)+systpar.steps;k++){
		
				float fixminimum_factor=1;
				//if((k-(systpar.steps+20)*bin)==3) fixminimum_factor=0.001;
				//else fixminimum_factor=1;

				float eq = EqualizeTemplate(tailFT[k]->Templ_P,tailFT[k]->Data);
				tailFT[k]->Templ_P=SimpleShiftHisto(tailFT[k]->Templ_P,eq);			
				tailFT[k]->RegularizeTemplateError();
				Do_TemplateFIT(tailFT[k],fits[Midbin-2+bin][0][5]->fitrangemin,4.5,constrainmin,constrainmax,IsFitNoise,false);
				if(fixminimum_factor*tailFT[k]->ChiSquare<chibest){
					chibest=fixminimum_factor*tailFT[k]->ChiSquare;
					sigmabest[bin]=sigmaindex;
				}
				sigmaindex++;
			}
		
			chitail.push_back(new TH1F(("mod_chi bin "+to_string(bin)).c_str(),("mod_chi bin "+to_string(bin)).c_str(),20,0,20));
			chibest=9e7;
			for(int k=0;k<20;k++) tailFT.push_back(tailFT[(systpar.steps+20)*bin +sigmabest[bin]]->Clone());
			for(int k=0;k<20;k++) {
			
				int n=(20+systpar.steps)*bin + (k+systpar.steps);

				for(int i=0;i<tailFT[n]->Templ_P->GetNbinsX();i++)
					if(tailFT[n]->Templ_P->GetBinCenter(i+1)>1.4) {
						float mod=((0.4+0.3*k)*mag)/(1+exp(-3*(tailFT[n]->Templ_P->GetBinCenter(i+1)-Mid)));
						//	float mod=(((0.5+0.2*k)*mag/2/1.57)*atan(Fast*(tailFT[n]->Templ_P->GetBinCenter(i+1)-Mid))+(0.5+0.2*k)*mag/2);
						tailFT[n]->Templ_P->SetBinContent(i+1,max((double)(1+mod)*tailFT[n]->Templ_P->GetBinContent(i+1),0.));
						tailFT[n]->Templ_P->SetBinError(i+1,max((double)(1+mod)*tailFT[n]->Templ_P->GetBinError(i+1),0.));
					
						} 		
				
				cout<<"*********************** TAIL FINE TUNING: "<<(0.4+0.3*k)<<" *********************************"<<endl;
				//tailFT[n]->RegularizeTemplateError();
				Do_TemplateFIT(tailFT[n],fits[Midbin-2+bin][0][5]->fitrangemin,fits[Midbin-2+bin][0][5]->fitrangemax,constrainmin,constrainmax,IsFitNoise,true);
				chitail[bin]->SetBinContent(k+1,tailFT[n]->ChiSquare);
				if(tailFT[n]->ChiSquare<chibest){
					chibest=tailFT[n]->ChiSquare;
					magbest[bin]=(0.4+0.3*k)*mag;
				}

			}	

		}	
	}

	for(int bin=0;bin<bins.size();bin++){

		float eq = 1;

		TH1F * DD = (TH1F*) fits[bin][0][5]->Data->Clone("Data");
		std::vector<TH1F*> TP;
		std::vector<TH1F*> TD;
		std::vector<TH1F*> THe;

		for(int sigma=0;sigma<systpar.steps;sigma++){
			TP.push_back((TH1F*)fits[bin][sigma][5]->Templ_P  ->Clone());
			TD.push_back((TH1F*)fits[bin][sigma][5]->Templ_D  ->Clone());
			THe.push_back((TH1F*)fits[bin][sigma][5]->Templ_He  ->Clone());
		}

	for(int sigma=0;sigma<systpar.steps;sigma++){
			for(int shift=0;shift<systpar.steps;shift++){
					//build all template combinations	
					fits[bin][sigma][shift]->Data=(TH1F*) DD->Clone();
					
					fits[bin][sigma][shift]->Templ_P =SimpleShiftHisto((TH1F*)TP[sigma]  ->Clone((("P_" +to_string(bin)+"_"+to_string(sigma)+"_"+to_string(shift)).c_str())),-systpar.shift+2*(systpar.shift/systpar.steps)*shift );
					fits[bin][sigma][shift]->Templ_D =SimpleShiftHisto((TH1F*)TD[sigma]  ->Clone((("D_" +to_string(bin)+"_"+to_string(sigma)+"_"+to_string(shift)).c_str())),-systpar.shift+2*(systpar.shift/systpar.steps)*shift );
					fits[bin][sigma][shift]->Templ_He=SimpleShiftHisto((TH1F*)THe[sigma] ->Clone((("He_"+to_string(bin)+"_"+to_string(sigma)+"_"+to_string(shift)).c_str())),-systpar.shift+2*(systpar.shift/systpar.steps)*shift );

				systpar.ShiftStart= -0.02;
				systpar.ShiftEnd=0.02;
				if(systpar.mode==1) {
					systpar.SigmaStart=-0.2;
					systpar.SigmaEnd=0.2;
				}
				if(systpar.mode==3) {
					systpar.SigmaStart=-0.2;
					systpar.SigmaEnd=0.2;
				}

					
				bool prioritizetail=false;
				// MAIN Template Fit Manipulations 
				if(adjousttail){

					float binf=bin;
					float bins=fits.size()/2.;
					float enmod=0;           //exp(-0.5*pow((binf-Midbin)/1,2));
					if(bin==Midbin-2) enmod=magbest[0];
					else if(bin==Midbin-1) enmod=magbest[1];
					else if(bin==Midbin) enmod=magbest[2];
					else if(bin==Midbin+1) enmod=magbest[3];
					else enmod=0;  
					cout<<"*********************** TUNED TAIL: "<<enmod<<" *********************************"<<endl;
						
					for(int i=0;i<fits[bin][sigma][shift]->Templ_P->GetNbinsX();i++)	
					if(fits[bin][sigma][shift]->Templ_P->GetBinCenter(i+1)>1.4) {
					float mod=enmod/(1+exp(-3*(fits[bin][sigma][shift]->Templ_P->GetBinCenter(i+1)-Mid)));
				//	float mod=enmod*((mag/2/1.57)*atan(Fast*(fits[bin][sigma][shift]->Templ_P->GetBinCenter(i+1)-Mid))+mag/2);
					fits[bin][sigma][shift]->Templ_P->SetBinContent(i+1,max((double)(1+mod)*fits[bin][sigma][shift]->Templ_P->GetBinContent(i+1),0.));
					fits[bin][sigma][shift]->Templ_P->SetBinError(i+1,max((double)(1+mod)*fits[bin][sigma][shift]->Templ_P->GetBinError(i+1),0.));
					//	prioritizetail=true;	
					}
				}

				if(adjoustfixedtail){
					float enmod=exp(-0.5*pow((bin-Midbin)/2,2));
					for(int i=0;i<fits[bin][sigma][shift]->Templ_P->GetNbinsX();i++)	
					if(1<2.1) {
                                        float mod=enmod*Mag/(1+exp(-4*(fits[bin][sigma][shift]->Templ_P->GetBinCenter(i+1)-Mid)));
					fits[bin][sigma][shift]->Templ_P->SetBinContent(i+1,max((double)(1+mod)*fits[bin][sigma][shift]->Templ_P->GetBinContent(i+1),0.));
                                      //fits[bin][sigma][shift]->Templ_P->SetBinError(i+1,1+mod);//max((double)(1+mod)*fits[bin][sigma][shift]->Templ_P->GetBinError(i+1),0.));
					}	
				}

				if(!fitDisabled) {
					cout<<endl;
					cout<<"Bin: "<<bin<<": "<<sigma<<" "<<shift<<endl;
				//	fits[bin][sigma][shift]->RegularizeTemplateError();	
					if(lowstatDmode) {
						fits[bin][sigma][shift]->Data->Rebin(2);
						fits[bin][sigma][shift]->Templ_P->Rebin(2);
						fits[bin][sigma][shift]->Templ_D->Rebin(2);
						fits[bin][sigma][shift]->Templ_He->Rebin(2);
					}
					Do_TemplateFIT(fits[bin][sigma][shift],fits[bin][0][5]->fitrangemin,fits[bin][0][5]->fitrangemax,constrainmin,constrainmax,IsFitNoise,prioritizetail);
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
						float fixminimum_factor=1;
		
			dcountsspread->SetBinContent(sigma+1,shift+1,fits[bin][sigma][shift]->DCounts);

				if(fits[bin][sigma][shift]->ChiSquare>0){
					
					tfitchisquare->SetBinContent(sigma+1,shift+1,fits[bin][sigma][shift]->ChiSquare);	
				}
				else  tfitchisquare->SetBinContent(sigma+1,shift+1,5000);
				tfitchisquare->SetBinError(sigma+1,shift+1,0.2);	
			}
		}

		DCountsSpread.push_back(dcountsspread);
		TFitChisquare.push_back(tfitchisquare);

	}

	Global_ChiSquare = (TH2F*) TFitChisquare[TFitChisquare.size()-1]->Clone("Global ChiSquare");
	for(int i=0;i<TFitChisquare.size()-2;i++) Global_ChiSquare->Add(TFitChisquare[i]);
	BestChiSquare = new BestChi();
	BestChiSquare->FindMinimum(Global_ChiSquare);



	for(int bin=0;bin<bins.size();bin++){
		BestChi * bestlocal = new BestChi();
		if(IsLocalConstrainedFit)  bestlocal->FindMinimum(TFitChisquare[bin],BestChiSquare->i);
		else bestlocal->FindMinimum(TFitChisquare[bin]);

		BestChi * besterror = new BestChi();
		besterror->FindMinimum(TFitChisquare[bin]);
		besterror->FindCentroid(TFitChisquare[bin]);


		if(IsLocalFit||IsLocalConstrainedFit) fits_best.push_back(fits[bin][bestlocal->i][bestlocal->j]->Clone());
		else fits_best.push_back(fits[bin][BestChiSquare->i][BestChiSquare->j]->Clone());

		//float Dmax= fits[bin][besterror->centroid_x][besterror->centroid_y]->DCounts;
		float Dmax = fits_best[bin]->DCounts;
		TH1F * weighteddcounts  = new TH1F(("WeightedDCountsBin " +to_string(bin)).c_str(),("WeightedDCountsBin " +to_string(bin)+";Counts;").c_str(),1000,0.3*Dmax,3*Dmax/*fits_best[bin]->DCounts*/);
		TH1F * weighteddcounts_unb  = new TH1F(("WeightedDCountsUnbBin " +to_string(bin)).c_str(),("WeightedDCountsUnbBin " +to_string(bin)+";Counts;").c_str(),1000,0.3*Dmax,3*Dmax/*fits_best[bin]->DCounts*/);
		float Pmax = fits_best[bin]->PCounts;
		TH1F * weightedpcounts  = new TH1F(("WeightedPCountsBin " +to_string(bin)).c_str(),("WeightedPCountsBin " +to_string(bin)+";Counts;").c_str(),1000,0.03*Dmax,3*Pmax/*fits_best[bin]->DCounts*/);
		TH1F * weightedpcounts_unb  = new TH1F(("WeightedPCountsUnbBin " +to_string(bin)).c_str(),("WeightedPCountsUnbBin " +to_string(bin)+";Counts;").c_str(),1000,0.3*Pmax,3*Pmax/*fits_best[bin]->DCounts*/);
		TF1* regfunc=new TF1("","gaus",-20,20);
		TH1F * weightedsigma  = new TH1F(("WeightedSigmaBin " +to_string(bin)).c_str(),("WeightedSigmaBin " +to_string(bin)+";Counts;").c_str(),100,-0.2,0.2);
		TH1F * weightedshift  = new TH1F(("WeightedShiftBin " +to_string(bin)).c_str(),("WeightedShiftBin " +to_string(bin)+";Counts;").c_str(),100,-0.02,0.02);
	

		regfunc->SetParameter(0,1);
		regfunc->SetParameter(1,0);
		regfunc->SetParameter(2,1);



		float chimin=fits[bin][bestlocal->i][bestlocal->j]->ChiSquare;	


		for(int sigma=0;sigma<systpar.steps;sigma++)
			for(int shift=0;shift<systpar.steps;shift++)
				if(fits[bin][sigma][shift]->DCounts&&fits[bin][0][5]->DCounts){
					float adjousttailchi=1;
					if(adjousttail){
						if(bin==Midbin-2||bin==Midbin-1||bin==Midbin||bin==Midbin+1){ 
							if(sigma==sigmabest[(int)(bin - Midbin+2)])  adjousttailchi=0.1;
							else adjousttailchi=1;
						}
					}
				
					float fixminimum_factor=1;

					/////************ WEIGHTING *******************//
					float ww = EvalFitProbability(fixminimum_factor*adjousttailchi*fits[bin][sigma][shift]->ChiSquare,sigma+1,shift+1,besterror->centroid_x,besterror->centroid_y,chimin);
					float ww_unb = EvalFitProbability(fixminimum_factor*adjousttailchi*fits[bin][sigma][shift]->ChiSquare,sigma+1,shift+1,besterror->centroid_x,besterror->centroid_y,chimin);
					
					/////************ REGULARIZATION *******************//
				
			//		ww*= std::max((1-0.2*fabs(bestlocal->i-sigma)),0.);
			//		ww*= std::max((1-0.2*fabs(bestlocal->j-shift)),0.);
				
					//ww*= std::max((1-0.025*fabs(sigma)),0.1);
					//ww*= std::max((1-0.025*fabs(shift)),0.1);



					/////////////////////////////////////////////////
	
					float totcounts=fits[bin][sigma][shift]->DCounts;
					//					weighteddcounts -> Fill(totcounts,EvalFitProbability(fits[bin][sigma][shift]->ChiSquare,TFitChisquare[bin]->GetBinContent(besterror->i+1,besterror->j+1) ));
					//					weighteddcounts -> Fill(totcounts,EvalFitProbability(fits[bin][sigma][shift]->ChiSquare,TFitChisquare[bin]->GetMaximum()));
					weighteddcounts -> Fill(fits[bin][sigma][shift]->DCounts,ww);
					weightedpcounts -> Fill(fits[bin][sigma][shift]->PCounts,ww);
					weighteddcounts_unb -> Fill(fits[bin][sigma][shift]->DCounts,ww_unb);
					weightedpcounts_unb -> Fill(fits[bin][sigma][shift]->PCounts,ww_unb);


					weightedsigma->Fill( systpar.GetSigmaStart()+(sigma*fabs(systpar.GetSigmaEnd()-systpar.GetSigmaStart())/systpar.steps) ,ww_unb);	
					weightedshift->Fill( systpar.GetShiftStart()+(shift*fabs(systpar.GetShiftEnd()-systpar.GetShiftStart())/systpar.steps) ,ww_unb);	


				}
		int nbin=2*(2*fits_best[bin]->DCounts-0.4*fits_best[bin]->DCounts)/weighteddcounts->GetStdDev();
		//weighteddcounts->Rebin(100/nbin);
		WeightedDCounts.push_back(weighteddcounts);
		WeightedPCounts.push_back(weightedpcounts);
		WeightedDCounts_unb.push_back(weighteddcounts_unb);
		WeightedPCounts_unb.push_back(weightedpcounts_unb);
		WeightedSigma.push_back(weightedsigma);
		WeightedShift.push_back(weightedshift);
	
	}


/////************ REGULARIZATION 2*******************//

	if(isreg){
		cout<<"*********************** REGULARIZATION ***************************"<<endl;
		TF1* regfunc2=new TF1("","gaus",-20,20);
		regfunc2->SetParameter(0,1);
		regfunc2->SetParameter(1,0);
		regfunc2->SetParameter(2,0.05);


		TH1F*	PC    = new TH1F("PC","PC",bins.size(),0,bins.size()) ;
		TH1F*	DC    = new TH1F("DC","DC",bins.size(),0,bins.size()) ;

		for(int bin=0;bin<bins.size();bin++){
			PC->SetBinContent(bin+1,WeightedPCounts[bin]->GetMean());
			DC->SetBinContent(bin+1,WeightedDCounts[bin]->GetMean());
		}

		DC->Divide(PC);
		for(int bin=0;bin<bins.size();bin++){
			WeightedDCounts[bin]->Reset();
			WeightedPCounts[bin]->Reset();
			BestChi * bestlocal = new BestChi();
			if(IsLocalConstrainedFit)  bestlocal->FindMinimum(TFitChisquare[bin],BestChiSquare->i);
			else bestlocal->FindMinimum(TFitChisquare[bin]);
			float chimin=fits[bin][bestlocal->i][bestlocal->j]->ChiSquare;	


	
			for(int sigma=0;sigma<systpar.steps;sigma++)
				for(int shift=0;shift<systpar.steps;shift++){
					float fixminimum_factor=1;
					float adjousttailchi=1;

					float ww = EvalFitProbability(fixminimum_factor*adjousttailchi*fits[bin][sigma][shift]->ChiSquare,sigma+1,shift+1,0,0,chimin);
					float r1 = DC->GetBinContent(bin);
					float r2 = DC->GetBinContent(bin+2);
					
					float ideal;
					if(bin!=0&&bin!=bins.size()-1) ideal = fabs(r2+r1)/2.;
					else if(bin==0) ideal = fabs(r2+regmin*DC->GetBinContent(bin+1))/2.;
					else ideal = fabs(r1+regmax*DC->GetBinContent(bin+1))/2.;
					
					float real = fits[bin][sigma][shift]->DCounts/fits[bin][sigma][shift]->PCounts;
					ww*=regfunc2->Eval(fabs(ideal-real)/ideal);
					WeightedPCounts[bin] -> Fill(fits[bin][sigma][shift]->PCounts,ww);
					WeightedDCounts[bin] -> Fill(fits[bin][sigma][shift]->DCounts,ww);
				}
			//      PC->SetBinContent(bin+1,WeightedPCounts[bin]->GetMean());
			//      DC->SetBinContent(bin+1,WeightedDCounts[bin]->GetMean());
		}
	}
	/////////////////////////////////////////////////////////

EvalFinalErrors();
CalculateFinalPDCounts();	
EvalFinalParameters();



/// TIME REGULARIZED FIT
for(int bin=0;bin<bins.size();bin++){
	BestChi * besterror = new BestChi();
	besterror->FindMinimum(TFitChisquare[bin]);
	besterror->FindCentroid(TFitChisquare[bin]);

	//if(!(adjousttail||adjoustfixedtail)) BuildRegularizedTemplates(bin,besterror->centroid_x,besterror->centroid_y,TFitChisquare[bin]);
}	
	

return;
}


void TemplateFIT::EvalFinalParameters(){
	BestSigma = (TH1F*) BestChiSquares->Clone("BestSigma");
	BestShift = (TH1F*) BestChiSquares->Clone("BestShift");
	BestSigma->Clear();
	BestShift->Clear();
	for(int bin=0;bin<bins.size();bin++){

		if(!IsLocalFit&&!IsLocalConstrainedFit) BestChiSquares     ->SetBinContent(bin+1,TFitChisquare[bin]->GetBinContent(BestChiSquare->i+1,BestChiSquare->j+1));
		else{
			BestChi * bestlocal = new BestChi();
			if(IsLocalConstrainedFit)  bestlocal->FindMinimum(TFitChisquare[bin],BestChiSquare->i);
			else bestlocal->FindMinimum(TFitChisquare[bin]);
			bestlocal->FindCentroid(TFitChisquare[bin]);
			BestChiSquares     ->SetBinContent(bin+1,TFitChisquare[bin]->GetBinContent(bestlocal->i+1,bestlocal->j+1));
	
			BestSigma->SetBinContent(bin+1,WeightedSigma[bin]->GetMean());
                        BestSigma->SetBinError(bin+1,WeightedSigma[bin]->GetStdDev());
			BestShift->SetBinContent(bin+1,WeightedShift[bin]->GetMean());
                        BestShift->SetBinError(bin+1,WeightedShift[bin]->GetStdDev());
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
	
			//ADD systematic if template mod
			double systmodfactor=1;
			if(adjousttail) systmodfactor*=1.2;
			if(adjoustpeak) systmodfactor*=1.2;
			if(sharpenfactor) systmodfactor*=1.2;
			//find best error
			float x,y;
               		 float BestED = 9999999;
  			for(int sigma=0;sigma<systpar.steps;sigma++)
                        	for(int shift=0;shift<systpar.steps;shift++)
                                	if(fits[bin][sigma][shift]->StatErrD<BestED){
                                        	if(fits[bin][sigma][shift] ->ContribP>constrainmin[0]+0.001&&fits[bin][sigma][shift] ->ContribP<constrainmax[0]-0.001)
				           	if(fits[bin][sigma][shift] ->ContribD>constrainmin[1]+0.001&&fits[bin][sigma][shift] ->ContribD<constrainmax[1]-0.001){
						BestED=fits[bin][sigma][shift]->StatErrD;
                        			x=sigma;y=shift;   
						}
				}	

			 float BestEP = 9999999;
  			for(int sigma=0;sigma<systpar.steps;sigma++)
                        	for(int shift=0;shift<systpar.steps;shift++)
                                	if(fits[bin][sigma][shift]->StatErrP<BestEP){
                                        	if(fits[bin][sigma][shift] ->ContribP>constrainmin[0]+0.001&&fits[bin][sigma][shift] ->ContribP<constrainmax[0]-0.001)
				           	if(fits[bin][sigma][shift] ->ContribD>constrainmin[1]+0.001&&fits[bin][sigma][shift] ->ContribD<constrainmax[1]-0.001){
					       	BestEP=fits[bin][sigma][shift]->StatErrP;
						}
                           	}	

			cout<<"STATERR BIN: "<<bin<<" "<<BestED<<" "<<x<<" "<<y<<endl;
			StatErrorP -> SetBinContent(bin+1,BestEP);
			StatErrorD -> SetBinContent(bin+1,BestED);//  max((double)BestED,0.015));
			StatErrorT -> SetBinContent(bin+1,1.2*fits[bin][bestlocal->i][bestlocal->j]->StatErrT);
		
			if(IsLocalConstrainedFit) SystError -> SetBinContent(bin+1,0.15*WeightedDCounts[bin]->GetStdDev()/fits[bin][bestlocal->i][bestlocal->j]->DCounts);
			else SystError -> SetBinContent(bin+1,systmodfactor*0.15*WeightedDCounts[bin]->GetStdDev()/fits[bin][bestlocal->i][bestlocal->j]->DCounts);

			if(IsLocalConstrainedFit) SystError_unb -> SetBinContent(bin+1,0.15*WeightedDCounts_unb[bin]->GetStdDev()/fits[bin][bestlocal->i][bestlocal->j]->DCounts);
			else SystError_unb -> SetBinContent(bin+1,systmodfactor*0.15*WeightedDCounts_unb[bin]->GetStdDev()/fits[bin][bestlocal->i][bestlocal->j]->DCounts);

			if(lowstatDmode) SystError -> SetBinContent(bin+1,SystError -> GetBinContent(bin+1)+0.01);	
			if(lowstatDmode) SystError_unb -> SetBinContent(bin+1,SystError_unb -> GetBinContent(bin+1)+0.01);	

			SystErrorP->SetBinContent(bin+1, SystError -> GetBinContent(bin+1)*fits[bin][bestlocal->i][bestlocal->j]->DCounts/fits[bin][bestlocal->i][bestlocal->j]->PCounts);
			SystErrorP_unb->SetBinContent(bin+1, SystError_unb -> GetBinContent(bin+1)*fits[bin][bestlocal->i][bestlocal->j]->DCounts/fits[bin][bestlocal->i][bestlocal->j]->PCounts);


			StatErrorP -> SetBinError(bin+1,0);
			StatErrorD -> SetBinError(bin+1,0);
			StatErrorT -> SetBinError(bin+1,0);
			SystError -> SetBinError(bin+1,0);
			SystErrorP -> SetBinError(bin+1,0);
			SystError_unb -> SetBinError(bin+1,0);
			SystErrorP_unb -> SetBinError(bin+1,0);


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

		float CountsP; 	
		if(!usebestonly) CountsP      = WeightedPCounts[bin]->GetMean()   - 8.5 * fits_best[bin]->TCounts;
		else CountsP=fits[bin][bestlocal->i][bestlocal->j]->PCounts - 8.5 * fits_best[bin]->TCounts; 

		float CountsP_unb; 	
		if(!usebestonly) CountsP_unb      = WeightedPCounts_unb[bin]->GetMean()   - 8.5 * fits_best[bin]->TCounts;
		else CountsP_unb=fits[bin][bestlocal->i][bestlocal->j]->PCounts - 8.5 * fits_best[bin]->TCounts; 
	
		float CountsD;	
				
		if(!lowstatDmode){
			if(!usebestonly) CountsD = WeightedDCounts[bin]->GetMean();
			else CountsD=fits[bin][bestlocal->i][bestlocal->j]->DCounts;
			if(fits_best[bin]->TCounts>0) CountsD=CountsD-0.035*CountsD; //Tritium contamination stable in time
		}
	
		else { 
			if(!usebestonly) CountsD = WeightedDCounts[bin]->GetMean() ;//- 1.5 * fits_best[bin]->TCounts;
			else CountsD=fits[bin][bestlocal->i][bestlocal->j]->DCounts ;//- 1.5*fits[bin][bestlocal->i][bestlocal->j]->TCounts;
			CountsD=CountsD-0.033*CountsD; //Tritium contamination stable in time
		}

		float CountsD_unb;	
				
		if(!lowstatDmode){
			if(!usebestonly) CountsD_unb = WeightedDCounts_unb[bin]->GetMean();
			else CountsD_unb=fits[bin][bestlocal->i][bestlocal->j]->DCounts;
			if(fits_best[bin]->TCounts>0) CountsD_unb=CountsD_unb-0.035*CountsD_unb; //Tritium contamination stable in time
		}
	
		else { 
			if(!usebestonly) CountsD_unb = WeightedDCounts_unb[bin]->GetMean() ;//- 1.5 * fits_best[bin]->TCounts;
			else CountsD_unb=fits[bin][bestlocal->i][bestlocal->j]->DCounts ;//- 1.5*fits[bin][bestlocal->i][bestlocal->j]->TCounts;
			CountsD_unb=CountsD_unb-0.033*CountsD; //Tritium contamination stable in time
		}
	






		float Countsnoise= fits[bin][bestlocal->i][bestlocal->j]->Templ_P->Integral(fits_best[bin]->Data->FindBin(1.55), fits_best[bin]->Data->GetNbinsX() ) ;
		ProtonTail->SetBinContent(bin+1,Countsnoise);
		
		CountsP=fits_best[bin]->Data->Integral()-(CountsD- 1.5 * fits_best[bin]->TCounts)- 8.5 * fits_best[bin]->TCounts;	
	

		float CountsT      = fits[bin][bestlocal->i][bestlocal->j]->TCounts;
		float CountsP_prim = fits_best[bin]->PCountsPrim - 8.5 * fits_best[bin]->TCounts;
		float CountsD_prim = fits_best[bin]->DCountsPrim - 1.0 * fits_best[bin]->TCounts;

		float conterrP= 0.1 * fits_best[bin]->TCounts;
		float conterrD= 0.1 * fits_best[bin]->TCounts;
		float conterrT= 0.1 * fits_best[bin]->TCounts;


		HeContError -> SetBinContent(bin+1,conterrD);

		float staterrP= StatErrorP -> GetBinContent(bin+1) * CountsP;
		float staterrD= StatErrorD -> GetBinContent(bin+1) * CountsD;
		float staterrT= StatErrorT -> GetBinContent(bin+1) * CountsT;
	
		float systerrP= SystError -> GetBinContent(bin+1) * CountsP;
		float systerrD= SystError -> GetBinContent(bin+1) * CountsD;
		float systerrT= SystError -> GetBinContent(bin+1) * CountsT;

		ProtonCounts->SetBinContent(bin+1,CountsP);
		ProtonCounts->SetBinError(bin+1,pow(pow(staterrP,2)+pow(systerrP,2)+pow(conterrP,2),0.5));

		ProtonCounts_unb->SetBinContent(bin+1,CountsP_unb);
		ProtonCounts_unb->SetBinError(bin+1,pow(pow(staterrP,2)+pow(systerrP,2)+pow(conterrP,2),0.5));

		if(CountsD>0) DeuteronCounts->SetBinContent(bin+1,CountsD);
		if(CountsD>0) DeuteronCounts->SetBinError(bin+1,pow(pow(staterrD,2)+pow(systerrD,2)+pow(conterrD,2),0.5));		

		if(CountsD_unb>0) DeuteronCounts_unb->SetBinContent(bin+1,CountsD_unb);
		if(CountsD_unb>0) DeuteronCounts_unb->SetBinError(bin+1,pow(pow(staterrD,2)+pow(systerrD,2)+pow(conterrD,2),0.5));		

		TritiumCounts->SetBinContent(bin+1,CountsT);
		TritiumCounts->SetBinError(bin+1,pow(pow(staterrT,2)+pow(systerrT,2)+pow(conterrT,2),0.5));		


		if(fits_best[bin]->PCountsPrim>0){
			ProtonCountsPrim->SetBinContent(bin+1,fits_best[bin]->PCountsPrim);
			ProtonCountsPrim->SetBinError(bin+1,pow(pow(staterrP,2)+pow(systerrP,2)+pow(conterrP,2),0.5) * ProtonCountsPrim->GetBinContent(bin+1)/ProtonCounts->GetBinContent(bin+1) );
		}
		if(fits_best[bin]->DCountsPrim>0){
			DeuteronCountsPrim->SetBinContent(bin+1,fits_best[bin]->DCountsPrim);
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
