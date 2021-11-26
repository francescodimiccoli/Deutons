#include "Flux.h"
#include "FluxRebin.h"
#include "histUtils.h"

TH1F * Twentisevenify_(TH1F* input){

	std::string name = input->GetName();
	TH1F * output = (TH1F*)input->Clone((name+"_27").c_str());
	for(int i=0;i<output->GetNbinsX();i++){
		output->SetBinContent(i+1,input->GetBinContent(i+1)*pow(input->GetBinCenter(i+1),2.7));
		output->SetBinError(i+1,input->GetBinError(i+1)*pow(input->GetBinCenter(i+1),2.7));
	}

	return output;
}

void Regularize_Counts(TH1 * flux){
	TF1* trend = new TF1("trend","[0]*x^[1]",0.1,25);
	flux->Fit("trend","","",0,10);
	for(int i=0;i<flux->GetNbinsX();i++){
		flux->SetBinContent(i+1,trend->Eval(flux->GetBinCenter(i+1)));
		flux->SetBinError(i+1,flux->GetBinError(i+1)/5.);
	}
	return;
}


void Flux::Unfold_Counts(float fit_min, float fit_max,int knots,float offset,bool regularize){
	if(fit_min==0 && fit_max==0) {
		Unfolding_factor = (TH1F*) Counts_density_unfolded->Clone("Unfolding_factor");
        	Unfolding_factor->Reset();
		return;
	}

	std::vector<TH1D *> random_meas;
        std::vector<TSpline3 *> spline_random;
        std::vector<TF1 *> funct_random;
        int iter =3;
        for(int i=0; i<iter;i++) {
 	
		Counts_density->Multiply(ExposureTime_plot_rig);
               	Counts_density->Multiply(Eff_Acceptance_plot_rig);
	
		TH1D * tmp =(TH1D*) Counts_density->Clone("");
		if(i!=0){
			tmp->Reset();
			tmp->FillRandom(Counts_density,Counts_density->Integral());
		}	

		for(int j=0;j<Counts_density->GetNbinsX();j++) tmp->SetBinError(j+1,Counts_density->GetBinError(j+1));

		Counts_density->Divide(ExposureTime_plot_rig);
		Counts_density->Divide(Eff_Acceptance_plot_rig);
		tmp->Divide(ExposureTime_plot_rig);
		tmp->Divide(Eff_Acceptance_plot_rig);
	/*	if(regularize) Regularize_Counts(tmp);
		else{
			tmp->Smooth();
			tmp->Smooth();
			tmp->Smooth();
		}*/
		random_meas.push_back(tmp);
	}

	TGraph * expo = new TGraph(ExposureTime_plot_rig);
	TGraph * acce = new TGraph(Eff_Acceptance_plot_rig);
	
	for(int i=0;i<iter;i++) {
		UnfoldRes Res = Unfold(migr_matr,random_meas[i], expo,acce, fit_min, fit_max,knots,offset);	
		spline_random.push_back(Res.spline);
		funct_random.push_back(Res.funct);
	}

	Unfolding_factor = (TH1F*) random_meas[0]->Clone("Unfolding_factor");
	Unfolding_factor->Sumw2();
	Unfolding_factor_raw = (TH1F*) random_meas[0]->Clone("Unfolding_factor_raw");
	Unfolding_factor_raw->Sumw2();
	
        for(int i=0;i< Unfolding_factor->GetNbinsX();i++) {
                float mean=0;
                for(int j=0;j<spline_random.size();j++) {
                        mean+=splineIntegral(spline_random[j],Unfolding_factor->GetBinLowEdge(i+1),Unfolding_factor->GetBinLowEdge(i+2))/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]);
                }
                float devst=0;
                for(int j=0;j<spline_random.size();j++) {
                        float value = splineIntegral(spline_random[j],Unfolding_factor->GetBinLowEdge(i+1),Unfolding_factor->GetBinLowEdge(i+2))/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]) ;
                        devst+=pow(value - mean/spline_random.size(),2);
                }
              //  devst=pow(devst/spline_random.size(),0.5);
		if(Counts_density->GetBinContent(i+1)>0){
			Unfolding_factor_raw->SetBinContent(i+1, (mean/spline_random.size()) /Unfolding_factor->GetBinContent(i+1));
	                Unfolding_factor_raw->SetBinError(i+1,0);

		}
		else{
                        Unfolding_factor_raw->SetBinContent(i+1,0);
                        Unfolding_factor_raw->SetBinError(i+1,0);
                }
        }


	TF1 * modelunf = new TF1("modelunf","pol2",0,30);
	Unfolding_factor_raw->Fit("modelunf","w0");
	
	for(int i=0;i< Unfolding_factor->GetNbinsX();i++) {
		Unfolding_factor->SetBinContent(i+1,modelunf->Eval(Unfolding_factor->GetBinCenter(i+1)));
		Unfolding_factor->SetBinError(i+1,modelunf->Eval(Unfolding_factor->GetBinCenter(i+1))*0.025);
	
	}	

	for(int i=0;i<Counts_density_unfolded->GetNbinsX();i++)
		{ 
			if(!forcefolded) Counts_density_unfolded->SetBinContent(i+1,Counts_density_unfolded->GetBinContent(i+1)*Unfolding_factor->GetBinContent(i+1));
			float error = pow(Counts_density->GetBinError(i+1)/Counts_density->GetBinContent(i+1),2);
			error+=pow(Unfolding_factor->GetBinError(i+1)/Unfolding_factor->GetBinContent(i+1),2);
			Counts_density_unfolded->SetBinError(i+1,pow(error,0.5)*Counts_density_unfolded->GetBinContent(i+1));
		}


	test = new TCanvas(("Test Unfolding"+basename).c_str());
	test->cd();
	TH1F * Counts_density_draw= (TH1F*)Counts_density->Clone("Counts density");
	Counts_density_draw->SetLineColor(2);
	Counts_density_draw->SetLineWidth(2);
	//TH2F * frame=new TH2F("frame","frame",1000, 0.5*Counts_density_draw->GetBinLowEdge(1),1.5*Counts_density_draw->GetBinLowEdge(Counts_density_draw->GetNbinsX()),1000,0.01*Counts_density_draw->GetBinContent(Counts_density_draw->GetMaximumBin()),1.5*Counts_density_draw->GetBinContent(Counts_density_draw->GetMaximumBin()) );
	Counts_density_draw->Draw();
	for(int i=0; i<iter;i++) random_meas[i]->Draw("same");
	for(int i=0;i<iter;i++) {spline_random[i]->SetMarkerStyle(8); spline_random[i]->Draw("PCsame");}
	for(int i=0;i<iter;i++) funct_random[i]->Draw("same");




	return;

}


void Flux::Eval_Flux(float corr_acc, float fit_min, float fit_max,int knots, float offset,bool regularize){


	// COUNTS
	cout<<endl;
	cout<<"********************************************"<<endl;
	if(Counts>0) {	
		cout<<"Counts "<<Counts->GetName()<<" "<<Counts->GetEntries()<<endl;	

		Counts_density 		= ConvertBinnedHisto(Counts,"Counts density",bins,false);
		Counts_density_unfolded = ConvertBinnedHisto(Counts,"Counts density Unfolded",bins,false);
		Counts_density_Ekin 	= ConvertBinnedHisto(Counts,"Counts density Ekin",bins,true);
		Counts_density_Ekin_unf	= ConvertBinnedHisto(Counts,"Counts density Ekin unfolded",bins,true);
		Counts_density 		-> Sumw2();
		Counts_density_unfolded -> Sumw2();
		Counts_density_Ekin 	-> Sumw2();
		Counts_density_Ekin_unf	-> Sumw2();
		Counts_density27= Twentisevenify_(Counts_density);
		DeltaE=(TH1F*)Counts_density->Clone("DeltaE"); 	
		DeltaE->Reset();

	}
	else { FluxEstim = new TH1F("dummy","dummy",10,0,10);
		return;
	}	

	// DeltaE
	for(int i=0;i<bins.size();i++){
		Counts_density->SetBinError(i+1,Counts->GetBinError(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));
                Counts_density->SetBinContent(i+1,Counts->GetBinContent(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));	
		Counts_density_unfolded->SetBinError(i+1,Counts->GetBinError(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));
                Counts_density_unfolded->SetBinContent(i+1,Counts->GetBinContent(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));	
		Counts_density_Ekin->SetBinError(i+1,Counts->GetBinError(i+1)/(bins.EkPerMasTOIBins()[i+1]-bins.EkPerMasTOIBins()[i]));
                Counts_density_Ekin->SetBinContent(i+1,Counts->GetBinContent(i+1)/(bins.EkPerMasTOIBins()[i+1]-bins.EkPerMasTOIBins()[i]));
		Counts_density_Ekin_unf->SetBinError(i+1,Counts->GetBinError(i+1)/(bins.EkPerMasTOIBins()[i+1]-bins.EkPerMasTOIBins()[i]));
                Counts_density_Ekin_unf->SetBinContent(i+1,Counts->GetBinContent(i+1)/(bins.EkPerMasTOIBins()[i+1]-bins.EkPerMasTOIBins()[i]));
	
		Counts_density27->SetBinError(i+1,Counts_density27->GetBinError(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));
                Counts_density27->SetBinContent(i+1,Counts_density27->GetBinContent(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));	
		DeltaE->SetBinContent(i+1,(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));
	}

	// EXPOSURE TIME	
	if(ExposureTime) {
		ExposureTime_plot_rig =  ConvertBinnedHisto(ExposureTime,"ExposureTime R",bins,false);
		ExposureTime_plot     =  ConvertBinnedHisto(ExposureTime,"ExposureTime Ekin",bins,true);
		Counts_density_Ekin -> Divide(ExposureTime_plot);
        	Counts_density_Ekin_unf -> Divide(ExposureTime_plot);
                Counts_density -> Divide(ExposureTime_plot_rig);
	        Counts_density_unfolded -> Divide(ExposureTime_plot_rig);
	}


		
	if(Eff_Acceptance&&ExposureTime) {

		Eff_Acceptance_plot = ConvertBinnedHisto(Eff_Acceptance,"Eff. Acceptance Ekin",bins,true);
		Eff_Acceptance_plot_rig = ConvertBinnedHisto(Eff_Acceptance,"Eff. Acceptance R",bins,false);

		Eff_Acceptance_plot->Scale(corr_acc);
		Eff_Acceptance_plot_rig->Scale(corr_acc);
		
		Counts_density -> Divide(Eff_Acceptance_plot_rig);
	        Counts_density_unfolded -> Divide(Eff_Acceptance_plot_rig);
	
		Unfold_Counts(fit_min,fit_max,knots,offset,regularize);	

		Counts_density -> Multiply(Eff_Acceptance_plot_rig);
	        Counts_density_unfolded -> Multiply(Eff_Acceptance_plot_rig);
	

		if(Unfolding_factor) for(int i=0; i<Counts_density_Ekin_unf->GetNbinsX();i++) { Counts_density_Ekin_unf->SetBinContent(i+1,Counts_density_Ekin->GetBinContent(i+1)*Unfolding_factor->GetBinContent(i+1));}


	}

	FluxEstim = (TH1F *) Counts_density_Ekin->Clone();
	FluxEstim -> SetName((basename+"_Flux").c_str());
	FluxEstim -> SetTitle((basename+"_Flux").c_str());
	FluxEstim -> Sumw2();

	FluxEstim_ekin_unf = (TH1F *) Counts_density_Ekin_unf->Clone();
	FluxEstim_ekin_unf -> SetName((basename+"_Flux_ekin_unf").c_str());
	FluxEstim_ekin_unf -> SetTitle((basename+"_Flux_ekin_unf").c_str());
	FluxEstim_ekin_unf -> Sumw2();


	FluxEstim_rig = (TH1F*) Counts_density ->Clone();
	FluxEstim_rig -> SetName((basename+"_Flux_rig").c_str());
	FluxEstim_rig -> SetTitle((basename+"_Flux_rig").c_str());
	FluxEstim_rig -> Sumw2();
	
	FluxEstim_unf = (TH1F*) Counts_density_unfolded ->Clone();
	FluxEstim_unf -> SetName((basename+"_Flux_unf").c_str());
	FluxEstim_unf -> SetTitle((basename+"_Flux_unf").c_str());
	FluxEstim_unf -> Sumw2();
	

	//ACCEPTANCE
	if(Eff_Acceptance) {

		/*Eff_Acceptance_plot = ConvertBinnedHisto(Eff_Acceptance,"Eff. Acceptance Ekin",bins,true);
		Eff_Acceptance_plot_rig = ConvertBinnedHisto(Eff_Acceptance,"Eff. Acceptance R",bins,false);

		Eff_Acceptance_plot->Scale(corr_acc);
		Eff_Acceptance_plot_rig->Scale(corr_acc);
		*/
	
		//ModelAcceptanceWithSpline(0);
		
		FluxEstim -> Divide(Eff_Acceptance_plot);
		FluxEstim_ekin_unf -> Divide(Eff_Acceptance_plot);
		FluxEstim_rig -> Divide(Eff_Acceptance_plot_rig);
		FluxEstim_unf -> Divide(Eff_Acceptance_plot_rig);

	/*for(int i=0;i<FluxEstim_rig->GetNbinsX();i++){
			FluxEstim ->SetBinContent(i+1,FluxEstim->GetBinContent(i+1)/AcceptanceModel->Eval(bins.ekpermassbincent_TOI[i]));
			FluxEstim ->SetBinError(i+1,FluxEstim->GetBinError(i+1)/AcceptanceModel->Eval(bins.ekpermassbincent_TOI[i]));
			FluxEstim_rig ->SetBinContent(i+1,FluxEstim_rig->GetBinContent(i+1)/AcceptanceModel->Eval(bins.RigTOIBinsCent()[i]));
			FluxEstim_rig ->SetBinError(i+1,FluxEstim_rig->GetBinError(i+1)/AcceptanceModel->Eval(bins.RigTOIBinsCent()[i]));
			FluxEstim_unf ->SetBinContent(i+1,FluxEstim_unf->GetBinContent(i+1)/AcceptanceModel->Eval(bins.RigTOIBinsCent()[i]));
			FluxEstim_unf ->SetBinError(i+1,FluxEstim_unf->GetBinError(i+1)/AcceptanceModel->Eval(bins.RigTOIBinsCent()[i]));
		}
	*/
	}

	Eval_Errors();

	return;		
}

void Flux::SaveResults(FileSaver finalhistos){
	std::cout<<"Results: "<<FluxEstim<<" "<<ExposureTime<<" "<<Eff_Acceptance<<std::endl;

	if(FluxEstim) finalhistos.Add(FluxEstim); 	
	if(FluxEstim_ekin_unf) finalhistos.Add(FluxEstim_ekin_unf); 	
	if(FluxEstim_rig) finalhistos.Add(FluxEstim_rig);
	if(FluxEstim_unf) finalhistos.Add(FluxEstim_unf);
	 	
	if(ExposureTime) { 
				finalhistos.Add(ExposureTime);
				finalhistos.Add(ExposureTime_plot);
				finalhistos.Add(ExposureTime_plot_rig);
	}

	if(Eff_Acceptance) { finalhistos.Add(Eff_Acceptance);
		finalhistos.Add(Eff_Acceptance_plot);	
		finalhistos.Add(Eff_Acceptance_plot_rig);	
	}

	if(Counts){ finalhistos.Add(Counts);
		    finalhistos.Add(Counts_density);
		    finalhistos.Add(Counts_density27);
		    finalhistos.Add(Counts_density_unfolded);
		    finalhistos.Add(Counts_density_Ekin);
	}
        if(DeltaE) finalhistos.Add(DeltaE);
	if(Unfolding_factor_raw) finalhistos.Add(Unfolding_factor_raw);
	if(Unfolding_factor) finalhistos.Add(Unfolding_factor);
	if(test) finalhistos.Add(test);
	if(migr_matr) finalhistos.Add(migr_matr);

	if(Counts_statErr) {
		Counts_statErr = ConvertBinnedHisto(Counts_statErr,"Counts_Stat_Error",bins,false);
		finalhistos.Add(Counts_statErr);
	}
	if(Counts_systErr) {
		Counts_systErr = ConvertBinnedHisto(Counts_systErr,"Counts_Syst_Error",bins,false);
		finalhistos.Add(Counts_systErr);

	}
	if(Acc_Err) {
		Acc_Err = ConvertBinnedHisto(Acc_Err,"Acceptance_Syst_Error",bins,false);
		finalhistos.Add(Acc_Err);

	}

	if(Unfolding_Err) {
		Unfolding_Err = ConvertBinnedHisto(Unfolding_Err,"Unfolding_Syst_Error",bins,false);
		finalhistos.Add(Unfolding_Err);

	}

	if(FluxEstim_unf){
		FluxEstim_unf_stat = (TH1F*) FluxEstim_unf->Clone();
	 	FluxEstim_unf_stat -> SetName((basename+"_Flux_unf_stat").c_str());
        	FluxEstim_unf_stat -> SetTitle((basename+"_Flux_unf_stat").c_str());
        	FluxEstim_unf_stat -> Sumw2();

		for(int i=0;i<FluxEstim_unf_stat->GetNbinsX();i++){
			float error=0;
			if(Counts_statErr) error+=pow(Counts_statErr->GetBinContent(i+1),2);
			//if(Counts_systErr) error+=pow(Counts_systErr->GetBinContent(i+1),2);
			//if(Acc_Err) 	   error+=pow( pow(pow(Acc_Err->GetBinContent(i+1),2)-pow(0.02,2),0.5),2);
			if(Acc_Err) 	   error+=pow(Acc_Err->GetBinContent(i+1),2);
			if(Unfolding_Err)  error+=pow(Unfolding_Err->GetBinContent(i+1),2);
			FluxEstim_unf_stat->SetBinError(i+1, pow(error,0.5)*FluxEstim_unf_stat->GetBinContent(i+1));
		}
		finalhistos.Add(FluxEstim_unf_stat);

	}
	
	if(FluxEstim_rig){
		FluxEstim_rig_stat = (TH1F*) FluxEstim_rig->Clone();
	 	FluxEstim_rig_stat -> SetName((basename+"_Flux_rig_stat").c_str());
        	FluxEstim_rig_stat -> SetTitle((basename+"_Flux_rig_stat").c_str());
        	FluxEstim_rig_stat -> Sumw2();

		for(int i=0;i<FluxEstim_rig_stat->GetNbinsX();i++){
			float error=0;
			if(Counts_statErr) error+=pow(Counts_statErr->GetBinContent(i+1),2);
			//if(Counts_systErr) error+=pow(Counts_systErr->GetBinContent(i+1),2);
			if(Acc_Err) 	   error+=pow(Acc_Err->GetBinContent(i+1),2);
			if(Unfolding_Err)  error+=pow(Unfolding_Err->GetBinContent(i+1),2);
			FluxEstim_rig_stat->SetBinError(i+1, pow(error,0.5)*FluxEstim_rig_stat->GetBinContent(i+1));
		}
		finalhistos.Add(FluxEstim_rig_stat);

	}
	


	finalhistos.writeObjsInFolder(("Fluxes/"+basename).c_str());
}


void Flux::Eval_ExposureTime(Variables * vars, TTree * treeDT,FileSaver finalhistos,bool refill){
		
	if(!ExposureTime||refill){
	
		ExposureTime = new TH1F(exposurename.c_str(),exposurename.c_str(),bins.size(),0,bins.size());

		int ActualTime=0;
		vars->ReadBranches(treeDT);
		for(int i=0;i<treeDT->GetEntries()/FRAC;i++){
			UpdateProgressBar(i, treeDT->GetEntries()/FRAC);
			treeDT->GetEvent(i);
			vars->Update();
			if((int)vars->U_time!=ActualTime) {
				UpdateZoneLivetime(vars->Livetime,vars->Rcutoff_IGRFRTI,ExposureTime,bins);
				ActualTime=vars->U_time;
			}
		}
		finalhistos.Add(ExposureTime); 	
		finalhistos.writeObjsInFolder(("Fluxes/"+basename).c_str());

	}
	return;
}


TH1F * Flux::Eval_FluxRatio(Flux * Denominator,std::string name){

        TH1F * Numerator = (TH1F*) FluxEstim -> Clone();
        Numerator->SetName(name.c_str());
        Numerator->SetTitle(name.c_str());
        Numerator->Sumw2();

        Numerator->Divide(Denominator->GetFlux());
        for(int i=0;i<Numerator->GetNbinsX();i++){
                float A = FluxEstim -> GetBinContent(i+1);
                float B = Denominator->GetFlux() -> GetBinContent(i+1);
                float sigA = FluxEstim -> GetBinError(i+1);
                float sigB = Denominator->GetFlux() -> GetBinError(i+1);
                Numerator->SetBinError(i+1,pow(pow(sigA/A,2)+pow(sigB/B,2)+2*sigA/A*sigB/B,0.5)*A/B);
        }

        return Numerator;
}

void Flux::Eval_Errors(){
	if(Counts>0){
		if(!Counts_statErr){
			Counts_statErr = (TH1F*) Counts->Clone("Counts_statErr");
			Counts_systErr = (TH1F*) Counts->Clone("Counts_systErr");
			Counts_statErr->Reset();	
			Counts_systErr->Reset();	

			for(int i=0;i<Counts->GetNbinsX();i++){
				if(Counts->GetBinContent(i+1)){	
					Counts_statErr->SetBinContent(i+1,Counts->GetBinError(i+1)/Counts->GetBinContent(i+1));
					Counts_systErr->SetBinContent(i+1,0);
				}
			}
		}

		else{
			Counts_systErr = (TH1F*) Counts->Clone("Counts_systErr");
			Counts_systErr->Reset();	

			for(int i=0;i<Counts->GetNbinsX();i++){
				Counts_systErr->SetBinContent(i+1,pow(pow(Counts->GetBinError(i+1)/Counts->GetBinContent(i+1),2) - pow(Counts_statErr->GetBinContent(i+1),2),0.5));
								
			}
		}
	}

	if(Eff_Acceptance>0){
		Acc_Err = (TH1F*) Eff_Acceptance->Clone("Acc_Err");	
		Acc_Err->Reset();
		for(int i=0;i<Eff_Acceptance->GetNbinsX();i++){
                                Acc_Err->SetBinContent(i+1,Eff_Acceptance->GetBinError(i+1)/Eff_Acceptance->GetBinContent(i+1));
                        }

	}
	if(Unfolding_factor>0){
		Unfolding_Err = (TH1F*) Unfolding_factor->Clone("Unfolding_Err");	
		Unfolding_Err->Reset();
		for(int i=0;i<Unfolding_factor->GetNbinsX();i++){
                                Unfolding_Err->SetBinContent(i+1,Unfolding_factor->GetBinError(i+1)/Unfolding_factor->GetBinContent(i+1));
                        }

	}


}

void Flux::ModelAcceptanceWithSpline(float shift){


	int nodes = bins.size()/4+2;
	double p[nodes];
	double x[nodes];

	for(int i=0;i<nodes;i++){
		x[i]=bins.RigTOIBins()[0]+((bins.RigTOIBins()[bins.size()] - bins.RigTOIBins()[0])/(nodes-1))*i;
		p[i]=Eff_Acceptance_plot_rig->GetBinContent(Eff_Acceptance_plot_rig->FindBin(x[i]));
	}
	x[nodes-1]=bins.RigTOIBins()[bins.size()];
	p[nodes-1]=Eff_Acceptance_plot_rig->GetBinContent(Eff_Acceptance_plot_rig->FindBin(x[nodes-1])-1);
	
	cout<<"******** Model Acc. "<<basename<<" ***********"<<endl;
	std::vector<float> spline_x;
	for(int i=0;i<nodes;i++) {
		spline_x.push_back(shift + x[i]);
		cout<<x[i]<<" "<<p[i]<<endl;
	}
	Model *M = new Model(spline_x);

	AcceptanceModel = new TF1(basename.c_str(),M,&Model::Function,shift,50,nodes,"Model","Function");
 
	for(int i=0;i<nodes;i++) {
		AcceptanceModel->SetParameter(i,p[i]);
		AcceptanceModel->SetParLimits(i,0.8*p[i],1.2*p[i]);
		if(i==nodes-1)	AcceptanceModel->SetParLimits(i,0.95*p[i],1.02*p[i]);
	
	}	

	Eff_Acceptance_plot_rig->Fit(basename.c_str(),"","M",x[0],x[nodes-1]);

}
