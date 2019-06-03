#include "Flux.h"
#include "FluxRebin.h"
#include "histUtils.h"



TH1F * EvalEffAcc(Efficiency* Eff, TH1F* ForAcc, Binning bins, MCPar param){
	SetUpUsualBinning();	
	cout<<"********** MC flux reweighting **********"<<endl;
	Variables * vars = new Variables();
	cout<<"********** ForAcc **************"<<endl;
	float normalization = ForAcc->GetEntries()/FRAC*(pow(param.Trigrate,-1));
	cout<<ForAcc->GetEntries()<<" "<<normalization<<" "<<param.tot_trig<<endl;
	//float normalization = param.tot_trig/FRAC;
	TH1F * Numerator = (TH1F *) (Eff->GetAfter())->Clone();	
	
	//no_reweighting
	/*for(int i=0;i<bins.size();i++){
		if(EffAcc -> GetBinContent(i+1)>0){
			cout<<"********** MC flux reweighting **********"<<endl;
			float range;
			range = log(param.Rmax)-log(param.Rmin);
			
			float gen_bins= normalization*(log(bins.RigBin(i+1))-log(bins.RigBin(i)))/range;  
			EffAcc -> SetBinContent(i+1,EffAcc -> GetBinContent(i+1)/gen_bins); 
			EffAcc -> SetBinError(i+1,0);//pow(GenSpectrum -> GetBinContent(i+1),0.5)/gen_bins); 			
			
		
		}
	}*/
	/////////////// 

	//with reweighting
	TH1F * Denominator = new TH1F ("Denominator","Denominator",bins.size(),0,bins.size());
	Histogram Spectrum = vars->reweighter.getTo();
	Histogram LogNorm  = vars->reweighter.getFrom();

	cout<<"********** Gen. Spectrum **********"<<endl;
	TH1F * FullGenSpec = new TH1F("FullGenSpec","FullGenSpec",ForAcceptance.size(),0,ForAcceptance.size());
	//total triggers in range
	for(int i=0;i<bins.size();i++){
		float bincontent = normalization*(log(bins.RigBin(i+1))-log(bins.RigBin(i)))/log(param.Rmax)-log(param.Rmin);
		float meanweight = Spectrum.integrate(bins.RigBin(i),bins.RigBin(i+1)) / LogNorm.integrate(bins.RigBin(i),bins.RigBin(i+1));
		Denominator->SetBinContent(i+1,param.art_ratio*bincontent*meanweight); //vars->reweighter.getWeight(bins.RigBinsCent()[i]));
		Denominator->SetBinError(i+1,0);
	}
	//Denominator->Scale(param.art_ratio*rangefactor*normalization/(Denominator->Integral()));

	cout<<"********** THE DIVISION **********"<<endl;
	cout<<Numerator->GetName()<<endl;
	TH1F * EffAcc= (TH1F*)Numerator->Clone();
	EffAcc->Divide(Denominator);
	/////
	
	EffAcc ->SetName("Eff. Acceptance");
	EffAcc ->SetTitle("Eff. Acceptance"); 
	EffAcc -> Sumw2();
	EffAcc -> Scale(47.78/param.gen_factor);


	cout<<"*********** Eff. Acceptance parameters ************"<<endl;
	cout<<"TOT. Ev. gen: "<<normalization<<endl;
	cout<<"Trig. rate: "<<param.Trigrate<<endl;
	return EffAcc;
}



void Flux::Eval_Flux(){
	// COUNTS
	cout<<endl;
	cout<<"********************************************"<<endl;
	cout<<"Counts "<<Counts<<" "<<Counts->GetName()<<endl;	
	if(Counts>0) {			FluxEstim = (TH1F *) Counts->Clone();
					Counts_Err  = (TH1F *) Counts->Clone();
					for(int i=0;i<Counts_Err->GetNbinsX();i++) {
						if(Counts_Err->GetBinContent(i+1)>0&&Counts_Err->GetBinError(i+1)>0)
							Counts_Err->SetBinContent(i+1,Counts_Err->GetBinError(i+1)/Counts_Err->GetBinContent(i+1));
						else Counts_Err->SetBinContent(i+1,0);	
						Counts_Err->SetBinError(i+1,0);
					}
					Counts_Err->SetName("Counts Error");
        				Counts_Err->SetTitle("Counts Error");

	}
	else { FluxEstim = new TH1F("dummy","dummy",10,0,10);
		return;
	}	
	FluxEstim -> SetName((basename+"_Flux").c_str());
	FluxEstim -> SetTitle((basename+"_Flux").c_str());
	FluxEstim -> Sumw2();

	////
	
	/*if(MCEfficiency->GetStat_Error()){	
		Acc_StatErr = (TH1F *) MCEfficiency->GetStat_Error()->Clone();
		Acc_SystErr = (TH1F *) MCEfficiency->GetSyst_Error()->Clone();
		Acc_StatErr->SetName("Acceptance Stat. Error");
		Acc_StatErr->SetTitle("Acceptance Stat. Error");	
		Acc_SystErr->SetName("Acceptance Syst. Error");
		Acc_SystErr->SetTitle("Acceptance Syst. Error");	
	}*/
	
	// EXPOSURE TIME	
	if(ExposureTime) FluxEstim -> Divide(ExposureTime);
	////
	//ACCEPTANCE
	if(MCEfficiency->GetAfter()){
	cout<<MCEfficiency->GetAfter()<<endl;	
	Eff_Acceptance = (TH1F *) MCEfficiency->GetAfter()->Clone();
		Acc_StatErr = (TH1F *)  MCEfficiency->GetAfter()->Clone();
		Acc_SystErr = (TH1F *)  MCEfficiency->GetAfter()->Clone();
	



	         Eff_Acceptance = EvalEffAcc(MCEfficiency,ForAcceptance,bins,param);	
		 for(int i=0;i<EfficiencyCorrections.size();i++) Eff_Acceptance -> Multiply(EfficiencyCorrections[i]->GetGlobCorrection()); 
		 for(int i=0;i<EfficiencyFromData.size();i++)    Eff_Acceptance -> Multiply(EfficiencyFromData[i]->GetGlobEfficiency()); 

		 if(EfficiencyCorrections.size()>0){

			 Acc_StatErr = (TH1F*) EfficiencyCorrections[0]->GetStat_Err()->Clone();
			 Acc_SystErr = (TH1F*) EfficiencyCorrections[0]->GetSyst_Err()->Clone();
			 Acc_StatErr->SetName("Acceptance Stat. Error");
			 Acc_StatErr->SetTitle("Acceptance Stat. Error");	
			 Acc_SystErr->SetName("Acceptance Syst. Error");
			 Acc_SystErr->SetTitle("Acceptance Syst. Error");	
			
			Acc_StatErr->Multiply(Acc_StatErr);
                        Acc_SystErr->Multiply(Acc_SystErr);
			
    			 for(int i=EfficiencyCorrections.size()-1;i>0;i--) {
				 TH1F * stat = (TH1F*) EfficiencyCorrections[i]->GetStat_Err()->Clone();
				 stat->Multiply(EfficiencyCorrections[i]->GetStat_Err());
				 Acc_StatErr->Add(stat);
			
				 TH1F * syst = (TH1F*) EfficiencyCorrections[i]->GetSyst_Err()->Clone();
				 syst->Multiply(EfficiencyCorrections[i]->GetSyst_Err());
				 Acc_SystErr->Add(syst);
			 }
			
			 for(int i=0;i<Acc_StatErr->GetNbinsX();i++) {
				 Acc_StatErr->SetBinContent(i+1,pow(Acc_StatErr->GetBinContent(i+1),0.5));
				 Acc_StatErr->SetBinError(i+1,0);
				 Acc_SystErr->SetBinContent(i+1,pow(Acc_SystErr->GetBinContent(i+1),0.5));
				 Acc_SystErr->SetBinError(i+1,0);
			 }

 
		}
		
		     //FluxEstim -> Divide(Eff_Acceptance);
	}
	///
	
	// DeltaE (CONV IN RIG)
	FluxEstim_rig = (TH1F*) FluxEstim ->Clone();
	FluxEstim_rig -> SetName((basename+"_Flux_rig").c_str());
	FluxEstim_rig -> SetTitle((basename+"_Flux_rig").c_str());
	FluxEstim_rig -> Sumw2();
	for(int i=0;i<bins.size();i++){
		FluxEstim->SetBinError(i+1,FluxEstim->GetBinError(i+1)/(bins.EkPerMasTOIBins()[i+1]-bins.EkPerMasTOIBins()[i]));
		FluxEstim->SetBinContent(i+1,FluxEstim->GetBinContent(i+1)/(bins.EkPerMasTOIBins()[i+1]-bins.EkPerMasTOIBins()[i]));
		FluxEstim_rig->SetBinError(i+1,FluxEstim_rig->GetBinError(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));
		FluxEstim_rig->SetBinContent(i+1,FluxEstim_rig->GetBinContent(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));
	}
	///
	return;		
}

void Flux::SaveResults(FileSaver finalhistos){
	std::cout<<"Results: "<<FluxEstim<<" "<<ExposureTime<<" "<<Eff_Acceptance<<std::endl;

	if(FluxEstim) finalhistos.Add(FluxEstim); 	
	if(FluxEstim_rig) finalhistos.Add(FluxEstim_rig);
	 	
	if(ExposureTime) finalhistos.Add(ExposureTime);
	if(Eff_Acceptance)finalhistos.Add(Eff_Acceptance);
	if(Acc_StatErr) finalhistos.Add(Acc_StatErr);
	if(Acc_SystErr) finalhistos.Add(Acc_SystErr);
	if(Counts_Err) finalhistos.Add(Counts_Err);
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
				UpdateZoneLivetime(vars->Livetime,vars->Rcutoff,ExposureTime,bins);
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


