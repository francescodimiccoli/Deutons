#include "Flux.h"
#include "FluxRebin.h"
#include "histUtils.h"


void Flux::Eval_Flux(){
	// COUNTS
	cout<<endl;
	cout<<"********************************************"<<endl;
	if(Counts>0) {	
					cout<<"Counts "<<Counts->GetName()<<" "<<Counts->GetEntries()<<endl;	
	
					FluxEstim = (TH1F *) Counts->Clone();
					Counts_Err  = (TH1F *) Counts->Clone();
					for(int i=0;i<Counts_Err->GetNbinsX();i++) {
						//FluxEstim->SetBinContent(i+1,1);
						FluxEstim->SetBinError(i+1,0);
						if(Counts_Err->GetBinContent(i+1)>0&&Counts_Err->GetBinError(i+1)>0)
							Counts_Err->SetBinContent(i+1,0);//Counts_Err->GetBinError(i+1)/Counts_Err->GetBinContent(i+1));
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

	// EXPOSURE TIME	
	if(ExposureTime) FluxEstim -> Divide(ExposureTime);
	//
	//ACCEPTANCE
	if(Eff_Acceptance) FluxEstim -> Divide(Eff_Acceptance);

	// DeltaE (CONV IN RIG)
	FluxEstim_rig = (TH1F*) FluxEstim ->Clone();
	FluxEstim_rig -> SetName((basename+"_Flux_rig").c_str());
	FluxEstim_rig -> SetTitle((basename+"_Flux_rig").c_str());
	FluxEstim_rig -> Sumw2();
	for(int i=0;i<bins.size();i++){
		if((bins.EkPerMasTOIBins()[i+1]-bins.EkPerMasTOIBins()[i])>0){
		FluxEstim->SetBinError(i+1,FluxEstim->GetBinError(i+1)/(bins.EkPerMasTOIBins()[i+1]-bins.EkPerMasTOIBins()[i]));
		FluxEstim->SetBinContent(i+1,FluxEstim->GetBinContent(i+1)/(bins.EkPerMasTOIBins()[i+1]-bins.EkPerMasTOIBins()[i]));
		FluxEstim_rig->SetBinError(i+1,FluxEstim_rig->GetBinError(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));
		FluxEstim_rig->SetBinContent(i+1,FluxEstim_rig->GetBinContent(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));
		}
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


