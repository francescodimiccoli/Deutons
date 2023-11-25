#include "ResultMerger.h"


TH1F * Twentisevenify(TH1F* input){

	std::string name = input->GetName();
	TH1F * output = (TH1F*)input->Clone((name+"_27").c_str());
	for(int i=0;i<output->GetNbinsX();i++){
		output->SetBinContent(i+1,input->GetBinContent(i+1)*pow(input->GetBinCenter(i+1),2.7));
		output->SetBinError(i+1,input->GetBinError(i+1)*pow(input->GetBinCenter(i+1),2.7));
	}

	return output;
}



ResultMerger::ResultMerger(FileSaver finalhistos, std::string name, RangeMerger Global, Flux * FluxTOF, Flux * FluxNaF, Flux * FluxAgl, Particle particle, bool nafpriority, ResultMerger * Fluxsum, ResultMerger * Fluxtosubtract,TH1F * forbinning, Flux * FluxR,bool Aver){

	Name=name;
	is_ave=Aver;
	if(!FluxTOF || !FluxNaF || !FluxAgl) return;


	if(particle.getA()/particle.getZ() <=1.5){

		if(FluxTOF->GetEffAcceptance()&&FluxNaF->GetEffAcceptance()&&FluxAgl->GetEffAcceptance())
			effAcc = Global.MergeSubDResult_P(FluxTOF->GetEffAcceptance(),FluxNaF->GetEffAcceptance(),FluxAgl->GetEffAcceptance(),nafpriority);

		if(FluxTOF->GetMCAcceptance()&&FluxNaF->GetMCAcceptance()&&FluxAgl->GetMCAcceptance())
			effAccMC = Global.MergeSubDResult_P(FluxTOF->GetMCAcceptance(),FluxNaF->GetMCAcceptance(),FluxAgl->GetMCAcceptance(),nafpriority);

		if(FluxTOF->GetMCAcceptance_raw()&&FluxNaF->GetMCAcceptance_raw()&&FluxAgl->GetMCAcceptance_raw())
			effAccMC_raw = Global.MergeSubDResult_P(FluxTOF->GetMCAcceptance_raw(),FluxNaF->GetMCAcceptance_raw(),FluxAgl->GetMCAcceptance_raw(),nafpriority);


		if(FluxTOF->GetFlux()&&FluxNaF->GetFlux()&&FluxAgl->GetFlux())
			flux_ekin = Global.MergeSubDResult_P(FluxTOF->GetFlux(),FluxNaF->GetFlux(),FluxAgl->GetFlux(),nafpriority);

		if(FluxTOF->GetFlux_ekin_unf()&&FluxNaF->GetFlux_ekin_unf()&&FluxAgl->GetFlux_ekin_unf())
			flux_ekin_unf = Global.MergeSubDResult_P(FluxTOF->GetFlux_ekin_unf(),FluxNaF->GetFlux_ekin_unf(),FluxAgl->GetFlux_ekin_unf(),nafpriority);

		if(FluxTOF->GetFlux_rig()&&FluxNaF->GetFlux_rig()&&FluxAgl->GetFlux_rig())
			flux = Global.MergeSubDResult_P(FluxTOF->GetFlux_rig(),FluxNaF->GetFlux_rig(),FluxAgl->GetFlux_rig(),nafpriority);

		if(FluxTOF->GetFlux_rig_stat()&&FluxNaF->GetFlux_rig_stat()&&FluxAgl->GetFlux_rig_stat())
			flux_stat = Global.MergeSubDResult_P(FluxTOF->GetFlux_rig_stat(),FluxNaF->GetFlux_rig_stat(),FluxAgl->GetFlux_rig_stat(),nafpriority);

		if(FluxTOF->GetFlux_unf()&&FluxNaF->GetFlux_unf()&&FluxAgl->GetFlux_unf())
			flux_unf= Global.MergeSubDResult_P(FluxTOF->GetFlux_unf(),FluxNaF->GetFlux_unf(),FluxAgl->GetFlux_unf(),nafpriority);
 
		if(FluxTOF->GetFlux_unf_stat()&&FluxNaF->GetFlux_unf_stat()&&FluxAgl->GetFlux_unf_stat())
			flux_unf_stat= Global.MergeSubDResult_P(FluxTOF->GetFlux_unf_stat(),FluxNaF->GetFlux_unf_stat(),FluxAgl->GetFlux_unf_stat()),nafpriority; 

		if(FluxTOF->GetCounts()&&FluxNaF->GetCounts()&&FluxAgl->GetCounts())
			counts = Global.MergeSubDResult_P(FluxTOF->GetCounts(),FluxNaF->GetCounts(),FluxAgl->GetCounts(),nafpriority);

		if(FluxTOF->GetUnfoldingFactor()&&FluxNaF->GetUnfoldingFactor()&&FluxAgl->GetUnfoldingFactor())
			unfolding = Global.MergeSubDResult_P(FluxTOF->GetUnfoldingFactor(),FluxNaF->GetUnfoldingFactor(),FluxAgl->GetUnfoldingFactor(),nafpriority);

		if(FluxTOF->GetRooUnfoldingFactor()&&FluxNaF->GetRooUnfoldingFactor()&&FluxAgl->GetRooUnfoldingFactor())
			roounfolding = Global.MergeSubDResult_P(FluxTOF->GetRooUnfoldingFactor(),FluxNaF->GetRooUnfoldingFactor(),FluxAgl->GetRooUnfoldingFactor(),nafpriority);

		if(FluxTOF->GetRooUnfoldingFactor_Raw()&&FluxNaF->GetRooUnfoldingFactor_Raw()&&FluxAgl->GetRooUnfoldingFactor_Raw())
			roounfolding_raw = Global.MergeSubDResult_P(FluxTOF->GetRooUnfoldingFactor_Raw(),FluxNaF->GetRooUnfoldingFactor_Raw(),FluxAgl->GetRooUnfoldingFactor_Raw(),nafpriority);

		if(FluxTOF->GetStatError()&&FluxNaF->GetStatError()&&FluxAgl->GetStatError())
			statErr = Global.MergeSubDResult_P(FluxTOF->GetStatError(),FluxNaF->GetStatError(),FluxAgl->GetStatError(),nafpriority);

		if(FluxTOF->GetSystError_unb()&&FluxNaF->GetSystError_unb()&&FluxAgl->GetSystError_unb())
			systErr= Global.MergeSubDResult_P(FluxTOF->GetSystError_unb(),FluxNaF->GetSystError_unb(),FluxAgl->GetSystError_unb(),nafpriority);

		if(FluxTOF->GetAccError()&&FluxNaF->GetAccError()&&FluxAgl->GetAccError())
			accErr= Global.MergeSubDResult_P(FluxTOF->GetAccError(),FluxNaF->GetAccError(),FluxAgl->GetAccError(),nafpriority);

		if(FluxTOF->GetUnfError()&&FluxNaF->GetUnfError()&&FluxAgl->GetUnfError())
			unfErr= Global.MergeSubDResult_P(FluxTOF->GetUnfError(),FluxNaF->GetUnfError(),FluxAgl->GetUnfError(),nafpriority);

		if(FluxTOF->GetRooUnfError()&&FluxNaF->GetRooUnfError()&&FluxAgl->GetRooUnfError())
			roounfErr= Global.MergeSubDResult_P(FluxTOF->GetRooUnfError(),FluxNaF->GetRooUnfError(),FluxAgl->GetRooUnfError(),nafpriority);


		if(	effAcc		) 	effAcc	      		= ConvertBinnedHisto( effAcc	          , (Name+"_EffAcc").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	effAccMC	) 	effAccMC      		= ConvertBinnedHisto( effAccMC          , (Name+"_EffAccMC").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	effAccMC_raw	) 	effAccMC_raw      	= ConvertBinnedHisto( effAccMC_raw       , (Name+"_EffAccMC_raw").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	flux		) 	flux	        	= ConvertBinnedHisto( flux	          , (Name+"").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	flux_ekin	) 	flux_ekin        	= ConvertBinnedHisto( flux_ekin           , (Name+"_ekin").c_str()     ,Global.GetGlobalPBins(),true); 
		if(	flux_ekin_unf	) 	flux_ekin_unf        	= ConvertBinnedHisto( flux_ekin_unf       , (Name+"_ekin_unf").c_str()     ,Global.GetGlobalPBins(),true); 
		if( 	flux_stat       )  	flux_stat     		= ConvertBinnedHisto( flux_stat         , (Name+"_rig_stat").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	flux_unf	) 	flux_unf      		= ConvertBinnedHisto( flux_unf          , (Name+"_unf").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	flux_unf_stat	) 	flux_unf_stat 		= ConvertBinnedHisto( flux_unf_stat     , (Name+"_unf_stat").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	counts		) 	counts	      		= ConvertBinnedHisto( counts	          , (Name+"_Counts").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	unfolding	) 	unfolding      		= ConvertBinnedHisto( unfolding	          , (Name+"_unfolding").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	roounfolding	) 	roounfolding   		= ConvertBinnedHisto( roounfolding	  , (Name+"_roounfolding").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	roounfolding_raw) 	roounfolding_raw 	= ConvertBinnedHisto( roounfolding_raw	  , (Name+"_roounfolding_raw").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	statErr		) 	statErr	      		= ConvertBinnedHisto( statErr	          , (Name+"_stat").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	systErr		) 	systErr	      		= ConvertBinnedHisto( systErr	          , (Name+"_syste").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	accErr		) 	accErr	      		= ConvertBinnedHisto( accErr	          , (Name+"_acce").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	unfErr	        )  	unfErr	      		= ConvertBinnedHisto( unfErr	          , (Name+"_unfe").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	roounfErr	)  	roounfErr	 	= ConvertBinnedHisto( roounfErr	          , (Name+"_roounfe").c_str()  ,Global.GetGlobalPBins(),false); 




	}	

	else{
		if(FluxTOF->GetEffAcceptance()&&FluxNaF->GetEffAcceptance()&&FluxAgl->GetEffAcceptance())
			effAcc = Global.MergeSubDResult_D(FluxTOF->GetEffAcceptance(),FluxNaF->GetEffAcceptance(),FluxAgl->GetEffAcceptance(),nafpriority);

		if(FluxTOF->GetMCAcceptance()&&FluxNaF->GetMCAcceptance()&&FluxAgl->GetMCAcceptance())
			effAccMC = Global.MergeSubDResult_D(FluxTOF->GetMCAcceptance(),FluxNaF->GetMCAcceptance(),FluxAgl->GetMCAcceptance(),nafpriority);

		if(FluxTOF->GetMCAcceptance_raw()&&FluxNaF->GetMCAcceptance_raw()&&FluxAgl->GetMCAcceptance_raw())
			effAccMC_raw = Global.MergeSubDResult_D(FluxTOF->GetMCAcceptance_raw(),FluxNaF->GetMCAcceptance_raw(),FluxAgl->GetMCAcceptance_raw(),nafpriority);


		if(FluxTOF->GetFlux()&&FluxNaF->GetFlux()&&FluxAgl->GetFlux())
			flux_ekin = Global.MergeSubDResult_D(FluxTOF->GetFlux(),FluxNaF->GetFlux(),FluxAgl->GetFlux(),nafpriority);

		if(FluxTOF->GetFlux_ekin_unf()&&FluxNaF->GetFlux_ekin_unf()&&FluxAgl->GetFlux_ekin_unf())
			flux_ekin_unf = Global.MergeSubDResult_D(FluxTOF->GetFlux_ekin_unf(),FluxNaF->GetFlux_ekin_unf(),FluxAgl->GetFlux_ekin_unf(),nafpriority);

		if(FluxTOF->GetFlux_rig()&&FluxNaF->GetFlux_rig()&&FluxAgl->GetFlux_rig())
			flux = Global.MergeSubDResult_D(FluxTOF->GetFlux_rig(),FluxNaF->GetFlux_rig(),FluxAgl->GetFlux_rig(),nafpriority);

		if(FluxTOF->GetFlux_rig_stat()&&FluxNaF->GetFlux_rig_stat()&&FluxAgl->GetFlux_rig_stat())
			flux_stat = Global.MergeSubDResult_D(FluxTOF->GetFlux_rig_stat(),FluxNaF->GetFlux_rig_stat(),FluxAgl->GetFlux_rig_stat(),nafpriority);

		if(FluxTOF->GetFlux_unf()&&FluxNaF->GetFlux_unf()&&FluxAgl->GetFlux_unf())
			flux_unf= Global.MergeSubDResult_D(FluxTOF->GetFlux_unf(),FluxNaF->GetFlux_unf(),FluxAgl->GetFlux_unf(),nafpriority);

		if(FluxTOF->GetFlux_unf_stat()&&FluxNaF->GetFlux_unf_stat()&&FluxAgl->GetFlux_unf_stat())
			flux_unf_stat= Global.MergeSubDResult_D(FluxTOF->GetFlux_unf_stat(),FluxNaF->GetFlux_unf_stat(),FluxAgl->GetFlux_unf_stat(),nafpriority);

		if(FluxTOF->GetCounts()&&FluxNaF->GetCounts()&&FluxAgl->GetCounts())
			counts = Global.MergeSubDResult_D(FluxTOF->GetCounts(),FluxNaF->GetCounts(),FluxAgl->GetCounts(),nafpriority);

		if(FluxTOF->GetUnfoldingFactor()&&FluxNaF->GetUnfoldingFactor()&&FluxAgl->GetUnfoldingFactor())
			unfolding = Global.MergeSubDResult_D(FluxTOF->GetUnfoldingFactor(),FluxNaF->GetUnfoldingFactor(),FluxAgl->GetUnfoldingFactor(),nafpriority);

		if(FluxTOF->GetRooUnfoldingFactor()&&FluxNaF->GetRooUnfoldingFactor()&&FluxAgl->GetRooUnfoldingFactor())
			roounfolding = Global.MergeSubDResult_D(FluxTOF->GetRooUnfoldingFactor(),FluxNaF->GetRooUnfoldingFactor(),FluxAgl->GetRooUnfoldingFactor(),nafpriority);

		if(FluxTOF->GetRooUnfoldingFactor_Raw()&&FluxNaF->GetRooUnfoldingFactor_Raw()&&FluxAgl->GetRooUnfoldingFactor_Raw())
			roounfolding_raw = Global.MergeSubDResult_D(FluxTOF->GetRooUnfoldingFactor_Raw(),FluxNaF->GetRooUnfoldingFactor_Raw(),FluxAgl->GetRooUnfoldingFactor_Raw(),nafpriority);


		if(FluxTOF->GetStatError()&&FluxNaF->GetStatError()&&FluxAgl->GetStatError())
			statErr = Global.MergeSubDResult_D(FluxTOF->GetStatError(),FluxNaF->GetStatError(),FluxAgl->GetStatError(),nafpriority);

		if(FluxTOF->GetSystError_unb()&&FluxNaF->GetSystError_unb()&&FluxAgl->GetSystError_unb())
			systErr= Global.MergeSubDResult_D(FluxTOF->GetSystError_unb(),FluxNaF->GetSystError_unb(),FluxAgl->GetSystError_unb(),nafpriority);

		if(FluxTOF->GetAccError()&&FluxNaF->GetAccError()&&FluxAgl->GetAccError())
			accErr= Global.MergeSubDResult_D(FluxTOF->GetAccError(),FluxNaF->GetAccError(),FluxAgl->GetAccError(),nafpriority);

		if(FluxTOF->GetUnfError()&&FluxNaF->GetUnfError()&&FluxAgl->GetUnfError())
			unfErr= Global.MergeSubDResult_D(FluxTOF->GetUnfError(),FluxNaF->GetUnfError(),FluxAgl->GetUnfError(),nafpriority);

		if(FluxTOF->GetRooUnfError()&&FluxNaF->GetRooUnfError()&&FluxAgl->GetRooUnfError())
			roounfErr= Global.MergeSubDResult_D(FluxTOF->GetRooUnfError(),FluxNaF->GetRooUnfError(),FluxAgl->GetRooUnfError(),nafpriority);


		if(	effAcc	        ) 	effAcc	           		= ConvertBinnedHisto( effAcc	          , (Name+"_EffAcc").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	effAccMC        ) 	effAccMC      		= ConvertBinnedHisto( effAccMC          , (Name+"_EffAccMC").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	effAccMC_raw        ) 	effAccMC_raw      	= ConvertBinnedHisto( effAccMC_raw          , (Name+"_EffAccMC_raw").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	flux	        ) 	flux	        	= ConvertBinnedHisto( flux	          , (Name+"").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	flux_ekin        ) 	flux_ekin        	= ConvertBinnedHisto( flux_ekin          , (Name+"_ekin").c_str()     ,Global.GetGlobalDBins(),true); 
		if(	flux_ekin_unf        ) 	flux_ekin_unf        	= ConvertBinnedHisto( flux_ekin_unf      , (Name+"_ekin_unf").c_str()     ,Global.GetGlobalDBins(),true); 
		if( 	flux_stat       )  	flux_stat     		= ConvertBinnedHisto( flux_stat         , (Name+"_rig_stat").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	flux_unf        ) 	flux_unf      		= ConvertBinnedHisto( flux_unf          , (Name+"_unf").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	flux_unf_stat	) 	flux_unf_stat 		= ConvertBinnedHisto( flux_unf_stat     , (Name+"_unf_stat").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	counts	        ) 	counts	           		= ConvertBinnedHisto( counts	          , (Name+"_Counts").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	unfolding	) 	unfolding      		= ConvertBinnedHisto( unfolding	          , (Name+"_unfolding").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	roounfolding	) 	roounfolding   		= ConvertBinnedHisto( roounfolding	  , (Name+"_roounfolding").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	roounfolding_raw) 	roounfolding_raw  	= ConvertBinnedHisto( roounfolding_raw	  , (Name+"_roounfolding_raw").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	statErr	        ) 	statErr	           	= ConvertBinnedHisto( statErr	          , (Name+"_stat").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	systErr	        ) 	systErr	           	= ConvertBinnedHisto( systErr	          , (Name+"_syste").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	accErr	        ) 	accErr	           	= ConvertBinnedHisto( accErr	          , (Name+"_acce").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	unfErr	        )  	unfErr	           	= ConvertBinnedHisto( unfErr	          , (Name+"_unfe").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	roounfErr	)  	roounfErr	 	= ConvertBinnedHisto( roounfErr	          , (Name+"_roounfe").c_str()  ,Global.GetGlobalDBins(),false); 


	}	

	if(Fluxsum&&Fluxtosubtract) {
		TGraphErrors * sumflux   = new TGraphErrors((TH1F*) Fluxsum->flux_unf->Clone());
		TGraphErrors * fluxtosub = new TGraphErrors((TH1F*) Fluxtosubtract->flux_unf->Clone());

		TH1F * totflux   = (TH1F*) forbinning->Clone((Name+"_unf").c_str());
		TH1F * Staterr   = (TH1F*) forbinning->Clone((Name+"_stat").c_str());
		TH1F * Systerr   = (TH1F*) forbinning->Clone((Name+"_syste").c_str());
		TH1F * Accerr    = (TH1F*) forbinning->Clone((Name+"_acce").c_str());
		
		totflux->Reset();
		Staterr->Reset();
		Systerr->Reset();
		Accerr ->Reset();

//		fluxtosub->Smooth(3);
//		sumflux->Smooth();

		for(int i=0;i<totflux->GetNbinsX();i++)
			if(flux_unf->GetBinLowEdge(flux_unf->FindBin(totflux->GetBinCenter(i+1)))<=flux_unf->GetBinLowEdge(flux_unf->GetNbinsX())){
				totflux->SetBinContent(i+1, flux_unf->GetBinContent(flux_unf->FindBin(totflux->GetBinCenter(i+1))));
				totflux->SetBinError(i+1, flux_unf->GetBinError(flux_unf->FindBin(totflux->GetBinCenter(i+1))));
				Staterr->SetBinContent(i+1,statErr->GetBinContent(i+1));
				Systerr->SetBinContent(i+1,systErr->GetBinContent(i+1));
				Accerr->SetBinContent(i+1,accErr->GetBinContent(i+1));
			}
			else{
				float bincenter= totflux->GetBinCenter(i+1);
		
				totflux->SetBinContent(i+1,sumflux->Eval(bincenter) - sumflux->Eval(bincenter)  );
		
				Staterr->SetBinContent(i+1,0.04);

				Systerr->SetBinContent(i+1,0.04);
				Accerr->SetBinContent(i+1,Fluxtosubtract->accErr->GetBinContent(2));
			

				totflux->SetBinError(i+1,0.05*totflux->GetBinContent(i+1));

			}
			
		flux_unf = (TH1F*) totflux->Clone((Name+"_unf").c_str());
		statErr = (TH1F*) Staterr->Clone((Name+"_stat").c_str());
		systErr = (TH1F*) Systerr->Clone((Name+"_syste").c_str());
		accErr = (TH1F*) Accerr->Clone((Name+"_acce").c_str());
	}
	
	if(FluxR){
		TH1F * fluxR_unf = (TH1F*) FluxR->GetFlux_rig()->Clone();
		TH1F * fluxR_state = (TH1F*) FluxR->GetStatError()->Clone();
		TH1F * fluxR_acce = (TH1F*) FluxR->GetAccError()->Clone();
		TH1F * fluxR_syste = (TH1F*) FluxR->GetSystError_unb()->Clone();
		TH1F * fluxR_unfe = (TH1F*) FluxR->GetRooUnfError()->Clone();





		//realligning fluxes
		if(fluxR_unf->GetBinLowEdge(1)<flux_unf->GetBinLowEdge(flux_unf->GetNbinsX())){
		float mean_new = (fluxR_unf->GetBinContent(2)+fluxR_unf->GetBinContent(3)+fluxR_unf->GetBinContent(4))/3.;
		float meanold = (flux_unf->GetBinContent(flux_unf->FindBin(fluxR_unf->GetBinCenter(2)))+
				flux_unf->GetBinContent(flux_unf->FindBin(fluxR_unf->GetBinCenter(3)))+
				flux_unf->GetBinContent(flux_unf->FindBin(fluxR_unf->GetBinCenter(4))))/3.;

		fluxR_unf->Scale(meanold/mean_new);
		}
		///
		//
		std::vector<float> newEdges;
		for (int i=0;i<flux->GetNbinsX();i++) newEdges.push_back(flux->GetBinLowEdge(i+1));
		for (int i=0;i<fluxR_unf->GetNbinsX()+1;i++) 
			if(fluxR_unf->GetBinLowEdge(i+1)>flux->GetBinLowEdge(flux->GetNbinsX())
				&& fluxR_unf->GetBinLowEdge(i+1) > fluxR_unf->GetBinLowEdge(i) //protection against strange bins
				) 

				newEdges.push_back(fluxR_unf->GetBinLowEdge(i+1));

		

		TH1F * newflux = new TH1F((Name+"_unf").c_str(),(Name+"_unf").c_str(),newEdges.size()-1,newEdges.data());
		TH1F * newflux_stat = new TH1F((Name+"_unf_stat").c_str(),(Name+"_unf_stat").c_str(),newEdges.size()-1,newEdges.data());
		TH1F * newsyst = new TH1F((Name+"_syste").c_str(),(Name+"_syste").c_str(),newEdges.size()-1,newEdges.data());
		TH1F * newacce = new TH1F((Name+"_acce").c_str(),(Name+"_acce").c_str(),newEdges.size()-1,newEdges.data());
		TH1F * newstate = new TH1F((Name+"_state").c_str(),(Name+"_state").c_str(),newEdges.size()-1,newEdges.data());
		TH1F * newunfe = new TH1F((Name+"_roounfe").c_str(),(Name+"_roounfe").c_str(),newEdges.size()-1,newEdges.data());
		for (int i=0;i<newflux->GetNbinsX();i++) {
			if(newflux->GetBinLowEdge(i+1)<fluxR_unf->GetBinLowEdge(1) ){
				newflux->SetBinContent(i+1,flux_unf->GetBinContent(i+1));
				newflux->SetBinError(i+1,flux_unf->GetBinError(i+1));
				newflux_stat->SetBinContent(i+1,flux_unf_stat->GetBinContent(i+1));
				newflux_stat->SetBinError(i+1,flux_unf_stat->GetBinError(i+1));
				newsyst->SetBinContent(i+1,systErr->GetBinContent(i+1));
				newsyst->SetBinError(i+1,systErr->GetBinError(i+1));
				newacce->SetBinContent(i+1,accErr->GetBinContent(i+1));
						newacce->SetBinError(i+1,accErr->GetBinError(i+1));
				newstate->SetBinContent(i+1,statErr->GetBinContent(i+1));
				newstate->SetBinError(i+1,statErr->GetBinError(i+1));
				newunfe->SetBinContent(i+1,roounfErr->GetBinContent(i+1));
				newunfe->SetBinError(i+1,roounfErr->GetBinError(i+1));
				}
			else if(newflux->GetBinLowEdge(i+1)<=flux->GetBinLowEdge(flux->GetNbinsX())     ){
				newflux->SetBinContent(i+1,((is_ave)*fluxR_unf->GetBinContent( fluxR_unf->FindBin(newflux->GetBinCenter(i+1)))+flux_unf->GetBinContent(i+1))/(1+is_ave)   );
				newflux->SetBinError(i+1,flux_unf->GetBinError(i+1));
				newflux_stat->SetBinContent(i+1,((is_ave)*fluxR_unf->GetBinContent(fluxR_unf->FindBin(newflux_stat->GetBinCenter(i+1)))+flux_unf->GetBinContent(i+1))/(1+is_ave)   );
				newflux_stat->SetBinError(i+1,flux_unf_stat->GetBinError(i+1));
				newsyst->SetBinContent(i+1,systErr->GetBinContent(i+1));
				newsyst->SetBinError(i+1,systErr->GetBinError(i+1));
				newacce->SetBinContent(i+1,accErr->GetBinContent(i+1));
				newacce->SetBinError(i+1,accErr->GetBinError(i+1));
				newstate->SetBinContent(i+1,statErr->GetBinContent(i+1));
				newstate->SetBinError(i+1,statErr->GetBinError(i+1));
				newunfe->SetBinContent(i+1,roounfErr->GetBinContent(i+1));
				newunfe->SetBinError(i+1,roounfErr->GetBinError(i+1));
			
			}
			else {
				newflux->SetBinContent(i+1,fluxR_unf->GetBinContent(fluxR_unf->FindBin(newflux->GetBinCenter(i+1))));
				newflux->SetBinError(i+1,fluxR_unf->GetBinError(fluxR_unf->FindBin(newflux->GetBinCenter(i+1)))+0.02*fluxR_unf->GetBinContent(fluxR_unf->FindBin(newflux_stat->GetBinCenter(i+1))));
				newflux_stat->SetBinContent(i+1,fluxR_unf->GetBinContent(fluxR_unf->FindBin(newflux_stat->GetBinCenter(i+1))));
				newflux_stat->SetBinError(i+1,fluxR_unf->GetBinError(fluxR_unf->FindBin(newflux_stat->GetBinCenter(i+1))));
				newsyst->SetBinContent(i+1,fluxR_syste->GetBinContent(fluxR_unf->FindBin(newflux_stat->GetBinCenter(i+1))));
				newsyst->SetBinError(i+1,0);
				newacce->SetBinContent(i+1,fluxR_acce->GetBinContent(fluxR_unf->FindBin(newflux_stat->GetBinCenter(i+1))));
				newacce->SetBinError(i+1,0);
				newstate->SetBinContent(i+1,fluxR_state->GetBinContent(fluxR_unf->FindBin(newflux_stat->GetBinCenter(i+1))));
				newstate->SetBinError(i+1,0);
				newunfe->SetBinContent(i+1,fluxR_unfe->GetBinContent(fluxR_unf->FindBin(newflux_stat->GetBinCenter(i+1))));
				newunfe->SetBinError(i+1,0);
			}
		}
	
		flux_unf = (TH1F*) newflux->Clone((Name+"_unf").c_str());
		flux_unf_stat = (TH1F*) newflux_stat->Clone((Name+"_unf_stat").c_str());
		statErr = (TH1F*) newstate->Clone((Name+"_stat").c_str());
		systErr = (TH1F*) newsyst->Clone((Name+"_syste").c_str());
		accErr = (TH1F*)  newacce->Clone((Name+"_acce").c_str());
		roounfErr= (TH1F*) newunfe->Clone((Name+"_roounfe").c_str());

		}



		if(	flux	        )  flux27 = Twentisevenify(flux);
                if( 	flux_stat       )  flux_stat27 = Twentisevenify(flux_stat);
                if(	flux_unf        )  flux_unf27 = Twentisevenify(flux_unf);
                if(	flux_unf_stat	)  flux_unf_stat27 = Twentisevenify(flux_unf_stat);

	return;
}




void ResultMerger::SaveResults(FileSaver finalhistos){

	if(	effAcc		)	finalhistos.Add(effAcc	     	); 
	if(	effAccMC	)	finalhistos.Add(effAccMC     	); 
	if(	effAccMC_raw	)	finalhistos.Add(effAccMC_raw   	); 
	
	if(	flux		)	finalhistos.Add(flux	     	); 
	if(	flux_ekin	)	finalhistos.Add(flux_ekin     	); 
	if(	flux_ekin_unf	)	finalhistos.Add(flux_ekin_unf  	); 
	if( 	flux_stat       )	finalhistos.Add(flux_stat    	); 
	if(	flux_unf	)	finalhistos.Add(flux_unf     	); 
	if(	flux_unf_stat	)	finalhistos.Add(flux_unf_stat	); 
	if(	flux27		)	finalhistos.Add(flux27	     	); 
	if( 	flux_stat27       )	finalhistos.Add(flux_stat27    	); 
	if(	flux_unf27	)	finalhistos.Add(flux_unf27     	); 
	if(	flux_unf_stat27	)	finalhistos.Add(flux_unf_stat27	); 
	if(	counts		)	finalhistos.Add(counts	     	); 
	if(	unfolding	)	finalhistos.Add(unfolding     	); 
	if(	roounfolding	)	finalhistos.Add(roounfolding  	); 
	if(	roounfolding_raw)	finalhistos.Add(roounfolding_raw); 
	if(	statErr		)	finalhistos.Add(statErr	     	); 
	if(	systErr		)	finalhistos.Add(systErr	     	); 
	if(	accErr		)	finalhistos.Add(accErr	     	); 
	if(	unfErr	        )	finalhistos.Add(unfErr	     	); 
	if(	roounfErr	)	finalhistos.Add(roounfErr	); 



	finalhistos.writeObjsInFolder(("MergedResults/"+Name).c_str());

}




void ResultMerger::ResidualCorrectionWithFluxes(TH1F * fluxekin1,TH1F * fluxekin2){

	for(int i=0;i<flux_ekin->GetNbinsX();i++){

		float ekin = flux_ekin->GetBinCenter(i+1);
		float ekinbin = fluxekin2->FindBin(ekin);
		float bin_this = flux_ekin->FindBin(ekin);
		float corr ;
		if(fluxekin1->GetBinContent(ekinbin)>0){
			corr= fluxekin2->GetBinContent(ekinbin)/fluxekin1->GetBinContent(ekinbin);
			flux_stat->SetBinContent(bin_this, flux_stat->GetBinContent(bin_this)*corr);
			flux->SetBinContent(bin_this, flux->GetBinContent(bin_this)*corr);
			flux_unf->SetBinContent(bin_this, flux_unf->GetBinContent(bin_this)*corr);
			flux_unf_stat->SetBinContent(bin_this, flux_unf_stat->GetBinContent(bin_this)*corr);
			flux27->SetBinContent(bin_this, flux27->GetBinContent(bin_this)*corr);
			flux_stat27->SetBinContent(bin_this, flux_stat27->GetBinContent(bin_this)*corr);
			flux_unf27->SetBinContent(bin_this, flux_unf27->GetBinContent(bin_this)*corr);
			flux_unf_stat27->SetBinContent(bin_this, flux_unf_stat27->GetBinContent(bin_this)*corr);
			effAcc->SetBinContent(bin_this, effAcc->GetBinContent(bin_this)*corr);
	
			flux_ekin->SetBinContent(i+1,flux->GetBinContent(i+1)
				/(flux->GetBinLowEdge(i+2)-flux->GetBinLowEdge(i+1)) 
				* (flux_ekin->GetBinLowEdge(i+2)-flux_ekin->GetBinLowEdge(i+1))  );
			flux_ekin_unf->SetBinContent(i+1,flux_unf->GetBinContent(i+1)
				/(flux_unf->GetBinLowEdge(i+2)-flux_unf->GetBinLowEdge(i+1)) 
				* (flux_ekin_unf->GetBinLowEdge(i+2)-flux_ekin_unf->GetBinLowEdge(i+1))  );
			

		}

	}
}
	

