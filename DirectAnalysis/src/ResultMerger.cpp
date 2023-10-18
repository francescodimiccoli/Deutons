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



ResultMerger::ResultMerger(FileSaver finalhistos, std::string name, RangeMerger Global, Flux * FluxTOF, Flux * FluxNaF, Flux * FluxAgl, Particle particle, bool nafpriority, Flux * Fluxsum, Flux * Fluxtosubtract,TH1F * forbinning){

	Name=name;
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

		if(FluxTOF->GetSystError()&&FluxNaF->GetSystError()&&FluxAgl->GetSystError())
			systErr= Global.MergeSubDResult_P(FluxTOF->GetSystError(),FluxNaF->GetSystError(),FluxAgl->GetSystError(),nafpriority);

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

		if(FluxTOF->GetSystError()&&FluxNaF->GetSystError()&&FluxAgl->GetSystError())
			systErr= Global.MergeSubDResult_D(FluxTOF->GetSystError(),FluxNaF->GetSystError(),FluxAgl->GetSystError(),nafpriority);

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
		TH1F * sumflux   = (TH1F*) Fluxsum->GetFlux_rig()->Clone();
		TH1F * fluxtosub = (TH1F*) Fluxtosubtract->GetFlux_unf()->Clone();

		TH1F * totflux   = (TH1F*) forbinning->Clone((Name+"_unf").c_str());
		TH1F * Staterr   = (TH1F*) forbinning->Clone((Name+"_stat").c_str());
		TH1F * Systerr   = (TH1F*) forbinning->Clone((Name+"_syste").c_str());
		TH1F * Accerr    = (TH1F*) forbinning->Clone((Name+"_acce").c_str());
		
		totflux->Reset();
		Staterr->Reset();
		Systerr->Reset();
		Accerr ->Reset();

		fluxtosub->Smooth(3);
		sumflux->Smooth();

		for(int i=0;i<totflux->GetNbinsX();i++)
			if(flux_unf->FindBin(totflux->GetBinCenter(i+1))<=flux_unf->GetNbinsX()){
				totflux->SetBinContent(i+1, flux_unf->GetBinContent(flux_unf->FindBin(totflux->GetBinCenter(i+1))));
				totflux->SetBinError(i+1, flux_unf->GetBinError(flux_unf->FindBin(totflux->GetBinCenter(i+1))));
				Staterr->SetBinContent(i+1,statErr->GetBinContent(i+1));
				Systerr->SetBinContent(i+1,systErr->GetBinContent(i+1));
				Accerr->SetBinContent(i+1,accErr->GetBinContent(i+1));
			}
			else{
				float bincenter= totflux->GetBinCenter(i+1);
		
				totflux->SetBinContent(i+1,sumflux->GetBinContent(sumflux->FindBin(bincenter)) - fluxtosub->GetBinContent(fluxtosub->FindBin(bincenter)));
		
				Staterr->SetBinContent(i+1,sqrt(pow(Fluxtosubtract->GetStatError()->GetBinContent(sumflux->FindBin(bincenter))*fluxtosub->GetBinContent(sumflux->FindBin(bincenter)),2)+ 
						     	   pow(Fluxsum->GetStatError()->GetBinContent(sumflux->FindBin(bincenter))*sumflux->GetBinContent(sumflux->FindBin(bincenter)),2))
							   /totflux->GetBinContent(sumflux->FindBin(bincenter)));

				Systerr->SetBinContent(i+1,Fluxtosubtract->GetSystError()->GetBinContent(sumflux->FindBin(bincenter))+0.05);
				Accerr->SetBinContent(i+1,Fluxtosubtract->GetAccError()->GetBinContent(sumflux->FindBin(bincenter)));
			
				float error1 = (fluxtosub->GetBinError(fluxtosub->FindBin(bincenter))/fluxtosub->GetBinContent(fluxtosub->FindBin(bincenter)));
				float error2 = (sumflux->GetBinError(sumflux->FindBin(bincenter))/sumflux->GetBinContent(sumflux->FindBin(bincenter)));
				float error3 = Systerr->GetBinContent(i+1);


				totflux->SetBinError(i+1,sqrt(pow(error1,2)+pow(error2,2)+pow(error3,3))*totflux->GetBinContent(i+1));

			}
			
		flux_unf = (TH1F*) totflux->Clone((Name+"_unf").c_str());
		statErr = (TH1F*) Staterr->Clone((Name+"_stat").c_str());
		systErr = (TH1F*) Systerr->Clone((Name+"_syste").c_str());
		accErr = (TH1F*) Accerr->Clone((Name+"_acce").c_str());
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






