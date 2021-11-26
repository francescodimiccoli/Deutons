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



ResultMerger::ResultMerger(FileSaver finalhistos, std::string name, RangeMerger Global, Flux * FluxTOF, Flux * FluxNaF, Flux * FluxAgl, Particle particle){

	Name=name;
	if(!FluxTOF || !FluxNaF || !FluxAgl) return;


	if(particle.getA()/particle.getZ() <=1.5){

		if(FluxTOF->GetFlux()&&FluxNaF->GetFlux()&&FluxAgl->GetFlux())
			flux_ekin = Global.MergeSubDResult_P(FluxTOF->GetFlux(),FluxNaF->GetFlux(),FluxAgl->GetFlux());

		if(FluxTOF->GetFlux_ekin_unf()&&FluxNaF->GetFlux_ekin_unf()&&FluxAgl->GetFlux_ekin_unf())
			flux_ekin_unf = Global.MergeSubDResult_P(FluxTOF->GetFlux_ekin_unf(),FluxNaF->GetFlux_ekin_unf(),FluxAgl->GetFlux_ekin_unf());

		if(FluxTOF->GetFlux_rig()&&FluxNaF->GetFlux_rig()&&FluxAgl->GetFlux_rig())
			flux = Global.MergeSubDResult_P(FluxTOF->GetFlux_rig(),FluxNaF->GetFlux_rig(),FluxAgl->GetFlux_rig());

		if(FluxTOF->GetFlux_rig_stat()&&FluxNaF->GetFlux_rig_stat()&&FluxAgl->GetFlux_rig_stat())
			flux_stat = Global.MergeSubDResult_P(FluxTOF->GetFlux_rig_stat(),FluxNaF->GetFlux_rig_stat(),FluxAgl->GetFlux_rig_stat());

		if(FluxTOF->GetFlux_unf()&&FluxNaF->GetFlux_unf()&&FluxAgl->GetFlux_unf())
			flux_unf= Global.MergeSubDResult_P(FluxTOF->GetFlux_unf(),FluxNaF->GetFlux_unf(),FluxAgl->GetFlux_unf());
 
		if(FluxTOF->GetFlux_unf_stat()&&FluxNaF->GetFlux_unf_stat()&&FluxAgl->GetFlux_unf_stat())
			flux_unf_stat= Global.MergeSubDResult_P(FluxTOF->GetFlux_unf_stat(),FluxNaF->GetFlux_unf_stat(),FluxAgl->GetFlux_unf_stat()); 

		if(FluxTOF->GetEffAcceptance()&&FluxNaF->GetEffAcceptance()&&FluxAgl->GetEffAcceptance())
			effAcc = Global.MergeSubDResult_P(FluxTOF->GetEffAcceptance(),FluxNaF->GetEffAcceptance(),FluxAgl->GetEffAcceptance());

		if(FluxTOF->GetMCAcceptance()&&FluxNaF->GetMCAcceptance()&&FluxAgl->GetMCAcceptance())
			effAccMC = Global.MergeSubDResult_P(FluxTOF->GetMCAcceptance(),FluxNaF->GetMCAcceptance(),FluxAgl->GetMCAcceptance());

		if(FluxTOF->GetCounts()&&FluxNaF->GetCounts()&&FluxAgl->GetCounts())
			counts = Global.MergeSubDResult_P(FluxTOF->GetCounts(),FluxNaF->GetCounts(),FluxAgl->GetCounts());

		if(FluxTOF->GetUnfoldingFactor()&&FluxNaF->GetUnfoldingFactor()&&FluxAgl->GetUnfoldingFactor())
			unfolding = Global.MergeSubDResult_P(FluxTOF->GetUnfoldingFactor(),FluxNaF->GetUnfoldingFactor(),FluxAgl->GetUnfoldingFactor());

		if(FluxTOF->GetStatError()&&FluxNaF->GetStatError()&&FluxAgl->GetStatError())
			statErr = Global.MergeSubDResult_P(FluxTOF->GetStatError(),FluxNaF->GetStatError(),FluxAgl->GetStatError());

		if(FluxTOF->GetSystError()&&FluxNaF->GetSystError()&&FluxAgl->GetSystError())
			systErr= Global.MergeSubDResult_P(FluxTOF->GetSystError(),FluxNaF->GetSystError(),FluxAgl->GetSystError());

		if(FluxTOF->GetAccError()&&FluxNaF->GetAccError()&&FluxAgl->GetAccError())
			accErr= Global.MergeSubDResult_P(FluxTOF->GetAccError(),FluxNaF->GetAccError(),FluxAgl->GetAccError());

		if(FluxTOF->GetUnfError()&&FluxNaF->GetUnfError()&&FluxAgl->GetUnfError())
			unfErr= Global.MergeSubDResult_P(FluxTOF->GetUnfError(),FluxNaF->GetUnfError(),FluxAgl->GetUnfError());




		if(	flux		) 	flux	        	= ConvertBinnedHisto( flux	          , (Name+"").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	flux_ekin	) 	flux_ekin        	= ConvertBinnedHisto( flux_ekin           , (Name+"_ekin").c_str()     ,Global.GetGlobalPBins(),true); 
		if(	flux_ekin_unf	) 	flux_ekin_unf        	= ConvertBinnedHisto( flux_ekin_unf       , (Name+"_ekin_unf").c_str()     ,Global.GetGlobalPBins(),true); 
		if( 	flux_stat       )  	flux_stat     		= ConvertBinnedHisto( flux_stat         , (Name+"_rig_stat").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	flux_unf	) 	flux_unf      		= ConvertBinnedHisto( flux_unf          , (Name+"_unf").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	flux_unf_stat	) 	flux_unf_stat 		= ConvertBinnedHisto( flux_unf_stat     , (Name+"_unf_stat").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	effAcc		) 	effAcc	      		= ConvertBinnedHisto( effAcc	          , (Name+"_EffAcc").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	effAccMC	) 	effAccMC      		= ConvertBinnedHisto( effAccMC          , (Name+"_EffAccMC").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	counts		) 	counts	      		= ConvertBinnedHisto( counts	          , (Name+"_Counts").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	unfolding	) 	unfolding      		= ConvertBinnedHisto( unfolding	          , (Name+"_unfolding").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	statErr		) 	statErr	      		= ConvertBinnedHisto( statErr	          , (Name+"_stat").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	systErr		) 	systErr	      		= ConvertBinnedHisto( systErr	          , (Name+"_syste").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	accErr		) 	accErr	      		= ConvertBinnedHisto( accErr	          , (Name+"_acce").c_str()     ,Global.GetGlobalPBins(),false); 
		if(	unfErr	        )  	unfErr	      		= ConvertBinnedHisto( unfErr	          , (Name+"_unfe").c_str()     ,Global.GetGlobalPBins(),false); 



	}	

	else{

		if(FluxTOF->GetFlux()&&FluxNaF->GetFlux()&&FluxAgl->GetFlux())
			flux_ekin = Global.MergeSubDResult_D(FluxTOF->GetFlux(),FluxNaF->GetFlux(),FluxAgl->GetFlux());

		if(FluxTOF->GetFlux_ekin_unf()&&FluxNaF->GetFlux_ekin_unf()&&FluxAgl->GetFlux_ekin_unf())
			flux_ekin_unf = Global.MergeSubDResult_D(FluxTOF->GetFlux_ekin_unf(),FluxNaF->GetFlux_ekin_unf(),FluxAgl->GetFlux_ekin_unf());

		if(FluxTOF->GetFlux_rig()&&FluxNaF->GetFlux_rig()&&FluxAgl->GetFlux_rig())
			flux = Global.MergeSubDResult_D(FluxTOF->GetFlux_rig(),FluxNaF->GetFlux_rig(),FluxAgl->GetFlux_rig());

		if(FluxTOF->GetFlux_rig_stat()&&FluxNaF->GetFlux_rig_stat()&&FluxAgl->GetFlux_rig_stat())
			flux_stat = Global.MergeSubDResult_D(FluxTOF->GetFlux_rig_stat(),FluxNaF->GetFlux_rig_stat(),FluxAgl->GetFlux_rig_stat());

		if(FluxTOF->GetFlux_unf()&&FluxNaF->GetFlux_unf()&&FluxAgl->GetFlux_unf())
			flux_unf= Global.MergeSubDResult_D(FluxTOF->GetFlux_unf(),FluxNaF->GetFlux_unf(),FluxAgl->GetFlux_unf());

		if(FluxTOF->GetFlux_unf_stat()&&FluxNaF->GetFlux_unf_stat()&&FluxAgl->GetFlux_unf_stat())
			flux_unf_stat= Global.MergeSubDResult_D(FluxTOF->GetFlux_unf_stat(),FluxNaF->GetFlux_unf_stat(),FluxAgl->GetFlux_unf_stat());

		if(FluxTOF->GetEffAcceptance()&&FluxNaF->GetEffAcceptance()&&FluxAgl->GetEffAcceptance())
			effAcc = Global.MergeSubDResult_D(FluxTOF->GetEffAcceptance(),FluxNaF->GetEffAcceptance(),FluxAgl->GetEffAcceptance());

		if(FluxTOF->GetMCAcceptance()&&FluxNaF->GetMCAcceptance()&&FluxAgl->GetMCAcceptance())
			effAccMC = Global.MergeSubDResult_D(FluxTOF->GetMCAcceptance(),FluxNaF->GetMCAcceptance(),FluxAgl->GetMCAcceptance());

		if(FluxTOF->GetCounts()&&FluxNaF->GetCounts()&&FluxAgl->GetCounts())
			counts = Global.MergeSubDResult_D(FluxTOF->GetCounts(),FluxNaF->GetCounts(),FluxAgl->GetCounts());

		if(FluxTOF->GetUnfoldingFactor()&&FluxNaF->GetUnfoldingFactor()&&FluxAgl->GetUnfoldingFactor())
			unfolding = Global.MergeSubDResult_D(FluxTOF->GetUnfoldingFactor(),FluxNaF->GetUnfoldingFactor(),FluxAgl->GetUnfoldingFactor());

		if(FluxTOF->GetStatError()&&FluxNaF->GetStatError()&&FluxAgl->GetStatError())
			statErr = Global.MergeSubDResult_D(FluxTOF->GetStatError(),FluxNaF->GetStatError(),FluxAgl->GetStatError());

		if(FluxTOF->GetSystError()&&FluxNaF->GetSystError()&&FluxAgl->GetSystError())
			systErr= Global.MergeSubDResult_D(FluxTOF->GetSystError(),FluxNaF->GetSystError(),FluxAgl->GetSystError());

		if(FluxTOF->GetAccError()&&FluxNaF->GetAccError()&&FluxAgl->GetAccError())
			accErr= Global.MergeSubDResult_D(FluxTOF->GetAccError(),FluxNaF->GetAccError(),FluxAgl->GetAccError());

		if(FluxTOF->GetUnfError()&&FluxNaF->GetUnfError()&&FluxAgl->GetUnfError())
			unfErr= Global.MergeSubDResult_D(FluxTOF->GetUnfError(),FluxNaF->GetUnfError(),FluxAgl->GetUnfError());



		if(	flux	        ) 	flux	        	= ConvertBinnedHisto( flux	          , (Name+"").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	flux_ekin        ) 	flux_ekin        	= ConvertBinnedHisto( flux_ekin          , (Name+"_ekin").c_str()     ,Global.GetGlobalDBins(),true); 
		if(	flux_ekin_unf        ) 	flux_ekin_unf        	= ConvertBinnedHisto( flux_ekin_unf      , (Name+"_ekin_unf").c_str()     ,Global.GetGlobalDBins(),true); 
		if( 	flux_stat       )  	flux_stat     		= ConvertBinnedHisto( flux_stat         , (Name+"_rig_stat").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	flux_unf        ) 	flux_unf      		= ConvertBinnedHisto( flux_unf          , (Name+"_unf").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	flux_unf_stat	) 	flux_unf_stat 		= ConvertBinnedHisto( flux_unf_stat     , (Name+"_unf_stat").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	effAcc	        ) 	effAcc	           		= ConvertBinnedHisto( effAcc	          , (Name+"_EffAcc").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	effAccMC        ) 	effAccMC      		= ConvertBinnedHisto( effAccMC          , (Name+"_EffAccMC").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	counts	        ) 	counts	           		= ConvertBinnedHisto( counts	          , (Name+"_Counts").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	unfolding	) 	unfolding      		= ConvertBinnedHisto( unfolding	          , (Name+"_unfolding").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	statErr	        ) 	statErr	           		= ConvertBinnedHisto( statErr	          , (Name+"_stat").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	systErr	        ) 	systErr	           		= ConvertBinnedHisto( systErr	          , (Name+"_syste").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	accErr	        ) 	accErr	           		= ConvertBinnedHisto( accErr	          , (Name+"_acce").c_str()     ,Global.GetGlobalDBins(),false); 
		if(	unfErr	        )  	unfErr	           		= ConvertBinnedHisto( unfErr	          , (Name+"_unfe").c_str()     ,Global.GetGlobalDBins(),false); 

	}	

		if(	flux	        )  flux27 = Twentisevenify(flux);
                if( 	flux_stat       )  flux_stat27 = Twentisevenify(flux_stat);
                if(	flux_unf        )  flux_unf27 = Twentisevenify(flux_unf);
                if(	flux_unf_stat	)  flux_unf_stat27 = Twentisevenify(flux_unf_stat);

	return;
}




void ResultMerger::SaveResults(FileSaver finalhistos){


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
	if(	effAcc		)	finalhistos.Add(effAcc	     	); 
	if(	effAccMC	)	finalhistos.Add(effAccMC     	); 
	if(	counts		)	finalhistos.Add(counts	     	); 
	if(	unfolding	)	finalhistos.Add(unfolding     	); 
	if(	statErr		)	finalhistos.Add(statErr	     	); 
	if(	systErr		)	finalhistos.Add(systErr	     	); 
	if(	accErr		)	finalhistos.Add(accErr	     	); 
	if(	unfErr	        )	finalhistos.Add(unfErr	     	); 


	finalhistos.writeObjsInFolder(("MergedResults/"+Name).c_str());

}






