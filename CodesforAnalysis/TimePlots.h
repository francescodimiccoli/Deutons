	TSpline3 *Corr_L1[40];
	TSpline3 *Corr_TOFU[40];
	TSpline3 *Corr_Track[40];
	TSpline3 *Corr_TOFD[40];
	TGraphErrors *CorrLATpre_Spl[3][40];
	TGraphErrors *CorrLAT_DistTOF_Spl[40];
	TGraphErrors *CorrLAT_LikTOF_Spl[40];	
	TGraphErrors *CorrLAT_DistNaF_Spl[40];
        TGraphErrors *CorrLAT_LikNaF_Spl[40];
	TGraphErrors *CorrLAT_DistAgl_Spl[40];
        TGraphErrors *CorrLAT_LikAgl_Spl[40];
	TGraphErrors *preDVSMC_P[3][40];
	TGraphErrors *LikDVSMC_P[40];
        TGraphErrors *DistDVSMC_P[40];

	TGraphErrors *preDVSMC_PTOF[3][40];
	TGraphErrors *LikDVSMC_PTOF[40];
        TGraphErrors *DistDVSMC_PTOF[40];
		
	TGraphErrors *preDVSMC_PNaF[3][40];
	TGraphErrors *LikDVSMC_PNaF[40];
        TGraphErrors *DistDVSMC_PNaF[40];
		
	TGraphErrors *preDVSMC_PAgl[3][40];
	TGraphErrors *LikDVSMC_PAgl[40];
        TGraphErrors *DistDVSMC_PAgl[40];
	
	TH1F * TrackerEff[40];
	TH1F * TriggerEff[40];
	TGraphErrors *P_Fluxes[40];
	TGraphErrors *P_Fluxesratio[40];
	TGraphErrors *D_FluxesTOF[40];
        TGraphErrors *D_FluxesratioTOF[40];
	TGraphErrors *D_FluxesNaF[40];
        TGraphErrors *D_FluxesratioNaF[40];
	TGraphErrors *D_FluxesAgl[40];
        TGraphErrors *D_FluxesratioAgl[40];	

	TGraphErrors *PD_FluxesTOF[40];
	TGraphErrors *PD_FluxesNaF[40];
	TGraphErrors *PD_FluxesAgl[40];

