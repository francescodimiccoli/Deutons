class MCMCEntry {
    std::string suffix;
    TH3D * MCMCPTemplates;
    TH3D * MCMCDTemplates;
    TH3D * MCMCHeTemplates;
    TH2D * MCMCData;
public:

    MCMCEntry( const std::string & suffix,
               int nBm, double minBm, double maxBm,
               int nRm, double minRm, double maxRm,
               int nRt, double minRt, double maxRt ) 
    {
        std::string nP = "MCMCPTemplates"  + suffix;
        std::string nD = "MCMCDTemplates"  + suffix;
        std::string nH = "MCMCHeTemplates" + suffix;
        std::string nT = "MCMCData" + suffix;

        MCMCPTemplates  = new TH3D(nP.c_str(), nP.c_str(), nBm, minBm, maxBm, nRm, minRm, maxRm, nRt, minRt, maxRt);
        MCMCDTemplates  = new TH3D(nD.c_str(), nD.c_str(), nBm, minBm, maxBm, nRm, minRm, maxRm, nRt, minRt, maxRt);
        MCMCHeTemplates = new TH3D(nH.c_str(), nH.c_str(), nBm, minBm, maxBm, nRm, minRm, maxRm, nRt, minRt, maxRt);
        MCMCData        = new TH2D(nT.c_str(), nT.c_str(), nBm, minBm, maxBm, nRm, minRm, maxRm );
    }

    void FillMC(double Bm, double Rm, double Rt)
    {
        if( Massa_gen < 1 && Massa_gen > 0.5) MCMCPTemplates->Fill(Bm, Rm, Rt);
        if( Massa_gen < 2 && Massa_gen > 1.5) MCMCDTemplates->Fill(Bm, Rm, Rt);
        if( Massa_gen < 4 && Massa_gen > 2.5) MCMCHeTemplates->Fill(Bm, Rm, Rt);
    }
    void FillData(double Bm , double Rm) { MCMCData->Fill(Bm, Rm); }

    void Write()
    {
        MCMCPTemplates->Write();
        MCMCDTemplates->Write(); 
        MCMCHeTemplates->Write();
        MCMCData->Write();
    }

};
                                        /*| B measured   | R measured     | R true    |*/
MCMCEntry * MCMC_TOF = new MCMCEntry("TOF", 50, 0.3,  1.2, 200, 0.0, 100.0, 300, 0, 30);
MCMCEntry * MCMC_NaF = new MCMCEntry("NaF", 50, 0.75, 1.2, 200, 0.0, 100.0, 300, 0, 30);
MCMCEntry * MCMC_Agl = new MCMCEntry("Agl", 50, 0.95, 1.1, 200, 0.0, 100.0, 300, 0, 30);

void MCMC_Fill(TNtuple *ntupla, int l)
{
    int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	if(Likcut&&Distcut){
		MCMC_TOF->FillMC( Beta, R, Momento_gen);
        if((((int)Cutmask)>>11)==512) MCMC_NaF->FillMC( BetaRICH, R, Momento_gen);
        if((((int)Cutmask)>>11)==0  ) MCMC_Agl->FillMC( BetaRICH, R, Momento_gen);
	}
}

void MCMCDATA_Fill(TNtuple *ntupla, int l){
    int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	if( Likcut && Distcut && R>1.2*Rcutoff)
    {
		MCMC_TOF->FillData( Beta, R );
        if((((int)Cutmask)>>11)==512) MCMC_NaF->FillData( BetaRICH, R );
        if((((int)Cutmask)>>11)==0  ) MCMC_Agl->FillData( BetaRICH, R );
    }
}

void MCMC_Write(){
    MCMC_TOF->Write();
    MCMC_NaF->Write();
    MCMC_Agl->Write();
}

void MCMC(TFile * file){
    std::cout << "That function should not be called.\n";
}
