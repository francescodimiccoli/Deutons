TH3D *MCMCPTemplatesTOF = new TH3D("MCMCPTemplatesTOF","MCMCPTemplatesTOF",nbinsbeta,0,nbinsbeta,500,0,1.3,500,0,30);
TH3D *MCMCPTemplatesNaF = new TH3D("MCMCPTemplatesNaF","MCMCPTemplatesNaF",nbinsbeta,0,nbinsbeta,500,0,1.3,500,0,30);
TH3D *MCMCPTemplatesAgl = new TH3D("MCMCPTemplatesAgl","MCMCPTemplatesAgl",nbinsbeta,0,nbinsbeta,500,0.8,1.3,500,0,30);

TH3D *MCMCDTemplatesTOF = new TH3D("MCMCDTemplatesTOF","MCMCDTemplatesTOF",nbinsbeta,0,nbinsbeta,500,0,1.3,500,0,30);
TH3D *MCMCDTemplatesNaF = new TH3D("MCMCDTemplatesNaF","MCMCDTemplatesNaF",nbinsbeta,0,nbinsbeta,500,0,1.3,500,0,30);
TH3D *MCMCDTemplatesAgl = new TH3D("MCMCDTemplatesAgl","MCMCDTemplatesAgl",nbinsbeta,0,nbinsbeta,500,0.8,1.3,500,0,30);

TH3D *MCMCHeTemplatesTOF = new TH3D("MCMCHeTemplatesTOF","MCMCHeTemplatesTOF",nbinsbeta,0,nbinsbeta,500,0,1.3,500,0,30);
TH3D *MCMCHeTemplatesNaF = new TH3D("MCMCHeTemplatesNaF","MCMCHeTemplatesNaF",nbinsbeta,0,nbinsbeta,500,0,1.3,500,0,30);
TH3D *MCMCHeTemplatesAgl = new TH3D("MCMCHeTemplatesAgl","MCMCHeTemplatesAgl",nbinsbeta,0,nbinsbeta,500,0.8,1.3,500,0,30);

TH2D *MCMCDataTOF = new TH2D("MCMCDataTOF","MCMCDataTOF",nbinsbeta,0,nbinsbeta,500,0,1.3);
TH2D *MCMCDataNaF = new TH2D("MCMCDataNaF","MCMCDataNaF",nbinsbeta,0,nbinsbeta,500,0,1.3);
TH2D *MCMCDataAgl = new TH2D("MCMCDataAgl","MCMCDataAgl",nbinsbeta,0,nbinsbeta,500,0.8,1.3);

void MCMC_Fill(TNtuple *ntupla, int l){
        int k = ntupla->GetEvent(l);
	
	if(Beta<=0||R<=0) return;
        float mm=0;
	if(Likcut&&Distcut){
		if(Massa_gen<1&&Massa_gen>0.5) 
			for(int m=0;m<nbinsbeta;m++)  if(Var>BetaD[m]&&Var<=BetaD[m+1]) MCMCPTemplatesTOF->Fill(m,Beta,Momento_gen);
		if(Massa_gen<2&&Massa_gen>1.5)
                        for(int m=0;m<nbinsbeta;m++)  if(Var>BetaD[m]&&Var<=BetaD[m+1]) MCMCDTemplatesTOF->Fill(m,Beta,Momento_gen);	
		 if(Massa_gen<4&&Massa_gen>2.5)
                        for(int m=0;m<nbinsbeta;m++)  if(Var>BetaD[m]&&Var<=BetaD[m+1]) MCMCHeTemplatesTOF->Fill(m,Beta,Momento_gen);
	}
	if(Likcut&&Distcut){
		if(Massa_gen<1&&Massa_gen>0.5)
			for(int m=0;m<nbinsbeta;m++) if((((int)Cutmask)>>11)==512&&Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) MCMCPTemplatesNaF->Fill(m,BetaRICH,Momento_gen);
		if(Massa_gen<2&&Massa_gen>1.5)
			for(int m=0;m<nbinsbeta;m++) if((((int)Cutmask)>>11)==512&&Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) MCMCDTemplatesNaF->Fill(m,BetaRICH,Momento_gen);
		if(Massa_gen<4&&Massa_gen>2.5)
			for(int m=0;m<nbinsbeta;m++) if((((int)Cutmask)>>11)==512&&Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) MCMCHeTemplatesNaF->Fill(m,BetaRICH,Momento_gen);
	}
	if(Likcut&&Distcut){
		if(Massa_gen<1&&Massa_gen>0.5)
			for(int m=0;m<nbinsbeta;m++) if((((int)Cutmask)>>11)==0&&Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) MCMCPTemplatesAgl->Fill(m,BetaRICH,Momento_gen);
		if(Massa_gen<2&&Massa_gen>1.5)
			for(int m=0;m<nbinsbeta;m++) if((((int)Cutmask)>>11)==0&&Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) MCMCDTemplatesAgl->Fill(m,BetaRICH,Momento_gen);
		if(Massa_gen<4&&Massa_gen>2.5)
			for(int m=0;m<nbinsbeta;m++) if((((int)Cutmask)>>11)==0&&Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) MCMCHeTemplatesAgl->Fill(m,BetaRICH,Momento_gen);
	}

	
        return;
}

void MCMCDATA_Fill(TNtuple *ntupla, int l){
        int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	
	float mm=0;
	if(Likcut&&Distcut)
		 for(int m=0;m<nbinsbeta;m++)  if(Var>BetaD[m]&&Var<=BetaD[m+1]) if(R>1.2*Rcutoff)  MCMCDataTOF->Fill(m,Beta);
						
					
	if(Likcut&&Distcut)
		for(int m=0;m<nbinsbeta;m++) 
				if((((int)Cutmask)>>11)==512) if(Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) if(R>1.2*Rcutoff) MCMCDataNaF->Fill(m,BetaRICH);

	if(Likcut&&Distcut)
                for(int m=0;m<nbinsbeta;m++) 
                                if((((int)Cutmask)>>11)==0) if(Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) if(R>1.2*Rcutoff) MCMCDataAgl->Fill(m,BetaRICH);
}

void MCMC_Write(){
        MCMCPTemplatesTOF->Write();
        MCMCDTemplatesTOF->Write();
        MCMCHeTemplatesTOF->Write();
        MCMCPTemplatesNaF->Write();
        MCMCDTemplatesNaF->Write();
        MCMCHeTemplatesNaF->Write();
        MCMCPTemplatesAgl->Write();
        MCMCDTemplatesAgl->Write();
        MCMCHeTemplatesAgl->Write();
        MCMCDataTOF->Write();
        MCMCDataNaF->Write();
        MCMCDataAgl->Write();
}

void MCMC(TFile * file){
	TH3D * MCMCPTemplatesTOF=(TH3D*)file->Get("MCMCPTemplatesTOF");
	TH3D * MCMCDTemplatesTOF=(TH3D*)file->Get("MCMCDTemplatesTOF");
	TH3D * MCMCHeTemplatesTOF=(TH3D*)file->Get("MCMCHeTemplatesTOF");
	TH3D * MCMCPTemplatesNaF=(TH3D*)file->Get("MCMCPTemplatesNaF");
	TH3D * MCMCDTemplatesNaF=(TH3D*)file->Get("MCMCDTemplatesNaF");
	TH3D * MCMCHeTemplatesNaF=(TH3D*)file->Get("MCMCHeTemplatesNaF");
	TH3D * MCMCPTemplatesAgl=(TH3D*)file->Get("MCMCPTemplatesAgl");
	TH3D * MCMCDTemplatesAgl=(TH3D*)file->Get("MCMCDTemplatesAgl");
	TH3D * MCMCHeTemplatesAgl=(TH3D*)file->Get("MCMCHeTemplatesAgl");
	TH2D * MCMCDataTOF=(TH2D*)file->Get("MCMCDataTOF");
	TH2D * MCMCDataNaF=(TH2D*)file->Get("MCMCDataNaF");
	TH2D * MCMCDataAgl=(TH2D*)file->Get("MCMCDataAgl");
}

