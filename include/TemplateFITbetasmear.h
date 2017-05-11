using namespace std;





struct TFit {
   TH1F * Templ_P ;
   TH1F * Templ_D ;
   TH1F * Templ_He;
   TH1F * Data;
   TH1F * DataPrim;	
   float wheightP,wheightD,wheightHe;
   TFractionFitter *Tfit;
   int Tfit_outcome;
   
   float DCounts=0;
   float PCounts=0;
   float ChiSquare=0;			
   float StatErr=0;	

   TFit(){}	
   TFit(TH1F * templ_P, TH1F * templ_D, TH1F * data, TH1F * dataPrim) { Templ_P= templ_P; Templ_D=templ_D; Data=data; DataPrim=dataPrim; }

};


struct Systpar{
	int steps;
	float sigma;
	float shift;
};

struct BestChi {
	int i=0;
	int j=0;
	float chimin=0;
	void FindMinimum(TH2F * Histo) {
		float Best = 9999999;
		for(int x=0;x<Histo->GetNbinsX();x++)
			for(int y=0;y<Histo->GetNbinsY();y++)
				if(Histo->GetBinContent(x+1,y+1)<Best){
					Best=Histo->GetBinContent(x+1,y+1);
				        i=x;
					j=y;
					chimin=Best;	
				}
	}
};

class TemplateFIT {

	private:
	std::vector<std::vector<std::vector<TFit *>>> fits;
	std::vector<BestChi *> BestChiSquare;	
	std::vector<TH1F *> TransferFunction;


	std::vector<TH2F *> DCountsSpread;
	std::vector<TH1F *> WeightedDCounts;
	std::vector<TH2F *> TFitChisquare;

	TH1F * StatError;
	TH1F * SystError;

	TH1F * ProtonCounts;
	TH1F * DeuteronCounts;

	TH1F * BestFitSigma;
	TH1F * BestFitShift;

	Binning bins;
        std::string var;
        std::string cut;
	std::string cutprimary;
	std::string discr_var;

	std::string basename;

	Systpar systpar;
	float fitrangemin;
	float fitrangemax;
	bool fitDisabled=false;

	public:	
	//standard constructor
	TemplateFIT(std::string Basename,Binning Bins, std::string Cut, int Nbins, float Xmin, float Xmax,int steps=11,float sigma=50,float shift=40){
		
		for(int bin=0;bin<Bins.size();bin++){
			fits.push_back(std::vector<std::vector<TFit *>>());
			for(int i=0;i<steps;i++){
				fits[bin].push_back(std::vector<TFit *>());
				for(int j=0;j<steps;j++){

					TFit * fit = new TFit;
					string named    =Basename + "_Data_" +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string namedprim=Basename + "_DataPrim_" +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameP    =Basename + "_MCP_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameD    =Basename + "_MCD_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);

					fit->Templ_P =  new TH1F(nameP.c_str(),nameP.c_str(),Nbins,Xmin,Xmax);
					fit->Templ_D =  new TH1F(nameD.c_str(),nameD.c_str(),Nbins,Xmin,Xmax);
					fit->Data    =  new TH1F(named.c_str(),named.c_str(),Nbins,Xmin,Xmax);
					fit->DataPrim=  new TH1F(namedprim.c_str(),namedprim.c_str(),Nbins,Xmin,Xmax);
					fits[bin][i].push_back(fit);
				}
			}
		}	
		basename=Basename;
		cut = Cut;
		cutprimary=Cut+"&IsPrimary";
		bins=Bins;
		
		StatError  = new TH1F("StatError","StatError",bins.size(),0,bins.size()) ;
        	SystError  = new TH1F("SystError","SystError",bins.size(),0,bins.size()) ;

		ProtonCounts    = new TH1F("Proton Counts","Proton Counts",bins.size(),0,bins.size()) ;
        	DeuteronCounts  = new TH1F("Deuteron Counts","Deuteron Counts",bins.size(),0,bins.size()) ;
	
		BestFitSigma  = new TH1F("Best FIt Sigma","Best FIt Sigma",bins.size(),0,bins.size()) ;
        	BestFitShift  = new TH1F("Best Fit Shift","Best FIt Shift",bins.size(),0,bins.size()) ;


	
		systpar.sigma=sigma;
		systpar.shift=shift;
		systpar.steps=steps;

		fitrangemin=0.6;
		fitrangemax=2.5;	
	}

	//reading constructor
	
	TemplateFIT(FileSaver  File, std::string Basename,Binning Bins,int steps=11,float sigma=50,float shift=40){

		TFile * file = File.GetFile();

		for(int bin=0;bin<Bins.size();bin++){
			fits.push_back(std::vector<std::vector<TFit *>>());
			for(int i=0;i<steps;i++){
				fits[bin].push_back(std::vector<TFit *>());
				for(int j=0;j<steps;j++){

					TFit * fit = new TFit;
					string named    =Basename + "/Bin "+ to_string(bin)+"/Data/" + Basename + "_Data_" +to_string(bin)+" "+to_string(0)+" "+to_string(5);
					string namedprim=Basename + "/Bin "+ to_string(bin)+"/Data/" + Basename + "_DataPrim_" +to_string(bin)+" "+to_string(0)+" "+to_string(5);
					string nameP    =Basename + "/Bin "+ to_string(bin)+"/TemplateP/" + Basename + "_MCP_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameD    =Basename + "/Bin "+ to_string(bin)+"/TemplateD/" + Basename + "_MCD_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);

					fit->Templ_P =  (TH1F *)file->Get(nameP.c_str());
					fit->Templ_D =  (TH1F *)file->Get(nameD.c_str());
					fit->Data    =  (TH1F *)file->Get(named.c_str());
					fit->DataPrim=  (TH1F *)file->Get(namedprim.c_str());
					fits[bin][i].push_back(fit);
				}
			}
		}	



			 basename=Basename;
			 bins=Bins; 

			 StatError  = new TH1F("StatError","StatError",bins.size(),0,bins.size()) ;
			 SystError  = new TH1F("SystError","SystError",bins.size(),0,bins.size()) ;

			 ProtonCounts    = new TH1F("Proton Counts","Proton Counts",bins.size(),0,bins.size()) ;
			 DeuteronCounts  = new TH1F("Deuteron Counts","Deuteron Counts",bins.size(),0,bins.size()) ;

			 BestFitSigma  = new TH1F("Best Fit Sigma","Best Fit Sigma",bins.size(),0,bins.size()) ;
        		 BestFitShift  = new TH1F("Best Fit Shift","Best Fit Shift",bins.size(),0,bins.size()) ;



			systpar.sigma=sigma;
			systpar.shift=shift;
			systpar.steps=steps;

			 fitrangemin=0.6;
			 fitrangemax=2.5;	



	}

	void Eval_TransferFunction();
	void Fill(TNtuple * treeMC,TNtuple * treeDT, Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars) );
	float SmearBeta(float Beta, float stepsigma, float stepshift);

	void FillEventByEventData(float var, float discr_var, bool CUT, bool CUTPRIM, float weight);
	void FillEventByEventMC(float var, float discr_var, bool CUTP, bool CUTD, float weight);

	void Do_TemplateFIT(TFit * Fit);
	void ExtractCounts(FileSaver finalhisto);
	void EvalFinalParameters();
	void CalculateFinalPDCounts();
	void Save(FileSaver finalhisto,bool recreate=false);	
	void SaveFitResults(FileSaver finalhisto);

	void SetSystematicParameters(int steps,float sigma,float shift){ systpar.steps=steps; systpar.shift=shift; systpar.sigma=sigma; return;};
	void SetFitRange(float min, float max){fitrangemin=min; fitrangemax=max; return;}

	void DisableFit(){fitDisabled=true;}

	std::string GetName(){return basename;}
	TH1F * GetStatError(){ return StatError;}
	TH1F * GetSystError() { return SystError;}	
	Binning  GetBinning() {return bins;}


	TH2F * GetDCountsSpread(int bin)	{ return DCountsSpread[bin];}
	TH2F * GetChiSquareSpread(int bin)      { return TFitChisquare[bin];}	
	TH1F * GetWeightedDCounts(int bin)     { return WeightedDCounts[bin];}

};

void TemplateFIT::Eval_TransferFunction(){
/*	for(int bin=0;bin<DataPrim.size();bin++){
		TH1F * transferfunction = (TH1F *) DataPrim[bin]->Clone();
		transferfunction->Sumw2();
		transferfunction->Divide(Data[bin]);
		TransferFunction.push_back(transferfunction);
	}*/
	return;
}

void TemplateFIT::Fill(TNtuple * treeMC,TNtuple * treeDT, Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars) ){

	cout<<basename.c_str()<<" Filling ... (Data)"<< endl;
	vars->ReadAnalysisBranches(treeDT);

	for(int i=0;i<treeDT->GetEntries();i++){
		vars->AnalysisVariablseReset();		
		UpdateProgressBar(i, treeDT->GetEntries());
		treeDT->GetEvent(i);
		FillEventByEventData(var(vars), discr_var(vars),ApplyCuts(cut,vars),ApplyCuts(cutprimary,vars),vars->mcweight);
	}
	
	cout<<basename.c_str()<<" Filling ... (MC Protons)"<< endl;
	vars->ReadAnalysisBranches(treeMC);

	for(int i=0;i<treeMC->GetEntries();i++){
		vars->AnalysisVariablseReset();		
		UpdateProgressBar(i, treeMC->GetEntries());
		treeMC->GetEvent(i);
		std::string cutP=cut+"&IsProtonMC";
		std::string cutD=cut+"&IsDeutonMC";
		FillEventByEventMC(vars->R, discr_var(vars),ApplyCuts(cutP,vars),ApplyCuts(cutD,vars),vars->mcweight);
	}

	return;
}


void TemplateFIT::FillEventByEventData(float var, float discr_var, bool CUT, bool CUTPRIM, float weight){

	int kbin;
	kbin = 	bins.GetBin(discr_var);
	if(CUT&&kbin>0){
		for(int i=0;i<systpar.steps;i++)
			for(int j=0;j<systpar.steps;j++){
				fits[kbin][i][j]->Data->Fill(var,weight);		
				if(CUTPRIM) fits[kbin][i][j]->DataPrim->Fill(var,weight);
				}
	}
	return;	

}

float TemplateFIT::SmearBeta(float Beta, float stepsigma, float stepshift){

	float time = 1.2/(Beta*3e-4);
	float shiftstart=-systpar.shift;
	time = time + (shiftstart+(2*systpar.shift/(float)systpar.steps)*stepshift) + Rand->Gaus(0,(float)((2*systpar.sigma/systpar.steps)*stepsigma));
	return 1.2/(time*3e-4);

}

void TemplateFIT::FillEventByEventMC(float var, float discr_var, bool CUTP, bool CUTD, float weight){

	if((CUTP||CUTD)){
		for(int i=0;i<systpar.steps;i++)
			for(int j=0;j<systpar.steps;j++){
				float betasmear= SmearBeta(discr_var,(float)i,(float)j);
				int kbin =  bins.GetBin(betasmear);
				float mass = var/betasmear * pow((1-pow(betasmear,2)),0.5);
				if(CUTP&&kbin>0) fits[kbin][i][j]->Templ_P->Fill(mass,weight);		
				if(CUTD&&kbin>0) fits[kbin][i][j]->Templ_D->Fill(mass,weight);
			}
	}
	return;	

}


void TemplateFIT::Save(FileSaver finalhisto,bool recreate){

	for(int bin=0;bin<bins.size();bin++){ 
		finalhisto.Add(fits[bin][0][5]->Data);
		finalhisto.Add(fits[bin][0][5]->DataPrim);

		finalhisto.writeObjsInFolder((basename + "/Bin "+ to_string(bin)+"/Data").c_str(),recreate);
	}

	for(int bin=0;bin<bins.size();bin++){ 
		for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){

				finalhisto.Add(fits[bin][i][j]->Templ_P);
		}
		finalhisto.writeObjsInFolder((basename + "/Bin "+ to_string(bin)+"/TemplateP").c_str(),recreate);
	}

	for(int bin=0;bin<bins.size();bin++){ 
		for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){

				finalhisto.Add(fits[bin][i][j]->Templ_D);
		}
		finalhisto.writeObjsInFolder((basename + "/Bin "+ to_string(bin)+"/TemplateD").c_str(),recreate);
	}



	return;
}

void TemplateFIT::SaveFitResults(FileSaver finalhisto){


	for(int bin=0;bin<bins.size();bin++){
	TH1F * OriginalP=(TH1F*)fits[bin][0][5]->Templ_P->Clone();
	TH1F * BestP=(TH1F*)fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_P->Clone();
	OriginalP->SetName("Original Proton MC ");
	BestP->SetName("Best #chi^{2} Mod. Proton MC ");
	finalhisto.Add(OriginalP);
	finalhisto.Add(BestP);
	finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesP/Bin"+to_string(bin)).c_str());
	
	TH1F * OriginalD=(TH1F*)fits[bin][0][5]->Templ_D->Clone();
	TH1F * BestD=(TH1F*)fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_D->Clone();
	OriginalD->SetName("Original Deuton MC ");
	BestD->SetName("Best #chi^{2} Mod. Deuton MC ");
	finalhisto.Add(OriginalD);
	finalhisto.Add(BestD);
	finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesD/Bin"+to_string(bin)).c_str());
	
	}
	
	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(!(i==0&&i==5))finalhisto.Add(fits[bin][i][j]->Templ_P);
			}	
	finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesP/Bin"+to_string(bin)).c_str());
	}

	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(!(i==0&&i==5))finalhisto.Add(fits[bin][i][j]->Templ_D);
				}
	finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesD/Bin"+to_string(bin)).c_str());
	}

/*	for(int i=0;i<bins.size();i++) 
                finalhisto.Add(TransferFunction[i]);
        finalhisto.writeObjsInFolder((basename + "/Data").c_str());
*/
	for(int i=0;i<bins.size();i++){
		finalhisto.Add(DCountsSpread[i]);
		finalhisto.writeObjsInFolder((basename+"/Fit Results/Spreads/DCounts").c_str());
	}
	for(int i=0;i<bins.size();i++){
		finalhisto.Add(TFitChisquare[i]);
		finalhisto.writeObjsInFolder((basename+"/Fit Results/Spreads/ChiSquare").c_str());	
	}
	for(int i=0;i<bins.size();i++){
		finalhisto.Add(WeightedDCounts[i]);
		finalhisto.writeObjsInFolder((basename+"/Fit Results/Spreads/Weighted D counts").c_str());	
	}


	finalhisto.Add(BestFitSigma);
	finalhisto.Add(BestFitShift);
	finalhisto.Add(StatError);
	finalhisto.Add(SystError);
	finalhisto.Add(ProtonCounts);
	finalhisto.Add(DeuteronCounts);

	finalhisto.writeObjsInFolder((basename+"/Fit Results/").c_str());

}


void TemplateFIT::Do_TemplateFIT(TFit * Fit){

	TObjArray *Tpl;
	Tpl = new TObjArray(2);
	Tpl -> Add( Fit ->  Templ_P );
	Tpl -> Add( Fit ->  Templ_D );
	Fit -> Tfit = new TFractionFitter(Fit -> Data, Tpl ,"q");
	float highP=1;
	float lowP=0.9;
	float highD=0.1;
	float lowD=0.00001;
	float min=fitrangemin;
	float max=fitrangemax;

	Fit -> Tfit -> SetRangeX(Fit -> Data -> FindBin(min), Fit -> Data -> FindBin(max));
	
	bool fitcondition = (Fit -> Data->Integral()>5000)&&(Fit -> Templ_P->Integral()>1000) &&(Fit -> Templ_D->Integral()>500);

	if(fitcondition) {  
		Fit -> Tfit_outcome = Fit -> Tfit -> Fit();

		for(int fit_attempt=0; fit_attempt<20; fit_attempt++) {
			cout<<fit_attempt<<endl;
			if(Fit -> Tfit_outcome == 0) break;
			else {
				cout<<fit_attempt<<endl;
				Fit -> Tfit -> Constrain(0, lowP+(float)fit_attempt/1000 ,highP );
				Fit -> Tfit -> Constrain(1, lowD-(float)fit_attempt/10000 ,highD-fit_attempt/10000 );
				Fit -> Tfit_outcome = Fit -> Tfit -> Fit();
			}
		}

		if(Fit -> Tfit_outcome==0){
			TH1F * Result = (TH1F *) Fit-> Tfit -> GetPlot();
			float itot= Result->Integral();
			double w1,e1 = 0;
			double w2,e2 = 0;
			Fit -> Tfit ->GetResult(0,w1,e1);
			Fit -> Tfit ->GetResult(1,w2,e2);
			float i1 = Fit-> Templ_P  ->Integral(Fit->Templ_P -> FindBin(min), Fit->Templ_P -> FindBin(max));
			float i2 = Fit-> Templ_D  ->Integral(Fit->Templ_D -> FindBin(min), Fit->Templ_D -> FindBin(max));
			Fit ->wheightP= w1*itot/i1;
			Fit ->wheightD= w2*itot/i2;
			cout<<w1<<" "<<w2<<endl;

			Fit ->  Templ_P  -> Scale(Fit ->wheightP);
			Fit ->  Templ_D  -> Scale(Fit ->wheightD);

			float Cov01=0;

			Cov01= Fit -> Tfit->GetFitter()->GetCovarianceMatrixElement(0,1);
			Fit -> StatErr = pow(pow(w2*e2,2)+pow(w1*e1,2)-2*Cov01*w1*w2,0.5);

			Fit -> ChiSquare = Fit -> Tfit -> GetChisquare()/(float) (Fit ->  Templ_P ->GetNbinsX()*1.7);
			Fit -> DCounts = Fit ->  Templ_D -> Integral();
			Fit -> PCounts = Fit ->  Templ_P -> Integral();
		}
		else{
			Fit ->wheightP= 0;
			Fit ->wheightD= 0;
		}
	}
	return;
}


void TemplateFIT::ExtractCounts(FileSaver finalhisto){
	
	Eval_TransferFunction();

	if(!fitDisabled){

	for(int bin=0;bin<bins.size();bin++){

		// FITS for systematic error
		TH2F * dcountsspread = new TH2F(("DCountsSpread Bin " +to_string(bin)).c_str(),("DCountsSpread Bin " +to_string(bin)).c_str(),systpar.steps,0,2*systpar.sigma,systpar.steps,-systpar.shift,systpar.shift);
        	TH2F * tfitchisquare = new TH2F(("ChiSquare Bin " +to_string(bin)).c_str(),("ChiSquare Bin " +to_string(bin)).c_str(),systpar.steps,0,2*systpar.sigma,systpar.steps,-systpar.shift,systpar.shift);
		TH1F * weighteddcounts  = new TH1F(("Weighted Counts Bin " +to_string(bin)).c_str(),("Weighted Counts Bin " +to_string(bin)).c_str(),25,0.4*fits[bin][0][5]->DCounts,1.7*fits[bin][0][5]->DCounts);
		BestChi * MinimumChi = new BestChi();
		

		for(int sigma=0;sigma<systpar.steps;sigma++){
			for(int shift=0;shift<systpar.steps;shift++){

				cout<<fits[bin][sigma][shift]<<endl;

				Do_TemplateFIT(fits[bin][sigma][shift]);

				
				dcountsspread->SetBinContent(sigma+1,shift+1,fits[bin][sigma][shift]->DCounts);
				if(fits[bin][sigma][shift]->ChiSquare>0)
						tfitchisquare->SetBinContent(sigma+1,shift+1,fits[bin][sigma][shift]->ChiSquare);
				else  tfitchisquare->SetBinContent(sigma+1,shift+1,500);
				weighteddcounts -> Fill(dcountsspread->GetBinContent(sigma+1,shift+1),exp(-fits[bin][sigma][shift]->ChiSquare));

			}

		}

		MinimumChi->FindMinimum(tfitchisquare);	
		

		BestChiSquare.push_back(MinimumChi);			
		DCountsSpread.push_back(dcountsspread);
		TFitChisquare.push_back(tfitchisquare);
		WeightedDCounts.push_back(weighteddcounts);


	}

	EvalFinalParameters();
	CalculateFinalPDCounts();	
	}
	
	return;
}

void TemplateFIT::EvalFinalParameters(){
	for(int bin=0;bin<bins.size();bin++){
		float BinSTD=0;
		for(int i=0;i<systpar.steps;i++){
			for(int j=0;j<systpar.steps;j++){
				BinSTD+=pow(TFitChisquare[bin]->GetBinContent(i+1,j+1)-BestChiSquare[bin]->chimin,2);
			}
		}
		BinSTD=pow(BinSTD,0.5)/(systpar.steps*systpar.steps);

		BestFitSigma->SetBinContent(bin+1,(float)((2*systpar.sigma/systpar.steps)*BestChiSquare[bin]->i));
		BestFitShift->SetBinContent(bin+1,-systpar.shift+(2*systpar.shift/(float)systpar.steps)*BestChiSquare[bin]->j);
		
		if(BinSTD>0){
		BestFitSigma->SetBinError(bin+1,4/BinSTD);
		BestFitShift->SetBinError(bin+1,4/BinSTD);
		}


		if(fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts>0){
			StatError -> SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->StatErr);
			SystError -> SetBinContent(bin+1,WeightedDCounts[bin]->GetStdDev()/fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts);
			StatError -> SetBinError(bin+1,0);
			SystError -> SetBinError(bin+1,0);
		}

	}
	return;

}


void TemplateFIT::CalculateFinalPDCounts(){
	for(int bin=0;bin<bins.size();bin++){

		StatError -> Smooth(1);
		SystError -> Smooth(1);
		
		float toterr= pow(pow(StatError -> GetBinContent(bin+1),2) + pow(SystError -> GetBinContent(bin+1),2) ,0.5);
		
		ProtonCounts->SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCounts);
		ProtonCounts->SetBinError(bin+1,toterr*fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCounts);

		DeuteronCounts->SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts);
                DeuteronCounts->SetBinError(bin+1,toterr*fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts);		
		
	}
		

}

