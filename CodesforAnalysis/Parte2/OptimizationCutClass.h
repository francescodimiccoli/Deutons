using namespace std;


class OptimizationCut{

	private: 

	 std::vector<TH1F *> OptimizationCurves;
	int bins;
	TH1F * cuts;
	TH1F * EffcutP;
	TH1F * EffcutD;

	public:

	TH1 * Distrib_P   ;
        TH1 * Distrib_D   ;
	
	
	//creation constructors
	//standard
		
	OptimizationCut(std::string basename ,int Nbins ,float val_min , float val_max)	{
		 Distrib_P =   new TH2F((basename + "_P"   ).c_str(),(basename + "_P"   ).c_str(),100,val_min,val_max,Nbins,0,Nbins);
         	 Distrib_D =   new TH2F((basename + "_D"   ).c_str(),(basename + "_D"   ).c_str(),100,val_min,val_max,Nbins,0,Nbins);	
	
	  bins= Nbins;
	}

	//reading constructors
	 OptimizationCut(TFile *file, std::string basename) {
		Distrib_P   =  (TH1 *)file->Get((basename   + "_P"     ).c_str());
                Distrib_D   =  (TH1 *)file->Get((basename   + "_D"     ).c_str());
		bins = Distrib_P->GetNbinsY();
		cuts = new TH1F((basename + "_cuts"   ).c_str(),(basename + "_cuts"   ).c_str(),bins,0,bins);
		EffcutP=new TH1F((basename + "_EffP"   ).c_str(),(basename + "_EffP"   ).c_str(),bins,0,bins);
		EffcutD=new TH1F((basename + "_EffD"   ).c_str(),(basename + "_EffD"   ).c_str(),bins,0,bins);
	}

	OptimizationCut(TFile *file, std::string basename, std::string dirname) {
                Distrib_P   =  (TH1 *)file->Get((basename   + "_P"     ).c_str());
                Distrib_D   =  (TH1 *)file->Get((basename   + "_D"     ).c_str());
                bins = Distrib_P->GetNbinsY();
                cuts = (TH1F*)file->Get(("/" + dirname +"/" + basename + "_cuts"   ).c_str());
        	EffcutP=new TH1F((basename + "_EffP"   ).c_str(),(basename + "_EffP"   ).c_str(),bins,0,bins);
                EffcutD=new TH1F((basename + "_EffD"   ).c_str(),(basename + "_EffD"   ).c_str(),bins,0,bins);
	}


	// Methods
	void Write();
	TH1* GetBinP(int nbin);
	TH1* GetBinD(int nbin);
	void OptimizationBin(int nbin);
	void NormalizeDistributions();
	
	void Optimization();
	TH1F * GetOptimizationCurve(int nbin) {return OptimizationCurves[nbin];}
	TH1F * Getcuts() {return cuts;}
	
	void EfficiencyCutBin(int nbin);
	void Eval_Efficiencies();
	TH1F * GetEfficiency_P() {return EffcutP;}
	TH1F * GetEfficiency_D() {return EffcutD;}
};


void  OptimizationCut::Write()
{

    Distrib_P -> Write();
    Distrib_D -> Write();
 
    return;
}

void OptimizationCut::NormalizeDistributions(){
	Distrib_P -> Scale(1/Distrib_P->GetEntries());
	Distrib_D -> Scale(1/Distrib_D->GetEntries());
}


TH1*  OptimizationCut::GetBinP(int nbin) {
	string nome = "D_" + to_string(nbin);
	TH1F * Slice =  ProjectionXtoTH1F( (TH2F* )Distrib_P,nome,nbin,nbin+1);
	return Slice;
}	

TH1*  OptimizationCut::GetBinD(int nbin) {
	string nome = "P_" + to_string(nbin);
        TH1F * Slice =  ProjectionXtoTH1F( (TH2F *)Distrib_D,nome,nbin,nbin+1);
        return Slice;
}


void OptimizationCut::OptimizationBin(int nbin){
	TH1F * HistoP = (TH1F*)GetBinP(nbin);
	TH1F * HistoD = (TH1F*)GetBinD(nbin);

	TH1F * Opt =(TH1F *) HistoP -> Clone();	
	float effP = 0;
	float effD = 0;

	for (int x = HistoP ->GetNbinsX(); x>0; x--) Opt -> SetBinContent(x,0);

	for (int x = HistoP ->GetNbinsX(); x>0; x--){
		effP = 1 - HistoP->Integral(0,x)/(float)HistoP->Integral();
		effD = 1 - HistoD->Integral(0,x)/(float)HistoD->Integral();
		
		if(effP>1e-5&&effD>0.001) 
		     Opt->SetBinContent(x+1,pow(effD,1.5)/effP);
		else Opt->SetBinContent(x+1,0);
	}
	OptimizationCurves.push_back(Opt);
	float cutvalue=Opt->GetBinCenter(Opt->FindBin(0.6)/*Opt->GetMaximumBin()i*/);
	if(cutvalue>0.35) cuts->SetBinContent(nbin+1,cutvalue);
	else cuts->SetBinContent(nbin+1,cutvalue);
	return;	
}



void OptimizationCut::Optimization(){
	for (int i =0;i<bins;i++) OptimizationCut::OptimizationBin(i);

	return;
}


void OptimizationCut::EfficiencyCutBin(int nbin){
	TH1F * HistoP = (TH1F*)GetBinP(nbin);
        TH1F * HistoD = (TH1F*)GetBinD(nbin);
	
	int x = HistoP->FindBin(cuts->GetBinContent(nbin+1));

	EffcutP->SetBinContent(nbin+1,1-HistoP->Integral(0,x)/(float)HistoP->Integral());	
	EffcutD->SetBinContent(nbin+1,1-HistoD->Integral(0,x)/(float)HistoD->Integral());

	return;
}


void OptimizationCut::Eval_Efficiencies(){
	for (int i =0;i<bins;i++) {OptimizationCut::EfficiencyCutBin(i); cout<<"bin "<<i<<endl;}
	
	return;
}
