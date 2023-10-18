#include "RatioMaker.h"

TH1F * DoRatio_(TH1F* numerator, TH1F * denominator,std::string name){

        TH1F * ratio;
        if(numerator->GetBinLowEdge(0)>denominator->GetBinLowEdge(0))
        ratio=(TH1F*) numerator->Clone(name.c_str());
        else
        ratio=(TH1F*) denominator->Clone(name.c_str());

        ratio->Reset();

        for(int i=0;i<ratio->GetNbinsX();i++){
                int bin_num=numerator->FindBin(ratio->GetBinCenter(i+1));
                int bin_den=denominator->FindBin(ratio->GetBinCenter(i+1));
                if(denominator->GetBinContent(bin_den)>0){
                        ratio->SetBinContent(i+1,numerator->GetBinContent(bin_num)/denominator->GetBinContent(bin_den));
                        float error = pow(denominator->GetBinError(bin_den)/denominator->GetBinContent(bin_den),2);
                        error += pow(numerator->GetBinError(bin_num)/numerator->GetBinContent(bin_num),2);
                        ratio->SetBinError(i+1,pow(error,0.5)*ratio->GetBinContent(i+1));
                }
        }
        return ratio;
}



void RatioMaker::DoRatioErrStat(float C){

	ratio_errstat = (TH1F*)ratio_unf->Clone((Name+ "err_stat").c_str());
	ratio_errstat->Reset();

	for(int i=0; i<ratio_errstat->GetNbinsX();i++){
		float x = ratio_errstat->GetBinCenter(i+1);
		float n = Num->flux_unf->GetBinContent(Num->flux_unf->FindBin(x));
		float d = Den->flux_unf->GetBinContent(Den->flux_unf->FindBin(x));
		float dn = Num->statErr->GetBinContent(Num->statErr->FindBin(x))*n;
		float dd = Den->statErr->GetBinContent(Den->statErr->FindBin(x))*d;
		float r = ratio_unf->GetBinContent(ratio->FindBin(x));
		ratio_errstat->SetBinError(i+1,0);
		ratio_errstat->SetBinContent(i+1, sqrt((dn*dn)/(d*d)+(n*n/pow(d,4))*(dd*dd) - 2*(n/(d*d*d))*C)/r);

	}

}

void RatioMaker::DoRatioErrSyst(float C){

	ratio_errsyst = (TH1F*)ratio_unf->Clone((Name+ "err_syst").c_str());
	ratio_errsyst->Reset();

	for(int i=0; i<ratio_errsyst->GetNbinsX();i++){
		float x = ratio_errsyst->GetBinCenter(i+1);
		float n = Num->flux_unf->GetBinContent(Num->flux_unf->FindBin(x));
		float d = Den->flux_unf->GetBinContent(Den->flux_unf->FindBin(x));
		float dn = sqrt(pow(Num->systErr->GetBinContent(Num->systErr->FindBin(x)),2)+pow(Num->unfErr->GetBinContent(Num->unfErr->FindBin(x)),2))*n;
		float dd = sqrt(pow(Den->systErr->GetBinContent(Den->systErr->FindBin(x)),2)+pow(Den->unfErr->GetBinContent(Den->unfErr->FindBin(x)),2))*d;
		float r = ratio_unf->GetBinContent(ratio->FindBin(x));
		if(dn/n>0.1) dn=0.1*n;
		ratio_errsyst->SetBinError(i+1,0);
		cout<<"RATIOERR "<<sqrt(dn*dn/(d*d)+(n*n/pow(d,4))*(dd*dd) - 2*n/(d*d*d)*C)/r<<" "<<dn<<" "<<dd<<" "<<n<<" "<<d<<endl;
		ratio_errsyst->SetBinContent(i+1, sqrt(dn*dn/(d*d)+(n*n/pow(d,4))*(dd*dd) - 2*n/(d*d*d)*C)/r);

	}

}



RatioMaker::RatioMaker(ResultMerger * num, ResultMerger * den, std::string name){
	
	Num = num;
	Den = den;
	Name = name;

	ratio 	  	  = DoRatio_(num->flux		,den->flux,	name.c_str());
	ratio_stat	  = DoRatio_(num->flux_stat	,den->flux_stat,(name + "_stat").c_str());
	ratio_unf 	  = DoRatio_(num->flux_unf	,den->flux_unf,	(name + "_unf").c_str());
	ratio_unf_stat 	  = DoRatio_(num->flux_unf_stat  ,den->flux_unf_stat,	(name + "_unf_stat").c_str());
	ratioacceptance   = DoRatio_(num->effAcc  ,den->effAcc,	(name + "_acc").c_str());
	ratiocounts       = DoRatio_(num->counts  ,den->counts,	(name + "_counts").c_str());
	ratiounfolding    = DoRatio_(num->roounfolding  ,den->roounfolding,	(name + "_roounf").c_str());



	DoRatioErrStat(0);
	DoRatioErrSyst(0);

}



void RatioMaker::SaveResults(FileSaver finalhistos){


	if(ratio 	 )     	finalhistos.Add(ratio 	 );
        if(ratio_stat    )	finalhistos.Add(ratio_stat    );
        if(ratio_unf 	 )      finalhistos.Add(ratio_unf 	 );
        if(ratio_unf_stat) 	finalhistos.Add(ratio_unf_stat);
	if(ratio_errstat )	finalhistos.Add(ratio_errstat );
	if(ratio_errsyst )	finalhistos.Add(ratio_errsyst );
        if(ratioacceptance    )	finalhistos.Add(ratioacceptance    );
        if(ratiocounts 	 )      finalhistos.Add(ratiocounts 	 );
        if(ratiounfolding) 	finalhistos.Add(ratiounfolding);
	

	
	finalhistos.writeObjsInFolder(("MergedRatios/"+Name+"_dir/").c_str());

}
