using namespace std;

class FitFunction {

	private: 
	TH1 * RoughValues;
	TH1 * FittedValues;
	int n_smooth=0;
	int mc_types=0;

	public:
	FitFunction( TH1 * values, int nsmooth=1){
		RoughValues = (TH1 *) values -> Clone();
		FittedValues= (TH1 *) values -> Clone(); 	
	
		n_smooth = nsmooth;
		mc_types = RoughValues -> GetNbinsY();	
	}	
		

	void FitValues();
	TH1 * ReturnRoughValues() {return RoughValues;}
	TH1 * ReturnFittedValues() {return FittedValues;}


};


float FitError(TH1F * fit,TH1F * values,float bins,int param){

        float fitX2=0;
        for(int i=1;i<bins;i++) {
                if((fit->GetBinContent(i+1)-values->GetBinContent(i+1))/values->GetBinError(i+1)<10)
			fitX2+=pow((fit->GetBinContent(i+1)-values->GetBinContent(i+1))/values->GetBinError(i+1),2);
                }
        fitX2=(fitX2)/(bins-param);

        float sigmamean=0;
        for(int i=1;i<bins;i++)
                sigmamean+=pow(values->GetBinError(i+1),2);
        sigmamean=pow(sigmamean,0.5)/(bins-param);

        float errorefit=fitX2*sigmamean;
        return errorefit;
}




void FitFunction::FitValues(){

	if(mc_types==1){

		FittedValues -> Smooth(n_smooth);

		float fiterror = FitError((TH1F*)FittedValues,(TH1F*)RoughValues,FittedValues ->GetNbinsX(),0);	

		for(int nbins = 0; nbins < FittedValues ->GetNbinsX(); nbins ++){
			FittedValues -> SetBinError(nbins +1, fiterror);
		}
	}
	
	return;

}










