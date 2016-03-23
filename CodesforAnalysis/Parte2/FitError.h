using namespace std;

float FitError(TH1 * fit,TH1 * values,float bins,int param){
	float fitX2=0;
        for(int i=1;i<bins;i++) fitX2+=pow((fit->GetBinContent(i+1)-values->GetBinContent(i+1))/values->GetBinError(i+1),2);
	fitX2=(fitX2)/(bins-param);
        float sigmamean=0;
        for(int i=1;i<bins;i++) sigmamean+=pow(values->GetBinError(i+1),2);
        sigmamean=pow(sigmamean,0.5)/(bins-param);
        float errorefit=fitX2*sigmamean;
	return errorefit;
}


