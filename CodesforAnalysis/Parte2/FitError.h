using namespace std;

float FitError(float fitvalues[],float values[],float sigmas[],float bins,int param){
	float fitX2=0;
        for(int i=1;i<bins;i++) fitX2+=pow((fitvalues[i]-values[i])/sigmas[i],2);
	fitX2=(fitX2)/(bins-param);
        float sigmamean=0;
        for(int i=1;i<bins;i++) sigmamean+=pow(sigmas[i],2);
        sigmamean=pow(sigmamean,0.5)/(bins-param);
        float errorefit=fitX2*sigmamean;
	return errorefit;
}


