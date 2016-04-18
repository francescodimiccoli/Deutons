
TemplateFIT * FitTOF_Dbins = new TemplateFIT("FitTOF_Dbins",nbinsToF,0,3,11);
TemplateFIT * FitNaF_Dbins = new TemplateFIT("FitNaF_Dbins",nbinsNaF,0,3,11);
TemplateFIT * FitAgl_Dbins = new TemplateFIT("FitAgl_Dbins",nbinsAgl,0,3,11);
TemplateFIT * FitTOF_Pbins = new TemplateFIT("FitTOF_Pbins",nbinsToF,0,3,11);
TemplateFIT * FitNaF_Pbins = new TemplateFIT("FitNaF_Pbins",nbinsNaF,0,3,11);
TemplateFIT * FitAgl_Pbins = new TemplateFIT("FitAgl_Pbins",nbinsAgl,0,3,11);


void DeutonsMC_Fill(TNtuple *ntupla, int l){
	int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	float mass = 0;
	if(!(Likcut&&Distcut)) return;
	for(int m=0;m<nbinsToF;m++){ //TOF
		mass = ((R/Beta)*pow((1-pow(Beta,2)),0.5));
		if(Var>BetaD[m]&&Var<=BetaD[m+1]){
			if(Massa_gen<1&&Massa_gen>0.5) FitTOF_Dbins -> TemplateP -> Fill(mass,m);
			if(Massa_gen<2&&Massa_gen>1.5) FitTOF_Dbins -> TemplateD -> Fill(mass,m);
			if(Massa_gen<4&&Massa_gen>2.5) FitTOF_Dbins -> TemplateHe-> Fill(mass,m);
		}
		if(Var>BetaP[m]&&Var<=BetaP[m+1]) {
			if(Massa_gen<1&&Massa_gen>0.5) FitTOF_Pbins -> TemplateP -> Fill(mass,m);
			if(Massa_gen<2&&Massa_gen>1.5) FitTOF_Pbins -> TemplateD -> Fill(mass,m);
			if(Massa_gen<4&&Massa_gen>2.5) FitTOF_Pbins -> TemplateHe-> Fill(mass,m);
		}
	}
	for(int m=0;m<nbinsNaF;m++) { //NaF
		if((((int)Cutmask)>>11)==512){
			mass = ((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5));
			if(Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5) FitNaF_Dbins -> TemplateP -> Fill(mass,m);
				if(Massa_gen<2&&Massa_gen>1.5) FitNaF_Dbins -> TemplateD -> Fill(mass,m);
				if(Massa_gen<4&&Massa_gen>2.5) FitNaF_Dbins -> TemplateHe-> Fill(mass,m);
			}
			if(Var2>BetaNaFP[m]&&Var2<=BetaNaFP[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5) FitNaF_Pbins -> TemplateP -> Fill(mass,m);
				if(Massa_gen<2&&Massa_gen>1.5) FitNaF_Pbins -> TemplateD -> Fill(mass,m);
				if(Massa_gen<4&&Massa_gen>2.5) FitNaF_Pbins -> TemplateHe-> Fill(mass,m);
			}
		}
	}
	for(int m=0;m<nbinsAgl;m++) { //Agl
		if((((int)Cutmask)>>11)==0){
			mass = ((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5));
			if(Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5) FitAgl_Dbins -> TemplateP -> Fill(mass,m);
				if(Massa_gen<2&&Massa_gen>1.5) FitAgl_Dbins -> TemplateD -> Fill(mass,m);
				if(Massa_gen<4&&Massa_gen>2.5) FitAgl_Dbins -> TemplateHe-> Fill(mass,m);
			}
			if(Var2>BetaAglP[m]&&Var2<=BetaAglP[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5) FitAgl_Pbins -> TemplateP -> Fill(mass,m);
				if(Massa_gen<2&&Massa_gen>1.5) FitAgl_Pbins -> TemplateD -> Fill(mass,m);
				if(Massa_gen<4&&Massa_gen>2.5) FitAgl_Pbins -> TemplateHe-> Fill(mass,m);
			}
		}
	}
} 


void DeutonsDATA_Fill(TNtuple *ntupla, int l,int zona){
	int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	float mass = 0;
	if(!(Likcut&&Distcut)) return;
	for(int m=0;m<nbinsToF;m++){ //TOF
		mass = ((R/Beta)*pow((1-pow(Beta,2)),0.5));
		if(Var>BetaD[m]&&Var<=BetaD[m+1]){
			if(R>1.2*Rcutoff) FitTOF_Dbins -> Data_Prim -> Fill(mass,m);
			((TH3*)FitTOF_Dbins -> Data_Geomag) -> Fill(mass,m,zona);
		}
		if(Var>BetaP[m]&&Var<=BetaP[m+1]) {
			if(R>1.2*Rcutoff) FitTOF_Pbins -> Data_Prim -> Fill(mass,m);
		}
	}
	for(int m=0;m<nbinsNaF;m++){//NaF
		if((((int)Cutmask)>>11)==512){
			mass = ((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5));
			if(Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) {
				if(R>1.2*Rcutoff) FitNaF_Dbins -> Data_Prim -> Fill(mass,m);
				((TH3*)FitNaF_Dbins -> Data_Geomag) -> Fill(mass,m,zona);
			}
			if(Var2>BetaNaFP[m]&&Var2<=BetaNaFP[m+1]) {
				if(R>1.2*Rcutoff) FitNaF_Pbins -> Data_Prim -> Fill(mass,m);
			}
		}
	}
	for(int m=0;m<nbinsAgl;m++){ //Agl
		if((((int)Cutmask)>>11)==0){
			mass = ((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5));
			if(Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) {
				if(R>1.2*Rcutoff) FitAgl_Dbins -> Data_Prim -> Fill(mass,m);
				((TH3*)FitAgl_Dbins -> Data_Geomag) -> Fill(mass,m,zona);
			}
			if(Var2>BetaAglP[m]&&Var2<=BetaAglP[m+1]) {
				if(R>1.2*Rcutoff) FitAgl_Pbins -> Data_Prim -> Fill(mass,m);
			}
		}
	}
	return;
}


void DeutonsMC_Write(){

	FitTOF_Dbins -> Write();
	FitNaF_Dbins -> Write();
	FitAgl_Dbins -> Write();
	
	FitTOF_Pbins -> Write();
	FitNaF_Pbins -> Write();
	FitAgl_Pbins -> Write();
	return;
}


void DeutonsTemplFits(){
	string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	TFile * file1 = TFile::Open(nomefile.c_str(),"READ");

	TemplateFIT * FitTOF_Dbins	= new TemplateFIT(file1,"FitTOF_Dbins",0,3,11);
	TemplateFIT * FitNaF_Dbins	= new TemplateFIT(file1,"FitNaF_Dbins",0,3,11);
	TemplateFIT * FitAgl_Dbins	= new TemplateFIT(file1,"FitAgl_Dbins",0,3,11);

	TemplateFIT * FitTOF_Pbins	= new TemplateFIT(file1,"FitTOF_Pbins",0,3,11);
	TemplateFIT * FitNaF_Pbins	= new TemplateFIT(file1,"FitNaF_Pbins",0,3,11);
	TemplateFIT * FitAgl_Pbins	= new TemplateFIT(file1,"FitAgl_Pbins",0,3,11);

	cout<<"******************** DEUTONS TEMPlATE FITS ************************"<<endl;

	FitTOF_Dbins -> TemplateFits(); 
	FitNaF_Dbins -> TemplateFits();
	FitAgl_Dbins -> TemplateFits();

	FitTOF_Pbins -> TemplateFits();
	FitNaF_Pbins -> TemplateFits();
	FitAgl_Pbins -> TemplateFits();

	cout<<"*** Updating P1 file ****"<<endl;

        nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
        file1 = TFile::Open(nomefile.c_str(),"UPDATE");

	FitTOF_Dbins -> DCounts -> Write ("D_FluxCounts_TOF");
	FitNaF_Dbins -> DCounts -> Write ("D_FluxCounts_NaF");
        FitAgl_Dbins -> DCounts -> Write ("D_FluxCounts_Agl");
                       
        FitTOF_Pbins -> PCounts -> Write ("P_FluxCounts_TOF");
	FitNaF_Pbins -> PCounts -> Write ("P_FluxCounts_NaF");
        FitAgl_Pbins -> PCounts -> Write ("P_FluxCounts_Agl");

	file1 -> Write();
	file1 -> Close(); 	
	
	TCanvas *c30_TOF[12][nbinsToF];
	TCanvas *c30_NaF[12][nbinsNaF];
	TCanvas *c30_Agl[12][nbinsAgl];
	
	for(int bin=0; bin <nbinsToF; bin++){
		c30_TOF[0][bin] = new TCanvas(("TOF bin:" + to_string(bin)).c_str());
		FitTOF_Dbins -> TemplateFitPlot(c30_TOF[0][bin],"Mass [GeV/C^2]",bin);
	}
	for(int bin=0; bin <nbinsNaF; bin++){
		c30_NaF[0][bin] = new TCanvas(("NaF bin:" + to_string(bin)).c_str());
		FitNaF_Dbins -> TemplateFitPlot(c30_NaF[0][bin],"Mass [GeV/C^2]",bin);
	}
	for(int bin=0; bin <nbinsAgl; bin++){
		c30_Agl[0][bin] = new TCanvas(("Agl bin:" + to_string(bin)).c_str());
		FitAgl_Dbins -> TemplateFitPlot(c30_Agl[0][bin],"Mass [GeV/C^2]",bin);
	}
	
	for(int bin=0; bin <nbinsToF; bin++){
		c30_TOF[11][bin] = new TCanvas(("TOF P bin:" + to_string(bin)).c_str());
		FitTOF_Pbins -> TemplateFitPlot(c30_TOF[11][bin],"Mass [GeV/C^2]",bin);
	}
	for(int bin=0; bin <nbinsNaF; bin++){
		c30_NaF[11][bin] = new TCanvas(("NaF P bin:" + to_string(bin)).c_str());
		FitNaF_Pbins -> TemplateFitPlot(c30_NaF[11][bin],"Mass [GeV/C^2]",bin);
	}
	for(int bin=0; bin <nbinsAgl; bin++){
		c30_Agl[11][bin] = new TCanvas(("Agl P bin:" + to_string(bin)).c_str());
		FitAgl_Pbins -> TemplateFitPlot(c30_Agl[11][bin],"Mass [GeV/C^2]",bin);
	}
	

	cout<<"*** Updating Results file ***"<<endl;
        nomefile="./Final_plots/"+mese+".root";
        TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
 
        f_out->mkdir("Mass Template Fits/TOF/Primaries/Dbins");
        f_out->cd("Mass Template Fits/TOF/Primaries/Dbins");
        for(int bin=0; bin <nbinsToF; bin++)
		c30_TOF[0][bin]->Write();
        f_out->mkdir("Mass Template Fits/NaF/Primaries/Dbins");
        f_out->cd("Mass Template Fits/NaF/Primaries/Dbins");
	for(int bin=0; bin <nbinsNaF; bin++)
		c30_NaF[0][bin]->Write();
	f_out->mkdir("Mass Template Fits/Agl/Primaries/Dbins");
        f_out->cd("Mass Template Fits/Agl/Primaries/Dbins");
	for(int bin=0; bin <nbinsAgl; bin++)
		c30_Agl[0][bin]->Write();

	f_out->mkdir("Mass Template Fits/TOF/Primaries/Pbins");
        f_out->cd("Mass Template Fits/TOF/Primaries/Pbins");
        for(int bin=0; bin <nbinsToF; bin++)
		c30_TOF[11][bin]->Write();
        f_out->mkdir("Mass Template Fits/NaF/Primaries/Pbins");
        f_out->cd("Mass Template Fits/NaF/Primaries/Pbins");
	for(int bin=0; bin <nbinsNaF; bin++)
		c30_NaF[11][bin]->Write();
	f_out->mkdir("Mass Template Fits/Agl/Primaries/Pbins");
        f_out->cd("Mass Template Fits/Agl/Primaries/Pbins");
	for(int bin=0; bin <nbinsAgl; bin++)
		c30_Agl[11][bin]->Write();




	f_out->Write();
        f_out->Close();



	return;


}
