TemplateFIT * FitTOF_Dbins = new TemplateFIT("FitTOF_Dbins",nbinsToF,0,3,6);
TemplateFIT * FitNaF_Dbins = new TemplateFIT("FitNaF_Dbins",nbinsNaF,0,3,6);
TemplateFIT * FitAgl_Dbins = new TemplateFIT("FitAgl_Dbins",nbinsAgl,0,3,6);

TemplateFIT * FitTOFgeo_Dbins = new TemplateFIT("FitTOFgeo_Dbins",nbinsToF,0,3,6,11);
TemplateFIT * FitNaFgeo_Dbins = new TemplateFIT("FitNaFgeo_Dbins",nbinsNaF,0,3,6,11);
TemplateFIT * FitAglgeo_Dbins = new TemplateFIT("FitAglgeo_Dbins",nbinsAgl,0,3,6,11);

TemplateFIT * FitTOF_Pbins = new TemplateFIT("FitTOF_Pbins",nbinsToF,0,3,6);
TemplateFIT * FitNaF_Pbins = new TemplateFIT("FitNaF_Pbins",nbinsNaF,0,3,6);
TemplateFIT * FitAgl_Pbins = new TemplateFIT("FitAgl_Pbins",nbinsAgl,0,3,6);


void DeutonsMC_Fill(TNtuple *ntupla, int l){
	 ntupla->GetEvent(l);
	if(Tup.Beta<=0||Tup.R<=0) return;
	float mass = 0;
	//cuts
	if(!(Likcut&&Distcut)) return;
	//
	for(int m=0;m<nbinsToF;m++){ //TOF
		mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));
		if(RUsed>ToFDB.MomBins()[m]&&RUsed<=ToFDB.MomBins()[m+1]){
			if(Massa_gen<1&&Massa_gen>0.5) FitTOF_Dbins -> TemplateP -> Fill(mass,m);
			if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitTOF_Dbins -> TemplateD) -> Fill(mass,m,ReturnMCGenType());
			if(Massa_gen<4&&Massa_gen>2.5) FitTOF_Dbins -> TemplateHe-> Fill(mass,m);
		}
		if(RUsed>ToFPB.MomBins()[m]&&RUsed<=ToFPB.MomBins()[m+1]) {
			if(Massa_gen<1&&Massa_gen>0.5) FitTOF_Pbins -> TemplateP -> Fill(mass,m);
			if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitTOF_Pbins -> TemplateD) -> Fill(mass,m,ReturnMCGenType());
			if(Massa_gen<4&&Massa_gen>2.5) FitTOF_Pbins -> TemplateHe-> Fill(mass,m);
		}
	}
	for(int m=0;m<nbinsNaF;m++) { //NaF
		if(cmask.isFromNaF()){
			mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
			if(RUsed>NaFDB.MomBins()[m]&&RUsed<=NaFDB.MomBins()[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5) FitNaF_Dbins -> TemplateP -> Fill(mass,m);
				if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitNaF_Dbins -> TemplateD) -> Fill(mass,m,ReturnMCGenType());
				if(Massa_gen<4&&Massa_gen>2.5) FitNaF_Dbins -> TemplateHe-> Fill(mass,m);
			}
			if(RUsed>NaFPB.MomBins()[m]&&RUsed<=NaFPB.MomBins()[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5) FitNaF_Pbins -> TemplateP -> Fill(mass,m);
				if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitNaF_Pbins -> TemplateD) -> Fill(mass,m,ReturnMCGenType());
				if(Massa_gen<4&&Massa_gen>2.5) FitNaF_Pbins -> TemplateHe-> Fill(mass,m);
			}
		}
	}
	for(int m=0;m<nbinsAgl;m++) { //Agl
		if(cmask.isFromAgl()){
			mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
			if(RUsed>AglDB.MomBins()[m]&&RUsed<=AglDB.MomBins()[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5) FitAgl_Dbins -> TemplateP -> Fill(mass,m);
				if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitAgl_Dbins -> TemplateD) -> Fill(mass,m,ReturnMCGenType());
				if(Massa_gen<4&&Massa_gen>2.5) FitAgl_Dbins -> TemplateHe-> Fill(mass,m);
			}
			if(RUsed>AglPB.MomBins()[m]&&RUsed<=AglPB.MomBins()[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5) FitAgl_Pbins -> TemplateP -> Fill(mass,m);
				if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitAgl_Pbins -> TemplateD) -> Fill(mass,m,ReturnMCGenType());
				if(Massa_gen<4&&Massa_gen>2.5) FitAgl_Pbins -> TemplateHe-> Fill(mass,m);
			}
		}
	}
} 


void DeutonsDATA_Fill(TNtuple *ntupla, int l,int zona){
	 ntupla->GetEvent(l);
	if(Tup.Beta<=0||Tup.R<=0) return;
	float mass = 0;
	//cuts
	if(!(Likcut&&Distcut)) return;
	//
	for(int m=0;m<nbinsToF;m++){ //TOF
		mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));
		if(RUsed>ToFDB.MomBins()[m]&&RUsed<=ToFDB.MomBins()[m+1]){
			if(Tup.R>1.2*Tup.Rcutoff) FitTOF_Dbins -> DATA -> Fill(mass,m);
			((TH3*)FitTOFgeo_Dbins -> DATA) -> Fill(mass,m,zona);
		}
		if(RUsed>ToFPB.MomBins()[m]&&RUsed<=ToFPB.MomBins()[m+1]) {
			if(Tup.R>1.2*Tup.Rcutoff) FitTOF_Pbins -> DATA -> Fill(mass,m);
		}
	}
	for(int m=0;m<nbinsNaF;m++){//NaF
		if(cmask.isFromNaF()){
			mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
			if(RUsed>NaFDB.MomBins()[m]&&RUsed<=NaFDB.MomBins()[m+1]) {
				if(Tup.R>1.2*Tup.Rcutoff) FitNaF_Dbins -> DATA -> Fill(mass,m);
				((TH3*)FitNaFgeo_Dbins -> DATA) -> Fill(mass,m,zona);
			}
			if(RUsed>NaFPB.MomBins()[m]&&RUsed<=NaFPB.MomBins()[m+1]) {
				if(Tup.R>1.2*Tup.Rcutoff) FitNaF_Pbins -> DATA -> Fill(mass,m);
			}
		}
	}
	for(int m=0;m<nbinsAgl;m++){ //Agl
		if(cmask.isFromAgl()){
			mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
			if(RUsed>AglDB.MomBins()[m]&&RUsed<=AglDB.MomBins()[m+1]) {
				if(Tup.R>1.2*Tup.Rcutoff) FitAgl_Dbins -> DATA -> Fill(mass,m);
				((TH3*)FitAglgeo_Dbins -> DATA) -> Fill(mass,m,zona);
			}
			if(RUsed>AglPB.MomBins()[m]&&RUsed<=AglPB.MomBins()[m+1]) {
				if(Tup.R>1.2*Tup.Rcutoff) FitAgl_Pbins -> DATA -> Fill(mass,m);
			}
		}
	}
	return;
}


void DeutonsMC_Write(){

	FitTOF_Dbins -> Write();
	FitNaF_Dbins -> Write();
	FitAgl_Dbins -> Write();

	FitTOFgeo_Dbins -> Write();
	FitNaFgeo_Dbins -> Write();
	FitAglgeo_Dbins -> Write();

	
	FitTOF_Pbins -> Write();
	FitNaF_Pbins -> Write();
	FitAgl_Pbins -> Write();
	return;
}


void DeutonsTemplFits(){
	string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	TFile * file1 = TFile::Open(nomefile.c_str(),"READ");

	TemplateFIT * FitTOF_Dbins	= new TemplateFIT(file1,"FitTOF_Dbins","FitTOF_Dbins");
	TemplateFIT * FitNaF_Dbins	= new TemplateFIT(file1,"FitNaF_Dbins","FitNaF_Dbins");
	TemplateFIT * FitAgl_Dbins	= new TemplateFIT(file1,"FitAgl_Dbins","FitAgl_Dbins");
                                                                                              
	TemplateFIT * FitTOFgeo_Dbins	= new TemplateFIT(file1,"FitTOF_Dbins","FitTOFgeo_Dbins",11);
	TemplateFIT * FitNaFgeo_Dbins	= new TemplateFIT(file1,"FitNaF_Dbins","FitNaFgeo_Dbins",11);
	TemplateFIT * FitAglgeo_Dbins	= new TemplateFIT(file1,"FitAgl_Dbins","FitAglgeo_Dbins",11);
                                                                                              
	TemplateFIT * FitTOF_Pbins	= new TemplateFIT(file1,"FitTOF_Pbins","FitTOF_Pbins");
	TemplateFIT * FitNaF_Pbins	= new TemplateFIT(file1,"FitNaF_Pbins","FitNaF_Pbins");
	TemplateFIT * FitAgl_Pbins	= new TemplateFIT(file1,"FitAgl_Pbins","FitAgl_Pbins");

	cout<<"******************** DEUTONS TEMPlATE FITS ************************"<<endl;

	FitTOF_Dbins 	-> SetFitConstraints(0.8,1,0.00,0.2,0.0001,0.0025); 
	FitNaF_Dbins 	-> SetFitConstraints(0.8,1,0.00,0.2,0.0001,0.0015);
	FitAgl_Dbins 	-> SetFitConstraints(0.8,1,0.00,0.2,0.0001,0.0005);

	FitTOFgeo_Dbins -> SetFitConstraints(0.8,1,0.00,0.2,0.0001,0.0025);
	FitNaFgeo_Dbins -> SetFitConstraints(0.8,1,0.00,0.2,0.0001,0.0015);
	FitAglgeo_Dbins -> SetFitConstraints(0.8,1,0.00,0.2,0.0001,0.0005);

	FitTOF_Pbins 	-> SetFitConstraints(0.8,1,0.00001,0.02,0.0001,0.0025); 
	FitNaF_Pbins 	-> SetFitConstraints(0.8,1,0.00001,0.02,0.0005,0.0015);
	FitAgl_Pbins 	-> SetFitConstraints(0.8,1,0.00001,0.02,0.0001,0.0005);

	cout<<"** TOF **"<<endl;
	FitTOF_Dbins 	-> TemplateFits(); 
	FitTOFgeo_Dbins -> TemplateFits();
	FitTOF_Pbins    -> TemplateFits();

	cout<<"** NaF **"<<endl;	
	FitNaF_Dbins 	-> TemplateFits();
	FitNaFgeo_Dbins -> TemplateFits();
	FitNaF_Pbins    -> TemplateFits();
	
	cout<<"** Agl **"<<endl;
	FitAgl_Dbins 	-> TemplateFits();
	FitAglgeo_Dbins -> TemplateFits();
	FitAgl_Pbins 	-> TemplateFits();

	cout<<"***** TemplateFits Outcome (Mass) ******"<<endl;
	cout<<"** TOF **"<<endl; 
	for(int bin =0; bin <nbinsToF;bin++)
		cout<<FitTOF_Dbins->GetFitOutcome(bin)<<" ";
	cout<<endl;
	for(int bin =0; bin <nbinsToF;bin++)
                cout<<FitTOF_Pbins->GetFitOutcome(bin)<<" ";
        cout<<endl;
	cout<<"** NaF **"<<endl;
        for(int bin =0; bin <nbinsNaF;bin++)
                cout<<FitNaF_Dbins->GetFitOutcome(bin)<<" ";
        cout<<endl;
	for(int bin =0; bin <nbinsNaF;bin++)
                cout<<FitNaF_Pbins->GetFitOutcome(bin)<<" ";
        cout<<endl;
	cout<<"** Agl **"<<endl;
        for(int bin =0; bin <nbinsAgl;bin++)
                cout<<FitAgl_Dbins->GetFitOutcome(bin)<<" ";
        cout<<endl;
	for(int bin =0; bin <nbinsAgl;bin++) 
		cout<<FitAgl_Pbins->GetFitOutcome(bin)<<" ";
        cout<<endl;

	cout<<"*** Updating P1 file ****"<<endl;

        nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
        file1 = TFile::Open(nomefile.c_str(),"UPDATE");

	FitTOF_Dbins -> DCounts -> Write ("D_FluxCounts_TOF");
	FitNaF_Dbins -> DCounts -> Write ("D_FluxCounts_NaF");
        FitAgl_Dbins -> DCounts -> Write ("D_FluxCounts_Agl");
        
	FitTOFgeo_Dbins -> DCounts -> Write ("D_Flux_geoCounts_TOF");
        FitNaFgeo_Dbins -> DCounts -> Write ("D_Flux_geoCounts_NaF");
        FitAglgeo_Dbins -> DCounts -> Write ("D_Flux_geoCounts_Agl");
               
        FitTOF_Pbins -> PCounts -> Write ("P_FluxCounts_TOF");
	FitNaF_Pbins -> PCounts -> Write ("P_FluxCounts_NaF");
        FitAgl_Pbins -> PCounts -> Write ("P_FluxCounts_Agl");

	file1 -> Write();
	file1 -> Close(); 	
	
	TCanvas *c30_TOF[2][nbinsToF];
	TCanvas *c30_NaF[2][nbinsNaF];
	TCanvas *c30_Agl[2][nbinsAgl];
	
	TCanvas *c30_TOFgeo[11];
	TCanvas *c30_NaFgeo[11];
        TCanvas *c30_Aglgeo[11];
	
	//Primaries
	for(int bin=0; bin <nbinsToF; bin++){
		c30_TOF[0][bin] = new TCanvas(("TOF bin:" + to_string(bin)).c_str());
		c30_TOF[0][bin]->cd();
		FitTOF_Dbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin);
	}
	for(int bin=0; bin <nbinsNaF; bin++){
		c30_NaF[0][bin] = new TCanvas(("NaF bin:" + to_string(bin)).c_str());
		c30_NaF[0][bin]->cd();
		FitNaF_Dbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin);
	}
	for(int bin=0; bin <nbinsAgl; bin++){
		c30_Agl[0][bin] = new TCanvas(("Agl bin:" + to_string(bin)).c_str());
		c30_Agl[0][bin]->cd();
		FitAgl_Dbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin);
	}
	
	for(int bin=0; bin <nbinsToF; bin++){
		c30_TOF[1][bin] = new TCanvas(("TOF P bin:" + to_string(bin)).c_str());
		c30_TOF[1][bin]->cd();
		FitTOF_Pbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin);
	}
	for(int bin=0; bin <nbinsNaF; bin++){
		c30_NaF[1][bin] = new TCanvas(("NaF P bin:" + to_string(bin)).c_str());
		c30_NaF[1][bin]->cd();
		FitNaF_Pbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin);
	}
	for(int bin=0; bin <nbinsAgl; bin++){
		c30_Agl[1][bin] = new TCanvas(("Agl P bin:" + to_string(bin)).c_str());
		c30_Agl[1][bin]->cd();
		FitAgl_Pbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin);
	}
	//Geo. Zones
	for(int lat=1;lat<11;lat++){
		

		c30_TOFgeo[lat] = new TCanvas(("TOF bins - Latitude: " + to_string(lat)).c_str());
		c30_NaFgeo[lat] = new TCanvas(("NaF bins - Latitude: " + to_string(lat)).c_str());
		c30_Aglgeo[lat] = new TCanvas(("Agl bins - Latitude: " + to_string(lat)).c_str());
	
		c30_TOFgeo[lat] -> Divide(6,3);
        	c30_NaFgeo[lat] -> Divide(6,3);
        	c30_Aglgeo[lat] -> Divide(6,3);	
		
	
		for(int bin=0; bin <nbinsToF; bin++){
			c30_TOFgeo[lat] -> cd (bin +1);
			FitTOFgeo_Dbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin,lat);
		}
		for(int bin=0; bin <nbinsNaF; bin++){
			c30_NaFgeo[lat] -> cd (bin +1);
			FitNaFgeo_Dbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin,lat);
		}
		for(int bin=0; bin <nbinsAgl; bin++){
			c30_Aglgeo[lat] -> cd (bin +1);
			FitAglgeo_Dbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin,lat);
		}
	}

	cout<<"*** Updating Results file ***"<<endl;
        nomefile="./Final_plots/"+mese+".root";
        TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
	//Primaries
	f_out->mkdir("Mass Template Fits/TOF/TOF Primaries/Dbins");
	f_out->cd("Mass Template Fits/TOF/TOF Primaries/Dbins");
	for(int bin=0; bin <nbinsToF; bin++)
		c30_TOF[0][bin]->Write();
	f_out->mkdir("Mass Template Fits/NaF/NaF Primaries/Dbins");
	f_out->cd("Mass Template Fits/NaF/NaF Primaries/Dbins");
	for(int bin=0; bin <nbinsNaF; bin++)
		c30_NaF[0][bin]->Write();
	f_out->mkdir("Mass Template Fits/Agl/Agl Primaries/Dbins");
	f_out->cd("Mass Template Fits/Agl/Agl Primaries/Dbins");
	for(int bin=0; bin <nbinsAgl; bin++)
		c30_Agl[0][bin]->Write();

	f_out->mkdir("Mass Template Fits/TOF/TOF Primaries/Pbins");
	f_out->cd("Mass Template Fits/TOF/TOF Primaries/Pbins");
	for(int bin=0; bin <nbinsToF; bin++)
		c30_TOF[1][bin]->Write();
	f_out->mkdir("Mass Template Fits/NaF/NaF Primaries/Pbins");
	f_out->cd("Mass Template Fits/NaF/NaF Primaries/Pbins");
	for(int bin=0; bin <nbinsNaF; bin++)
		c30_NaF[1][bin]->Write();
	f_out->mkdir("Mass Template Fits/Agl/Agl Primaries/Pbins");
	f_out->cd("Mass Template Fits/Agl/Agl Primaries/Pbins");
	for(int bin=0; bin <nbinsAgl; bin++)
		c30_Agl[1][bin]->Write();
	//Geom. Zones
	for(int lat=1;lat<11;lat++){
		f_out->mkdir("Mass Template Fits/TOF/TOF Geom. Zones");
		f_out->cd("Mass Template Fits/TOF/TOF Geom. Zones");
		c30_TOFgeo[lat] -> Write(("TOF Geo. Zone: " + to_string(lat)).c_str());
		f_out->mkdir("Mass Template Fits/NaF/NaF Geom. Zones ");
		f_out->cd("Mass Template Fits/NaF/NaF Geom. Zones ");
		c30_NaFgeo[lat] -> Write(("NaF Geo. Zone: " + to_string(lat)).c_str());
		f_out->mkdir("Mass Template Fits/Agl/Agl Geom. Zones  ");
		f_out->cd("Mass Template Fits/Agl/Agl Geom. Zones  ");
		c30_Aglgeo[lat] -> Write(("Agl Geo. Zone: " + to_string(lat)).c_str());
	}
	f_out->Write();
	f_out->Close();


	return;


}
