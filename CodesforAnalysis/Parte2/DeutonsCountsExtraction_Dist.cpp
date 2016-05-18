TemplateFIT * FitTOF_Dbins_Dist  = new TemplateFIT("FitTOF_Dbins_Dist ",nbinsToF,-1.2,1.2,6);
TemplateFIT * FitNaF_Dbins_Dist  = new TemplateFIT("FitNaF_Dbins_Dist ",nbinsNaF,-1.2,1.2,6);
TemplateFIT * FitAgl_Dbins_Dist  = new TemplateFIT("FitAgl_Dbins_Dist ",nbinsAgl,-1.2,1.2,6);

TemplateFIT * FitTOFgeo_Dbins_Dist  = new TemplateFIT("FitTOFgeo_Dbins_Dist ",nbinsToF,-1.2,1.2,6,11);
TemplateFIT * FitNaFgeo_Dbins_Dist  = new TemplateFIT("FitNaFgeo_Dbins_Dist ",nbinsNaF,-1.2,1.2,6,11);
TemplateFIT * FitAglgeo_Dbins_Dist  = new TemplateFIT("FitAglgeo_Dbins_Dist ",nbinsAgl,-1.2,1.2,6,11);

TemplateFIT * FitTOF_Pbins_Dist  = new TemplateFIT("FitTOF_Pbins_Dist ",nbinsToF,-1.2,1.2,6);
TemplateFIT * FitNaF_Pbins_Dist  = new TemplateFIT("FitNaF_Pbins_Dist ",nbinsNaF,-1.2,1.2,6);
TemplateFIT * FitAgl_Pbins_Dist  = new TemplateFIT("FitAgl_Pbins_Dist ",nbinsAgl,-1.2,1.2,6);


void DeutonsMC_Dist_Fill(TNtuple *ntupla, int l){
	int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	float Distance_Discr = 0;
	if(!(Likcut&&Distcut)) return;
	for(int m=0;m<nbinsToF;m++){ //TOF
		Distance_Discr = ((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
		if(Var>BetaD[m]&&Var<=BetaD[m+1]){
			if(Massa_gen<1&&Massa_gen>0.5) FitTOF_Dbins_Dist -> TemplateP -> Fill(Distance_Discr,m);
			if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitTOF_Dbins_Dist -> TemplateD) -> Fill(Distance_Discr,m,ReturnMCGenType());
			if(Massa_gen<4&&Massa_gen>2.5) FitTOF_Dbins_Dist -> TemplateHe-> Fill(Distance_Discr,m);
		}
		if(Var>BetaP[m]&&Var<=BetaP[m+1]) {
			if(Massa_gen<1&&Massa_gen>0.5) FitTOF_Pbins_Dist -> TemplateP -> Fill(Distance_Discr,m);
			if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitTOF_Pbins_Dist -> TemplateD) -> Fill(Distance_Discr,m,ReturnMCGenType());
			if(Massa_gen<4&&Massa_gen>2.5) FitTOF_Pbins_Dist -> TemplateHe-> Fill(Distance_Discr,m);
		}
	}
	for(int m=0;m<nbinsNaF;m++) { //NaF
		if((((int)Cutmask)>>11)==512){
			Distance_Discr = ((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
			if(Var2>NaFDB.MomBins()[m]&&Var2<=NaFDB.MomBins()[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5) FitNaF_Dbins_Dist -> TemplateP -> Fill(Distance_Discr,m);
				if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitNaF_Dbins_Dist -> TemplateD) -> Fill(Distance_Discr,m,ReturnMCGenType());
				if(Massa_gen<4&&Massa_gen>2.5) FitNaF_Dbins_Dist -> TemplateHe-> Fill(Distance_Discr,m);
			}
			if(Var2>BetaNaFP[m]&&Var2<=BetaNaFP[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5) FitNaF_Pbins_Dist -> TemplateP -> Fill(Distance_Discr,m);
				if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitNaF_Pbins_Dist -> TemplateD) -> Fill(Distance_Discr,m,ReturnMCGenType());
				if(Massa_gen<4&&Massa_gen>2.5) FitNaF_Pbins_Dist -> TemplateHe-> Fill(Distance_Discr,m);
			}
		}
	}
	for(int m=0;m<nbinsAgl;m++) { //Agl
		if((((int)Cutmask)>>11)==0){
			Distance_Discr = ((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
			if(Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5) FitAgl_Dbins_Dist -> TemplateP -> Fill(Distance_Discr,m);
				if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitAgl_Dbins_Dist -> TemplateD) -> Fill(Distance_Discr,m,ReturnMCGenType());
				if(Massa_gen<4&&Massa_gen>2.5) FitAgl_Dbins_Dist -> TemplateHe-> Fill(Distance_Discr,m);
			}
			if(Var2>BetaAglP[m]&&Var2<=BetaAglP[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5) FitAgl_Pbins_Dist -> TemplateP -> Fill(Distance_Discr,m);
				if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitAgl_Pbins_Dist -> TemplateD) -> Fill(Distance_Discr,m,ReturnMCGenType());
				if(Massa_gen<4&&Massa_gen>2.5) FitAgl_Pbins_Dist -> TemplateHe-> Fill(Distance_Discr,m);
			}
		}
	}
} 


void DeutonsDATA_Dist_Fill(TNtuple *ntupla, int l,int zona){
	int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	float Distance_Discr = 0;
	if(!(Likcut&&Distcut)) return;
	for(int m=0;m<nbinsToF;m++){ //TOF
		Distance_Discr = ((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
		if(Var>BetaD[m]&&Var<=BetaD[m+1]){
			if(R>1.2*Rcutoff) FitTOF_Dbins_Dist -> DATA -> Fill(Distance_Discr,m);
			((TH3*)FitTOFgeo_Dbins_Dist -> DATA) -> Fill(Distance_Discr,m,zona);
		}
		if(Var>BetaP[m]&&Var<=BetaP[m+1]) {
			if(R>1.2*Rcutoff) FitTOF_Pbins_Dist -> DATA -> Fill(Distance_Discr,m);
		}
	}
	for(int m=0;m<nbinsNaF;m++){//NaF
		if((((int)Cutmask)>>11)==512){
			Distance_Discr =  ((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
			if(Var2>NaFDB.MomBins()[m]&&Var2<=NaFDB.MomBins()[m+1]) {
				if(R>1.2*Rcutoff) FitNaF_Dbins_Dist -> DATA -> Fill(Distance_Discr,m);
				((TH3*)FitNaFgeo_Dbins_Dist -> DATA) -> Fill(Distance_Discr,m,zona);
			}
			if(Var2>BetaNaFP[m]&&Var2<=BetaNaFP[m+1]) {
				if(R>1.2*Rcutoff) FitNaF_Pbins_Dist -> DATA -> Fill(Distance_Discr,m);
			}
		}
	}
	for(int m=0;m<nbinsAgl;m++){ //Agl
		if((((int)Cutmask)>>11)==0){
			Distance_Discr =  ((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
			if(Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) {
				if(R>1.2*Rcutoff) FitAgl_Dbins_Dist -> DATA -> Fill(Distance_Discr,m);
				((TH3*)FitAglgeo_Dbins_Dist -> DATA) -> Fill(Distance_Discr,m,zona);
			}
			if(Var2>BetaAglP[m]&&Var2<=BetaAglP[m+1]) {
				if(R>1.2*Rcutoff) FitAgl_Pbins_Dist -> DATA -> Fill(Distance_Discr,m);
			}
		}
	}
	return;
}


void DeutonsMC_Dist_Write(){

	FitTOF_Dbins_Dist -> Write();
	FitNaF_Dbins_Dist -> Write();
	FitAgl_Dbins_Dist -> Write();

	FitTOFgeo_Dbins_Dist -> Write();
	FitNaFgeo_Dbins_Dist -> Write();
	FitAglgeo_Dbins_Dist -> Write();

	
	FitTOF_Pbins_Dist -> Write();
	FitNaF_Pbins_Dist -> Write();
	FitAgl_Pbins_Dist -> Write();
	return;
}


void DeutonsTemplFits_Dist(){
	string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	TFile * file1 = TFile::Open(nomefile.c_str(),"READ");
	
	TemplateFIT * FitTOF_Dbins_Dist	= new TemplateFIT(file1,"FitTOF_Dbins_Dist ","FitTOF_Dbins_Dist ",0,3);
	TemplateFIT * FitNaF_Dbins_Dist	= new TemplateFIT(file1,"FitNaF_Dbins_Dist ","FitNaF_Dbins_Dist ",0,3);
	TemplateFIT * FitAgl_Dbins_Dist	= new TemplateFIT(file1,"FitAgl_Dbins_Dist ","FitAgl_Dbins_Dist ",0,3);
                                                                                              
	TemplateFIT * FitTOFgeo_Dbins_Dist	= new TemplateFIT(file1,"FitTOF_Dbins_Dist ","FitTOFgeo_Dbins_Dist ",0,3,11);
	TemplateFIT * FitNaFgeo_Dbins_Dist	= new TemplateFIT(file1,"FitNaF_Dbins_Dist ","FitNaFgeo_Dbins_Dist ",0,3,11);
	TemplateFIT * FitAglgeo_Dbins_Dist	= new TemplateFIT(file1,"FitAgl_Dbins_Dist ","FitAglgeo_Dbins_Dist ",0,3,11);
                                                                                              
	TemplateFIT * FitTOF_Pbins_Dist	= new TemplateFIT(file1,"FitTOF_Pbins_Dist ","FitTOF_Pbins_Dist ",0,3);
	TemplateFIT * FitNaF_Pbins_Dist	= new TemplateFIT(file1,"FitNaF_Pbins_Dist ","FitNaF_Pbins_Dist ",0,3);
	TemplateFIT * FitAgl_Pbins_Dist	= new TemplateFIT(file1,"FitAgl_Pbins_Dist ","FitAgl_Pbins_Dist ",0,3);

	cout<<"******************** DEUTONS DISTANCE TEMPlATE FITS ************************"<<endl;

	FitTOF_Dbins_Dist 	-> SetFitConstraints(0.5,1,0.005,0.1,0.0001,0.025); 
	FitNaF_Dbins_Dist 	-> SetFitConstraints(0.8,1,0.0001,0.1,0.0005,0.0015);
	FitAgl_Dbins_Dist 	-> SetFitConstraints(0.8,1,0.0001,0.1,0.0001,0.0005);
                                
	FitTOFgeo_Dbins_Dist 	-> DisableFit();
	FitNaFgeo_Dbins_Dist 	-> DisableFit();
	FitAglgeo_Dbins_Dist 	-> DisableFit();
                                
	FitTOF_Pbins_Dist 	-> SetFitConstraints(0.5,1,0.005,0.1,0.0001,0.025); 
	FitNaF_Pbins_Dist 	-> SetFitConstraints(0.8,1,0.0001,0.1,0.0005,0.0015);	
	FitAgl_Pbins_Dist 	-> SetFitConstraints(0.8,1,0.0001,0.1,0.0001,0.0005);

	cout<<"** TOF **"<<endl;
        FitTOF_Dbins_Dist    -> TemplateFits();
        FitTOFgeo_Dbins_Dist -> TemplateFits();
        FitTOF_Pbins_Dist    -> TemplateFits();

        cout<<"** NaF **"<<endl;
        FitNaF_Dbins_Dist    -> TemplateFits();
        FitNaFgeo_Dbins_Dist -> TemplateFits();
        FitNaF_Pbins_Dist    -> TemplateFits();

        cout<<"** Agl **"<<endl;
        FitAgl_Dbins_Dist    -> TemplateFits();
        FitAglgeo_Dbins_Dist -> TemplateFits();
        FitAgl_Pbins_Dist    -> TemplateFits();

	cout<<"***** TemplateFits Outcome (Dist) ******"<<endl;
        cout<<"** TOF **"<<endl;
        for(int bin =0; bin <nbinsToF;bin++)
                cout<<FitTOF_Dbins_Dist->GetFitOutcome(bin)<<" ";
        cout<<endl;
	for(int bin =0; bin <nbinsToF;bin++)
                cout<<FitTOF_Pbins_Dist->GetFitOutcome(bin)<<" ";
        cout<<endl;

	cout<<"** NaF **"<<endl;
        for(int bin =0; bin <nbinsNaF;bin++)
                cout<<FitNaF_Dbins_Dist->GetFitOutcome(bin)<<" ";
        cout<<endl;
	 for(int bin =0; bin <nbinsNaF;bin++)
                cout<<FitNaF_Pbins_Dist->GetFitOutcome(bin)<<" ";
        cout<<endl;
	cout<<"** Agl **"<<endl;
        for(int bin =0; bin <nbinsAgl;bin++)
                cout<<FitAgl_Dbins_Dist->GetFitOutcome(bin)<<" ";
        cout<<endl;
	for(int bin =0; bin <nbinsAgl;bin++)
                cout<<FitAgl_Pbins_Dist->GetFitOutcome(bin)<<" ";
        cout<<endl;

	cout<<"*** Updating P1 file ****"<<endl;

        nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
        file1 = TFile::Open(nomefile.c_str(),"UPDATE");

	FitTOF_Dbins_Dist -> DCounts -> Write ("D_Flux_DistCounts_TOF");
	FitNaF_Dbins_Dist -> DCounts -> Write ("D_Flux_DistCounts_NaF");
        FitAgl_Dbins_Dist -> DCounts -> Write ("D_Flux_DistCounts_Agl");
                       
        FitTOFgeo_Dbins_Dist -> DCounts -> Write ("D_Flux_geo_DistCounts_TOF");
        FitNaFgeo_Dbins_Dist -> DCounts -> Write ("D_Flux_geo_DistCounts_NaF");
        FitAglgeo_Dbins_Dist -> DCounts -> Write ("D_Flux_geo_DistCounts_Agl");
	
	FitTOF_Pbins_Dist -> PCounts -> Write ("P_Flux_DistCounts_TOF");
	FitNaF_Pbins_Dist -> PCounts -> Write ("P_Flux_DistCounts_NaF");
        FitAgl_Pbins_Dist -> PCounts -> Write ("P_Flux_DistCounts_Agl");

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
		c30_TOF[0][bin] = new TCanvas(("TOF bin:" + to_string(bin) + "(Dist.)").c_str());
		c30_TOF[0][bin]->cd();
		FitTOF_Dbins_Dist -> TemplateFitPlot(gPad,"Distance Discr.",bin);
	}
	for(int bin=0; bin <nbinsNaF; bin++){
		c30_NaF[0][bin] = new TCanvas(("NaF bin:" + to_string(bin) + "(Dist.)").c_str());
		c30_NaF[0][bin]->cd();
		FitNaF_Dbins_Dist -> TemplateFitPlot(gPad,"Distance Discr.",bin);
	}
	for(int bin=0; bin <nbinsAgl; bin++){
		c30_Agl[0][bin] = new TCanvas(("Agl bin:" + to_string(bin) + "(Dist.)").c_str());
		c30_Agl[0][bin]->cd();
		FitAgl_Dbins_Dist -> TemplateFitPlot(gPad,"Distance Discr.",bin);
	}
	
	for(int bin=0; bin <nbinsToF; bin++){
		c30_TOF[1][bin] = new TCanvas(("TOF P bin:" + to_string(bin) + "(Dist.)").c_str());
		c30_TOF[1][bin]->cd();
		FitTOF_Pbins_Dist -> TemplateFitPlot(gPad,"Distance Discr.",bin);
	}
	for(int bin=0; bin <nbinsNaF; bin++){
		c30_NaF[1][bin] = new TCanvas(("NaF P bin:" + to_string(bin) + "(Dist.)").c_str());
		c30_NaF[1][bin]->cd();
		FitNaF_Pbins_Dist -> TemplateFitPlot(gPad,"Distance Discr.",bin);
	}
	for(int bin=0; bin <nbinsAgl; bin++){
		c30_Agl[1][bin] = new TCanvas(("Agl P bin:" + to_string(bin) + "(Dist.)").c_str());
		c30_Agl[1][bin]->cd();
		FitAgl_Pbins_Dist -> TemplateFitPlot(gPad,"Distance Discr.",bin);
	}
	//Geo. Zones
	for(int lat=1;lat<11;lat++){
		

		c30_TOFgeo[lat] = new TCanvas(("TOF bins - Latitude: " + to_string(lat) + "(Dist.)").c_str());
		c30_NaFgeo[lat] = new TCanvas(("NaF bins - Latitude: " + to_string(lat) + "(Dist.)").c_str());
		c30_Aglgeo[lat] = new TCanvas(("Agl bins - Latitude: " + to_string(lat) + "(Dist.)").c_str());
	
		c30_TOFgeo[lat] -> Divide(6,3);
        	c30_NaFgeo[lat] -> Divide(6,3);
        	c30_Aglgeo[lat] -> Divide(6,3);	
		
	
		for(int bin=0; bin <nbinsToF; bin++){
			c30_TOFgeo[lat] -> cd (bin +1);
			FitTOFgeo_Dbins_Dist -> TemplateFitPlot(gPad,"Distance Discr.",bin,lat);
		}
		for(int bin=0; bin <nbinsNaF; bin++){
			c30_NaFgeo[lat] -> cd (bin +1);
			FitNaFgeo_Dbins_Dist -> TemplateFitPlot(gPad,"Distance Discr.",bin,lat);
		}
		for(int bin=0; bin <nbinsAgl; bin++){
			c30_Aglgeo[lat] -> cd (bin +1);
			FitAglgeo_Dbins_Dist -> TemplateFitPlot(gPad,"Distance Discr.",bin,lat);
		}
	}

	cout<<"*** Updating Results file ***"<<endl;
        nomefile="./Final_plots/"+mese+".root";
        TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
	//Primaries
	f_out->mkdir("Distance Template Fits/TOF/TOF Primaries/Dbins");
	f_out->cd("Distance Template Fits/TOF/TOF Primaries/Dbins");
	for(int bin=0; bin <nbinsToF; bin++)
		c30_TOF[0][bin]->Write();
	f_out->mkdir("Distance Template Fits/NaF/NaF Primaries/Dbins");
	f_out->cd("Distance Template Fits/NaF/NaF Primaries/Dbins");
	for(int bin=0; bin <nbinsNaF; bin++)
		c30_NaF[0][bin]->Write();
	f_out->mkdir("Distance Template Fits/Agl/Agl Primaries/Dbins");
	f_out->cd("Distance Template Fits/Agl/Agl Primaries/Dbins");
	for(int bin=0; bin <nbinsAgl; bin++)
		c30_Agl[0][bin]->Write();

	f_out->mkdir("Distance Template Fits/TOF/TOF Primaries/Pbins");
	f_out->cd("Distance Template Fits/TOF/TOF Primaries/Pbins");
	for(int bin=0; bin <nbinsToF; bin++)
		c30_TOF[1][bin]->Write();
	f_out->mkdir("Distance Template Fits/NaF/NaF Primaries/Pbins");
	f_out->cd("Distance Template Fits/NaF/NaF Primaries/Pbins");
	for(int bin=0; bin <nbinsNaF; bin++)
		c30_NaF[1][bin]->Write();
	f_out->mkdir("Distance Template Fits/Agl/Agl Primaries/Pbins");
	f_out->cd("Distance Template Fits/Agl/Agl Primaries/Pbins");
	for(int bin=0; bin <nbinsAgl; bin++)
		c30_Agl[1][bin]->Write();
	//Geom. Zones
	for(int lat=1;lat<11;lat++){
		f_out->mkdir("Distance Template Fits/TOF/TOF Geom. Zones");
		f_out->cd("Distance Template Fits/TOF/TOF Geom. Zones");
		c30_TOFgeo[lat] -> Write(("TOF Geo. Zone: " + to_string(lat)).c_str());
		f_out->mkdir("Distance Template Fits/NaF/NaF Geom. Zones ");
		f_out->cd("Distance Template Fits/NaF/NaF Geom. Zones ");
		c30_NaFgeo[lat] -> Write(("NaF Geo. Zone: " + to_string(lat)).c_str());
		f_out->mkdir("Distance Template Fits/Agl/Agl Geom. Zones  ");
		f_out->cd("Distance Template Fits/Agl/Agl Geom. Zones  ");
		c30_Aglgeo[lat] -> Write(("Agl Geo. Zone: " + to_string(lat)).c_str());
	}
	f_out->Write();
	f_out->Close();


	return;


}
