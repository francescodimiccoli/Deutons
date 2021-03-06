

TemplateFIT * FitTOF_Dbins_Dist  = new TemplateFIT("FitTOF_Dbins_Dist ",nbinsToF,-1.2,1.2,6);
TemplateFIT * FitNaF_Dbins_Dist  = new TemplateFIT("FitNaF_Dbins_Dist ",nbinsNaF,-1.2,1.2,6);
TemplateFIT * FitAgl_Dbins_Dist  = new TemplateFIT("FitAgl_Dbins_Dist ",nbinsAgl,-1.2,1.2,6);

TemplateFIT * FitTOFgeo_Dbins_Dist  = new TemplateFIT("FitTOFgeo_Dbins_Dist ",nbinsToF,-1.2,1.2,6,11);
TemplateFIT * FitNaFgeo_Dbins_Dist  = new TemplateFIT("FitNaFgeo_Dbins_Dist ",nbinsNaF,-1.2,1.2,6,11);
TemplateFIT * FitAglgeo_Dbins_Dist  = new TemplateFIT("FitAglgeo_Dbins_Dist ",nbinsAgl,-1.2,1.2,6,11);

TemplateFIT * FitTOF_Pbins_Dist  = new TemplateFIT("FitTOF_Pbins_Dist ",nbinsToF,-1.2,1.2,6);
TemplateFIT * FitNaF_Pbins_Dist  = new TemplateFIT("FitNaF_Pbins_Dist ",nbinsNaF,-1.2,1.2,6);
TemplateFIT * FitAgl_Pbins_Dist  = new TemplateFIT("FitAgl_Pbins_Dist ",nbinsAgl,-1.2,1.2,6);


void DeutonsMC_Dist_Fill()
{

	float Distance_Discr = 0;
	int Kbin;
	if(!(Likcut&&Distcut)) return;
		Distance_Discr = ((Tup.Dist5D_P-Tup.Dist5D)/(Tup.Dist5D_P+Tup.Dist5D));
		Kbin=ToFDB.GetBin(RUsed);
			if(Massa_gen<1&&Massa_gen>0.5) ((TH2*)FitTOF_Dbins_Dist -> TemplateP) -> Fill(Distance_Discr,Kbin,Tup.mcweight);
			if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitTOF_Dbins_Dist -> TemplateD) -> Fill(Distance_Discr,Kbin,ReturnMCGenType(),Tup.mcweight);
			if(Massa_gen<4&&Massa_gen>2.5) ((TH2*)FitTOF_Dbins_Dist -> TemplateHe)-> Fill(Distance_Discr,Kbin,Tup.mcweight);
		Kbin=ToFPB.GetBin(RUsed);	
			if(Massa_gen<1&&Massa_gen>0.5) ((TH2*)FitTOF_Pbins_Dist -> TemplateP) -> Fill(Distance_Discr,Kbin,Tup.mcweight);
			if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitTOF_Pbins_Dist -> TemplateD) -> Fill(Distance_Discr,Kbin,ReturnMCGenType(),Tup.mcweight);
			if(Massa_gen<4&&Massa_gen>2.5) ((TH2*)FitTOF_Pbins_Dist -> TemplateHe)-> Fill(Distance_Discr,Kbin,Tup.mcweight);
	
		if(cmask.isFromNaF()) {
			Distance_Discr = ((Tup.Dist5D_P-Tup.Dist5D)/(Tup.Dist5D_P+Tup.Dist5D));
			Kbin=NaFDB.GetBin(RUsed);	
				if(Massa_gen<1&&Massa_gen>0.5) ((TH2*)FitNaF_Dbins_Dist -> TemplateP) -> Fill(Distance_Discr,Kbin,Tup.mcweight);
				if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitNaF_Dbins_Dist -> TemplateD) -> Fill(Distance_Discr,Kbin,ReturnMCGenType(),Tup.mcweight);
				if(Massa_gen<4&&Massa_gen>2.5) ((TH2*)FitNaF_Dbins_Dist -> TemplateHe)-> Fill(Distance_Discr,Kbin,Tup.mcweight);
			Kbin=NaFPB.GetBin(RUsed);		
				if(Massa_gen<1&&Massa_gen>0.5) ((TH2*)FitNaF_Pbins_Dist -> TemplateP) -> Fill(Distance_Discr,Kbin,Tup.mcweight);
				if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitNaF_Pbins_Dist -> TemplateD) -> Fill(Distance_Discr,Kbin,ReturnMCGenType(),Tup.mcweight);
				if(Massa_gen<4&&Massa_gen>2.5) ((TH2*)FitNaF_Pbins_Dist -> TemplateHe)-> Fill(Distance_Discr,Kbin,Tup.mcweight);
		}
		if(cmask.isFromAgl()) {
			Distance_Discr = ((Tup.Dist5D_P-Tup.Dist5D)/(Tup.Dist5D_P+Tup.Dist5D));
			Kbin=AglDB.GetBin(RUsed);
				if(Massa_gen<1&&Massa_gen>0.5) ((TH2*)FitAgl_Dbins_Dist -> TemplateP) -> Fill(Distance_Discr,Kbin,Tup.mcweight);
				if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitAgl_Dbins_Dist -> TemplateD) -> Fill(Distance_Discr,Kbin,ReturnMCGenType(),Tup.mcweight);
				if(Massa_gen<4&&Massa_gen>2.5) ((TH2*)FitAgl_Dbins_Dist -> TemplateHe)-> Fill(Distance_Discr,Kbin,Tup.mcweight);
			Kbin=AglPB.GetBin(RUsed);	
				if(Massa_gen<1&&Massa_gen>0.5) ((TH2*)FitAgl_Pbins_Dist -> TemplateP) -> Fill(Distance_Discr,Kbin,Tup.mcweight);
				if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitAgl_Pbins_Dist -> TemplateD) -> Fill(Distance_Discr,Kbin,ReturnMCGenType(),Tup.mcweight);
				if(Massa_gen<4&&Massa_gen>2.5) ((TH2*)FitAgl_Pbins_Dist -> TemplateHe)-> Fill(Distance_Discr,Kbin,Tup.mcweight);
			}
}


void DeutonsDATA_Dist_Fill(int zona)
{
	float Distance_Discr = 0;
	int Kbin;
	if(!(Likcut&&Distcut)) return;
		Distance_Discr = ((Tup.Dist5D_P-Tup.Dist5D)/(Tup.Dist5D_P+Tup.Dist5D));
		Kbin=ToFDB.GetBin(RUsed);	
			if(Tup.R>SF*Tup.Rcutoff) FitTOF_Dbins_Dist -> DATA -> Fill(Distance_Discr,Kbin);
			((TH3*)FitTOFgeo_Dbins_Dist -> DATA) -> Fill(Distance_Discr,Kbin,zona);
		Kbin=ToFPB.GetBin(RUsed);	
			if(Tup.R>SF*Tup.Rcutoff) FitTOF_Pbins_Dist -> DATA -> Fill(Distance_Discr,Kbin);
		if(cmask.isFromNaF()) {
			Distance_Discr =  ((Tup.Dist5D_P-Tup.Dist5D)/(Tup.Dist5D_P+Tup.Dist5D));
			Kbin=NaFDB.GetBin(RUsed);
				if(Tup.R>SF*Tup.Rcutoff) FitNaF_Dbins_Dist -> DATA -> Fill(Distance_Discr,Kbin);
				((TH3*)FitNaFgeo_Dbins_Dist -> DATA) -> Fill(Distance_Discr,Kbin,zona);
			Kbin=NaFPB.GetBin(RUsed);	
				if(Tup.R>SF*Tup.Rcutoff) FitNaF_Pbins_Dist -> DATA -> Fill(Distance_Discr,Kbin);
		}
		if(cmask.isFromAgl()) {
			Distance_Discr =  ((Tup.Dist5D_P-Tup.Dist5D)/(Tup.Dist5D_P+Tup.Dist5D));
			Kbin=AglDB.GetBin(RUsed);	
				if(Tup.R>SF*Tup.Rcutoff) FitAgl_Dbins_Dist -> DATA -> Fill(Distance_Discr,Kbin);
				((TH3*)FitAglgeo_Dbins_Dist -> DATA) -> Fill(Distance_Discr,Kbin,zona);
			Kbin=AglPB.GetBin(RUsed);	
				if(Tup.R>SF*Tup.Rcutoff) FitAgl_Pbins_Dist -> DATA -> Fill(Distance_Discr,Kbin);
			}
	return;
}


void DeutonsMC_Dist_Write()
{

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


void DeutonsTemplFits_Dist(string filename,string frac)
{

	cout<<"******************** DEUTONS TEMPlATE FITS ************************"<<endl;
	cout<<"*** Reading  P1 file ****"<<endl;
	TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");



	TemplateFIT * FitTOF_Dbins_Dist	= new TemplateFIT(inputHistoFile,"FitTOF_Dbins_Dist ","FitTOF_Dbins_Dist ");
	TemplateFIT * FitNaF_Dbins_Dist	= new TemplateFIT(inputHistoFile,"FitNaF_Dbins_Dist ","FitNaF_Dbins_Dist ");
	TemplateFIT * FitAgl_Dbins_Dist	= new TemplateFIT(inputHistoFile,"FitAgl_Dbins_Dist ","FitAgl_Dbins_Dist ");

	TemplateFIT * FitTOFgeo_Dbins_Dist	= new TemplateFIT(inputHistoFile,"FitTOF_Dbins_Dist ","FitTOFgeo_Dbins_Dist ",11);
	TemplateFIT * FitNaFgeo_Dbins_Dist	= new TemplateFIT(inputHistoFile,"FitNaF_Dbins_Dist ","FitNaFgeo_Dbins_Dist ",11);
	TemplateFIT * FitAglgeo_Dbins_Dist	= new TemplateFIT(inputHistoFile,"FitAgl_Dbins_Dist ","FitAglgeo_Dbins_Dist ",11);

	TemplateFIT * FitTOF_Pbins_Dist	= new TemplateFIT(inputHistoFile,"FitTOF_Pbins_Dist ","FitTOF_Pbins_Dist ");
	TemplateFIT * FitNaF_Pbins_Dist	= new TemplateFIT(inputHistoFile,"FitNaF_Pbins_Dist ","FitNaF_Pbins_Dist ");
	TemplateFIT * FitAgl_Pbins_Dist	= new TemplateFIT(inputHistoFile,"FitAgl_Pbins_Dist ","FitAgl_Pbins_Dist ");

	cout<<"******************** DEUTONS DISTANCE TEMPlATE FITS ************************"<<endl;

	FitTOF_Dbins_Dist 	-> DisableFit();//SetFitConstraints(0.5,1,0.0000,0.01,0.0001,0.0025);
	FitNaF_Dbins_Dist 	-> DisableFit();//SetFitConstraints(0.8,1,0.0000,0.01,0.0005,0.0015);
	FitAgl_Dbins_Dist 	-> DisableFit();//SetFitConstraints(0.8,1,0.0000,0.01,0.0001,0.0005);

	FitTOFgeo_Dbins_Dist	-> DisableFit();
	FitNaFgeo_Dbins_Dist	-> DisableFit();
	FitAglgeo_Dbins_Dist	-> DisableFit();

	FitTOF_Pbins_Dist 	-> DisableFit();// SetFitConstraints(0.5,1,0.0000,0.01,0.0001,0.0025);
	FitNaF_Pbins_Dist 	-> DisableFit();// SetFitConstraints(0.8,1,0.0000,0.01,0.0005,0.0015);
	FitAgl_Pbins_Dist 	-> DisableFit();// SetFitConstraints(0.8,1,0.0000,0.01,0.0001,0.0005);

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
	for(int bin =0; bin <nbinsToF; bin++)
		cout<<FitTOF_Dbins_Dist->GetFitOutcome(bin)<<" ";
	cout<<endl;
	for(int bin =0; bin <nbinsToF; bin++)
		cout<<FitTOF_Pbins_Dist->GetFitOutcome(bin)<<" ";
	cout<<endl;

	cout<<"** NaF **"<<endl;
	for(int bin =0; bin <nbinsNaF; bin++)
		cout<<FitNaF_Dbins_Dist->GetFitOutcome(bin)<<" ";
	cout<<endl;
	for(int bin =0; bin <nbinsNaF; bin++)
		cout<<FitNaF_Pbins_Dist->GetFitOutcome(bin)<<" ";
	cout<<endl;
	cout<<"** Agl **"<<endl;
	for(int bin =0; bin <nbinsAgl; bin++)
		cout<<FitAgl_Dbins_Dist->GetFitOutcome(bin)<<" ";
	cout<<endl;
	for(int bin =0; bin <nbinsAgl; bin++)
		cout<<FitAgl_Pbins_Dist->GetFitOutcome(bin)<<" ";
	cout<<endl;

	cout<<"*** Updating P1 file ****"<<endl;

	inputHistoFile->ReOpen("UPDATE");

	FitTOF_Dbins_Dist -> DCounts -> SetName ("D_Flux_DistCounts_TOF");
	FitNaF_Dbins_Dist -> DCounts -> SetName ("D_Flux_DistCounts_NaF");
	FitAgl_Dbins_Dist -> DCounts -> SetName ("D_Flux_DistCounts_Agl");

	FitTOFgeo_Dbins_Dist -> DCounts -> SetName ("D_Flux_geo_DistCounts_TOF");
	FitNaFgeo_Dbins_Dist -> DCounts -> SetName ("D_Flux_geo_DistCounts_NaF");
	FitAglgeo_Dbins_Dist -> DCounts -> SetName ("D_Flux_geo_DistCounts_Agl");

	FitTOF_Pbins_Dist -> PCounts -> SetName ("P_Flux_DistCounts_TOF");
	FitNaF_Pbins_Dist -> PCounts -> SetName ("P_Flux_DistCounts_NaF");
	FitAgl_Pbins_Dist -> PCounts -> SetName ("P_Flux_DistCounts_Agl");

	finalHistos.Add(FitTOF_Dbins_Dist -> DCounts         );
	finalHistos.Add(FitNaF_Dbins_Dist -> DCounts         );
	finalHistos.Add(FitAgl_Dbins_Dist -> DCounts         );

	finalHistos.Add(FitTOFgeo_Dbins_Dist -> DCounts   );
	finalHistos.Add(FitNaFgeo_Dbins_Dist -> DCounts   );
	finalHistos.Add(FitAglgeo_Dbins_Dist -> DCounts   );

	finalHistos.Add(FitTOF_Pbins_Dist -> PCounts         );
	finalHistos.Add(FitNaF_Pbins_Dist -> PCounts         );
	finalHistos.Add(FitAgl_Pbins_Dist -> PCounts         );
	finalHistos.writeObjsInFolder("Results");

	cout<<"*** Plotting ...  ****"<<endl;

	string varname = " Dist. ";
	DeutonsTemplateFits_Plot(FitTOF_Dbins_Dist,
			FitNaF_Dbins_Dist ,
			FitAgl_Dbins_Dist,

			FitTOFgeo_Dbins_Dist,
			FitNaFgeo_Dbins_Dist,
			FitAglgeo_Dbins_Dist,

			FitTOF_Pbins_Dist,
			FitNaF_Pbins_Dist ,
			FitAgl_Pbins_Dist,
			varname
			);


	return;
}


