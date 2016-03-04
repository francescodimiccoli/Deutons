using namespace std;

void FillIstogram(int INDX,string frac,string mese)

{


	//string percorso="/home/francesco/PhD/LocalCNAF/";
	string percorso="/storage/gpfs_ams/ams/users/fdimicco/Deutons";
	TH1F * Esposizione[11];
	string numero[11]={"0","1","2","3","4","5","6","7","8","9","10"};
	string tagli[10]={"Trigger","3of4 TOF","TRD Segments","Rigidity exists","Chi^2 R","Matching TOF","Matching TRD","In TRD Accept.","1 Particle","1 Tr. Track"};
	string nome;
	float DXWind[43];
	for(int i=0;i<11;i++) {
		nome="Esposizione"+numero[i];
		Esposizione[i]= new TH1F(nome.c_str(),nome.c_str(),43,0,43);
	}


	float fraz=1;
	float Q=0;
	float temp=0;
	float Zona=0;
	int UnbiasPre=9;
	float tempi[11]={0};
	float Esposizionegeo[43][11]={{0}};
	double geomag[12]={0,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.3};
	float DistTOF,DistTrack,DistTRD=0;
	
	cout<<"*********************** CALIB. READING *********************"<<endl;

	string nomecal=("/storage/gpfs_ams/ams/users/fdimicco/Deutons/CodesforAnalysis/CALIBRAZIONI/"+mese+".root");
        TFile *calib = TFile::Open(nomecal.c_str());
        cout<<"calibrazione: "<<calib<<endl;	
	Rig = (TSpline3 *) calib->Get("Fit Results/Splines/Rig");
        beta = (TSpline3 *) calib->Get("Fit Results/Splines/beta");
        betaNaF = (TF1 *) calib->Get("Fit Results/Splines/SigmaInvBetaNaF_spl");
        betaAgl = (TF1 *) calib->Get("Fit Results/Splines/SigmaInvBetaAgl_spl");
        eL1 = (TSpline3 *) calib->Get("Fit Results/Splines/eL1");
        etofu =  (TSpline3 *) calib->Get("Fit Results/Splines/etofu");
        etrack =  (TSpline3 *) calib->Get("Fit Results/Splines/etrack");
        etofd =  (TSpline3 *) calib->Get("Fit Results/Splines/etofd");
        EdepL1beta =  (TSpline3 *) calib->Get("Fit Results/Splines/EdepL1beta");
        EdepTOFbeta =  (TSpline3 *) calib->Get("Fit Results/Splines/EdepTOFbeta");
        EdepTrackbeta =  (TSpline3 *) calib->Get("Fit Results/Splines/EdepTrackbeta");
        EdepTOFDbeta =  (TSpline3 *) calib->Get("Fit Results/Splines/EdepTOFDbeta");
        Corr_L1 =  (TSpline3 *) calib->Get("Fit Results/Splines/Corr_L1");
        Corr_TOFU =  (TSpline3 *) calib->Get("Fit Results/Splines/Corr_TOFU");
        Corr_Track =  (TSpline3 *) calib->Get("Fit Results/Splines/Corr_Track");
        Corr_TOFD =  (TSpline3 *) calib->Get("Fit Results/Splines/Corr_TOFD");
        cout<<Rig<<" "<<beta<<" "<<" "<<betaNaF<<" "<<betaAgl<<" "<<eL1<<" "<<etofu<<" "<<etrack<<" "<<etofd<<" "<<EdepL1beta<<" "<<EdepTOFbeta<<" "<<EdepTrackbeta<<" "<<EdepTOFDbeta<<" "<<Corr_L1<<" "<<Corr_TOFU<<" "<<Corr_Track<<" "<<Corr_TOFD<<endl;
	cout<<"******************************"<<endl;

	string nomefile=percorso + "/Risultati/risultati/"+mese+"_"+frac+"_P1.root";
 	TFile *file =TFile::Open(nomefile.c_str());	
	if(!file){
	nomefile=percorso + "/Risultati/"+mese+"/"+mese+"_"+frac+"_P1.root";
	file =TFile::Open(nomefile.c_str());
	}
	TFile *file1;
	TNtuple *ntupla1;
	TNtuple *ntupla0;
	TFile *file2;
        TNtuple *ntupla3;
        TNtuple *ntupla2;
	if(!file) {cout<<"## Histograms file not detected: rebuilding from trigger ##"<<endl; INDX=0;}
	if(INDX!=2){	
	nomefile=percorso+"/Risultati/risultati/RisultatiMC_"+frac+".root";	
	file1 =TFile::Open(nomefile.c_str());
	if(!file1)  {nomefile=percorso + "/Risultati/"+mese+"/RisultatiMC_"+frac+".root";
		    file1 =TFile::Open(nomefile.c_str());}
	ntupla1=(TNtuple*)file1->Get("grandezzesepd");
	ntupla0=(TNtuple*)file1->Get("trig");
	nomefile=percorso+"/Risultati/risultati/RisultatiDATI_"+frac+".root";	
	file2 =TFile::Open(nomefile.c_str());
	if(!file2)  {nomefile=percorso+"/Risultati/"+mese+"/RisultatiDATI_"+frac+".root";
                    file2 =TFile::Open(nomefile.c_str());}
        ntupla3=(TNtuple*)file2->Get("grandezzesepd");
        ntupla2=(TNtuple*)file2->Get("trig");
	
	ntupla0->SetBranchAddress("Momento_gen",&Momento_gen);
	ntupla0->SetBranchAddress("Ev_Num",&Ev_Num);
	ntupla0->SetBranchAddress("Trig_Num",&Trig_Num);	
	ntupla0->SetBranchAddress("R_pre",&R_pre);
	ntupla0->SetBranchAddress("Beta_pre",&Beta_pre);
	ntupla0->SetBranchAddress("Cutmask",&CUTMASK);	
	ntupla0->SetBranchAddress("Massagen",&Massa_gen);
	ntupla0->SetBranchAddress("EdepL1",&EdepL1);	
	ntupla0->SetBranchAddress("EdepTOFU",&EdepTOFU);	
	ntupla0->SetBranchAddress("EdepTOFD",&EdepTOFD);
	ntupla0->SetBranchAddress("EdepTrack",&EdepTrack);
	ntupla0->SetBranchAddress("EdepECAL",&EdepECAL);
	ntupla0->SetBranchAddress("BetaRICH",&BetaRICH);
	ntupla0->SetBranchAddress("Unbias",&Unbias);

	ntupla1->SetBranchAddress("Momentogen",&Momento_gen);
	ntupla1->SetBranchAddress("R",&R);
        ntupla1->SetBranchAddress("Beta",&Beta);
	ntupla1->SetBranchAddress("BetaRICH_new",&BetaRICH);
        ntupla1->SetBranchAddress("EdepL1",&EdepL1);
        ntupla1->SetBranchAddress("Rmin",&Rmin);
        ntupla1->SetBranchAddress("EdepTOF",&EdepTOFU);
        ntupla1->SetBranchAddress("EdepTrack",&EdepTrack);
        ntupla1->SetBranchAddress("EdepTOFD",&EdepTOFD);
	ntupla1->SetBranchAddress("Massagen",&Massa_gen);
        ntupla1->SetBranchAddress("LDiscriminant",&LDiscriminant);
        ntupla1->SetBranchAddress("BDT_response",&BDT_response);
        ntupla1->SetBranchAddress("Cutmask",&CUTMASK);
	ntupla1->SetBranchAddress("Dist5D",&Dist5D);
	ntupla1->SetBranchAddress("Dist5D_P",&Dist5D_P);	
	
	ntupla2->SetBranchAddress("Rcutoff",&Rcutoff);
        ntupla2->SetBranchAddress("Ev_Num",&Ev_Num);
        ntupla2->SetBranchAddress("Trig_Num",&Trig_Num);
        ntupla2->SetBranchAddress("R_pre",&R_pre);
        ntupla2->SetBranchAddress("Beta_pre",&Beta_pre);
	ntupla2->SetBranchAddress("Cutmask",&CUTMASK);
        ntupla2->SetBranchAddress("Latitude",&Latitude);
	ntupla2->SetBranchAddress("EdepECAL",&EdepECAL);
	ntupla2->SetBranchAddress("EdepL1",&EdepL1);
        ntupla2->SetBranchAddress("EdepTOFU",&EdepTOFU);
        ntupla2->SetBranchAddress("EdepTOFD",&EdepTOFD);
        ntupla2->SetBranchAddress("EdepTrack",&EdepTrack);
        ntupla2->SetBranchAddress("BetaRICH",&BetaRICH);
        ntupla2->SetBranchAddress("Unbias",&Unbias);


        ntupla3->SetBranchAddress("R",&R);
        ntupla3->SetBranchAddress("Beta",&Beta);
        ntupla3->SetBranchAddress("BetaRICH_new",&BetaRICH);
	ntupla3->SetBranchAddress("EdepL1",&EdepL1);
        ntupla3->SetBranchAddress("Rcutoff",&Rcutoff);
	ntupla3->SetBranchAddress("Rmin",&Rmin);
        ntupla3->SetBranchAddress("EdepTOFU",&EdepTOFU);
        ntupla3->SetBranchAddress("EdepTrack",&EdepTrack);
        ntupla3->SetBranchAddress("EdepTOFD",&EdepTOFD);
        ntupla3->SetBranchAddress("Latitude",&Latitude);
	ntupla3->SetBranchAddress("Massagen",&Massa_gen);
        ntupla3->SetBranchAddress("LDiscriminant",&LDiscriminant);
        ntupla3->SetBranchAddress("BDT_response",&BDT_response);
        ntupla3->SetBranchAddress("Cutmask",&CUTMASK);
        ntupla3->SetBranchAddress("Dist5D",&Dist5D);
        ntupla3->SetBranchAddress("Dist5D_P",&Dist5D_P);
	}
	cout<<"*********************** MC READING *********************"<<endl;
	if(INDX==0)
	for(int i=0; i<ntupla0->GetEntries()/fraz;i++) {
		int k = ntupla0->GetEvent(i);
		Cutmask=CUTMASK;
		Var3=Momento_gen;
                Var=R_pre;
                Var2=R_pre;
                /*Var=Beta_pre;
                  Var2=BetaRICH;  
                  Var3=Beta_gen;*/
		if(100*(i/(float)(ntupla0->GetEntries()/fraz))>avanzamento) {cout<<avanzamento<<endl;avanzamento=(int)(100*(i/(float)(ntupla0->GetEntries()/fraz)))+1;}
		Beta_gen=(pow(pow(Momento_gen/Massa_gen,2)/(1+pow(Momento_gen/Massa_gen,2)),0.5));
		MCpreseff_Fill(ntupla0,i);
		MCUnbiaseff_Fill(ntupla0,i);
		MCTrackeff_Fill(ntupla0,i);
		MCpreSeleff_Fill(ntupla0,i);
		MCpreCheck_Fill(ntupla0,i);
		MigrationMatrix_Fill(ntupla0,i);
		Correlazione_Preselezioni(ntupla0,i);
		DVSMCpreSeleff_Fill(ntupla0,i);
		DVSMCTrackeff_Fill(ntupla0,i);
	}
	if(INDX==0||INDX==1){
	avanzamento=0;
	for(int i=0; i<ntupla1->GetEntries();i++) {
		int k = ntupla1->GetEvent(i);
		Cutmask=CUTMASK;
        	if(100*(i/(float)(ntupla1->GetEntries()))>avanzamento) {cout<<avanzamento<<endl;avanzamento=(int)(100*(i/(float)(ntupla1->GetEntries())))+1;}
		Cuts();
		Var3=Momento_gen;
 	        Var=R;
  		Var2=R;
        	/*Var=Beta;
         	Var2=BetaRICH;
                Var3=Beta_gen;*/
		HecutMC_Fill(ntupla1,i);
		MCQualeff_Fill(ntupla1,i);
		DVSMCQualeff2_Fill(ntupla1,i);	
		DeutonsMC_Fill(ntupla1,i);
		DeutonsMC_Dist_Fill(ntupla1,i);
		MCMC_Fill(ntupla1,i);
	}
	}
	if(INDX==1||INDX==2){
		 MCpreeff_Copy(file);
		 MCUnbiaseff_Copy(file); 
		 MCTrackeff_Copy(file);
		 MCpreSeleff_Copy(file);
		 MCpreCheck_Copy(file);
		 MigrationMatrix_Copy(file);
		 Correlazione_Preselezioni(file);
		 DVSMCpreSeleff_Copy(file);
		 DVSMCTrackeff_Copy(file);	
        }
        if(INDX==2){
                 HecutMC_Copy(file);
		 MCQualeff_Copy(file);
		 DeutonsMC_Copy(file);
		 DeutonsMC_Dist_Copy(file);
        }
	
	cout<<"*********************** DATA READING *********************"<<endl;
        if(INDX!=2){
	Tempi = (TH1F *)file2->Get("Tempi");
	esposizionegeo = (TH2F *)file2->Get("esposizionegeo");
	esposizionepgeo = (TH2F*)file2->Get("esposizionepgeo");
        esposizionepgeoNaF = (TH2F*)file2->Get("esposizionepgeoNaF");
        esposizionepgeoAgl = (TH2F*)file2->Get("esposizionepgeoAgl");
        esposizionedgeo = (TH2F*)file2->Get("esposizionedgeo");
        esposizionedgeoNaF = (TH2F*)file2->Get("esposizionedgeoNaF");
        esposizionedgeoAgl = (TH2F*)file2->Get("esposizionedgeoAgl");
	}
	else{
	Tempi = (TH1F *)file->Get("Tempi");
        esposizionegeo = (TH2F *)file->Get("esposizionegeo");
	esposizionepgeo = (TH2F*)file->Get("esposizionepgeo");
        esposizionepgeoNaF = (TH2F*)file->Get("esposizionepgeoNaF");
        esposizionepgeoAgl = (TH2F*)file->Get("esposizionepgeoAgl");
        esposizionedgeo = (TH2F*)file->Get("esposizionedgeo");
        esposizionedgeoNaF = (TH2F*)file->Get("esposizionedgeoNaF");
        esposizionedgeoAgl = (TH2F*)file->Get("esposizionedgeoAgl");
	}
	avanzamento=0;
	if(INDX==0)
	for(int i=0; i<ntupla2->GetEntries()/fraz;i++) {
                int k = ntupla2->GetEvent(i);
		Cutmask=CUTMASK;
		Massa=pow(fabs(pow(fabs(R_pre)*pow((1-pow(Beta_pre,2)),0.5)/Beta_pre,2)),0.5);
	        for(int z=0;z<12;z++){
                        double geo= geomag[z]  ;
                        double geo2=geomag[z+1];
                        if(Latitude>geo && Latitude<geo2) Zona=z;
        	}
		Var3=Momento_gen;
                Var=R_pre;
                Var2=R_pre;
                /*Var=Beta_pre;
                  Var2=BetaRICH;  
                  Var3=Beta_gen;*/
		if(100*(i/(float)(ntupla2->GetEntries()/fraz))>avanzamento) {cout<<avanzamento<<endl;avanzamento=(int)(100*(i/(float)(ntupla2->GetEntries()/fraz)))+1;}
                
		Beta_gen=(pow(pow(Momento_gen/Massa_gen,2)/(1+pow(Momento_gen/Massa_gen,2)),0.5));    
		DATAUnbiaseff_Fill(ntupla2,i);
		DATApreSeleff_Fill(ntupla2,i,Zona);
		DVSMCpreSeleff_D_Fill(ntupla2,i,Zona);
		DVSMCTrackeff_D_Fill(ntupla0,i);		
	}
        if(INDX==0||INDX==1){
        avanzamento=0;
        for(int i=0; i<ntupla3->GetEntries();i++) {
                int k = ntupla3->GetEvent(i);
                Cutmask=CUTMASK;
                for(int z=0;z<12;z++){
                        double geo= geomag[z]  ;
                        double geo2=geomag[z+1];
                        if(Latitude>geo && Latitude<geo2) Zona=z;
                }
                Cuts();
		Var3=Momento_gen;
                Var=R;
                Var2=R;
                /*Var=Beta;
                  Var2=BetaRICH;  
                  Var3=Beta_gen;*/
		if(100*(i/(float)(ntupla3->GetEntries()))>avanzamento) {cout<<avanzamento<<endl;avanzamento=(int)(100*(i/(float)(ntupla3->GetEntries())))+1;}
		HecutD_Fill(ntupla3,i);	
		DATAQualeff_Fill(ntupla3,i,Zona);
		DVSMCQualeff2_D_Fill(ntupla3,i,Zona);
		ProtonFlux_Fill(ntupla3,i,Zona);	
		DeutonsDATA_Fill(ntupla3,i,Zona);
		DeutonsDATA_Dist_Fill(ntupla3,i,Zona);
		MCMCDATA_Fill(ntupla3,i);
		}
        }
        if(INDX==1||INDX==2){
                 DATAUnbiaseff_Copy(file);
		 DATApreSeleff_Copy(file);
        }
        if(INDX==2){
		DATAQualeff_Copy(file);
		DVSMCQualeff2_Copy(file);
		ProtonFlux_Copy(file);
		MCMC_Copy(file);	
        }
	
	cout<<"************************ SAVING DATA ************************"<<endl;
	nomefile=percorso+"/Risultati/risultati/RisultatiMC_"+frac+".root";
        file1 =TFile::Open(nomefile.c_str());
	if(file1)
	nomefile=percorso+"/Risultati/risultati/"+mese+"_"+frac+"_P1.root";
	else
	nomefile=percorso+"/Risultati/"+mese+"/"+mese+"_"+frac+"_P1.root";
	TFile *f_out=new TFile(nomefile.c_str(), "RECREATE");
	
	DATAQualeff_Write();
	DATApreSeleff_Write();
	Correlazione_Preselezioni_Write();
	HecutMC_Write();
	DATAUnbiaseff_Write();
	DeutonsMC_Write();
	DeutonsMC_Dist_Write();
	DVSMCpreSeleff_Write();
	DVSMCQualeff2_Write();
	DVSMCTrackeff_Write();
	MCMC_Write();
	MCpreeff_Write();
	MCpreCheck_Write();
	MCpreSeleff_Write();
	MCQualeff_Write();
	MCTrackeff_Write();	
	MCUnbiaseff_Write();
	MigrationMatrix_Write();
	ProtonFlux_Write();

	Tempi->Write();
	esposizionegeo->Write();
	esposizionepgeo->Write();
	esposizionepgeoNaF->Write();
	esposizionepgeoAgl->Write();
	esposizionedgeo->Write();
	esposizionedgeoNaF->Write();
	esposizionedgeoAgl->Write();

	/*MigrMatrix->Write();
	EffpreselMCP1->Write();
	EffpreselMCP2->Write();
	EffpreselMCP1NaF->Write();
        EffpreselMCP2NaF->Write();
	EffpreselMCP1Agl->Write();
        EffpreselMCP2Agl->Write();
	EffpreselMCP1_R->Write();
	EffpreselMCP2_R->Write();
	EffpreselMCD1->Write();
	EffpreselMCD2->Write();
	EffpreselMCD1NaF->Write();
        EffpreselMCD2NaF->Write();
	EffpreselMCD1Agl->Write();
        EffpreselMCD2Agl->Write();
	EffpreselMCD1_R->Write();
	EffpreselMCD2_R->Write();
	EffUnbiasMCP1->Write();
        EffUnbiasMCP2->Write();
        EffUnbiasMCP1_R->Write();
        EffUnbiasMCP2_R->Write();
        EffUnbiasMCD1->Write();
        EffUnbiasMCD2->Write();
        EffUnbiasMCD1_R->Write();
        EffUnbiasMCD2_R->Write();
	EffpreSelMCP1->Write();
        EffpreSelMCP2->Write();
        EffpreSelMCP1_R->Write();
        EffpreSelMCP2_R->Write();
        EffpreSelMCD1->Write();
        EffpreSelMCD2->Write();
        EffpreSelMCD1_R->Write();
        EffpreSelMCD2_R->Write();
	EffpreCheckP1_R->Write();
        EffpreCheckP2_R->Write();
	EffTriggMCP1->Write();
        EffTriggMCP2->Write();
        EffTriggMCP1_R->Write();
        EffTriggMCP2_R->Write();
        EffTriggMCD1->Write();
        EffTriggMCD2->Write();
        EffTriggMCD1_R->Write();
        EffTriggMCD2_R->Write();
	EffQTOFMCP1->Write();
        EffQTOFMCP2->Write();
        EffQTOFMCP1_R->Write();
        EffQTOFMCP2_R->Write();
        EffQTOFMCD1->Write();
        EffQTOFMCD2->Write();
        EffQTOFMCD1_R->Write();
        EffQTOFMCD2_R->Write();
	EffTrackMCP1->Write();
        EffTrackMCP2->Write();
        EffTrackMCP1_R->Write();
        EffTrackMCP2_R->Write();
        EffTrackMCD1->Write();
        EffTrackMCD2->Write();
        EffTrackMCD1_R->Write();
        EffTrackMCD2_R->Write();
	EffTOFMCP1->Write();
        EffTOFMCP2->Write();
        EffTOFMCP1_R->Write();
        EffTOFMCP2_R->Write();
        EffTOFMCD1->Write();
        EffTOFMCD2->Write();
        EffTOFMCD1_R->Write();
        EffTOFMCD2_R->Write();
	EffMassMCP1_R->Write();
        EffMassMCP2_R->Write();
	EffQualCheckMCP1->Write();
	EffQualCheckMCP2->Write();
	EffTOFUCheckMCP1->Write();
        EffTOFUCheckMCP2->Write();
	EffTrackCheckMCP1->Write();
        EffTrackCheckMCP2->Write();
	EffTOFDCheckMCP1->Write();
        EffTOFDCheckMCP2->Write();
	EffLikCheckMCP1->Write();
        EffLikCheckMCP2->Write();
	EffDistCheckMCP1->Write();
        EffDistCheckMCP2->Write();
        EffLik2CheckMCP1->Write(); 
        EffLik2CheckMCP2->Write();
	EffQualMCP->Write(); 
	EffQualMCP_Beta->Write();
	EffQualMCP_BetaNaF->Write();
	EffQualMCP_BetaAgl->Write(); 
	EffQualMCD->Write(); 
	EffQualMCD_Beta->Write();
	EffQualMCD_BetaNaF->Write();
	EffQualMCD_BetaAgl->Write();
	EffDistMCP->Write();
	EffDistMCP_Beta->Write();
	EffDistMCP_BetaNaF->Write();
	EffDistMCP_BetaAgl->Write();
	EffLikMCP->Write();
	EffLikMCP_Beta->Write();
	EffLikMCP_BetaNaF->Write();
	EffLikMCP_BetaAgl->Write();
	EffDistMCD->Write();
	EffDistMCD_Beta->Write();
	EffDistMCD_BetaNaF->Write();
	EffDistMCD_BetaAgl->Write();
        EffLikMCD->Write();
        EffLikMCD_Beta->Write();
	EffLikMCD_BetaNaF->Write();
	EffLikMCD_BetaAgl->Write();
	EffUnbiasDATA1->Write();
        EffUnbiasDATA2->Write();
        EffUnbiasDATA1_R->Write();
        EffUnbiasDATA2_R->Write();
        EffpreSelDATA1_R ->Write();
        EffpreSelDATA2_R ->Write();
	CorrelazionePreselezioni->Write();
	EffpreSelMCvsD1->Write();
	EffpreSelMCvsD2->Write();
	EffpreSelMCvsD1_R->Write();
	EffpreSelMCvsD2_R->Write();
	EffpreSelMCvsD1_D->Write();
        EffpreSelMCvsD2_D->Write();
        EffpreSelMCvsD1_R_D->Write();
        EffpreSelMCvsD2_R_D->Write();
	EffDistMCvsDP1->Write();
        EffDistMCvsDP2->Write();
        EffLik2MCvsDP1->Write();
        EffLik2MCvsDP2->Write();
	ECALvsR_D->Write();
	ECALvsR_MC->Write();
	PCountsgeo->Write();
	PCounts->Write();
	PCountsgeo_prim->Write();
	PCounts_pre->Write();
	PCounts_sel->Write();
	DTemplatesTOF_Dist->Write();
        PTemplatesTOF_Dist->Write();
        HeTemplatesTOF_Dist->Write();
	DTemplatesNaF->Write();
        PTemplatesNaF->Write();
        HeTemplatesNaF->Write();
	DTemplatesAgl->Write();
        PTemplatesAgl->Write();
        HeTemplatesAgl->Write();
	DTemplatesTOF_Dist2->Write();
        PTemplatesTOF_Dist2->Write();
        HeTemplatesTOF_Dist2->Write();
        DTemplatesNaF2->Write();
        PTemplatesNaF2->Write();
        HeTemplatesNaF2->Write();
        DTemplatesAgl2->Write();
        PTemplatesAgl2->Write();
        HeTemplatesAgl2->Write();
	DhistosgeoTOF_Dist->Write();
	DhistosgeoNaF->Write();
	DhistosgeoAgl->Write();
	DhistosTOF_Dist->Write();
        DhistosNaF->Write();
        DhistosAgl->Write();
	DhistosTOF_Dist2->Write();
        DhistosNaF2->Write();
        DhistosAgl2->Write();
	MCMCPTemplatesTOF->Write();
        MCMCDTemplatesTOF->Write();
        MCMCHeTemplatesTOF->Write();
        MCMCPTemplatesNaF->Write();
        MCMCDTemplatesNaF->Write();
        MCMCHeTemplatesNaF->Write();
        MCMCPTemplatesAgl->Write();
        MCMCDTemplatesAgl->Write();
        MCMCHeTemplatesAgl->Write();
        MCMCDataTOF->Write();
        MCMCDataNaF->Write();
        MCMCDataAgl->Write();*/
	f_out->Write();
	f_out->Close();


	return;
}
