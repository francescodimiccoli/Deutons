using namespace std;


void SetRisultatiBranchAddresses (TNtuple* ntupMCSepD, TNtuple* ntupMCTrig, TNtuple* ntupDataSepD, TNtuple* ntupDataTrig);
void LoopOnMCTrig(TNtuple*  ntupMCTrig);
void LoopOnMCSepD(TNtuple* ntupMCSepD);
void LoopOnDataTrig(TNtuple* ntupDataTrig);
void LoopOnDataSepD(TNtuple* ntupDataSepD);
void UpdateProgressBar(int currentevent, int totalentries);
float getGeoZone(float latitude);

void FillIstogramAndDoAnalysis(mode INDX,string frac,string mese, string outputpath)
{

   cout<<"*********************** CALIB. READING *********************"<<endl;

   string inputpath="/storage/gpfs_ams/ams/users/fdimicco/Deutons";
   //string inputpath="/home/francesco/PhD/Deutons"; 
   string nomecal=inputpath + "/CodesforAnalysis/CALIBRAZIONI/"+mese+".root";
   TFile *calib = TFile::Open(nomecal.c_str());
   if(calib) cout<<"MC calibration for month "<<mese<<" ... ok"<<endl;
   else cout<<"ERROR: MC calibration not found"<<endl;
   Rig           = (TSpline3 *) calib->Get("Fit Results/Splines/Rig");
   beta          = (TSpline3 *) calib->Get("Fit Results/Splines/beta");
   etofu         = (TSpline3 *) calib->Get("Fit Results/Splines/etofu");
   etrack        = (TSpline3 *) calib->Get("Fit Results/Splines/etrack");
   EdepL1beta    = (TSpline3 *) calib->Get("Fit Results/Splines/EdepL1beta");
   EdepTOFbeta   = (TSpline3 *) calib->Get("Fit Results/Splines/EdepTOFbeta");
   EdepTrackbeta = (TSpline3 *) calib->Get("Fit Results/Splines/EdepTrackbeta");

   betaNaF = (TF1 *) calib->Get("Fit Results/Splines/SigmaInvBetaNaF_spl");
   betaAgl = (TF1 *) calib->Get("Fit Results/Splines/SigmaInvBetaAgl_spl");
   cout<<"******************************"<<endl;



   TFile *fileMC;
   TNtuple *ntupMCSepD;
   TNtuple *ntupMCTrig;
   TFile *fileData;
   TNtuple *ntupDataSepD;
   TNtuple *ntupDataTrig;

   string filename=outputpath+"Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
   inputHistoFile=new TFile(filename.c_str(),"READ");
   if (inputHistoFile->IsZombie()) {
      cout<<"## Histograms file not detected: rebuilding from trigger ##"<<endl;
      INDX=BUILDALL;
      cout<<"Running in Mode 0 ..."<<endl;
   }


   if(INDX!=READ) {
      string nomefile=inputpath + "/Risultati/"+mese+"/RisultatiMC_"+frac+".root";
      fileMC =TFile::Open(nomefile.c_str());
      nomefile=inputpath+"/Risultati/"+mese+"/RisultatiDATI_"+frac+".root";
      fileData =TFile::Open(nomefile.c_str(), "READ");
      ntupMCSepD=(TNtuple*)fileMC->Get("grandezzesepd");
      ntupMCTrig=(TNtuple*)fileMC->Get("trig");
      ntupDataSepD=(TNtuple*)fileData->Get("grandezzesepd");
      ntupDataTrig=(TNtuple*)fileData->Get("trig");
      SetRisultatiBranchAddresses(ntupMCSepD, ntupMCTrig, ntupDataSepD, ntupDataTrig);
   }

   cout<<"*********************** MC READING *********************"<<endl;

   if(INDX==BUILDALL) {
      LoopOnMCTrig(ntupMCTrig);
   }
   if(INDX==BUILDALL||INDX==BUILDSEPD) {
      LoopOnMCSepD(ntupMCSepD);
   }

   cout<<endl<<"*********************** DATA READING *********************"<<endl;
   TFile *usedfile=(INDX==READ?inputHistoFile:fileData);

   if(INDX==BUILDALL) {
      LoopOnDataTrig(ntupDataTrig);

   }
   if(INDX==BUILDALL||INDX==BUILDSEPD) {
      LoopOnDataSepD(ntupDataSepD);
   }

   

    cout<<endl<<"************************ SAVING DATA ************************"<<endl;
  


    if(INDX==BUILDALL||INDX==BUILDSEPD) {

      TFile * outputHistoFile=TFile::Open(filename.c_str(),"RECREATE");
	
  
      ExposureTime_Write();	  
      DATAQualeff_Write();
      DATARICHeff_Write();
      DATApreSeleff_Write();
      FluxFactorizationtest_Write();
      Correlazione_Preselezioni_Write();
      HecutMC_Write();
      SlidesforPlot_Write();
      DistanceCut_Write();
      AntiDCutOptimization_Write();
      DATAUnbiaseff_Write();
      DATAUnbiaseffQ_Write();
      DeutonsMC_Write();
      DeutonsMC_Dist_Write();
      DVSMCPreSeleff_Write();
      DVSMCPreSeleffD_Write();
      DVSMCQualeff2_Write();
      DVSMCRICHeff_Write();
      DVSMCQualeffD_Write();
      DVSMCTrackeff_Write();
      MCMC_Write();
      MCUnbiaseff_Write();
      MCpreeff_Write();
      MCControlsamplecuteff_Write();
      MCQualeff_Write();
      MCTrackeff_Write();
      MigrationMatrix_Write();
      ProtonFlux_Write();

   outputHistoFile -> Write();
   outputHistoFile -> Close();	
   }
   
   cout<<"************************* ANALYSIS **********************************************************************"<<endl;
   

   string finalfilename="./Final_plots/"+mese+".root";
   finalPlots.setName(finalfilename);
   finalHistos.setName(filename);  
   INDX = READ; 	

   if(INDX==READ) {
      if(frac=="tot") Hecut(filename);
      SlidesforPlot(filename);
      //DistanceCut(filename);
      Correlazione_Preselezioni(filename);

      MCpreeff(filename);
      MCUnbiaseff(filename);
      MCControlsamplecuteff(filename);
      MCQualeff(filename);
      FluxFactorizationtest(filename);
      MCTrackeff(filename);
      AntiDCutOptimization(filename);	
      AntiDEfficiencies(filename);
      MCFullseteff(filename);
      MigrationMatrix(filename);
      DVSMCTrackeff(filename);
      DATAUnbiaseff(filename);
      DATApreSeleff(filename);
      DATAQualeff(filename);
      DATARICHeff(filename);
      DeutonsTemplFits(filename);
      DeutonsTemplFits_Dist(filename);

      CorrLAT(filename);
      DVSMCPreSeleff(filename);
      DVSMCPreSeleffD(filename);
      DVSMCRICHeff(filename);
      DVSMCQualeff2(filename);
      DVSMCQualeffD(filename);
      Acceptance(filename);
      AntiDpredictions(filename);
      ProtonFlux(filename);
      DeutonFlux(filename);
      OtherExperimentsComparison(filename);
   }



   return;
}




void SetRisultatiBranchAddresses(TNtuple* ntupMCSepD, TNtuple* ntupMCTrig, TNtuple* ntupDataSepD, TNtuple* ntupDataTrig)
{

   ntupMCTrig->SetBranchAddress("Momento_gen",&Tup.Momento_gen);
   ntupMCTrig->SetBranchAddress("Ev_Num",&Tup.Ev_Num);
   ntupMCTrig->SetBranchAddress("Trig_Num",&Tup.Trig_Num);
   ntupMCTrig->SetBranchAddress("R_pre",&Tup.R_pre);
   ntupMCTrig->SetBranchAddress("Beta_pre",&Tup.Beta_pre);
   ntupMCTrig->SetBranchAddress("Cutmask",&Tup.Cutmask);
   ntupMCTrig->SetBranchAddress("MC_type",&Tup.MC_type);
   ntupMCTrig->SetBranchAddress("EdepL1",&Tup.EdepL1);
   ntupMCTrig->SetBranchAddress("EdepTOFU",&Tup.EdepTOFU);
   ntupMCTrig->SetBranchAddress("EdepTOFD",&Tup.EdepTOFD);
   ntupMCTrig->SetBranchAddress("EdepTrack",&Tup.EdepTrack);
   ntupMCTrig->SetBranchAddress("EdepECAL",&Tup.EdepECAL);
   ntupMCTrig->SetBranchAddress("BetaRICH",&Tup.BetaRICH);
   ntupMCTrig->SetBranchAddress("PhysBPatt",&Tup.PhysBPatt);
   ntupMCTrig->SetBranchAddress("mcweight",&Tup.mcweight);

   ntupMCSepD->SetBranchAddress("Momentogen",&Tup.Momento_gen);
   ntupMCSepD->SetBranchAddress("R",&Tup.R);
   ntupMCSepD->SetBranchAddress("Beta",&Tup.Beta);
   ntupMCSepD->SetBranchAddress("BetaRICH_new",&Tup.BetaRICH);
   ntupMCSepD->SetBranchAddress("EdepL1",&Tup.EdepL1);
   ntupMCSepD->SetBranchAddress("PhysBPatt",&Tup.PhysBPatt);
   ntupMCSepD->SetBranchAddress("EdepTOF",&Tup.EdepTOFU);
   ntupMCSepD->SetBranchAddress("EdepTrack",&Tup.EdepTrack);
   ntupMCSepD->SetBranchAddress("EdepTOFD",&Tup.EdepTOFD);
   ntupMCSepD->SetBranchAddress("MC_type",&Tup.MC_type);
   ntupMCSepD->SetBranchAddress("LDiscriminant",&Tup.LDiscriminant);
   ntupMCSepD->SetBranchAddress("BDT_response",&Tup.BDT_response);
   ntupMCSepD->SetBranchAddress("Cutmask",&Tup.Cutmask);
   ntupMCSepD->SetBranchAddress("Dist5D",&Tup.Dist5D);
   ntupMCSepD->SetBranchAddress("Dist5D_P",&Tup.Dist5D_P);
   ntupMCSepD->SetBranchAddress("mcweight",&Tup.mcweight);
   

   ntupDataTrig->SetBranchAddress("Rcutoff",&Tup.Rcutoff);
   ntupDataTrig->SetBranchAddress("R_pre",&Tup.R_pre);
   ntupDataTrig->SetBranchAddress("Beta_pre",&Tup.Beta_pre);
   ntupDataTrig->SetBranchAddress("Cutmask",&Tup.Cutmask);
   ntupDataTrig->SetBranchAddress("Latitude",&Tup.Latitude);
   ntupDataTrig->SetBranchAddress("EdepECAL",&Tup.EdepECAL);
   ntupDataTrig->SetBranchAddress("EdepL1",&Tup.EdepL1);
   ntupDataTrig->SetBranchAddress("EdepTOFU",&Tup.EdepTOFU);
   ntupDataTrig->SetBranchAddress("EdepTOFD",&Tup.EdepTOFD);
   ntupDataTrig->SetBranchAddress("EdepTrack",&Tup.EdepTrack);
   ntupDataTrig->SetBranchAddress("BetaRICH",&Tup.BetaRICH);
   ntupDataTrig->SetBranchAddress("PhysBPatt",&Tup.PhysBPatt);
   ntupDataTrig->SetBranchAddress("Seconds",&Tup.U_time);	
   ntupDataTrig->SetBranchAddress("Livetime",&Tup.Livetime);	


   ntupDataSepD->SetBranchAddress("R",&Tup.R);
   ntupDataSepD->SetBranchAddress("Beta",&Tup.Beta);
   ntupDataSepD->SetBranchAddress("BetaRICH_new",&Tup.BetaRICH);
   ntupDataSepD->SetBranchAddress("EdepL1",&Tup.EdepL1);
   ntupDataSepD->SetBranchAddress("Rcutoff",&Tup.Rcutoff);
   ntupDataSepD->SetBranchAddress("PhysBPatt",&Tup.PhysBPatt);
   ntupDataSepD->SetBranchAddress("EdepTOFU",&Tup.EdepTOFU);
   ntupDataSepD->SetBranchAddress("EdepTrack",&Tup.EdepTrack);
   ntupDataSepD->SetBranchAddress("EdepTOFD",&Tup.EdepTOFD);
   ntupDataSepD->SetBranchAddress("Latitude",&Tup.Latitude);
   ntupDataSepD->SetBranchAddress("LDiscriminant",&Tup.LDiscriminant);
   ntupDataSepD->SetBranchAddress("BDT_response",&Tup.BDT_response);
   ntupDataSepD->SetBranchAddress("Cutmask",&Tup.Cutmask);
   ntupDataSepD->SetBranchAddress("Dist5D",&Tup.Dist5D);
   ntupDataSepD->SetBranchAddress("Dist5D_P",&Tup.Dist5D_P);


}



void LoopOnMCTrig(TNtuple*  ntupMCTrig)
{
   int nentries=ntupMCTrig->GetEntries();
   EffUnbiasMCP = new Efficiency("EffUnbiasMCP");
   EffUnbiasMCD = new Efficiency("EffUnbiasMCD");
   for(int i=0; i<ntupMCTrig->GetEntries(); i++) {
      ntupMCTrig->GetEvent(i);
      cmask.setMask(Tup.Cutmask);
      trgpatt.SetTriggPatt(Tup.PhysBPatt); 
      Cuts_Pre();
      Massa_gen = ReturnMass_Gen();
      RUsed=Tup.R_pre;
      UpdateProgressBar(i, nentries);
//      Disable_MCreweighting();

      MCpreseff_Fill();
      MCUnbiaseff_Fill();
      MCTrackeff_Fill();
      MigrationMatrix_Fill();
      Correlazione_Preselezioni();
      FluxFactorizationtest_Pre_Fill();
      DVSMCTrackeff_Fill();
      DVSMCPreSeleff_Fill();
      DVSMCPreSeleffD_Fill();
   }
   cout << endl;
   return;
}




void LoopOnMCSepD(TNtuple* ntupMCSepD)
{
   int nentries=ntupMCSepD->GetEntries();
   for(int i=0; i<ntupMCSepD->GetEntries(); i++) {
      ntupMCSepD->GetEvent(i);
      if(Tup.Beta<=0 || Tup.R<=0) continue;
      cmask.setMask(Tup.Cutmask);
      trgpatt.SetTriggPatt(Tup.PhysBPatt);
      if(!cmask.isPreselected()) continue;
      Massa_gen = ReturnMass_Gen();
      UpdateProgressBar(i, nentries);
      Cuts();
      RUsed=Tup.R;
//      Disable_MCreweighting();

      HecutMC_Fill();
      SlidesforPlot_Fill();
      FluxFactorizationtest_Qual_Fill();
      DistanceCut_Fill();
      AntiDCutOptimization_Fill();
      MCControlsamplecuteff_Fill();
      MCQualeff_Fill();
      DVSMCQualeff2_Fill();
      DVSMCQualeffD_Fill();
      DVSMCRICHeff_Fill();
      DeutonsMC_Fill();
      DeutonsMC_Dist_Fill();
      MCMC_Fill();

   }
   cout << endl;
   return;
}


void LoopOnDataTrig(TNtuple* ntupDataTrig)
{
	int nentries=ntupDataTrig->GetEntries();
	for(int i=0; i<ntupDataTrig->GetEntries(); i++) {
		ntupDataTrig->GetEvent(i);
		cmask.setMask(Tup.Cutmask);
		trgpatt.SetTriggPatt(Tup.PhysBPatt);
		Cuts_Pre();
		RUsed=Tup.R_pre;
		if(i==0) {
			ActualTime = (int)Tup.U_time;	
			cout<<"Starting time: "<<ActualTime<<endl;	
		}

		UpdateProgressBar(i, nentries); 	

		float Zona=getGeoZone(Tup.Latitude);	     
		ExposureTime_Fill(Zona);

		if((cmask.isFromAgl()||cmask.isFromNaF())&&Tup.BetaRICH<0) continue;
                if(Tup.Beta_pre<=0) continue;


		DATApreSeleff_Fill(Zona);
		DVSMCTrackeff_D_Fill(); // < Check if needs the ones before
		DVSMCPreSeleff_D_Fill(Zona);
		DVSMCPreSeleffD_D_Fill(Zona);
	}
	cout << endl;
	return;
}



void LoopOnDataSepD(TNtuple* ntupDataSepD)
{
   //InitProtonFlux();
   int nentries=ntupDataSepD->GetEntries();
   for(int i=0; i<nentries; i++) {
      ntupDataSepD->GetEvent(i);
      cmask.setMask(Tup.Cutmask);
      trgpatt.SetTriggPatt(Tup.PhysBPatt);
      if((cmask.isFromAgl()||cmask.isFromNaF())&&Tup.BetaRICH<0) continue;
      if(Tup.Beta<=0||Tup.R<=0) continue;
      if(!cmask.isPreselected()) continue;

      Cuts();
      float Zona=getGeoZone(Tup.Latitude);
      RUsed=Tup.R;
      UpdateProgressBar(i, nentries);

      HecutD_Fill();
      SlidesforPlot_D_Fill();
      DATAUnbiaseffQ_Fill ();
      DATAQualeff_Fill(Zona);
      DATARICHeff_Fill(Zona);
      ProtonFlux_Fill(Zona);
      DVSMCQualeff2_D_Fill(Zona);
      DVSMCQualeffD_D_Fill(Zona);
      DVSMCRICHeff_D_Fill(Zona);
      DeutonsDATA_Fill(Zona);
      DeutonsDATA_Dist_Fill(Zona);
      MCMCDATA_Fill();
   }
   cout << endl;
   return;
}

void UpdateProgressBar(int currentevent, int totalentries)
{
   int newratio =(int)100*(currentevent/    (float)totalentries);
   int oldratio =(int)100*((currentevent-1)/(float)totalentries);
   if(newratio>oldratio)
      cout<<'\r' << "Progress : "<< newratio+1 << " %"<< flush; //+1 pour finir a 100%
}

float getGeoZone(float latitude)
{
   float zone=-1;
   double geomag[12]= {0,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.3};
   for(int z=0; z<12; z++) {
      if(latitude>geomag[z] && latitude<geomag[z+1])
         zone=z;
   }
   return zone;
}






