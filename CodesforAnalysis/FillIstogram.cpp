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
      

   filename=outputpath+"Final_plots/"+mese+".root";
   finalPlots.setName(filename);
   fileFinalPlots=new TFile(filename.c_str(),"UPDATE");


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
   Tempi = (TH1F *)usedfile->Get("Tempi");
   TH2F* esposizionegeo = (TH2F *)usedfile->Get("esposizionegeo");
   TH2F* esposizionepgeo = (TH2F*)usedfile->Get("esposizionepgeo");
   TH2F* esposizionepgeoNaF = (TH2F*)usedfile->Get("esposizionepgeoNaF");
   TH2F* esposizionepgeoAgl = (TH2F*)usedfile->Get("esposizionepgeoAgl");
   TH2F* esposizionedgeo = (TH2F*)usedfile->Get("esposizionedgeo");
   TH2F* esposizionedgeoNaF = (TH2F*)usedfile->Get("esposizionedgeoNaF");
   TH2F* esposizionedgeoAgl = (TH2F*)usedfile->Get("esposizionedgeoAgl");



   if(INDX==BUILDALL) {
      LoopOnDataTrig(ntupDataTrig);

   }
   if(INDX==BUILDALL||INDX==BUILDSEPD) {
      LoopOnDataSepD(ntupDataSepD);

      cout<<endl<<"************************ SAVING DATA ************************"<<endl;

      inputHistoFile->Flush();
      inputHistoFile->Close();
      delete inputHistoFile;
      inputHistoFile =new TFile(filename.c_str(), "RECREATE");
      inputHistoFile->cd();

      DATAQualeff_Write();
      DATARICHeff_Write();
      DATApreSeleff_Write();
      FluxFactorizationtest_Write();
      Correlazione_Preselezioni_Write();
      HecutMC_Write();
      SlidesforPlot_Write();
      DistanceCut_Write();
      DATAUnbiaseff_Write();
      DeutonsMC_Write();
      DeutonsMC_Dist_Write();
      DVSMCPreSeleff_Write();
      DVSMCPreSeleffD_Write();
      DVSMCQualeff2_Write();
      DVSMCRICHeff_Write();
      DVSMCQualeffD_Write();
      DVSMCTrackeff_Write();
      MCMC_Write();
      MCpreeff_Write();
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

      inputHistoFile->Write();
      inputHistoFile->Flush();
   }

   inputHistoFile->ReOpen("READ");

   cout<<"************************* ANALYSIS **********************************************************************"<<endl;
   if(INDX!=1) {
      if(frac=="tot") Hecut(inputHistoFile);
      SlidesforPlot(inputHistoFile);
      DistanceCut(inputHistoFile);
      Correlazione_Preselezioni(inputHistoFile);

      MCpreeff(inputHistoFile);
      MCUnbiaseff(inputHistoFile);
      MCQualeff(inputHistoFile);
      FluxFactorizationtest(inputHistoFile);
      MCTrackeff(inputHistoFile);
      MCFullseteff(inputHistoFile);
      MigrationMatrix(inputHistoFile);
      DVSMCTrackeff(inputHistoFile);
      DATAUnbiaseff(inputHistoFile);
      DATApreSeleff(inputHistoFile);
      DATAQualeff(inputHistoFile);
      DATARICHeff(inputHistoFile);
      if(frac=="tot") DeutonsTemplFits();
      if(frac=="tot") DeutonsTemplFits_Dist();

      CorrLAT();
      DVSMCPreSeleff();
      DVSMCPreSeleffD();
      DVSMCRICHeff();
      DVSMCQualeff2();
      DVSMCQualeffD();
      Acceptance();
      ProtonFlux();
      if(frac=="tot") DeutonFlux();
      if(frac=="tot") OtherExperimentsComparison();
   }

   fileFinalPlots->Close();
   inputHistoFile->Close();


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
   ntupMCTrig->SetBranchAddress("Unbias",&Tup.Unbias);

   ntupMCSepD->SetBranchAddress("Momentogen",&Tup.Momento_gen);
   ntupMCSepD->SetBranchAddress("R",&Tup.R);
   ntupMCSepD->SetBranchAddress("Beta",&Tup.Beta);
   ntupMCSepD->SetBranchAddress("BetaRICH_new",&Tup.BetaRICH);
   ntupMCSepD->SetBranchAddress("EdepL1",&Tup.EdepL1);
   ntupMCSepD->SetBranchAddress("Rmin",&Tup.Rmin);
   ntupMCSepD->SetBranchAddress("EdepTOF",&Tup.EdepTOFU);
   ntupMCSepD->SetBranchAddress("EdepTrack",&Tup.EdepTrack);
   ntupMCSepD->SetBranchAddress("EdepTOFD",&Tup.EdepTOFD);
   ntupMCSepD->SetBranchAddress("MC_type",&Tup.MC_type);
   ntupMCSepD->SetBranchAddress("LDiscriminant",&Tup.LDiscriminant);
   ntupMCSepD->SetBranchAddress("BDT_response",&Tup.BDT_response);
   ntupMCSepD->SetBranchAddress("Cutmask",&Tup.Cutmask);
   ntupMCSepD->SetBranchAddress("Dist5D",&Tup.Dist5D);
   ntupMCSepD->SetBranchAddress("Dist5D_P",&Tup.Dist5D_P);

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
   ntupDataTrig->SetBranchAddress("Unbias",&Tup.Unbias);


   ntupDataSepD->SetBranchAddress("R",&Tup.R);
   ntupDataSepD->SetBranchAddress("Beta",&Tup.Beta);
   ntupDataSepD->SetBranchAddress("BetaRICH_new",&Tup.BetaRICH);
   ntupDataSepD->SetBranchAddress("EdepL1",&Tup.EdepL1);
   ntupDataSepD->SetBranchAddress("Rcutoff",&Tup.Rcutoff);
   ntupDataSepD->SetBranchAddress("Rmin",&Tup.Rmin);
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
      Cuts_Pre();
      Massa_gen = ReturnMass_Gen();
      RUsed=Tup.R_pre;
      UpdateProgressBar(i, nentries);
      
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
      Massa_gen = ReturnMass_Gen();
      UpdateProgressBar(i, nentries);
      Cuts();
      RUsed=Tup.R;

      HecutMC_Fill();
      SlidesforPlot_Fill();
      FluxFactorizationtest_Qual_Fill();
      DistanceCut_Fill();
      MCQualeff_Fill();
      if(!(Tup.R<1.2*Tup.Rcutoff||Tup.Beta>protons->Eval(Tup.R)+0.1||Tup.Beta<protons->Eval(Tup.R)-0.1) && Herejcut) {
         DVSMCQualeff2_Fill();
         DVSMCQualeffD_Fill();
         DVSMCRICHeff_Fill();
      }
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
      if((cmask.isFromAgl()||cmask.isFromNaF())&&Tup.BetaRICH<0) continue;
      if(Tup.Beta_pre<=0) continue;
      Cuts_Pre();

      // Temporary Betarich check
      RUsed=Tup.R_pre;
      UpdateProgressBar(i, nentries);
      DATAUnbiaseff_Fill();
      float Zona=getGeoZone(Tup.Latitude);
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
   InitProtonFlux();
   int nentries=ntupDataSepD->GetEntries();
   for(int i=0; i<nentries; i++) {
      ntupDataSepD->GetEvent(i);
      cmask.setMask(Tup.Cutmask);
      if((cmask.isFromAgl()||cmask.isFromNaF())&&Tup.BetaRICH<0) continue;
      if(Tup.Beta<=0||Tup.R<=0) continue;
      
      Cuts();
      float Zona=getGeoZone(Tup.Latitude);
      RUsed=Tup.R;
      UpdateProgressBar(i, nentries);

      HecutD_Fill();
      SlidesforPlot_D_Fill();
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
