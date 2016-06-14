using namespace std;


void SetRisultatiBranchAddresses (TNtuple* ntupMCSepD, TNtuple* ntupMCTrig, TNtuple* ntupDataSepD, TNtuple* ntupDataTrig);
void LoopOnMCTrig(TNtuple*  ntupMCTrig);
void LoopOnMCSepD(TNtuple* ntupMCSepD);
void LoopOnDataTrig(TNtuple* ntupDataTrig);
void LoopOnDataSepD(TNtuple* ntupDataSepD);
void UpdateProgressBar(TNtuple* ntuple, int currentevent);


void FillIstogram(int INDX,string frac,string mese)
{

   cout<<"*********************** CALIB. READING *********************"<<endl;

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

   string nomefile=outputpath + "/Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
   TFile * file =TFile::Open(nomefile.c_str());

   TFile *fileMC;
   TNtuple *ntupMCSepD;
   TNtuple *ntupMCTrig;
   TFile *fileData;
   TNtuple *ntupDataSepD;
   TNtuple *ntupDataTrig;
   bool histonotexist = false;
   if(!file) {
      cout<<"## Histograms file not detected: rebuilding from trigger ##"<<endl;
      INDX=0;
      histonotexist = true;
      cout<<"Running in Mode 0 ..."<<endl;
   }

   if(INDX!=2) {
      nomefile=inputpath + "/Risultati/"+mese+"/RisultatiMC_"+frac+".root";
      fileMC =TFile::Open(nomefile.c_str());
      nomefile=inputpath+"/Risultati/"+mese+"/RisultatiDATI_"+frac+".root";
      fileData =TFile::Open(nomefile.c_str());
      ntupMCSepD=(TNtuple*)fileMC->Get("grandezzesepd");
      ntupMCTrig=(TNtuple*)fileMC->Get("trig");
      ntupDataSepD=(TNtuple*)fileData->Get("grandezzesepd");
      ntupDataTrig=(TNtuple*)fileData->Get("trig");
      SetRisultatiBranchAddresses(ntupMCSepD, ntupMCTrig, ntupDataSepD, ntupDataTrig);
   }

   cout<<"*********************** MC READING *********************"<<endl;

   if(INDX==0) {
      LoopOnMCTrig(ntupMCTrig);
   }
   if(INDX==0||INDX==1) {
      LoopOnMCSepD(ntupMCSepD);
   }

   cout<<endl<<"*********************** DATA READING *********************"<<endl;
   TFile *usedfile=(INDX==2?file:fileData);
   Tempi = (TH1F *)usedfile->Get("Tempi");
   TH2F* esposizionegeo = (TH2F *)usedfile->Get("esposizionegeo");
   TH2F* esposizionepgeo = (TH2F*)usedfile->Get("esposizionepgeo");
   TH2F* esposizionepgeoNaF = (TH2F*)usedfile->Get("esposizionepgeoNaF");
   TH2F* esposizionepgeoAgl = (TH2F*)usedfile->Get("esposizionepgeoAgl");
   TH2F* esposizionedgeo = (TH2F*)usedfile->Get("esposizionedgeo");
   TH2F* esposizionedgeoNaF = (TH2F*)usedfile->Get("esposizionedgeoNaF");
   TH2F* esposizionedgeoAgl = (TH2F*)usedfile->Get("esposizionedgeoAgl");



   if(INDX==0) {
      LoopOnDataTrig(ntupDataTrig);

   }
   if(INDX==0||INDX==1) {
      LoopOnDataSepD(ntupDataSepD);


      cout<<endl<<"************************ SAVING DATA ************************"<<endl;

      nomefile= outputpath + "Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
      TFile *f_out=new TFile(nomefile.c_str(), "RECREATE");

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
      //DVSMCTrackeff_Write();
      MCMC_Write();
      MCpreeff_Write();
      //MCpreCheck_Write();
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

      f_out->Write();
      f_out->Close();
   }

   if(histonotexist) INDX=2;
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
   for(int i=0; i<ntupMCTrig->GetEntries(); i++) {
      ntupMCTrig->GetEvent(i);
      Cuts_Pre();
      Massa_gen = ReturnMass_Gen();
      RUsed=Tup.R_pre;
      UpdateProgressBar(ntupMCTrig, i);
      
      MCpreseff_Fill(ntupMCTrig,i);
      MCUnbiaseff_Fill(ntupMCTrig,i);
      MCTrackeff_Fill(ntupMCTrig,i);
      MigrationMatrix_Fill(ntupMCTrig,i);
      Correlazione_Preselezioni(ntupMCTrig,i);
      FluxFactorizationtest_Pre_Fill(ntupMCTrig,i);
      DVSMCPreSeleff_Fill(ntupMCTrig,i);
      DVSMCPreSeleffD_Fill(ntupMCTrig,i);
   }
   return;
}




void LoopOnMCSepD(TNtuple* ntupMCSepD)
{
   for(int i=0; i<ntupMCSepD->GetEntries(); i++) {
      ntupMCSepD->GetEvent(i);
      Massa_gen = ReturnMass_Gen();
      UpdateProgressBar(ntupMCSepD, i);
      Cuts();
      RUsed=Tup.R;

      HecutMC_Fill(ntupMCSepD,i);
      SlidesforPlot_Fill(ntupMCSepD,i);
      FluxFactorizationtest_Qual_Fill(ntupMCSepD,i);
      DistanceCut_Fill(ntupMCSepD,i);
      MCQualeff_Fill(ntupMCSepD,i);
      DVSMCQualeff2_Fill(ntupMCSepD,i);
      DVSMCQualeffD_Fill(ntupMCSepD,i);
      DVSMCRICHeff_Fill(ntupMCSepD,i);
      DeutonsMC_Fill(ntupMCSepD,i);
      DeutonsMC_Dist_Fill(ntupMCSepD,i);
      MCMC_Fill(ntupMCSepD,i);
   }
   return;
}


void LoopOnDataTrig(TNtuple* ntupDataTrig)
{
   for(int i=0; i<ntupDataTrig->GetEntries(); i++) {
      ntupDataTrig->GetEvent(i);
      Cuts_Pre();
      double geomag[12]= {0,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.3};
      for(int z=0; z<12; z++) {
         double geo= geomag[z]  ;
         double geo2=geomag[z+1];
         if(Tup.Latitude>geo && Tup.Latitude<geo2) Zona=z;
      }
      // Temporary Betarich check
      if((((int)Tup.Cutmask>>11)==0||((int)Tup.Cutmask>>11)==512)&&Tup.BetaRICH<0) continue;
      RUsed=Tup.R_pre;
      UpdateProgressBar(ntupDataTrig, i);

      DATAUnbiaseff_Fill(ntupDataTrig,i);
      DATApreSeleff_Fill(ntupDataTrig,i,Zona);
      DVSMCPreSeleff_D_Fill(ntupDataTrig,i,Zona);
      DVSMCPreSeleffD_D_Fill(ntupDataTrig,i,Zona);
   }
   return;
}



void LoopOnDataSepD(TNtuple* ntupDataSepD)
{
   for(int i=0; i<ntupDataSepD->GetEntries(); i++) {
      ntupDataSepD->GetEvent(i);
      double geomag[12]= {0,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.3};
      for(int z=0; z<12; z++) {
         double geo= geomag[z]  ;
         double geo2=geomag[z+1];
         if(Tup.Latitude>geo && Tup.Latitude<geo2) Zona=z;
      }
      // Temporary Betarich check
      if((((int)Tup.Cutmask>>11)==0||((int)Tup.Cutmask>>11)==512)&&Tup.BetaRICH<0) continue;
      Cuts();
      RUsed=Tup.R;
      UpdateProgressBar(ntupDataSepD, i);


      HecutD_Fill(ntupDataSepD,i);
      SlidesforPlot_D_Fill(ntupDataSepD,i);
      DATAQualeff_Fill(ntupDataSepD,i,Zona);
      DATARICHeff_Fill(ntupDataSepD,i,Zona);
      ProtonFlux_Fill(ntupDataSepD,i,Zona);
      DVSMCQualeff2_D_Fill(ntupDataSepD,i,Zona);
      DVSMCQualeffD_D_Fill(ntupDataSepD,i,Zona);
      DVSMCRICHeff_D_Fill(ntupDataSepD,i,Zona);
      DeutonsDATA_Fill(ntupDataSepD,i,Zona);
      DeutonsDATA_Dist_Fill(ntupDataSepD,i,Zona);
      MCMCDATA_Fill(ntupDataSepD,i);
   }
   return;
}

void UpdateProgressBar(TNtuple* ntuple, int currentevent) {
   long nentries=ntuple->GetEntries();
   int newratio =(int)100*(currentevent/    (float)nentries);
   int oldratio =(int)100*((currentevent-1)/(float)nentries);
   if(newratio>oldratio) 
         cout<<'\r' << "Progress : "<< newratio << " %"<< flush;
}
