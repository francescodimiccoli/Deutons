using namespace std;

void FillIstogram(int INDX,string frac,string mese)
{


   float fraz=1;
   float Zona=0;
   double geomag[12]= {0,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.3};
   int progress=0;

   cout<<"*********************** CALIB. READING *********************"<<endl;

   string nomecal=inputpath + "/CodesforAnalysis/CALIBRAZIONI/"+mese+".root";
   TFile *calib = TFile::Open(nomecal.c_str());
   if(calib) cout<<"MC calibration for month "<<mese<<" ... ok"<<endl;
   else cout<<"ERROR: MC calibration not found"<<endl;
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
      ntupMCSepD=(TNtuple*)fileMC->Get("grandezzesepd");
      ntupMCTrig=(TNtuple*)fileMC->Get("trig");

      nomefile=inputpath+"/Risultati/"+mese+"/RisultatiDATI_"+frac+".root";
      fileData =TFile::Open(nomefile.c_str());
      ntupDataSepD=(TNtuple*)fileData->Get("grandezzesepd");
      ntupDataTrig=(TNtuple*)fileData->Get("trig");

      ntupMCTrig->SetBranchAddress("Momento_gen",&Momento_gen);
      ntupMCTrig->SetBranchAddress("Ev_Num",&Ev_Num);
      ntupMCTrig->SetBranchAddress("Trig_Num",&Trig_Num);
      ntupMCTrig->SetBranchAddress("R_pre",&R_pre);
      ntupMCTrig->SetBranchAddress("Beta_pre",&Beta_pre);
      ntupMCTrig->SetBranchAddress("Cutmask",&CUTMASK);
      ntupMCTrig->SetBranchAddress("MC_type",&MC_type);
      ntupMCTrig->SetBranchAddress("EdepL1",&EdepL1);
      ntupMCTrig->SetBranchAddress("EdepTOFU",&EdepTOFU);
      ntupMCTrig->SetBranchAddress("EdepTOFD",&EdepTOFD);
      ntupMCTrig->SetBranchAddress("EdepTrack",&EdepTrack);
      ntupMCTrig->SetBranchAddress("EdepECAL",&EdepECAL);
      ntupMCTrig->SetBranchAddress("BetaRICH",&BetaRICH);
      ntupMCTrig->SetBranchAddress("Unbias",&Unbias);

      ntupMCSepD->SetBranchAddress("Momentogen",&Momento_gen);
      ntupMCSepD->SetBranchAddress("R",&R);
      ntupMCSepD->SetBranchAddress("Beta",&Beta);
      ntupMCSepD->SetBranchAddress("BetaRICH_new",&BetaRICH);
      ntupMCSepD->SetBranchAddress("EdepL1",&EdepL1);
      ntupMCSepD->SetBranchAddress("Rmin",&Rmin);
      ntupMCSepD->SetBranchAddress("EdepTOF",&EdepTOFU);
      ntupMCSepD->SetBranchAddress("EdepTrack",&EdepTrack);
      ntupMCSepD->SetBranchAddress("EdepTOFD",&EdepTOFD);
      ntupMCSepD->SetBranchAddress("MC_type",&MC_type);
      ntupMCSepD->SetBranchAddress("LDiscriminant",&LDiscriminant);
      ntupMCSepD->SetBranchAddress("BDT_response",&BDT_response);
      ntupMCSepD->SetBranchAddress("Cutmask",&CUTMASK);
      ntupMCSepD->SetBranchAddress("Dist5D",&Dist5D);
      ntupMCSepD->SetBranchAddress("Dist5D_P",&Dist5D_P);

      ntupDataTrig->SetBranchAddress("Rcutoff",&Rcutoff);
      ntupDataTrig->SetBranchAddress("R_pre",&R_pre);
      ntupDataTrig->SetBranchAddress("Beta_pre",&Beta_pre);
      ntupDataTrig->SetBranchAddress("Cutmask",&CUTMASK);
      ntupDataTrig->SetBranchAddress("Latitude",&Latitude);
      ntupDataTrig->SetBranchAddress("EdepECAL",&EdepECAL);
      ntupDataTrig->SetBranchAddress("EdepL1",&EdepL1);
      ntupDataTrig->SetBranchAddress("EdepTOFU",&EdepTOFU);
      ntupDataTrig->SetBranchAddress("EdepTOFD",&EdepTOFD);
      ntupDataTrig->SetBranchAddress("EdepTrack",&EdepTrack);
      ntupDataTrig->SetBranchAddress("BetaRICH",&BetaRICH);
      ntupDataTrig->SetBranchAddress("Unbias",&Unbias);


      ntupDataSepD->SetBranchAddress("R",&R);
      ntupDataSepD->SetBranchAddress("Beta",&Beta);
      ntupDataSepD->SetBranchAddress("BetaRICH_new",&BetaRICH);
      ntupDataSepD->SetBranchAddress("EdepL1",&EdepL1);
      ntupDataSepD->SetBranchAddress("Rcutoff",&Rcutoff);
      ntupDataSepD->SetBranchAddress("Rmin",&Rmin);
      ntupDataSepD->SetBranchAddress("EdepTOFU",&EdepTOFU);
      ntupDataSepD->SetBranchAddress("EdepTrack",&EdepTrack);
      ntupDataSepD->SetBranchAddress("EdepTOFD",&EdepTOFD);
      ntupDataSepD->SetBranchAddress("Latitude",&Latitude);
      ntupDataSepD->SetBranchAddress("LDiscriminant",&LDiscriminant);
      ntupDataSepD->SetBranchAddress("BDT_response",&BDT_response);
      ntupDataSepD->SetBranchAddress("Cutmask",&CUTMASK);
      ntupDataSepD->SetBranchAddress("Dist5D",&Dist5D);
      ntupDataSepD->SetBranchAddress("Dist5D_P",&Dist5D_P);
   }
   cout<<"*********************** MC READING *********************"<<endl;
   if(INDX==0)
      for(int i=0; i<ntupMCTrig->GetEntries()/fraz; i++) {
         ntupMCTrig->GetEvent(i);
         Cutmask=CUTMASK;
         Cuts_Pre();
         Massa_gen = ReturnMass_Gen();
         Momento_gen=Momento_gen;
         Var=R_pre;
         Var2=R_pre;

         if(100*(i/(float)(ntupMCTrig->GetEntries()/fraz))>progress) {
            cout<<'\r' << "Progress : "<<progress << " %"<< flush;
            progress=(int)(100*(i/(float)(ntupMCTrig->GetEntries()/fraz)))+1;
         }
         Beta_gen=(pow(pow(Momento_gen/Massa_gen,2)/(1+pow(Momento_gen/Massa_gen,2)),0.5));
         MCpreseff_Fill(ntupMCTrig,i);
         MCUnbiaseff_Fill(ntupMCTrig,i);
         MCTrackeff_Fill(ntupMCTrig,i);
         MigrationMatrix_Fill(ntupMCTrig,i);
         Correlazione_Preselezioni(ntupMCTrig,i);
         FluxFactorizationtest_Pre_Fill(ntupMCTrig,i);
         DVSMCPreSeleff_Fill(ntupMCTrig,i);
         DVSMCPreSeleffD_Fill(ntupMCTrig,i);
         //DVSMCTrackeff_Fill(ntupMCTrig,i);*/
      }
   if(INDX==0||INDX==1) {
      progress=0;
      for(int i=0; i<ntupMCSepD->GetEntries(); i++) {
         ntupMCSepD->GetEvent(i);
         Cutmask=CUTMASK;
         Massa_gen = ReturnMass_Gen();
         if(100*(i/(float)(ntupMCSepD->GetEntries()))>progress) {
            cout<<'\r' << "Progress : "<<progress << " %"<< flush;
            progress=(int)(100*(i/(float)(ntupMCSepD->GetEntries()/fraz)))+1;
         }
         Cuts();
         Var=R;
         Var2=R;

         HecutMC_Fill(ntupMCSepD,i);
         SlidesforPlot_Fill(ntupMCSepD,i);
         FluxFactorizationtest_Qual_Fill(ntupMCTrig,i);
         DistanceCut_Fill(ntupMCSepD,i);
         MCQualeff_Fill(ntupMCSepD,i);
         DVSMCQualeff2_Fill(ntupMCSepD,i);
         DVSMCQualeffD_Fill(ntupMCSepD,i);
         DVSMCRICHeff_Fill(ntupMCSepD,i);
         DeutonsMC_Fill(ntupMCSepD,i);
         DeutonsMC_Dist_Fill(ntupMCSepD,i);
         MCMC_Fill(ntupMCSepD,i);
      }
   }
   cout<<endl<<"*********************** DATA READING *********************"<<endl;
   if(INDX!=2) {
      Tempi = (TH1F *)fileData->Get("Tempi");
      esposizionegeo = (TH2F *)fileData->Get("esposizionegeo");
      esposizionepgeo = (TH2F*)fileData->Get("esposizionepgeo");
      esposizionepgeoNaF = (TH2F*)fileData->Get("esposizionepgeoNaF");
      esposizionepgeoAgl = (TH2F*)fileData->Get("esposizionepgeoAgl");
      esposizionedgeo = (TH2F*)fileData->Get("esposizionedgeo");
      esposizionedgeoNaF = (TH2F*)fileData->Get("esposizionedgeoNaF");
      esposizionedgeoAgl = (TH2F*)fileData->Get("esposizionedgeoAgl");
   } else {
      Tempi = (TH1F *)file->Get("Tempi");
      esposizionegeo = (TH2F *)file->Get("esposizionegeo");
      esposizionepgeo = (TH2F*)file->Get("esposizionepgeo");
      esposizionepgeoNaF = (TH2F*)file->Get("esposizionepgeoNaF");
      esposizionepgeoAgl = (TH2F*)file->Get("esposizionepgeoAgl");
      esposizionedgeo = (TH2F*)file->Get("esposizionedgeo");
      esposizionedgeoNaF = (TH2F*)file->Get("esposizionedgeoNaF");
      esposizionedgeoAgl = (TH2F*)file->Get("esposizionedgeoAgl");
   }
   progress=0;
   if(INDX==0)
      for(int i=0; i<ntupDataTrig->GetEntries()/fraz; i++) {
         ntupDataTrig->GetEvent(i);
         Cutmask=CUTMASK;
         Cuts_Pre();
         Massa=pow(fabs(pow(fabs(R_pre)*pow((1-pow(Beta_pre,2)),0.5)/Beta_pre,2)),0.5);
         for(int z=0; z<12; z++) {
            double geo= geomag[z]  ;
            double geo2=geomag[z+1];
            if(Latitude>geo && Latitude<geo2) Zona=z;
         }
         // Temporary Betarich check
         if(((Cutmask>>11)==0||(Cutmask>>11)==512)&&BetaRICH<0) continue;
         Var=R_pre;
         Var2=R_pre;

         if(100*(i/(float)(ntupDataTrig->GetEntries()/fraz))>progress) {
            cout<<'\r' << "Progress : "<<progress << " %"<< flush;
            progress=(int)(100*(i/(float)(ntupDataTrig->GetEntries()/fraz)))+1;
         }

         Beta_gen=(pow(pow(Momento_gen/Massa_gen,2)/(1+pow(Momento_gen/Massa_gen,2)),0.5));
         DATAUnbiaseff_Fill(ntupDataTrig,i);
         DATApreSeleff_Fill(ntupDataTrig,i,Zona);
         DVSMCPreSeleff_D_Fill(ntupDataTrig,i,Zona);
         DVSMCPreSeleffD_D_Fill(ntupDataTrig,i,Zona);
         //DVSMCTrackeff_D_Fill(ntupMCTrig,i);
      }
   if(INDX==0||INDX==1) {
      progress=0;
      for(int i=0; i<ntupDataSepD->GetEntries(); i++) {
         ntupDataSepD->GetEvent(i);
         Cutmask=CUTMASK;
         for(int z=0; z<12; z++) {
            double geo= geomag[z]  ;
            double geo2=geomag[z+1];
            if(Latitude>geo && Latitude<geo2) Zona=z;
         }
         // Temporary Betarich check
         if(((Cutmask>>11)==0||(Cutmask>>11)==512)&&BetaRICH<0) continue;
         Cuts();
         Var=R;
         Var2=R;

         if(100*(i/(float)(ntupDataSepD->GetEntries()))>progress) {
            cout<<'\r' << "Progress : "<<progress << " %"<< flush;
            progress=(int)(100*(i/(float)(ntupDataSepD->GetEntries()/fraz)))+1;
         }
         HecutD_Fill(ntupDataSepD,i);
         SlidesforPlot_D_Fill(ntupMCSepD,i);
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
