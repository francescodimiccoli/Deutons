using namespace std;


LATcorr *LATpreSelDATA = new LATcorr("LATpreSelDATA",3);


void DATApreSeleff_Fill(int zona)
{

   if(Tup.Unbias!=0||Tup.R_pre<=0||Tup.Beta_pre>protons->Eval(Tup.R_pre)+0.1||Tup.Beta_pre<protons->Eval(Tup.R_pre)-0.1) return;
   if(!(Tup.EdepL1>0&&Tup.EdepL1<EdepL1beta->Eval(Tup.Beta_pre)+0.1&&Tup.EdepL1>EdepL1beta->Eval(Tup.Beta_pre)-0.1)) return;
   if(!Herejcut) return;  
 
   for(int S=0; S<3; S++) {
      int Kbin=PRB.GetRBin(fabs(Tup.R_pre));
      if(cmask.notPassed(S)) ((TH3 *)LATpreSelDATA->beforeR)->Fill(Kbin,zona,S);
      if(cmask.passed(S))    ((TH3 *)LATpreSelDATA->afterR )->Fill(Kbin,zona,S);
   }
   return;
}


void DATApreSeleff_Write()
{
   LATpreSelDATA->Write();
   return;
}


void DATApreSeleff(string histoName)
{
   inputHistoFile->cd();
   LATcorr * LATpreSelDATA = new LATcorr(inputHistoFile,"LATpreSelDATA",3);

   string tagli[3]= {"Matching TOF","Chi^2 R","1 Tr. Track"};
   string nome;

   cout<<"********************** DATA PRESELECTIONS EFFICIENCIES ****************************************"<<endl;

   LATpreSelDATA -> Eval_Efficiency();
   TH3F *LATpreSelDATA_R = (TH3F *) LATpreSelDATA     -> effR -> Clone();
   cout<<"********************** LAT. Eff. CORRECTION **************************************************"<<endl;

   LATpreSelDATA -> Eval_LATcorr(3);

   TH2F * preSelLATcorr = (TH2F *) LATpreSelDATA   -> LATcorrR -> Clone();

   TH2F * preSelLATcorr_fit = (TH2F *) LATpreSelDATA   -> LATcorrR_fit -> Clone();

   cout<<"*** Updating P1 file ****"<<endl;

   //inputHistoFile->Close();
   //delete inputHistoFile;
   inputHistoFile=TFile::Open(histoName.data(),  "UPDATE");

   inputHistoFile->cd("Results");
//inputHistoFile->ls();
   cout << LATpreSelDATA_R ->GetName() << " : " <<  LATpreSelDATA_R ->GetEntries() << endl;
   LATpreSelDATA_R  -> Write();
   preSelLATcorr    -> Write();
   preSelLATcorr_fit-> Write();
   inputHistoFile->Write();
	inputHistoFile->Close();


   TGraphErrors * Eff_preSelLAT[3][11];
   TCanvas *c14[4];
   for(int S=0; S<3; S++) {
      nome="Latitude Efficiency: "+tagli[S];
      c14[S]=new TCanvas(nome.c_str());
      c14[S]->Divide(2,1);
      c14[S]->cd(1);
      gPad->SetLogx();
      gPad->SetGridx();
      gPad->SetGridy();
      for(int l=1; l<11; l++) {
         Eff_preSelLAT[S][l]=new TGraphErrors();
         int point=0;
         for(int i=1; i<nbinsr; i++) {
            if(LATpreSelDATA_R->GetBinContent(i+1,l+1,S+1)>0) {
               Eff_preSelLAT[S][l]->SetPoint(point,PRB.RigBinCent(i),LATpreSelDATA_R->GetBinContent(i+1,l+1,S+1));
               Eff_preSelLAT[S][l]->SetPointError(point,0,LATpreSelDATA_R->GetBinError(i+1,l+1,S+1));
               point++;
            }
         }
      }
      Eff_preSelLAT[S][10]->SetMarkerColor(1);
      Eff_preSelLAT[S][10]->SetLineColor(1);
      Eff_preSelLAT[S][10]->SetMarkerStyle(8);
      Eff_preSelLAT[S][10]->GetXaxis()->SetTitle("R [GV]");
      Eff_preSelLAT[S][10]->GetYaxis()->SetTitle("Efficiency");
      Eff_preSelLAT[S][10]->GetYaxis()->SetRangeUser(0.1,1);
      Eff_preSelLAT[S][10]->Draw("AP");
      for(int l=1; l<10; l++) {
         Eff_preSelLAT[S][l]->SetMarkerColor(l);
         Eff_preSelLAT[S][l]->SetLineColor(l);
         Eff_preSelLAT[S][l]->SetMarkerStyle(8);
         Eff_preSelLAT[S][l]->Draw("Psame");
      }
   }



   TGraphErrors *CorrLATpre[3];
   TGraphErrors *CorrLATpre_Spl[3];

   for(int S=0; S<3; S++) {
      c14[S]->cd(2);
      gPad->SetGridy();
      gPad->SetGridx();
      nome="Latitude Efficiency: "+tagli[S];
      CorrLATpre[S]=new TGraphErrors();
      CorrLATpre[S]->SetTitle("Latitude Efficiency Corr.");
      CorrLATpre[S]->GetXaxis()->SetTitle("Latitude");
      CorrLATpre[S]->GetYaxis()->SetTitle("Eff. Corr. Factor");
      CorrLATpre[S]->SetMarkerStyle(8);
      for(int i=1; i<11; i++) {
         CorrLATpre[S]->SetPoint(i-1,geomagC[i],preSelLATcorr->GetBinContent(i+1,S+1));
         CorrLATpre[S]->SetPointError(i-1,0,preSelLATcorr->GetBinError(i+1,S+1));
      }
      CorrLATpre[S]->Fit(nome.c_str());
      CorrLATpre[S]->Draw("AP");

      nome="CorrLATpre_spl"+tagli[S];
      CorrLATpre_Spl[S]=new TGraphErrors(11);
      CorrLATpre_Spl[S]->SetName(tagli[S].c_str());
      int j=0;
      for(int i=1; i<11; i++) {
         CorrLATpre_Spl[S]->SetPoint(j,geomagC[i],preSelLATcorr_fit->GetBinContent(i+1,S+1));
         CorrLATpre_Spl[S]->SetPointError(j,0,preSelLATcorr_fit->GetBinError(i+1,S+1));
         j++;
      }
      CorrLATpre_Spl[S]->SetLineColor(2);
      CorrLATpre_Spl[S]->SetMarkerColor(2);
      CorrLATpre_Spl[S]->SetFillColor(2);
      CorrLATpre_Spl[S]->SetFillStyle(3001);
      CorrLATpre_Spl[S]->SetLineWidth(2);
      CorrLATpre_Spl[S]->Draw("PCsame");
   }

   
   for(int S=0; S<3; S++) finalPlots.Add(c14[S]);
   finalPlots.writeObjsInFolder("DATA-driven Results/Latitude effect/Clean-event Selections");
   for(int S=0; S<3; S++) finalPlots.Add(CorrLATpre_Spl[S]);
   finalPlots.writeObjsInFolder("Export");
   

}
