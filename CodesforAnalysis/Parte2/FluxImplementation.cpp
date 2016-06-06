using namespace std;


void FormatGraphPDRatios(TGraphErrors* pdratio);
TGraph* GetGalpropRatioFromFile(string filename);
TGraph* GetPowerWeightedGalpropRatioFromFile(string filename, float power, float yweight=1e7);
void FormatGraphGalprop(TGraph* graph);
void FormatGraphGalpropRatio(TGraph* graph);
void FormatGraphGalpropXaxis(TGraph* graph);
void FormatPadAndPlotGalpropRatios(TGraph* graph1, TGraph* graph2);
void SetDeutonsExposureGraphStyle(TGraphErrors* graph);
void SetProtonsExposureGraphStyle(TGraphErrors* graph);
void FormatDeutonsExposureGraphTitle(TGraphErrors* graph);
void FormatDFluxGeoWithDetZone(TGraphErrors* graph, int detector, int geozone) ;
TCanvas* MakeLogGridCanvas(string title);
TGraphErrors* FillGraphErrorFromBinningAndHisto(Binning bin, TH1* histo);


const int NGEOBINS=11;
const int NSUBDETECTORS=3;


void TGraphErrorsFormat(TGraphErrors* graph)
{
   graph->SetMarkerStyle(8);
   graph->GetXaxis()->SetTitleSize(0.045);
   graph->GetYaxis()->SetTitleSize(0.045);
   return;
}



Flux * P_Flux         = new Flux("P_Flux"    	  , "", RB  );
Flux * P_Flux_geo     = new Flux("P_Flux_geo"     , "", RB, NGEOBINS);
Flux * P_Flux_geo_prim= new Flux("P_Flux_geo_prim", "", RB, NGEOBINS);

//DvsMC check
Flux * P_Flux_pre = new Flux("P_Flux_pre" , "", RB);
Flux * P_Flux_sel = new Flux("P_Flux_sel" , "", RB);

void ProtonFlux_Fill(TNtuple *ntupla, int l,int zona)
{
   ntupla->GetEvent(l);
   if(Tup.Beta<=0||Tup.R<=0) return;
   int Kbin=RB.GetRBin(Tup.R);

   if(Tup.Dist5D_P<6 && Likcut) {
      P_Flux_geo-> Counts -> Fill(Kbin,zona);
      if(Tup.R>1.2*Tup.Rcutoff) {
         P_Flux -> Counts-> Fill(Kbin);
         P_Flux_geo_prim -> Counts -> Fill(Kbin,zona);
      }
   }

   if(Herejcut && Tup.R>1.2*Tup.Rcutoff) {
      P_Flux_pre -> Counts -> Fill(Kbin);
      if(Tup.Dist5D_P<6&&Likcut)  P_Flux_sel -> Counts -> Fill(Kbin);
   }

   return;
}


void ProtonFlux_Write()
{
   P_Flux      	->Write();
   P_Flux_geo  	->Write();
   P_Flux_geo_prim ->Write();
   P_Flux_pre  	->Write();
   P_Flux_sel  	->Write();

   return;
}




class PFluxComputing {
   public:
      PFluxComputing(string objname, bool geobin=1, string acceptname="Corr_AcceptanceP");
      TH1* ComputeFluxAndGetHisto();
   private:
      int ngeobins;
      Flux* PFlux;
      TH2F* H2ExpoGeo;
};

PFluxComputing::PFluxComputing(string objname, bool geobin, string acceptname)
{
   string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
   TFile * file1 = TFile::Open(nomefile.c_str(),"READ");
   ngeobins=1;
   if (geobin) {
      ngeobins=NGEOBINS;
      acceptname="Geomag_AcceptanceP";
      H2ExpoGeo=(TH2F*)file1->Get("esposizionegeo")->Clone();
   }
   PFlux=new Flux(file1, objname, "","Results", acceptname, ngeobins);
}


TH1* PFluxComputing::ComputeFluxAndGetHisto()
{
   if (ngeobins)  PFlux->Set_Exposure_Time (H2ExpoGeo);
   else           PFlux->Set_Exposure_Time (Tempi); // whatever it may be
   enum {protons, deutons};
   PFlux->Eval_Flux(ngeobins, protons);
   return PFlux->Fluxes;
}


void ProtonFlux()
{


   PFluxComputing pfBase("P_Flux",     0);
   PFluxComputing pfSel ("P_Flux_sel", 0);
   PFluxComputing pfPre ("P_Flux_pre", 0, "Corr_AcceptancePreP");
   PFluxComputing pfGeo ("P_Flux_geo");
   PFluxComputing pfPrim("P_Flux_geo_prim");



   cout<<"*************** PROTONS FLUXES CALCULATION *******************"<<endl;


   TH1F * ProtonsPrimaryFlux = (TH1F *)  pfBase.ComputeFluxAndGetHisto();
   TH2F * ProtonsGeomagFlux  = (TH2F *)  pfGeo .ComputeFluxAndGetHisto();
   TH1F * P_pre_PrimaryFlux  = (TH1F *)  pfPre .ComputeFluxAndGetHisto(); //P flux with pre-selections alone
   TH1F * P_sel_PrimaryFlux  = (TH1F *)  pfSel .ComputeFluxAndGetHisto(); //P flux with full-set selections


   cout<<"*** Updating P1 file ****"<<endl;

   string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
   TFile* file1 = TFile::Open(nomefile.c_str(),"UPDATE");
   file1->mkdir("Results/Fluxes");
   file1->cd("Results/Fluxes");

   ProtonsPrimaryFlux  ->Write("ProtonsPrimaryFlux" 	);
   ProtonsGeomagFlux   ->Write("ProtonsGeomagFlux"  	);
   P_pre_PrimaryFlux   ->Write("P_pre_PrimaryFlux"  	);
   P_sel_PrimaryFlux   ->Write("P_sel_PrimaryFlux"  	);

   file1->Write();
   file1->Close();


   TGraphErrors * P_Fluxgeo[NGEOBINS];
   TGraphErrors * PFlux;
   TGraphErrors * PFlux_pre;
   float potenza=0;



   PBinning PRB;
   PRB.Setbins(nbinsr, 0.5, 100, 2); // RB did not have Ek


   for(int j=0; j<NGEOBINS; j++) {
      P_Fluxgeo[j]=new TGraphErrors();
      TGraphErrorsFormat(P_Fluxgeo[j]);
      string nome="Protons Flux: Geo. Zone "+to_string(j);
      P_Fluxgeo[j]->SetName(nome.c_str());
      for(int i=0; i<PRB.EkBinsCent().size(); i++) {
         float ekin=PRB.EkBinCent(i);
         P_Fluxgeo[j]->SetPoint(i,ekin,ProtonsGeomagFlux->GetBinContent(i+1,j+1)*pow(ekin,potenza));
         P_Fluxgeo[j]->SetPointError(i,0,ProtonsGeomagFlux->GetBinError(i+1,j+1)*pow(ekin,potenza));
      }
      P_Fluxgeo[j]->SetMarkerColor(j-1);
      P_Fluxgeo[j]->SetLineColor(j-1);
      P_Fluxgeo[j]->SetLineWidth(2);
   }
   P_Fluxgeo[10]->SetTitle("Protons Flux: Geo. Zones;Kin. En./nucl. [GeV/nucl.];Flux [(m^2 sr GeV/nucl.)^-1]");
   P_Fluxgeo[10]->GetYaxis()->SetRangeUser(1e-2,1e4);

   TCanvas * c23 = MakeLogGridCanvas("Protons Flux: Geo. Zones");
   P_Fluxgeo[10]->Draw("AP");
   for(int j=0; j<10; j++) P_Fluxgeo[j]->Draw("Psame");



   PFlux=new TGraphErrors();
   for(int i=0; i<PRB.EkBinsCent().size(); i++) {
      float ekin=PRB.EkBinCent(i);
      PFlux->SetPoint(i,ekin,ProtonsPrimaryFlux->GetBinContent(i+1)*pow(ekin,potenza));
      PFlux->SetPointError(i,0,ProtonsPrimaryFlux->GetBinError(i+1)*pow(ekin,potenza));
   }
   PFlux->SetName("Protons Primary Flux");
   TGraphErrorsFormat(PFlux);
   PFlux->SetMarkerColor(2);
   PFlux->SetTitle("Primary Protons Flux;Kin. En./nucl. [GeV/nucl.];Flux [(m^2 sr GeV/nucl.)^-1]");
   PFlux->GetYaxis()->SetRangeUser(1e-2,1e4);

   TCanvas * c24 = MakeLogGridCanvas("Primary Protons Flux");
   PFlux->Draw("AP");
   P_Fluxgeo[10]->Draw("Psame");


   TGraph* galprop3P=GetPowerWeightedGalpropRatioFromFile("./Galprop/Trotta2011/Def/new_P200.txt",   potenza);
   TGraph* galprop3P2=GetPowerWeightedGalpropRatioFromFile("./Galprop/Trotta2011/Def/new_P1250.txt", potenza);
   galprop3P->Draw("sameC");
   galprop3P2->Draw("sameC");


   TGraphErrors * PFluxpre=new TGraphErrors();
   TGraphErrors * PFlux_=new TGraphErrors();
   PFlux_pre=new TGraphErrors();
   PFlux_pre->SetName("Protons Primary Flux (only pres.)");



   for(int ip=0; ip<PRB.EkBinsCent().size()-1; ip++) {
      PFlux_->   SetPoint(ip,PRB.RigBinCent(ip), 1);
      PFlux_pre->SetPoint(ip,PRB.EkBinCent(ip) , 1);
      float err=(P_sel_PrimaryFlux->GetBinError(ip+1,2)+P_pre_PrimaryFlux->GetBinError(ip+1,2))/P_pre_PrimaryFlux->GetBinContent(ip+1,1);
      PFlux_->SetPointError(ip,0,err);
      PFluxpre->SetPoint(ip,PRB.RigBinCent(ip+1),
                         P_sel_PrimaryFlux->GetBinContent(ip+1)/P_pre_PrimaryFlux->GetBinContent(ip+1,1));

   }

   PFlux_->SetMarkerStyle(4);
   PFlux_->SetFillStyle(3002);
   PFlux_->SetFillColor(4);
   PFlux_->SetMarkerColor(2);

   TGraphErrorsFormat(PFluxpre);
   PFluxpre->SetMarkerColor(2);
   PFluxpre->SetTitle("Primary Protons Flux;R [GV];Flux [(m^2 sr GeV/nucl.)^-1]");
   PFluxpre->GetYaxis()->SetRangeUser(0.8,1.2);




   TCanvas * c25 = new  TCanvas("Protons Flux: Pre vs Qual");
   c25->cd();
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   PFluxpre->Draw("AP");
   PFlux_->Draw("P4same");



   cout<<"*** Updating Results file ***"<<endl;
   nomefile="./Final_plots/"+mese+".root";
   TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
   f_out->mkdir("P Fluxes");
   f_out->cd("P Fluxes");
   c23->Write();
   c24->Write();
   c25->Write();
   f_out->cd("Export");
   PFlux->Write("Protons Primary Flux");
   f_out->Write();
   f_out->Close();

   return;
}


void DeutonFlux()
{
   string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
   TFile * file1 = TFile::Open(nomefile.c_str(),"READ");




   cout<<"*************** DEUTONS FLUXES CALCULATION *******************"<<endl;

   //Evaluation additional Fit Error


   std::string suf[4] = {"ToF", "NaF", "Agl", "R"};
   Binning binP[4] = {ToFPB, NaFPB, AglPB, RB};
   Binning binD[4] = {ToFDB, NaFDB, AglDB, RB};
   string expodnames[4] = {"esposizionedgeo", "esposizionedgeoNaF", "esposizionedgeoAgl", "esposizionegeo" };
   string expopnames[4] = {"esposizionepgeo", "esposizionepgeoNaF", "esposizionepgeoAgl", "esposizionegeo" };
   enum {protons, deutons};

   std::vector <TH1F*> vD_Flux (NSUBDETECTORS);
   std::vector <TH1F*> vD_Dist (NSUBDETECTORS);
   std::vector <TH1F*> vP_Flux (NSUBDETECTORS);
   std::vector <TH1F*> vP_Dist (NSUBDETECTORS);
   std::vector <TH2F*> vG_Flux (NSUBDETECTORS);
   std::vector <TH2F*> vG_Dist (NSUBDETECTORS);
   std::vector <TH1F*> vD_Exp  (NSUBDETECTORS);
   std::vector <TH1F*> vP_Exp  (NSUBDETECTORS);



   for (int ifx=0; ifx<4; ifx++) { // index flux

      //Deutons
      Flux * D_Flux       = new Flux(file1, "D_Flux"       ,"Results", suf[ifx], "Corr_AcceptanceD",binD[ifx]);
      Flux * D_Flux_Dist  = new Flux(file1, "D_Flux_Dist"  ,"Results", suf[ifx], "Corr_AcceptanceD",binD[ifx]);
      if(D_Flux -> Counts) {
         TH1F * Syst=(TH1F*) D_Flux -> Counts ->Clone();
         Syst -> Add ((TH1F*) D_Flux_Dist -> Counts ->Clone(), -1);
         // Checking systematics
         for(int i=0; i<Syst->GetNbinsX(); i++) {
            if ( ((TH1F*)D_Flux_Dist-> Counts )-> GetBinContent(i)==0
                  || ((TH1F*)D_Flux     -> Counts )-> GetBinContent(i)==0 )
               Syst -> SetBinContent (i,0);
         }
         D_Flux -> Add_SystFitError(1,Syst);
      }
      TH2F * expogeo  = (TH2F*)file1->Get(expodnames[ifx].data());
      D_Flux         -> Set_Exposure_Time (expogeo);
      D_Flux         -> Eval_Flux(1 , deutons, 2 );
      D_Flux_Dist    -> Set_Exposure_Time (expogeo);
      D_Flux_Dist    -> Eval_Flux(1 , deutons, 2 );

      // Protons
      Flux * P_Flux         = new Flux(file1, "P_Flux"        ,"Results", suf[ifx], "Corr_AcceptanceP",binD[ifx]);
      Flux * P_Flux_Dist    = new Flux(file1, "P_Flux_Dist"   ,"Results", suf[ifx], "Corr_AcceptanceP",binD[ifx]);
      expogeo  = (TH2F*)file1->Get(expopnames[ifx].data());
      P_Flux         -> Set_Exposure_Time (expogeo);
      P_Flux_Dist    -> Set_Exposure_Time (expogeo);
      P_Flux         -> Eval_Flux(1 , protons );
      P_Flux_Dist    -> Eval_Flux(1 , protons );

      Flux * D_Flux_geo     = new Flux(file1, "D_Flux_geo"           ,"Results", suf[ifx],"Geomag_AcceptanceD",binD[ifx]);
      Flux * D_Flux_geo_Dist     = new Flux(file1, "D_Flux_geo_Dist" ,"Results", suf[ifx],"Geomag_AcceptanceD",binD[ifx]);

      // Flux geo
      D_Flux_geo     -> Set_Exposure_Time (Tempi);
      D_Flux_geo_Dist-> Set_Exposure_Time (Tempi);
      D_Flux_geo     -> Eval_Flux(NGEOBINS, deutons, 2 );
      D_Flux_geo_Dist-> Eval_Flux(NGEOBINS, deutons, 2 );

      if (ifx<NSUBDETECTORS) { // No detector for Rigidity
         vD_Flux .push_back( (TH1F*)  D_Flux        ->Fluxes);
         vD_Dist.push_back( (TH1F*)  D_Flux_Dist    ->Fluxes);
         vP_Flux .push_back( (TH1F*)  P_Flux        ->Fluxes);
         vP_Dist.push_back( (TH1F*)  P_Flux_Dist    ->Fluxes);
         vG_Flux .push_back( (TH2F*) D_Flux_geo     ->Fluxes);
         vG_Dist.push_back( (TH2F*)  D_Flux_geo_Dist->Fluxes);
         vD_Exp.push_back( (TH1F*) D_Flux -> Exposure);
         vP_Exp.push_back( (TH1F*) D_Flux -> Exposure);


      }

   }


   //// Writing the files

   cout<<"*** Updating P1 file ****"<<endl;

   nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
   file1 = TFile::Open(nomefile.c_str(),"UPDATE");
   file1->cd("Results/Fluxes");


   TH1F * DP_ratio[NSUBDETECTORS];
   TH1F * DP_ratio_Dist[NSUBDETECTORS];

   for (int ifx=0; ifx<NSUBDETECTORS; ifx++) {
      //Fit on Mass
      vD_Flux[ifx] 	->Write(("DeutonsPrimaryFlux_"+suf[ifx]).data());
      vG_Flux[ifx]  	->Write(("DeutonsGeomagFlux_"+suf[ifx]).data());
      //Fit on Distance
      vD_Dist[ifx]	->Write(("DeutonsPrimaryFlux_Dist_"+suf[ifx]).data());
      vG_Dist[ifx] 	->Write(("DeutonsGeomagFlux_Dist_"+suf[ifx]).data());
      //for D/P ratio: P Flux
      vP_Flux[ifx] 	->Write(("ProtonsPrimaryFlux_"+suf[ifx]).data());
      vP_Dist[ifx]	->Write(("Protons_PrimaryFlux_Dist_"+suf[ifx]).data());

      //D/P ratio
      //Fit on Mass
      DP_ratio[ifx] = (TH1F*)vD_Flux[ifx] -> Clone();
      DP_ratio[ifx] ->  Divide(vP_Flux[ifx]			);
      DP_ratio[ifx] -> Write(("DP_ratio_"+suf[ifx+1]		).data());
      //Fit on Distance
      DP_ratio_Dist[ifx] = (TH1F*)vD_Flux[ifx] -> Clone();
      DP_ratio_Dist[ifx] ->  Divide(vP_Dist[ifx] );
      DP_ratio_Dist[ifx] -> Write(("DP_ratio_Dist_"+suf[ifx+1]	).data());
   }


   file1->Write();
   file1->Close();



   //// Drawing the stuff


   float potenza=0; ///< Constant, but could be changed ie to -2.7 to plot the flux with that power


   TGraph* galprop3P= GetPowerWeightedGalpropRatioFromFile("./Galprop/Trotta2011/Def/new_D500.txt",  potenza);
   TGraph* galprop3P2=GetPowerWeightedGalpropRatioFromFile("./Galprop/Trotta2011/Def/new_D1250.txt", potenza);
   FormatGraphGalprop(galprop3P);

   TCanvas * c32 = new TCanvas("Deutons Flux: Geo. Zones");
   c32-> Divide(2,1);
   c32->cd(1);
   FormatPadAndPlotGalpropRatios(galprop3P, galprop3P2);
   c32->cd(2);
   gPad->SetTitle("Deutons Flux: Geo. Zone (Distance FIT)");
   FormatPadAndPlotGalpropRatios(galprop3P, galprop3P2);





   /// Flux Geo

   TGraphErrors * D_Fluxgeo[NSUBDETECTORS][NGEOBINS];
   TGraphErrors * D_FluxgeoDist[NSUBDETECTORS][NGEOBINS];


   for (int idet=0; idet<NSUBDETECTORS; idet++) {
      for(int izone=0; izone<NGEOBINS; izone++) {
         string nome="Deutons Flux: Geo. Zone "+to_string(izone);
         TH1D* SliceOfFluxGeo=vG_Flux[idet]->ProjectionX("", izone+1, izone+1);
         D_Fluxgeo[idet][izone]=FillGraphErrorFromBinningAndHisto(binP[idet], SliceOfFluxGeo);
         D_Fluxgeo[idet][izone]->SetName(nome.c_str());

         TH1D* SliceOfFluxGeoDist=vG_Dist[idet]->ProjectionX("", izone+1, izone+1);
         D_FluxgeoDist[idet][izone]=FillGraphErrorFromBinningAndHisto(binP[idet], SliceOfFluxGeoDist);
         D_FluxgeoDist[idet][izone]->SetName((nome + "Distance Fit").c_str());

         FormatDFluxGeoWithDetZone(D_Fluxgeo    [idet][izone], idet, izone);
         FormatDFluxGeoWithDetZone(D_FluxgeoDist[idet][izone], idet, izone);

      }
      c32->cd(1);
      D_Fluxgeo[idet][10]->Draw("Psame");
      for(int j=0; j<NGEOBINS; j++) D_Fluxgeo[idet][j]->Draw("Psame");

      c32->cd(2);
      D_FluxgeoDist[idet][10]->Draw("Psame");
      for(int j=0; j<NGEOBINS; j++) D_FluxgeoDist[idet][j]->Draw("Psame");

   }

   /// Plot exposure time


   TCanvas * c33 = new TCanvas("Exposure Time");
   c33->Divide(1,NSUBDETECTORS);

   string name[NSUBDETECTORS]= {"ToF", "RICH NaF", "RICH Agl"};
   TGraphErrors * ged[NSUBDETECTORS];
   TGraphErrors * gep[NSUBDETECTORS];

   for (int i=0; i<NSUBDETECTORS; i++) {
      ged[i]=new TGraphErrors();
      gep[i]=new TGraphErrors();
      for(int m=0; m<binD[i].EkBinsCent().size(); m++) {
         ged[i]->SetPoint(m,binD[i].EkBinCent(m), vD_Exp[i]->GetBinContent(m+1));
         gep[i]->SetPoint(m,binP[i].EkBinCent(m), vP_Exp[i]->GetBinContent(m+1));
      }
      c33->cd(i+1);
      gPad->SetLogy();
      gPad->SetGridx();
      gPad->SetGridy();

      SetProtonsExposureGraphStyle(gep[i]);
      SetDeutonsExposureGraphStyle(ged[i]);
      FormatDeutonsExposureGraphTitle(ged[i]);
      ged[i]->SetTitle((name[i]+" range").data());
      ged[i]->Draw("APC");
      gep[i]->Draw("PCsame");

   }

   /// Flux Detectors


   TCanvas * c34 = new TCanvas("Deutons Flux: Primaries");
   c34-> Divide(2,1);



   // Loop on detectors


   TGraphErrors * D_Flux    [NSUBDETECTORS];
   TGraphErrors * D_FluxDist[NSUBDETECTORS];

   for (int idet=0; idet<NSUBDETECTORS; idet++) {

      string nome="Deutons Flux: Primaries" ;
      D_Flux[idet]=    FillGraphErrorFromBinningAndHisto(binD[idet], vD_Flux[idet]);
      D_FluxDist[idet]=FillGraphErrorFromBinningAndHisto(binD[idet], vD_Dist[idet]);
      D_Flux[idet]    ->SetName(nome.c_str());
      D_FluxDist[idet]->SetName((nome + "Distance Fit").c_str());


      FormatDFluxGeoWithDetZone(D_Flux    [idet],      idet, 0);
      FormatDFluxGeoWithDetZone(D_FluxDist[idet],  idet, 0);

      c34->cd(1);
      FormatPadAndPlotGalpropRatios(galprop3P, galprop3P2);
      D_Flux[idet]->Draw("Psame");

      c34->cd(2);
      FormatPadAndPlotGalpropRatios(galprop3P, galprop3P2);
      D_FluxDist[idet]->Draw("Psame");

   }



   // Some galprop init stuff of Canvas 35

   TCanvas * c35 = new TCanvas("D/P ratio");
   c35->Divide(2,1);

   TGraph* galpropRatio500= GetGalpropRatioFromFile("./Galprop/Trotta2011/PDratio/500.dat");
   TGraph* galpropRatio1000=GetGalpropRatioFromFile("./Galprop/Trotta2011/PDratio/1000.dat");

   galpropRatio500->GetXaxis()->SetRangeUser(0.1,10);
   galpropRatio500->GetYaxis()->SetRangeUser(1e-3,1e-1);
   galpropRatio500->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
   galpropRatio500->GetYaxis()->SetTitle("Flux ratio");



   // PD_Ratio

   TGraphErrors * PD_ratio     [NSUBDETECTORS];
   TGraphErrors * PD_ratio_Dist[NSUBDETECTORS];

   for (int idet=0; idet<NSUBDETECTORS; idet++) {

      PD_ratio[idet]     =FillGraphErrorFromBinningAndHisto(binP[idet], DP_ratio[idet]     );
      PD_ratio_Dist[idet]=FillGraphErrorFromBinningAndHisto(binP[idet], DP_ratio_Dist[idet]);
      FormatGraphPDRatios(PD_ratio[idet]);
      FormatGraphPDRatios(PD_ratio_Dist[idet]);

      c35->cd(1);
      FormatPadAndPlotGalpropRatios(galpropRatio500, galpropRatio1000);
      PD_ratio[idet]->Draw("Psame");

      c35->cd(2);
      FormatPadAndPlotGalpropRatios(galpropRatio500, galpropRatio1000);
      PD_ratio_Dist[idet]->Draw("Psame");

   }


   /// Updating results file

   cout<<"*** Updating Results file ***"<<endl;
   nomefile="./Final_plots/"+mese+".root";
   TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
   f_out->mkdir("D Fluxes");
   f_out->cd("D Fluxes");
   c33-> Write();
   c32->Write();
   c34-> Write();
   c35->Write();
   f_out->cd("Export");
   D_Flux[0]->Write("Deutons Primary Flux: TOF");
   D_Flux[1]->Write("Deutons Primary Flux: NaF");
   D_Flux[2]->Write("Deutons Primary Flux: Agl");
   f_out->Write();
   f_out->Close();


   return;
}

void FormatGraphPDRatios(TGraphErrors* pdratio)
{
   pdratio->SetMarkerStyle(8);
   pdratio->SetMarkerSize(1.5);
   pdratio->SetMarkerColor(4);
   pdratio->SetLineColor(4);
   pdratio->SetLineWidth(2);
   return;
}


TGraph* GetGalpropRatioFromFile(string filename)
{
   return GetPowerWeightedGalpropRatioFromFile(filename, 0, 1);
}


TGraph* GetPowerWeightedGalpropRatioFromFile(string filename, float power, float yweight)
{
   TGraph* graph=new TGraph();
   float x, y=0;
   int j=0;
   cout<<filename<<endl;
   ifstream fp(filename.c_str());
   while (!fp.eof()) {
      fp>>x>>y;
      if(x/1e3>0.05&&x/1e3<=100)
         graph->SetPoint(j++,x/1e3, y*yweight*pow(x/1e3,power));
   }
   return graph;
}

void FormatGraphGalprop(TGraph* graph)
{
   FormatGraphGalpropXaxis(graph);
   graph->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
   graph->GetXaxis()->SetTitleSize(0.045);
   graph->GetYaxis()->SetTitleSize(0.045);
   graph->GetYaxis()->SetRangeUser(1e-2,1e4);
   graph->SetTitle("Deutons Flux: Geo. Zones");
   return;
}

void FormatGraphGalpropRatio(TGraph* graph)
{
   FormatGraphGalpropXaxis(graph);
   graph->GetYaxis()->SetRangeUser(1e-3,1e-1);
   graph->GetYaxis()->SetTitle("Flux ratio");
   return;
}

void FormatGraphGalpropXaxis(TGraph* graph)
{
   graph->GetXaxis()->SetRangeUser(0.1,10);
   graph->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
   return;
}


void FormatPadAndPlotGalpropRatios(TGraph* graph1, TGraph* graph2)
{


   gPad->SetLogx(); // This makes me sad :(
   gPad->SetLogy();
   gPad->SetGridx();
   gPad->SetGridy();

   graph1->Draw("AC");
   graph2->Draw("sameC");

   return;

}


void SetDeutonsExposureGraphStyle(TGraphErrors* graph)
{
   graph->SetMarkerStyle(8);
   graph->SetMarkerColor(4);
   graph->SetLineColor(4);
   return;
}

void SetProtonsExposureGraphStyle(TGraphErrors* graph)
{
   graph->SetMarkerStyle(8);
   graph->SetMarkerColor(2);
   graph->SetLineColor(2);
   return;
}


void FormatDeutonsExposureGraphTitle(TGraphErrors* graph)
{
   graph->GetXaxis()->SetTitle("kin. En./nucl.");
   graph->GetYaxis()->SetTitle("Exposure Time");
   graph->GetXaxis()->SetTitleSize(0.05);
   graph->GetYaxis()->SetTitleSize(0.05);
   return;
}


void FormatDFluxGeoWithDetZone(TGraphErrors* graph, int detector, int geozone)
{

   int style[NSUBDETECTORS]= {8, 4, 3}; ///< for SetMarkerStyle

   graph->SetMarkerStyle(style[detector]);
   graph->SetMarkerSize(1.5);
   graph->SetMarkerColor(geozone);
   graph->SetLineColor(geozone);
   graph->SetLineWidth(2);

   return;

}


TCanvas* MakeLogGridCanvas(string title)
{
   TCanvas* canv=new TCanvas(title.data());
   canv->cd();
   gPad->SetLogx(); // This makes me sad :(
   gPad->SetLogy();
   gPad->SetGridx();
   gPad->SetGridy();
   return canv;
}



TGraphErrors* FillGraphErrorFromBinningAndHisto(Binning bin, TH1* histo)
{
   TGraphErrors* graph=new TGraphErrors();
   for(int ibin=0; ibin<bin.EkBinsCent().size(); ibin++) {
      graph->SetPoint     (ibin, bin.EkBinCent(ibin),  histo->GetBinContent(ibin));
      graph->SetPointError(ibin, 0,                 histo->GetBinError(ibin)  );
   }
   return graph;
}



