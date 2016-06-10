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

   std::vector <TH1F*> vD_Flux (3);
   std::vector <TH1F*> vD_Dist (3);
   std::vector <TH1F*> vP_Flux (3);
   std::vector <TH1F*> vP_Dist (3);
   std::vector <TH2F*> vG_Flux (3);
   std::vector <TH2F*> vG_Dist (3);
   std::vector <TH1F*> vD_Exp  (3);
   std::vector <TH1F*> vP_Exp  (3);



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
      D_Flux_geo     -> Eval_Flux(11, deutons, 2 );
      D_Flux_geo_Dist-> Eval_Flux(11, deutons, 2 );

      if (ifx<3) { // No detector for Rigidity
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


   TH1F * DP_ratio[3];
   TH1F * DP_ratio_Dist[3];

   for (int ifx=0; ifx<3; ifx++) {
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

   TGraphErrors * D_Fluxgeo[3][11];
   TGraphErrors * D_FluxgeoDist[3][11];


   for (int idet=0; idet<3; idet++) {
      for(int izone=0; izone<11; izone++) {
         string nome="Deutons Flux: Geo. Zone "+to_string(izone);
         D_Fluxgeo[idet][izone]=new TGraphErrors();
         D_Fluxgeo[idet][izone]->SetName(nome.c_str());
         D_FluxgeoDist[idet][izone]=new TGraphErrors();
         D_FluxgeoDist[idet][izone]->SetName((nome + "Distance Fit").c_str());

         for(int m=1; m<binP[idet+1].size(); m++) {
            D_Fluxgeo    [idet][izone]->SetPoint     (m-1, binP[idet].EkBinCent(m), vG_Flux[idet]->GetBinContent(m+1, izone+1));
            D_FluxgeoDist[idet][izone]->SetPoint     (m-1, binP[idet].EkBinCent(m), vG_Dist[idet]->GetBinContent(m+1, izone+1));
            D_Fluxgeo    [idet][izone]->SetPointError(m-1, 0,                    vG_Flux[idet]->GetBinError  (m+1, izone+1));
            D_FluxgeoDist[idet][izone]->SetPointError(m-1, 0,                    vG_Dist[idet]->GetBinError  (m+1, izone+1));
         }

         FormatDFluxGeoWithDetZone(D_Fluxgeo    [idet][izone], idet, izone);
         FormatDFluxGeoWithDetZone(D_FluxgeoDist[idet][izone], idet, izone);

      }
      c32->cd(1);
      D_Fluxgeo[idet][10]->Draw("Psame");
      for(int j=0; j<11; j++) D_Fluxgeo[idet][j]->Draw("Psame");

      c32->cd(2);
      D_FluxgeoDist[idet][10]->Draw("Psame");
      for(int j=0; j<11; j++) D_FluxgeoDist[idet][j]->Draw("Psame");

   }

   /// Plot exposure time


   TCanvas * c33 = new TCanvas("Exposure Time");
   c33->Divide(1,3);

   string name[3]= {"ToF", "RICH NaF", "RICH Agl"};
   TGraphErrors * ged[3];
   TGraphErrors * gep[3];

   for (int i=0; i<3; i++) {
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
   c34->cd(1);
   FormatPadAndPlotGalpropRatios(galprop3P, galprop3P2);
   c34->cd(2);
   FormatPadAndPlotGalpropRatios(galprop3P, galprop3P2);


   // Loop on detectors


   TGraphErrors * D_Flux    [3];
   TGraphErrors * D_FluxDist[3];

   for (int idet=0; idet<3; idet++) {

      string nome="Deutons Flux: Primaries" ;
      D_Flux[idet]=new TGraphErrors();
      D_Flux[idet]->SetName(nome.c_str());
      D_FluxDist[idet]=new TGraphErrors();
      D_FluxDist[idet]->SetName((nome + "Distance Fit").c_str());

      for(int m=0; m< binD[idet].EkBinsCent().size(); m++) {
         D_Flux[idet]    ->SetPoint     (m, binD[idet].EkBinCent(m), vD_Flux[idet]->GetBinContent(m+2));
         D_FluxDist[idet]->SetPoint     (m, binD[idet].EkBinCent(m), vD_Dist[idet]->GetBinContent(m+2));
         D_Flux[idet]    ->SetPointError(m, 0,                    vD_Flux[idet]->GetBinError  (m+2));
         D_FluxDist[idet]->SetPointError(m, 0,                    vD_Dist[idet]->GetBinError  (m+2));
      }

      FormatDFluxGeoWithDetZone(D_Flux    [idet],      idet, 0);
      FormatDFluxGeoWithDetZone(D_FluxDist[idet],  idet, 0);



      c34->cd(1);
      D_Flux[idet]->Draw("Psame");

      c34->cd(2);
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

   TGraphErrors * PD_ratio     [3];
   TGraphErrors * PD_ratio_Dist[3];

   for (int i=0; i<3; i++) {

      PD_ratio[i]     =new TGraphErrors();
      PD_ratio_Dist[i]=new TGraphErrors();
      for(int m=0; m<binP[i].EkBinsCent().size(); m++) {
         PD_ratio[i]     ->SetPoint     (m, binP[i].EkBinCent(m), DP_ratio[i]->     GetBinContent(m+2));
         PD_ratio_Dist[i]->SetPoint     (m, binP[i].EkBinCent(m), DP_ratio_Dist[i]->GetBinContent(m+2));
         PD_ratio[i]     ->SetPointError(m, 0,                    DP_ratio[i]->GetBinError(m+2));
         PD_ratio_Dist[i]->SetPointError(m, 0,                    DP_ratio_Dist[i]->GetBinError(m+2));
      }
      FormatGraphPDRatios(PD_ratio[i]);
      FormatGraphPDRatios(PD_ratio_Dist[i]);


      c35->cd(1);
      FormatPadAndPlotGalpropRatios(galpropRatio500, galpropRatio1000);
      PD_ratio[i]->Draw("Psame");

      c35->cd(2);
      FormatPadAndPlotGalpropRatios(galpropRatio500, galpropRatio1000);
      PD_ratio_Dist[i]->Draw("Psame");

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

void FormatGraphGalprop(TGraph* graph) {
   FormatGraphGalpropXaxis(graph);
   graph->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
   graph->GetXaxis()->SetTitleSize(0.045);
   graph->GetYaxis()->SetTitleSize(0.045);
   graph->GetYaxis()->SetRangeUser(1e-2,1e4);
   graph->SetTitle("Deutons Flux: Geo. Zones");
   return;
}

void FormatGraphGalpropRatio(TGraph* graph) {
   FormatGraphGalpropXaxis(graph);
   graph->GetYaxis()->SetRangeUser(1e-3,1e-1);
   graph->GetYaxis()->SetTitle("Flux ratio");
   return;
}

void FormatGraphGalpropXaxis(TGraph* graph) {
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


void FormatDFluxGeoWithDetZone(TGraphErrors* graph, int detector, int geozone) {

    int style[3]= {8, 4, 3}; ///< for SetMarkerStyle

    graph->SetMarkerStyle(style[detector]);
    graph->SetMarkerSize(1.5);
    graph->SetMarkerColor(geozone);
    graph->SetLineColor(geozone);
    graph->SetLineWidth(2);

         return;

}


