using namespace std;


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
      vD_Flux[ifx] 	   ->Write(("DeutonsPrimaryFlux_"+suf[ifx]).data());
      vG_Flux[ifx]  	->Write(("DeutonsGeomagFlux_"+suf[ifx]).data());
      //Fit on Distance
      vD_Dist[ifx]	->Write(("DeutonsPrimaryFlux_Dist_"+suf[ifx]).data());
      vG_Dist[ifx] 	->Write(("DeutonsGeomagFlux_Dist_"+suf[ifx]).data());
      //for D/P ratio: P Flux
      vP_Flux[ifx] 	  	->Write(("ProtonsPrimaryFlux_"+suf[ifx]).data());
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



   TCanvas * c34 = new TCanvas("Deutons Flux: Primaries");
   TCanvas * c35 = new TCanvas("D/P ratio");



   float potenza=0; ///< Stays 0, to be checked !


   TGraph* galprop3P=new TGraph();
   TGraph* galprop3P2=new TGraph();
   float x,y=0;
   int j=0;
   {
      string nomefile="./Galprop/Trotta2011/Def/new_D500.txt";
      cout<<nomefile<<endl;
      ifstream fp(nomefile.c_str());
      while (!fp.eof()) {
         fp>>x>>y;
         if(x/1e3>0.05&&x/1e3<=100)
            galprop3P->SetPoint(j,x/1e3,y*1e7*pow(x/1e3,potenza));
         j++;
      }
   }

   j=0;
   {
      string nomefile="./Galprop/Trotta2011/Def/new_D1250.txt";
      cout<<nomefile<<endl;
      ifstream fp(nomefile.c_str());
      while (!fp.eof()) {
         fp>>x>>y;
         if(x/1e3>0.05&&x/1e3<=100)
            galprop3P2->SetPoint(j,x/1e3,y*1e7*pow(x/1e3,potenza));
         j++;
      }
   }
   galprop3P->GetXaxis()->SetRangeUser(0.1,10);
   galprop3P->GetYaxis()->SetRangeUser(1e-3,1e3);

   TCanvas * c32 = new TCanvas("Deutons Flux: Geo. Zones");
   c32-> Divide(2,1);

   c32->cd(1);
   gPad->SetLogx();
   gPad->SetLogy();
   gPad->SetGridx();
   gPad->SetGridy();

   galprop3P->SetTitle("Deutons Flux: Geo. Zones");
   galprop3P->GetXaxis()->SetTitle("Kin.En./nucl. [GeV/nucl.]");
   galprop3P ->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
   galprop3P ->GetXaxis()->SetTitleSize(0.045);
   galprop3P->GetYaxis()->SetTitleSize(0.045);
   galprop3P ->GetYaxis()->SetRangeUser(1e-2,1e4);

   galprop3P->Draw("AC");
   galprop3P2->Draw("sameC");

   c32->cd(2);
   gPad->SetLogx();
   gPad->SetLogy();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetTitle("Deutons Flux: Geo. Zone (Distance FIT)");

   galprop3P->Draw("AC");
   galprop3P2->Draw("sameC");



   int style[3]= {8, 4, 3}; ///< for SetMarkerStyle


   /// Flux Geo

   TGraphErrors * D_Fluxgeo[3][11];
   TGraphErrors * D_FluxgeoDist[3][11];


   for (int i=0; i<3; i++) {
      for(int j=0; j<11; j++) {
         string nome="Deutons Flux: Geo. Zone "+to_string(j);
         D_Fluxgeo[i][j]=new TGraphErrors();
         D_Fluxgeo[i][j]->SetName(nome.c_str());
         D_FluxgeoDist[i][j]=new TGraphErrors();
         D_FluxgeoDist[i][j]->SetName((nome + "Distance Fit").c_str());

         for(int m=1; m<binP[i+1].size(); m++) {
            D_Fluxgeo    [i][j]->SetPoint     (m-1, binP[i].EkBinCent(m), vG_Flux[i]->GetBinContent(m+1, j+1));
            D_FluxgeoDist[i][j]->SetPoint     (m-1, binP[i].EkBinCent(m), vG_Dist[i]->GetBinContent(m+1, j+1));
            D_Fluxgeo    [i][j]->SetPointError(m-1, 0,                    vG_Flux[i]->GetBinError  (m+1, j+1));
            D_FluxgeoDist[i][j]->SetPointError(m-1, 0,                    vG_Dist[i]->GetBinError  (m+1, j+1));
         }

         D_Fluxgeo[i][j]->SetMarkerStyle(style[i]);
         D_Fluxgeo[i][j]->SetMarkerSize(1.5);
         D_Fluxgeo[i][j]->SetMarkerColor(j-1);
         D_Fluxgeo[i][j]->SetLineColor(j-1);
         D_Fluxgeo[i][j]->SetLineWidth(2);

         D_FluxgeoDist[i][j]->SetMarkerStyle(style[i]);
         D_FluxgeoDist[i][j]->SetMarkerSize(1.5);
         D_FluxgeoDist[i][j]->SetMarkerColor(j-1);
         D_FluxgeoDist[i][j]->SetLineColor(j-1);
         D_FluxgeoDist[i][j]->SetLineWidth(2);

      }
      c32->cd(1);
      D_Fluxgeo[i][10]->Draw("Psame");
      for(int j=0; j<11; j++) D_Fluxgeo[i][j]->Draw("Psame");

      c32->cd(2);
      D_FluxgeoDist[i][10]->Draw("Psame");
      for(int j=0; j<11; j++) D_FluxgeoDist[i][j]->Draw("Psame");

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

      gep[i]->SetMarkerStyle(8);
      gep[i]->SetMarkerColor(2);
      gep[i]->SetLineColor(2);
      ged[i]->SetMarkerStyle(8);
      ged[i]->SetMarkerColor(4);
      ged[i]->SetLineColor(4);
      ged[i]->GetXaxis()->SetTitle("kin. En./nucl.");
      ged[i]->GetYaxis()->SetTitle("Exposure Time");
      ged[i]->GetXaxis()->SetTitleSize(0.05);
      ged[i]->GetYaxis()->SetTitleSize(0.05);
      ged[i]->SetTitle((name[i]+" range").data());
      ged[i]->Draw("APC");
      gep[i]->Draw("PCsame");

   }

   /// Flux Detectors

   // Plot init

   c34-> Divide(2,1);

   c34->cd(1);
   gPad->SetLogx();
   gPad->SetLogy();
   gPad->SetGridx();
   gPad->SetGridy();
   galprop3P->Draw("AC");
   galprop3P2->Draw("sameC");
   c34->cd(2);
   gPad->SetLogx();
   gPad->SetLogy();
   gPad->SetGridx();
   gPad->SetGridy();
   galprop3P->Draw("AC");
   galprop3P2->Draw("sameC");

   // Loop on detectors


   TGraphErrors * D_Flux    [3];
   TGraphErrors * D_FluxDist[3];

   for (int i=0; i<3; i++) {

      string nome="Deutons Flux: Primaries" ;
      D_Flux[i]=new TGraphErrors();
      D_Flux[i]->SetName(nome.c_str());
      D_FluxDist[i]=new TGraphErrors();
      D_FluxDist[i]->SetName((nome + "Distance Fit").c_str());

      for(int m=0; m< binD[i].EkBinsCent().size(); m++) {
         D_Flux[i]    ->SetPoint     (m, binD[i].EkBinCent(m), vD_Flux[i]->GetBinContent(m+2));
         D_FluxDist[i]->SetPoint     (m, binD[i].EkBinCent(m), vD_Dist[i]->GetBinContent(m+2));
         D_Flux[i]    ->SetPointError(m, 0,                    vD_Flux[i]->GetBinError  (m+2));
         D_FluxDist[i]->SetPointError(m, 0,                    vD_Dist[i]->GetBinError  (m+2));
      }

      D_Flux[i]->SetMarkerStyle(style[i]);
      D_Flux[i]->SetMarkerSize(1.5);
      D_Flux[i]->SetMarkerColor(j-1);
      D_Flux[i]->SetLineColor(j-1);
      D_Flux[i]->SetLineWidth(2);

      D_FluxDist[i]->SetMarkerStyle(style[i]);
      D_FluxDist[i]->SetMarkerSize(1.5);
      D_FluxDist[i]->SetMarkerColor(j-1);
      D_FluxDist[i]->SetLineColor(j-1);
      D_FluxDist[i]->SetLineWidth(2);

      c34->cd(1);
      D_Flux[i]->Draw("Psame");

      c34->cd(2);
      D_FluxDist[i]->Draw("Psame");

   }



   // Some galprop init stuff of Canvas 35

   c35->Divide(2,1);
   c35->cd(1);
   gPad->SetLogx();
   gPad->SetLogy();
   gPad->SetGridx();
   gPad->SetGridy();
   c35->cd(2);
   gPad->SetLogx();
   gPad->SetLogy();
   gPad->SetGridx();
   gPad->SetGridy();

   TGraph* galpropratio1=new TGraph();
   TGraph* galpropratio2=new TGraph();
   j=0;
   {
      string nomefile="./Galprop/Trotta2011/PDratio/500.dat";
      cout<<nomefile<<endl;
      ifstream fp(nomefile.c_str());
      while (!fp.eof()) {
         fp>>x>>y;
         if(x/1e3>0.05&&x/1e3<=100)
            galpropratio1->SetPoint(j,x/1e3,y);
         j++;
      }
   }

   j=0;
   {
      string nomefile="./Galprop/Trotta2011/PDratio/1000.dat";
      cout<<nomefile<<endl;
      ifstream fp(nomefile.c_str());
      while (!fp.eof()) {
         fp>>x>>y;
         if(x/1e3>0.05&&x/1e3<=100)
            galpropratio2->SetPoint(j,x/1e3,y);
         j++;
      }
   }
   galpropratio1->GetXaxis()->SetRangeUser(0.1,10);
   galpropratio1->GetYaxis()->SetRangeUser(1e-3,1e-1);
   galpropratio1->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
   galpropratio1->GetYaxis()->SetTitle("Flux ratio");
   c35->cd(1);
   galpropratio1->Draw("AC");
   galpropratio2->Draw("sameC");
   c35->cd(2);
   galpropratio1->Draw("AC");
   galpropratio2->Draw("sameC");


   /// PD_Ratio

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
      PD_ratio[i]->SetMarkerStyle(8);
      PD_ratio[i]->SetMarkerSize(1.5);
      PD_ratio[i]->SetMarkerColor(4);
      PD_ratio[i]->SetLineColor(4);
      PD_ratio[i]->SetLineWidth(2);

      PD_ratio_Dist[i]->SetMarkerStyle(8);
      PD_ratio_Dist[i]->SetMarkerSize(1.5);
      PD_ratio_Dist[i]->SetMarkerColor(4);
      PD_ratio_Dist[i]->SetLineColor(4);
      PD_ratio_Dist[i]->SetLineWidth(2);

      c35->cd(1);
      PD_ratio[i]->Draw("Psame");

      c35->cd(2);
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
