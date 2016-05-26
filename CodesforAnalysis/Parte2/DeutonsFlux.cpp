using namespace std;

void CheckSyst(TH1F * Counts1,TH1F * Counts2, TH1F * Syst);

void DeutonFlux()
{
   string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
   TFile * file1 = TFile::Open(nomefile.c_str(),"READ");








   cout<<"*************** DEUTONS FLUXES CALCULATION *******************"<<endl;

   //Evaluation additional Fit Error


   std::string suf[4] = {"R", "ToF", "NaF", "Agl"};
   Binning binP[4] = {RB, TofPB, NaFPB, AglPB};
   Binning binD[4] = {RB, TofDB, NaFDB, AglDB};
   string expodnames[4] = ["esposizionegeo","esposizionedgeo", "esposizionedgeoNaF", "esposizionedgeoAgl" ];
   string expopnames[4] = ["esposizionegeo","esposizionepgeo", "esposizionepgeoNaF", "esposizionepgeoAgl" ];
   enum {protons, deutons};

   std::vector <TH1F*> vD_Flux (4);
   std::vector <TH1F*> vD_Dist (4);
   std::vector <TH1F*> vP_Flux (4);
   std::vector <TH1F*> vP_Dist (4);
   std::vector <TH2F*> vG_Flux (4);
   std::vector <TH2F*> vG_Dist (4);




   for (int ifx=0; ifx<4; ifx++) { // index flux

      //Deutons
      Flux * D_Flux       = new Flux(file1, "D_Flux"       ,"Results", suf[ifx], "Corr_AcceptanceD",binD[ifx]);
      Flux * D_Flux_Dist  = new Flux(file1, "D_Flux_Dist"  ,"Results", suf[ifx], "Corr_AcceptanceD",binD[ifx]);
      if(D_Flux -> Counts) {
         TH1F * Syst=(TH1F*) D_Flux -> Counts ->Clone();
         Syst -> Add ((TH1F*) D_Flux_Dist -> Counts ->Clone(), -1);
         CheckSyst( (TH1F*)D_Flux_Dist -> Counts , (TH1F*)D_Flux -> Counts,Syst);
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
      P_Flux         -> Set_Exposure_Time (esposizionegeo);
      P_Flux_Dist    -> Set_Exposure_Time (esposizionegeo);
      P_Flux         -> Eval_Flux(1 , protons );
      P_Flux_Dist    -> Eval_Flux(1 , protons );

      Flux * D_Flux_geo     = new Flux(file1, "D_Flux_geo"           ,"Results", suf[ifx],"Geomag_AcceptanceD",binD[ifx]);
      Flux * D_Flux_geo_Dist     = new Flux(file1, "D_Flux_geo_Dist" ,"Results", suf[ifx],"Geomag_AcceptanceD",binD[ifx]);
      // Flux geo

      D_Flux_geo     -> Set_Exposure_Time (Tempi);
      D_Flux_geo_Dist-> Set_Exposure_Time (Tempi);
      D_Flux_geo     -> Eval_Flux(11, deutons, 2 );
      D_Flux_geo_Dist-> Eval_Flux(11, deutons, 2 );


      if (ifx>=1) {
         vD_Flux .push_back( (TH1F*)  D_Flux        ->Fluxes);
         vD_Dist.push_back( (TH1F*)  D_Flux_Dist    ->Fluxes);
         vP_Flux .push_back( (TH1F*)  P_Flux        ->Fluxes);
         vP_Dist.push_back( (TH1F*)  P_Flux_Dist    ->Fluxes);
         vG_Flux .push_back( (TH2F*) P_Flux_geo     ->Fluxes);
         vG_Dist.push_back( (TH2F*)  P_Flux_geo_Dist->Fluxes);

      }

   }


   //// Writing the files

   cout<<"*** Updating P1 file ****"<<endl;

   nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
   file1 = TFile::Open(nomefile.c_str(),"UPDATE");
   file1->cd("Results/Fluxes");




   for (int ifx=0; ifx<3; ifx++) {
      //Fit on Mass
      vD_Flux[i] 	   ->Write(("DeutonsPrimaryFlux_"+suf[i+1]).data());
      vG_Flux[i]  	->Write(("DeutonsGeomagFlux_"+suf[i+1]).data());
      //Fit on Distance
      vD_Dist[i]	->Write(("DeutonsPrimaryFlux_Dist_"+suf[i+1]).data());
      vG_Dist[i] 	->Write(("DeutonsGeomagFlux_Dist_"+suf[i+1]).data());
      //for D/P ratio: P Flux
      vP_Flux[i] 	  	->Write(("ProtonsPrimaryFlux_"+suf[i+1]).data());
      vP_Dist[i]	->Write(("Protons_PrimaryFlux_Dist_"+suf[i+1]).data());

      //D/P ratio
      //Fit on Mass
      TH1F * DP_ratio = vD_Flux[i] -> Clone();
      DP_ratio ->  Divide(vP_Flux[i]			);
      DP_ratio 	 -> Write(("DP_ratio_"+suf[i+1]		).data());
      //Fit on Distance
      TH1F * DP_ratio_Dist = vD_Flux[i] -> Clone();
      DP_ratio_Dist ->  Divide(vP_Dist[i] );
      DP_ratioTOF_Dist -> Write(("DP_ratio_Dist_"+suf[i+1]	).data());
   }


   file1->Write();
   file1->Close();



   //// Drawing the stuff

   TCanvas * c32 = new TCanvas("Deutons Flux: Geo. Zones");
   TCanvas * c33 = new TCanvas("Exposure Time");
   c33->Divide(1,3);
   TCanvas * c34 = new TCanvas("Deutons Flux: Primaries");
   TCanvas * c35 = new TCanvas("D/P ratio");



   float potenza=0;



   c32-> Divide(2,1);
   c32->cd(1);

   c32->cd(2);
   gPad->SetLogx();
   gPad->SetLogy();
   gPad->SetGridx();
   gPad->SetGridy();
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
   gPad->SetTitle("Deutons Flux: Geo. Zone (Distance FIT)");
   galprop3P->Draw("AC");
   galprop3P2->Draw("sameC");


   TGraphErrors * D_FluxgeoTOF[11];
   TGraphErrors * D_FluxgeoNaF[11];
   TGraphErrors * D_FluxgeoAgl[11];
   TGraphErrors * D_FluxgeoDistTOF[11];
   TGraphErrors * D_FluxgeoDistNaF[11];
   TGraphErrors * D_FluxgeoDistAgl[11];


   string nome;
   int p=0;
   for(int j=0; j<11; j++) {
      nome="Deutons Flux: Geo. Zone "+to_string(j);
      D_FluxgeoTOF[j]=new TGraphErrors();
      D_FluxgeoTOF[j]->SetName(nome.c_str());
      D_FluxgeoDistTOF[j]=new TGraphErrors();
      D_FluxgeoDistTOF[j]->SetName((nome + "Distance Fit").c_str());

      p=0;
      for(int m=1; m<nbinsToF; m++) {
         D_FluxgeoTOF[j]->SetPoint(p,ToFPB.EkBinCent(m),vG_Flux[0]->GetBinContent(m+1,j+1));
         D_FluxgeoTOF[j]->SetPointError(p,0,vG_Flux[0]->GetBinError(m+1,j+1));
         D_FluxgeoDistTOF[j]->SetPoint(p,ToFPB.EkBinCent(m),vG_Dist[i]->GetBinContent(m+1,j+1));
         D_FluxgeoDistTOF[j]->SetPointError(p,0,vG_Dist[i]->GetBinError(m+1,j+1));
         p++;
      }
      D_FluxgeoTOF[j]->SetMarkerStyle(8);
      D_FluxgeoTOF[j]->SetMarkerSize(1.5);
      D_FluxgeoTOF[j]->SetMarkerColor(j-1);
      D_FluxgeoTOF[j]->SetLineColor(j-1);
      D_FluxgeoTOF[j]->SetLineWidth(2);

      D_FluxgeoDistTOF[j]->SetMarkerStyle(8);
      D_FluxgeoDistTOF[j]->SetMarkerSize(1.5);
      D_FluxgeoDistTOF[j]->SetMarkerColor(j-1);
      D_FluxgeoDistTOF[j]->SetLineColor(j-1);
      D_FluxgeoDistTOF[j]->SetLineWidth(2);

   }
   c32->cd(1);
   D_FluxgeoTOF[10]->Draw("Psame");
   for(int j=0; j<11; j++) D_FluxgeoTOF[j]->Draw("Psame");

   c32->cd(2);
   D_FluxgeoDistTOF[10]->Draw("Psame");
   for(int j=0; j<11; j++) D_FluxgeoDistTOF[j]->Draw("Psame");

   p=0;
   for(int j=0; j<11; j++) {
      nome="Deutons Flux: Geo. Zone "+to_string(j);
      D_FluxgeoNaF[j]=new TGraphErrors();
      D_FluxgeoNaF[j]->SetName(nome.c_str());
      D_FluxgeoDistNaF[j]=new TGraphErrors();
      D_FluxgeoDistNaF[j]->SetName((nome + "Distance Fit").c_str());

      p=0;
      for(int m=1; m<nbinsToF; m++) {
         D_FluxgeoNaF[j]->SetPoint(p,NaFPB.EkBinCent(m),vG_Flux[1]->GetBinContent(m+1,j+1));
         D_FluxgeoNaF[j]->SetPointError(p,0,vG_Flux[1]->GetBinError(m+1,j+1));
         D_FluxgeoDistNaF[j]->SetPoint(p,NaFPB.EkBinCent(m),vG_Dist[i]->GetBinContent(m+1,j+1));
         D_FluxgeoDistNaF[j]->SetPointError(p,0,vG_Dist[i]->GetBinError(m+1,j+1));
         p++;
      }
      D_FluxgeoNaF[j]->SetMarkerStyle(4);
      D_FluxgeoNaF[j]->SetMarkerSize(1.5);
      D_FluxgeoNaF[j]->SetMarkerColor(j-1);
      D_FluxgeoNaF[j]->SetLineColor(j-1);
      D_FluxgeoNaF[j]->SetLineWidth(2);

      D_FluxgeoDistNaF[j]->SetMarkerStyle(4);
      D_FluxgeoDistNaF[j]->SetMarkerSize(1.5);
      D_FluxgeoDistNaF[j]->SetMarkerColor(j-1);
      D_FluxgeoDistNaF[j]->SetLineColor(j-1);
      D_FluxgeoDistNaF[j]->SetLineWidth(2);

   }
   c32->cd(1);
   D_FluxgeoNaF[10]->Draw("Psame");
   for(int j=0; j<11; j++) D_FluxgeoNaF[j]->Draw("Psame");

   c32->cd(2);
   D_FluxgeoDistNaF[10]->Draw("Psame");
   for(int j=0; j<11; j++) D_FluxgeoDistNaF[j]->Draw("Psame");

   p=0;
   for(int j=0; j<11; j++) {
      nome="Deutons Flux: Geo. Zone "+to_string(j);
      D_FluxgeoAgl[j]=new TGraphErrors();
      D_FluxgeoAgl[j]->SetName(nome.c_str());
      D_FluxgeoDistAgl[j]=new TGraphErrors();
      D_FluxgeoDistAgl[j]->SetName((nome + "Distance Fit").c_str());

      p=0;
      for(int m=1; m<nbinsToF; m++) {
         D_FluxgeoAgl[j]->SetPoint(p,AglPB.EkBinCent(m),vG_Flux[2]->GetBinContent(m+1,j+1));
         D_FluxgeoAgl[j]->SetPointError(p,0,vG_Flux[2]->GetBinError(m+1,j+1));
         D_FluxgeoDistAgl[j]->SetPoint(p,AglPB.EkBinCent(m),vG_Dist[i]->GetBinContent(m+1,j+1));
         D_FluxgeoDistAgl[j]->SetPointError(p,0,vG_Dist[i]->GetBinError(m+1,j+1));
         p++;
      }
      D_FluxgeoAgl[j]->SetMarkerStyle(3);
      D_FluxgeoAgl[j]->SetMarkerSize(1.5);
      D_FluxgeoAgl[j]->SetMarkerColor(j-1);
      D_FluxgeoAgl[j]->SetLineColor(j-1);
      D_FluxgeoAgl[j]->SetLineWidth(2);

      D_FluxgeoDistAgl[j]->SetMarkerStyle(3);
      D_FluxgeoDistAgl[j]->SetMarkerSize(1.5);
      D_FluxgeoDistAgl[j]->SetMarkerColor(j-1);
      D_FluxgeoDistAgl[j]->SetLineColor(j-1);
      D_FluxgeoDistAgl[j]->SetLineWidth(2);

   }
   c32->cd(1);
   D_FluxgeoAgl[10]->Draw("Psame");
   for(int j=0; j<11; j++) D_FluxgeoAgl[j]->Draw("Psame");

   c32->cd(2);
   D_FluxgeoDistAgl[10]->Draw("Psame");
   for(int j=0; j<11; j++) D_FluxgeoDistAgl[j]->Draw("Psame");

   TGraphErrors * D_FluxTOF;
   TGraphErrors * D_FluxNaF;
   TGraphErrors * D_FluxAgl;
   TGraphErrors * D_FluxDistTOF;
   TGraphErrors * D_FluxDistNaF;
   TGraphErrors * D_FluxDistAgl;

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

   nome="Deutons Flux: Primaries" ;
   D_FluxTOF=new TGraphErrors();
   D_FluxTOF->SetName(nome.c_str());
   D_FluxDistTOF=new TGraphErrors();
   D_FluxDistTOF->SetName((nome + "Distance Fit").c_str());

   p=0;
   for(int m=1; m<nbinsToF; m++) {
      D_FluxTOF->SetPoint(p,ToFPB.EkBinCent(m),vD_Flux[0]->GetBinContent(m+1));
      D_FluxTOF->SetPointError(p,0,vD_Flux[0]->GetBinError(m+1));
      D_FluxDistTOF->SetPoint(p,ToFPB.EkBinCent(m),vD_Dist[i]->GetBinContent(m+1));
      D_FluxDistTOF->SetPointError(p,0,vD_Dist[i]->GetBinError(m+1));
      p++;
   }
   D_FluxTOF->SetMarkerStyle(8);
   D_FluxTOF->SetMarkerSize(1.5);
   D_FluxTOF->SetMarkerColor(j-1);
   D_FluxTOF->SetLineColor(j-1);
   D_FluxTOF->SetLineWidth(2);

   D_FluxDistTOF->SetMarkerStyle(8);
   D_FluxDistTOF->SetMarkerSize(1.5);
   D_FluxDistTOF->SetMarkerColor(j-1);
   D_FluxDistTOF->SetLineColor(j-1);
   D_FluxDistTOF->SetLineWidth(2);

   c34->cd(1);
   D_FluxTOF->Draw("Psame");

   c34->cd(2);
   D_FluxDistTOF->Draw("Psame");

   nome="Deutons Flux: Primaries" ;
   D_FluxNaF=new TGraphErrors();
   D_FluxNaF->SetName(nome.c_str());
   D_FluxDistNaF=new TGraphErrors();
   D_FluxDistNaF->SetName((nome + "Distance Fit").c_str());

   p=0;
   for(int m=1; m<nbinsToF; m++) {
      D_FluxNaF->SetPoint(p,NaFPB.EkBinCent(m),vD_Flux[1]->GetBinContent(m+1));
      D_FluxNaF->SetPointError(p,0,vD_Flux[1]->GetBinError(m+1));
      D_FluxDistNaF->SetPoint(p,NaFPB.EkBinCent(m),vD_Dist[i]->GetBinContent(m+1));
      D_FluxDistNaF->SetPointError(p,0,vD_Dist[i]->GetBinError(m+1));
      p++;
   }
   D_FluxNaF->SetMarkerStyle(4);
   D_FluxNaF->SetMarkerSize(1.5);
   D_FluxNaF->SetMarkerColor(j-1);
   D_FluxNaF->SetLineColor(j-1);
   D_FluxNaF->SetLineWidth(2);

   D_FluxDistNaF->SetMarkerStyle(4);
   D_FluxDistNaF->SetMarkerSize(1.5);
   D_FluxDistNaF->SetMarkerColor(j-1);
   D_FluxDistNaF->SetLineColor(j-1);
   D_FluxDistNaF->SetLineWidth(2);

   c34->cd(1);
   D_FluxNaF->Draw("Psame");

   c34->cd(2);
   D_FluxDistNaF->Draw("Psame");

   nome="Deutons Flux: Primaries" ;
   D_FluxAgl=new TGraphErrors();
   D_FluxAgl->SetName(nome.c_str());
   D_FluxDistAgl=new TGraphErrors();
   D_FluxDistAgl->SetName((nome + "Distance Fit").c_str());

   p=0;
   for(int m=1; m<nbinsToF; m++) {
      D_FluxAgl->SetPoint(p,AglPB.EkBinCent(m),vD_Flux[2]->GetBinContent(m+1));
      D_FluxAgl->SetPointError(p,0,vD_Flux[2]->GetBinError(m+1));
      D_FluxDistAgl->SetPoint(p,AglPB.EkBinCent(m),vD_Dist[i]->GetBinContent(m+1));
      D_FluxDistAgl->SetPointError(p,0,vD_Dist[i]->GetBinError(m+1));
      p++;
   }
   D_FluxAgl->SetMarkerStyle(3);
   D_FluxAgl->SetMarkerSize(1.5);
   D_FluxAgl->SetMarkerColor(j-1);
   D_FluxAgl->SetLineColor(j-1);
   D_FluxAgl->SetLineWidth(2);

   D_FluxDistAgl->SetMarkerStyle(3);
   D_FluxDistAgl->SetMarkerSize(1.5);
   D_FluxDistAgl->SetMarkerColor(j-1);
   D_FluxDistAgl->SetLineColor(j-1);
   D_FluxDistAgl->SetLineWidth(2);

   c34->cd(1);
   D_FluxAgl->Draw("Psame");

   c34->cd(2);
   D_FluxDistAgl->Draw("Psame");


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

   TGraphErrors * PD_ratioTOF=new TGraphErrors();
   TGraphErrors * PD_ratioTOF_Dist=new TGraphErrors();
   p=0;
   for(int m=1; m<nbinsToF; m++) {
      PD_ratioTOF->SetPoint(p,ToFPB.EkBinCent(m),DP_ratioTOF->GetBinContent(m+1));
      PD_ratioTOF_Dist->SetPoint(p,ToFPB.EkBinCent(m),DP_ratioTOF_Dist->GetBinContent(m+1));
      PD_ratioTOF->SetPointError(p,0,DP_ratioTOF->GetBinError(m+1));
      PD_ratioTOF_Dist->SetPointError(p,0,DP_ratioTOF_Dist->GetBinError(m+1));
      p++;
   }
   PD_ratioTOF->SetMarkerStyle(8);
   PD_ratioTOF->SetMarkerSize(1.5);
   PD_ratioTOF->SetMarkerColor(4);
   PD_ratioTOF->SetLineColor(4);
   PD_ratioTOF->SetLineWidth(2);

   PD_ratioTOF_Dist->SetMarkerStyle(8);
   PD_ratioTOF_Dist->SetMarkerSize(1.5);
   PD_ratioTOF_Dist->SetMarkerColor(4);
   PD_ratioTOF_Dist->SetLineColor(4);
   PD_ratioTOF_Dist->SetLineWidth(2);

   c35->cd(1);
   PD_ratioTOF->Draw("Psame");

   c35->cd(2);
   PD_ratioTOF_Dist->Draw("Psame");


   TGraphErrors * PD_ratioNaF=new TGraphErrors();
   TGraphErrors * PD_ratioNaF_Dist=new TGraphErrors();
   p=0;
   for(int m=1; m<nbinsToF; m++) {
      PD_ratioNaF->SetPoint(p,NaFPB.EkBinCent(m),DP_ratioNaF->GetBinContent(m+1));
      PD_ratioNaF_Dist->SetPoint(p,NaFPB.EkBinCent(m),DP_ratioNaF_Dist->GetBinContent(m+1));
      PD_ratioNaF->SetPointError(p,0,DP_ratioNaF->GetBinError(m+1));
      PD_ratioNaF_Dist->SetPointError(p,0,DP_ratioNaF_Dist->GetBinError(m+1));
      p++;
   }
   PD_ratioNaF->SetMarkerStyle(4);
   PD_ratioNaF->SetMarkerSize(1.5);
   PD_ratioNaF->SetMarkerColor(4);
   PD_ratioNaF->SetLineColor(4);
   PD_ratioNaF->SetLineWidth(2);

   PD_ratioNaF_Dist->SetMarkerStyle(4);
   PD_ratioNaF_Dist->SetMarkerSize(1.5);
   PD_ratioNaF_Dist->SetMarkerColor(4);
   PD_ratioNaF_Dist->SetLineColor(4);
   PD_ratioNaF_Dist->SetLineWidth(2);

   c35->cd(1);
   PD_ratioNaF->Draw("Psame");

   c35->cd(2);
   PD_ratioNaF_Dist->Draw("Psame");

   TGraphErrors * PD_ratioAgl=new TGraphErrors();
   TGraphErrors * PD_ratioAgl_Dist=new TGraphErrors();
   p=0;
   for(int m=1; m<nbinsToF; m++) {
      PD_ratioAgl->SetPoint(p,AglPB.EkBinCent(m),DP_ratioAgl->GetBinContent(m+1));
      PD_ratioAgl_Dist->SetPoint(p,AglPB.EkBinCent(m),DP_ratioAgl_Dist->GetBinContent(m+1));
      PD_ratioAgl->SetPointError(p,0,DP_ratioAgl->GetBinError(m+1));
      PD_ratioAgl_Dist->SetPointError(p,0,DP_ratioAgl_Dist->GetBinError(m+1));
      p++;
   }
   PD_ratioAgl->SetMarkerStyle(3);
   PD_ratioAgl->SetMarkerSize(1.5);
   PD_ratioAgl->SetMarkerColor(4);
   PD_ratioAgl->SetLineColor(4);
   PD_ratioAgl->SetLineWidth(2);

   PD_ratioAgl_Dist->SetMarkerStyle(3);
   PD_ratioAgl_Dist->SetMarkerSize(1.5);
   PD_ratioAgl_Dist->SetMarkerColor(4);
   PD_ratioAgl_Dist->SetLineColor(4);
   PD_ratioAgl_Dist->SetLineWidth(2);

   c35->cd(1);
   PD_ratioAgl->Draw("Psame");

   c35->cd(2);
   PD_ratioAgl_Dist->Draw("Psame");

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
   D_FluxTOF->Write("Deutons Primary Flux: TOF");
   D_FluxNaF->Write("Deutons Primary Flux: NaF");
   D_FluxAgl->Write("Deutons Primary Flux: Agl");
   f_out->Write();
   f_out->Close();



   return;
}


void CheckSyst(TH1F * Counts1,TH1F * Counts2, TH1F * Syst)
{
   for(int i=0; i<Syst->GetNbinsX(); i++) {
      if(Counts1 -> GetBinContent(i)==0 || Counts2 -> GetBinContent(i)==0 )
         Syst -> SetBinContent (i,0);
   }
   return;
}




void PlotExpo (TH1F* hdexpo, TH1F* hpexpo, int idx)
{
   string name[3]= {"ToF", "RICH NaF", "RICH Agl"}
                   TGraphErrors * ged=new TGraphErrors();
   TGraphErrors * gep=new TGraphErrors();
   for(int m=0; m<nbinsbeta; m++) {
      ged->SetPoint(m,ToFDB.EkBinCent(m), hdexpo-> GetBinContent(m+1));
      gep->SetPoint(m,ToFPB.EkBinCent(m), hpexpo-> GetBinContent(m+1));
   }
   c33->cd(idx);
   gPad->SetLogy();
   gPad->SetGridx();
   gPad->SetGridy();

   gep->SetMarkerStyle(8);
   gep->SetMarkerColor(2);
   gep->SetLineColor(2);
   ged->SetMarkerStyle(8);
   ged->SetMarkerColor(4);
   ged->SetLineColor(4);
   ged->GetXaxis()->SetTitle("kin. En./nucl.");
   ged->GetYaxis()->SetTitle("Exposure Time");
   ged->GetXaxis()->SetTitleSize(0.05);
   ged->GetYaxis()->SetTitleSize(0.05);
   ged->SetTitle((name[idx]+" range").data());
   ged->Draw("APC");
   gep->Draw("PCsame");
}
