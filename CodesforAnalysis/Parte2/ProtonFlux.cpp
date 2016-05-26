using namespace std;

Flux * P_Flux         = new Flux("P_Flux"    	  );
Flux * P_Flux_geo     = new Flux("P_Flux_geo"     ,11);
Flux * P_Flux_geo_prim= new Flux("P_Flux_geo_prim",11);

//DvsMC check
Flux * P_Flux_pre = new Flux("P_Flux_pre" );
Flux * P_Flux_sel = new Flux("P_Flux_sel" );

void ProtonFlux_Fill(TNtuple *ntupla, int l,int zona) {
    ntupla->GetEvent(l);
    if(Tup.Beta<=0||Tup.R<=0) return;
    int Kbin=RB.GetRBin(Tup.R);

    if(Tup.Dist5D_P<6 && Likcut) {
        P_Flux_geo-> Counts_R -> Fill(Kbin,zona);
        if(Tup.R>1.2*Tup.Rcutoff) {
            P_Flux -> Counts_R-> Fill(Kbin);
            P_Flux_geo_prim -> Counts_R -> Fill(Kbin,zona);
        }
    }

    if(Herejcut && Tup.R>1.2*Tup.Rcutoff) {
        P_Flux_pre -> Counts_R -> Fill(Kbin);
        if(Tup.Dist5D_P<6&&Likcut)  P_Flux_sel -> Counts_R -> Fill(Kbin);
    }

    return;
}


void ProtonFlux_Write() {
	P_Flux      	->Write();
	P_Flux_geo  	->Write();
	P_Flux_geo_prim ->Write();
	P_Flux_pre  	->Write();
	P_Flux_sel  	->Write();

	return;
}


void ProtonFlux() {
   string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
   TFile * file1 = TFile::Open(nomefile.c_str(),"READ");

   TH2F * esposizionegeo_R    = (TH2F*)file1->Get( "esposizionegeo"        );
   TH2F * esposizionepgeoTOF  = (TH2F*)file1->Get(	"esposizionepgeo"	);
   TH2F * esposizionepgeoNaF  = (TH2F*)file1->Get(	"esposizionepgeoNaF"	);
   TH2F * esposizionepgeoAgl  = (TH2F*)file1->Get(	"esposizionepgeoAgl"	);


   Flux * P_Flux         = new Flux(file1, "P_Flux"    	 ,"Results","Corr_AcceptanceP",1);
   Flux * P_Flux_geo     = new Flux(file1, "P_Flux_geo"     ,"Results","Geomag_AcceptanceP",11);
   Flux * P_Flux_geo_prim= new Flux(file1, "P_Flux_geo_prim","Results","Geomag_AcceptanceP",11);

   Flux * P_Flux_pre     = new Flux(file1, "P_Flux_pre" 	 ,"Results","Corr_AcceptancePreP",1);
   Flux * P_Flux_sel     = new Flux(file1, "P_Flux_sel"     ,"Results","Corr_AcceptanceP",1);

   cout<<"*************** PROTONS FLUXES CALCULATION *******************"<<endl;

   P_Flux         -> Set_Exposure_Time (esposizionegeo_R,esposizionepgeoTOF,esposizionepgeoNaF,esposizionepgeoAgl);
   P_Flux_geo     -> Set_Exposure_Time (Tempi);
   P_Flux_geo_prim-> Set_Exposure_Time (Tempi);

   P_Flux_pre     -> Set_Exposure_Time (esposizionegeo_R,esposizionepgeoTOF,esposizionepgeoNaF,esposizionepgeoAgl);
   P_Flux_sel     -> Set_Exposure_Time (esposizionegeo_R,esposizionepgeoTOF,esposizionepgeoNaF,esposizionepgeoAgl);

   enum {protons, deutons};

   P_Flux         -> Eval_Flux(1 , protons);
   P_Flux_geo     -> Eval_Flux(11, protons);
   P_Flux_geo_prim-> Eval_Flux(11, protons);

   P_Flux_pre     -> Eval_Flux(1 , protons);
   P_Flux_sel     -> Eval_Flux(1 , protons);

   TH1F * ProtonsPrimaryFlux = (TH1F *)P_Flux     -> Flux_R ;
   TH2F * ProtonsGeomagFlux  = (TH2F *)P_Flux_geo -> Flux_R ;
   TH1F * P_pre_PrimaryFlux  = (TH1F *)P_Flux_pre -> Flux_R ; //P flux with pre-selections alone
   TH1F * P_sel_PrimaryFlux  = (TH1F *)P_Flux_sel -> Flux_R ; //P flux with full-set selections

   cout<<"*** Updating P1 file ****"<<endl;

   nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
   file1 = TFile::Open(nomefile.c_str(),"UPDATE");
   file1->mkdir("Results/Fluxes");
   file1->cd("Results/Fluxes");
	
   ProtonsPrimaryFlux  ->Write("ProtonsPrimaryFlux" 	);
   ProtonsGeomagFlux   ->Write("ProtonsGeomagFlux"  	);
   P_pre_PrimaryFlux   ->Write("P_pre_PrimaryFlux"  	);
   P_sel_PrimaryFlux   ->Write("P_sel_PrimaryFlux"  	);

   file1->Write();
   file1->Close();



   TCanvas * c23 = new TCanvas("Protons Flux: Geo. Zones");
   TCanvas * c24 = new TCanvas("Primary Protons Flux");
   TCanvas * c25 = new TCanvas("Protons Flux: Pre vs Qual");


   TGraphErrors * P_Fluxgeo[11];
   TGraphErrors * PFlux;
   TGraphErrors * PFlux_pre;
   float potenza=0;

   c23->cd();
   gPad->SetLogx();
   gPad->SetLogy();
   gPad->SetGridx();
   gPad->SetGridy();
   string nome;
   for(int j=0; j<11; j++) {
      nome="Protons Flux: Geo. Zone "+to_string(j);
      P_Fluxgeo[j]=new TGraphErrors();
      P_Fluxgeo[j]->SetName(nome.c_str());
      for(int i=0; i<nbinsr; i++) {
         P_Fluxgeo[j]->SetPoint(i,encinprot[i],ProtonsGeomagFlux->GetBinContent(i+1,j+1)*pow(encinprot[i],potenza));
         P_Fluxgeo[j]->SetPointError(i,0,ProtonsGeomagFlux->GetBinError(i+1,j+1)*pow(encinprot[i],potenza));
      }
      P_Fluxgeo[j]->SetMarkerStyle(8);
      P_Fluxgeo[j]->SetMarkerColor(j-1);
      P_Fluxgeo[j]->SetLineColor(j-1);
      P_Fluxgeo[j]->SetLineWidth(2);
   }
   P_Fluxgeo[10]->SetTitle("Protons Flux: Geo. Zones");
   P_Fluxgeo[10]->GetXaxis()->SetTitle("Kin. En./nucl. [GeV/nucl.]");
   P_Fluxgeo[10]->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
   P_Fluxgeo[10]->GetXaxis()->SetTitleSize(0.045);
   P_Fluxgeo[10]->GetYaxis()->SetTitleSize(0.045);
   P_Fluxgeo[10]->GetYaxis()->SetRangeUser(1e-2,1e4);
   P_Fluxgeo[10]->Draw("AP");
   for(int j=0; j<11; j++) P_Fluxgeo[j]->Draw("Psame");



   c24->cd();
   gPad->SetLogx();
   gPad->SetLogy();
   gPad->SetGridx();
   gPad->SetGridy();
   PFlux=new TGraphErrors();
   for(int i=0; i<nbinsr; i++) {
      PFlux->SetPoint(i,encinprot[i],ProtonsPrimaryFlux->GetBinContent(i+1)*pow(encinprot[i],potenza));
      PFlux->SetPointError(i,0,ProtonsPrimaryFlux->GetBinError(i+1)*pow(encinprot[i],potenza));
   }
   PFlux->SetName("Protons Primary Flux");
   PFlux->SetMarkerStyle(8);
   PFlux->SetMarkerColor(2);
   PFlux->SetTitle("Primary Protons Flux");
   PFlux->GetXaxis()->SetTitle("Kin. En./nucl. [GeV/nucl.]");
   PFlux->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
   PFlux->GetXaxis()->SetTitleSize(0.045);
   PFlux->GetYaxis()->SetTitleSize(0.045);
   PFlux->GetYaxis()->SetRangeUser(1e-2,1e4);
   PFlux->Draw("AP");
   P_Fluxgeo[10]->Draw("Psame");
   TGraph* galprop3P=new TGraph();
   TGraph* galprop3P2=new TGraph();
   float x,y=0;
   int j=0;
   {
      string nomefile="./Galprop/Trotta2011/Def/new_P200.txt";
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
      string nomefile="./Galprop/Trotta2011/Def/new_P1250.txt";
      ifstream fp(nomefile.c_str());
      while (!fp.eof()) {
         fp>>x>>y;
         if(x/1e3>0.05&&x/1e3<=100)
            galprop3P2->SetPoint(j,x/1e3,y*1e7*pow(x/1e3,potenza));
         j++;
      }
   }

   galprop3P->Draw("sameC");
   galprop3P2->Draw("sameC");

   c25->cd();
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraphErrors * PFluxpre=new TGraphErrors();
   TGraphErrors * PFlux_=new TGraphErrors();
   PFlux_pre=new TGraphErrors();
   PFlux_pre->SetName("Protons Primary Flux (only pres.)");
   int p=0;
   for(int i=1; i<nbinsr; i++) {
      PFlux_->SetPoint(p,RB.RigBinCent(i),1);
      PFlux_pre->SetPoint(p,encinprot[i],1);
      PFlux_->SetPointError(p,0,(P_sel_PrimaryFlux->GetBinError(i+1,2)+P_pre_PrimaryFlux->GetBinError(i+1,2))/P_pre_PrimaryFlux->GetBinContent(i+1,1));
      p++;
   }
   p=0;
   for(int i=1; i<nbinsr; i++) {
      PFluxpre->SetPoint(p,RB.RigBinCent(i),P_sel_PrimaryFlux->GetBinContent(i+1)/P_pre_PrimaryFlux->GetBinContent(i+1,1));
      p++;
   }
   PFluxpre->SetMarkerStyle(8);
   PFluxpre->SetMarkerColor(2);
   PFlux_->SetMarkerStyle(4);
   PFlux_->SetFillStyle(3002);
   PFlux_->SetFillColor(4);
   PFlux_->SetMarkerColor(2);
   PFluxpre->SetTitle("Primary Protons Flux");
   PFluxpre->GetXaxis()->SetTitle("R [GV]");
   PFluxpre->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
   PFluxpre->GetXaxis()->SetTitleSize(0.045);
   PFluxpre->GetYaxis()->SetTitleSize(0.045);
   PFluxpre->GetYaxis()->SetRangeUser(0.8,1.2);
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
