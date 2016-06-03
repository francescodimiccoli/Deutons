using namespace std;

const int NGEOBINS=11;

void TGraphErrorsFormat(TGraphErrors* graph) {
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

PFluxComputing::PFluxComputing(string objname, bool geobin, string acceptname) {
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


TH1* PFluxComputing::ComputeFluxAndGetHisto() {
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


   string nome;

   PBinning PRB; PRB.Setbins(nbinsr, 0.5, 100, 2); // RB did not have Ek
   
   
   for(int j=0; j<NGEOBINS; j++) {
      TGraphErrorsFormat(P_Fluxgeo[j]);
      nome="Protons Flux: Geo. Zone "+to_string(j);
      P_Fluxgeo[j]=new TGraphErrors();
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

   TCanvas * c23 = new TCanvas("Protons Flux: Geo. Zones");
   c23->cd();
   gPad->SetLogx();
   gPad->SetLogy();
   gPad->SetGridx();
   gPad->SetGridy();
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

   TCanvas * c24 = new TCanvas("Primary Protons Flux");
   c24->cd();
   gPad->SetLogx();
   gPad->SetLogy();
   gPad->SetGridx();
   gPad->SetGridy();
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



   TGraphErrors * PFluxpre=new TGraphErrors();
   TGraphErrors * PFlux_=new TGraphErrors();
   PFlux_pre=new TGraphErrors();
   PFlux_pre->SetName("Protons Primary Flux (only pres.)");


   
   for(int ip=0; ip<PRB.EkBinsCent().size()-1; ip++) {
      PFlux_->   SetPoint(ip,PRB.RigBinCent(ip+1), 1);
      PFlux_pre->SetPoint(ip,PRB.EkBinCent(ip+1) , 1);
      float err=(P_sel_PrimaryFlux->GetBinError(ip+2,2)+P_pre_PrimaryFlux->GetBinError(ip+2,2))/P_pre_PrimaryFlux->GetBinContent(ip+2,1);
      PFlux_->SetPointError(ip,0,err);
      PFluxpre->SetPoint(ip,PRB.RigBinCent(ip+1),
         P_sel_PrimaryFlux->GetBinContent(ip+1+1)/P_pre_PrimaryFlux->GetBinContent(ip+1+1,1));

   }

   PFlux_->SetMarkerStyle(4);
   PFlux_->SetFillStyle(3002);
   PFlux_->SetFillColor(4);
   PFlux_->SetMarkerColor(2);

   TGraphErrorsFormat(PFluxpre);
   PFluxpre->SetMarkerColor(2);
   PFluxpre->SetTitle("Primary Protons Flux;R [GV];Flux [(m^2 sr GeV/nucl.)^-1]");
   PFluxpre->GetYaxis()->SetRangeUser(0.8,1.2);




   TCanvas * c25 = new TCanvas("Protons Flux: Pre vs Qual");
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
