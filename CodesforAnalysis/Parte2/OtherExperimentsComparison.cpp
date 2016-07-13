using namespace std;


void OtherExperimentsComparison(){
   inputHistoFile->ReOpen("READ");
	
	string filename2="../database_P.root";
        TFile * file2 = TFile::Open(filename2.c_str(),"READ");
	
	string filename3="../database_D.root";
        TFile * file3 = TFile::Open(filename3.c_str(),"READ");

	cout<<"******************** OTHER EXPERIMENTS COMP. ********************"<<endl;

	cout<<"*** Protons ***"<<endl;
	TList *Experiments = file2->GetListOfKeys();
	TIter next(Experiments);
	TKey * key;
	TObject * obj;	

	std::vector<TGraphAsymmErrors *> P_Graphs;	
	
	while((key = (TKey*)next())){
		obj = file2->Get(key->GetName());
		if(obj->InheritsFrom("TGraphAsymmErrors")) P_Graphs.push_back((TGraphAsymmErrors *)obj); 
	}
	
	TH1F * ThisWork_P = (TH1F*) inputHistoFile -> Get("Results/Fluxes/ProtonsPrimaryFlux");

	cout<<"*** Deutons ***"<<endl;

	std::vector<TGraphAsymmErrors *> D_Graphs;
	
	TList *ExperimentsD = file3->GetListOfKeys();
        TIter nextD(ExperimentsD);
        TKey * keyD;
	
        while((keyD = (TKey*)nextD())){
                obj = file3->Get(keyD->GetName());
		if(obj->InheritsFrom("TGraphAsymmErrors")) D_Graphs.push_back((TGraphAsymmErrors *)obj);
        }
	

	TH1F * ThisWork_DTOF = (TH1F*) inputHistoFile -> Get("Results/Fluxes/DeutonsPrimaryFlux_TOF");
	TH1F * ThisWork_DNaF = (TH1F*) inputHistoFile -> Get("Results/Fluxes/DeutonsPrimaryFlux_NaF"); 	
	TH1F * ThisWork_DAgl = (TH1F*) inputHistoFile -> Get("Results/Fluxes/DeutonsPrimaryFlux_Agl");

	TCanvas * c1 = new TCanvas ("Proton Flux");
	TCanvas * c2 = new TCanvas ("Deuton Flux");
	
	c1 -> cd();
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();	
	int potenza = 0;
	TGraphErrors * ThisWorkP = new TGraphErrors;
	for(int i=0; i<nbinsr; i++) {
                ThisWorkP->SetPoint(i,PRB.EkPerMassBinCent(i),ThisWork_P->GetBinContent(i+1)*pow(PRB.EkPerMassBinCent(i),potenza));
                ThisWorkP->SetPointError(i,0,ThisWork_P->GetBinError(i+1)*pow(PRB.EkPerMassBinCent(i),potenza));
        }
        ThisWorkP->SetName("Protons Primary Flux");
        ThisWorkP->SetMarkerStyle(8);
        ThisWorkP->SetMarkerColor(4);
	ThisWorkP->SetMarkerSize(1.4);
        ThisWorkP->SetTitle("Primary Protons Flux");
        ThisWorkP->SetName(("This Work (" + mese + ")").c_str());
	ThisWorkP->GetXaxis()->SetTitle("Kin. En./nucl. [GeV/nucl.]");
        ThisWorkP->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
        ThisWorkP->GetXaxis()->SetTitleSize(0.045);
        ThisWorkP->GetYaxis()->SetTitleSize(0.045);
        ThisWorkP->GetYaxis()->SetRangeUser(1e-2,1e4);
        ThisWorkP->Draw("AP");
	TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(ThisWorkP,("This Work (" + mese + ")").c_str(), "ep");
	
	for(uint n=0;n<P_Graphs.size();n++){	
		P_Graphs[n] ->Draw("Psame");
		 leg->AddEntry(P_Graphs[n],P_Graphs[n]->GetTitle(),"ep");
	}
	leg -> Draw("same");
	
	c2 -> cd();
	gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();
	potenza = 0;
        TGraphErrors * ThisWorkDTOF = new TGraphErrors;
        for(int i=0; i<nbinsToF; i++) {
                ThisWorkDTOF->SetPoint(i,ToFDB.EkPerMassBinCent(i),ThisWork_DTOF->GetBinContent(i+1)*pow(ToFDB.EkPerMassBinCent(i),potenza));
                ThisWorkDTOF->SetPointError(i,0,ThisWork_DTOF->GetBinError(i+1)*pow(ToFDB.EkPerMassBinCent(i),potenza));
        }
	ThisWorkDTOF->SetPoint(nbinsToF,30,0.000001);
        ThisWorkDTOF->SetPoint(nbinsToF+1,0.01,0.000001);
	ThisWorkDTOF->SetName("Protons Primary Flux");
        ThisWorkDTOF->SetMarkerStyle(8);
        ThisWorkDTOF->SetMarkerColor(4);
        ThisWorkDTOF->SetMarkerSize(1.4);
        ThisWorkDTOF->SetTitle("Primary Deutons Flux");
        ThisWorkDTOF->SetName(("This Work (" + mese + ")").c_str());
        ThisWorkDTOF->GetXaxis()->SetTitle("Kin. En./nucl. [GeV/nucl.]");
        ThisWorkDTOF->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
        ThisWorkDTOF->GetXaxis()->SetTitleSize(0.045);
        ThisWorkDTOF->GetYaxis()->SetTitleSize(0.045);
        ThisWorkDTOF->GetYaxis()->SetRangeUser(1e-3,1e3);
	ThisWorkDTOF->GetXaxis()->SetRangeUser(1e-2,40);

	TGraphErrors * ThisWorkDNaF = new TGraphErrors;
        for(int i=0; i<nbinsNaF; i++) {
                ThisWorkDNaF->SetPoint(i,NaFDB.EkPerMassBinCent(i),ThisWork_DNaF->GetBinContent(i+1)*pow(NaFDB.EkPerMassBinCent(i),potenza));
                ThisWorkDNaF->SetPointError(i,0,ThisWork_DNaF->GetBinError(i+1)*pow(NaFDB.EkPerMassBinCent(i),potenza));
        }
        ThisWorkDNaF->SetName("Protons Primary Flux");
        ThisWorkDNaF->SetMarkerStyle(4);
        ThisWorkDNaF->SetMarkerColor(4);
        ThisWorkDNaF->SetMarkerSize(1.4);
	
	TGraphErrors * ThisWorkDAgl = new TGraphErrors;
        for(int i=0; i<nbinsAgl; i++) {
                ThisWorkDAgl->SetPoint(i,AglDB.EkPerMassBinCent(i),ThisWork_DAgl->GetBinContent(i+1)*pow(AglDB.EkPerMassBinCent(i),potenza));
                ThisWorkDAgl->SetPointError(i,0,ThisWork_DAgl->GetBinError(i+1)*pow(AglDB.EkPerMassBinCent(i),potenza));
        }
        ThisWorkDAgl->SetName("Protons Primary Flux");
        ThisWorkDAgl->SetMarkerStyle(3);
        ThisWorkDAgl->SetMarkerColor(4);
        ThisWorkDAgl->SetMarkerSize(1.4);

	TLegend * legD =new TLegend(0.4, 0.7,0.95,0.95);
        legD->AddEntry(ThisWorkDTOF,("This Work (" + mese + ") - TOF").c_str(), "ep");
	legD->AddEntry(ThisWorkDNaF,("This Work (" + mese + ") - NaF").c_str(), "ep");
	legD->AddEntry(ThisWorkDAgl,("This Work (" + mese + ") - Agl").c_str(), "ep");

	ThisWorkDTOF->Draw("AP");
	ThisWorkDNaF->Draw("Psame");
	ThisWorkDAgl->Draw("Psame");
	for(uint n=0;n<D_Graphs.size();n++){
                D_Graphs[n] ->Draw("Psame");
                 legD->AddEntry(D_Graphs[n],D_Graphs[n]->GetTitle(),"ep");
        }
	legD->Draw("same");


        finalPlots.Add(c1);
	finalPlots.Add(c2);
        finalPlots.writeObjsInFolder("Comparison with others");

	return;

}
