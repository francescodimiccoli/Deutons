using namespace std;


void OtherExperimentsComparison(string filename){

	 cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");
	
	string filename2="../database_P.root";
        TFile * file2 = TFile::Open(filename2.c_str(),"READ");
	
	string filename3="../database_D.root";
        TFile * file3 = TFile::Open(filename3.c_str(),"READ");

	string filename4="../database_PD.root";
        TFile * file4 = TFile::Open(filename4.c_str(),"READ");

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
	
	TH1F * ThisWork_P = (TH1F*) inputHistoFile -> Get("Results/ProtonsPrimaryFlux");

	cout<<"*** Deutons ***"<<endl;

	std::vector<TGraphAsymmErrors *> D_Graphs;
	
	TList *ExperimentsD = file3->GetListOfKeys();
        TIter nextD(ExperimentsD);
        TKey * keyD;
	
        while((keyD = (TKey*)nextD())){
                obj = file3->Get(keyD->GetName());
		if(obj->InheritsFrom("TGraphAsymmErrors")) D_Graphs.push_back((TGraphAsymmErrors *)obj);
        }

	
	cout<<"*** D / P ratio ***"<<endl;

        std::vector<TGraphAsymmErrors *> PD_Graphs;

        TList *ExperimentsPD = file4->GetListOfKeys();
        TIter nextPD(ExperimentsPD);
        TKey * keyPD;

        while((keyPD = (TKey*)nextPD())){
                obj = file4->Get(keyPD->GetName());
                if(obj->InheritsFrom("TGraphAsymmErrors")) PD_Graphs.push_back((TGraphAsymmErrors *)obj);
        }

	cout<<"*** reading database completed ***"<<endl;	

	TH1F * ThisWork_DTOF = (TH1F*) inputHistoFile -> Get("Results/DeutonsPrimaryFlux_TOF");
	TH1F * ThisWork_DNaF = (TH1F*) inputHistoFile -> Get("Results/DeutonsPrimaryFlux_NaF"); 	
	TH1F * ThisWork_DAgl = (TH1F*) inputHistoFile -> Get("Results/DeutonsPrimaryFlux_Agl");

	TH1F * ThisWork_PDTOF = (TH1F*) inputHistoFile -> Get("Results/DP_ratioTOF");
	TH1F * ThisWork_PDNaF = (TH1F*) inputHistoFile -> Get("Results/DP_ratioNaF"); 	
	TH1F * ThisWork_PDAgl = (TH1F*) inputHistoFile -> Get("Results/DP_ratioAgl");


	TCanvas * c1 = new TCanvas ("Proton Flux");
	TCanvas * c2 = new TCanvas ("Deuton Flux");
	TCanvas * c3 = new TCanvas ("D / P ratio");
	
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
	ThisWorkP->SetMarkerSize(2);
        ThisWorkP->SetTitle("Primary Protons Flux");
        ThisWorkP->SetName(("This Work (" + mese + ")").c_str());
	ThisWorkP->GetXaxis()->SetTitle("Kin. En./nucl. [GeV/nucl.]");
        ThisWorkP->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
        ThisWorkP->GetXaxis()->SetTitleSize(0.045);
        ThisWorkP->GetYaxis()->SetTitleSize(0.045);
        ThisWorkP->GetYaxis()->SetRangeUser(1e-2,1e4);
        ThisWorkP->Draw("AP");
	TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(ThisWorkP,("This Work (" + mese + ")").c_str(), "p");
	
	for(uint n=0;n<P_Graphs.size();n++){	
		P_Graphs[n]->SetMarkerSize(2);
		P_Graphs[n] ->Draw("Psame");
		leg->AddEntry(P_Graphs[n],P_Graphs[n]->GetTitle(),"p");
	}
	TGraph* galpropP=new TGraph();
        TGraph* galpropP2=new TGraph();
        float x,y=0;
        int j=0;
        {
                string filename="./Galprop/Tom/prot_1500.dat";
                cout<<filename<<endl;
                ifstream fp(filename.c_str());
                while (!fp.eof()){
                        fp>>x>>y;
                        if(x/1e3>0.05&&x/1e3<=100)
                                galpropP->SetPoint(j,x/1e3,y*1e7*pow(x/1e3,potenza));
                        j++;
                }
        }

        j=0;
        {
                string filename="./Galprop/Tom/prot_100.dat";
                cout<<filename<<endl;
                ifstream fp(filename.c_str());
                while (!fp.eof()){
                        fp>>x>>y;
                        if(x/1e3>0.05&&x/1e3<=100)
                                galpropP2->SetPoint(j,x/1e3,y*1e7*pow(x/1e3,potenza));
                        j++;
                }
        }

        galpropP->SetTitle("GALPROP (#Phi = 400-1500 MV)");

        galpropP->SetLineColor(2);
        galpropP2->SetLineColor(2);
        galpropP->SetLineWidth(4);
        galpropP2->SetLineWidth(4);
        galpropP->SetLineStyle(4);
        galpropP2->SetLineStyle(4);

        galpropP->Draw("sameC");
        galpropP2->Draw("sameC");
	leg->AddEntry(galpropP,galpropP->GetTitle(),"l");
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
        ThisWorkDTOF->SetLineColor(4);
        ThisWorkDTOF->SetMarkerSize(2);
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
        ThisWorkDNaF->SetLineColor(4);
	ThisWorkDNaF->SetMarkerSize(2);
	
	TGraphErrors * ThisWorkDAgl = new TGraphErrors;
        for(int i=0; i<nbinsAgl; i++) {
                ThisWorkDAgl->SetPoint(i,AglDB.EkPerMassBinCent(i),ThisWork_DAgl->GetBinContent(i+1)*pow(AglDB.EkPerMassBinCent(i),potenza));
                ThisWorkDAgl->SetPointError(i,0,ThisWork_DAgl->GetBinError(i+1)*pow(AglDB.EkPerMassBinCent(i),potenza));
        }
        ThisWorkDAgl->SetName("Protons Primary Flux");
        ThisWorkDAgl->SetMarkerStyle(3);
        ThisWorkDAgl->SetMarkerColor(4);
        ThisWorkDAgl->SetLineColor(4);
	ThisWorkDAgl->SetMarkerSize(2);

	TLegend * legD =new TLegend(0.8, 0.1,0.95,0.95);
        legD->AddEntry(ThisWorkDTOF,("This Work (" + mese + ") - TOF").c_str(), "ep");
	legD->AddEntry(ThisWorkDNaF,("This Work (" + mese + ") - NaF").c_str(), "ep");
	legD->AddEntry(ThisWorkDAgl,("This Work (" + mese + ") - Agl").c_str(), "ep");

	ThisWorkDTOF->Draw("AP");
	ThisWorkDNaF->Draw("Psame");
	ThisWorkDAgl->Draw("Psame");
	for(uint n=0;n<D_Graphs.size();n++){
                D_Graphs[n] ->Draw("Psame");
                D_Graphs[n]->SetMarkerSize(2); 
		legD->AddEntry(D_Graphs[n],D_Graphs[n]->GetTitle(),"p");
        }
	
	
	TGraph* galprop3P=new TGraph();
        TGraph* galprop3P2=new TGraph();
        x,y=0;
        j=0;
        {
                string filename="./Galprop/Tom/deut_1500.dat";
                cout<<filename<<endl;
                ifstream fp(filename.c_str());
                while (!fp.eof()){
                        fp>>x>>y;
                        if(x/1e3>0.05&&x/1e3<=100)
                                galprop3P->SetPoint(j,x/1e3,y*1e7*pow(x/1e3,potenza));
                        j++;
                }
        }

        j=0;
        {
                string filename="./Galprop/Tom/deut_100.dat";
                cout<<filename<<endl;
                ifstream fp(filename.c_str());
                while (!fp.eof()){
                        fp>>x>>y;
                        if(x/1e3>0.05&&x/1e3<=100)
                                galprop3P2->SetPoint(j,x/1e3,y*1e7*pow(x/1e3,potenza));
                        j++;
                }
        }

        galprop3P->SetTitle("GALPROP (#Phi = 400-1500 MV)");

	galprop3P->SetLineColor(4);
	galprop3P2->SetLineColor(4);
	galprop3P->SetLineWidth(4);
        galprop3P2->SetLineWidth(4);
	galprop3P->SetLineStyle(4);
        galprop3P2->SetLineStyle(4);

        galprop3P->Draw("sameC");
        galprop3P2->Draw("sameC");

	legD->AddEntry(galprop3P,galprop3P->GetTitle(),"l");
	legD->Draw("same");


	c3 -> cd();
	gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
	potenza = 0;
        TGraphErrors * ThisWorkPDTOF = new TGraphErrors;
        for(int i=0; i<nbinsToF; i++) {
                ThisWorkPDTOF->SetPoint(i,ToFDB.EkPerMassBinCent(i),ThisWork_PDTOF->GetBinContent(i+1));
                ThisWorkPDTOF->SetPointError(i,0,ThisWork_PDTOF->GetBinError(i+1));
        }
	ThisWorkPDTOF->SetPoint(nbinsToF,30,0.000001);
        ThisWorkPDTOF->SetPoint(nbinsToF+1,0.01,0.000001);
	ThisWorkPDTOF->SetName("Protons Primary Flux");
        ThisWorkPDTOF->SetMarkerStyle(8);
        ThisWorkPDTOF->SetMarkerColor(2);
        ThisWorkPDTOF->SetLineColor(2);
	ThisWorkPDTOF->SetMarkerSize(1.4);
        ThisWorkPDTOF->SetTitle("D/P Fluxes ratio");
        ThisWorkPDTOF->SetName(("This Work (" + mese + ")").c_str());
        ThisWorkPDTOF->GetXaxis()->SetTitle("Kin. En./nucl. [GeV/nucl.]");
        ThisWorkPDTOF->GetYaxis()->SetTitle("Flux ratio");
        ThisWorkPDTOF->GetXaxis()->SetTitleSize(0.045);
        ThisWorkPDTOF->GetYaxis()->SetTitleSize(0.045);
        ThisWorkPDTOF->GetYaxis()->SetRangeUser(0.001,0.1);
	ThisWorkPDTOF->GetXaxis()->SetRangeUser(1e-2,40);

	TGraphErrors * ThisWorkPDNaF = new TGraphErrors;
        for(int i=0; i<nbinsNaF; i++) {
                ThisWorkPDNaF->SetPoint(i,NaFDB.EkPerMassBinCent(i),ThisWork_PDNaF->GetBinContent(i+1));
                ThisWorkPDNaF->SetPointError(i,0,ThisWork_PDNaF->GetBinError(i+1));
        }
        
        ThisWorkPDNaF->SetMarkerStyle(4);
        ThisWorkPDNaF->SetMarkerColor(2);
        ThisWorkPDNaF->SetLineColor(2);
	ThisWorkPDNaF->SetMarkerSize(1.4);
	
	TGraphErrors * ThisWorkPDAgl = new TGraphErrors;
        for(int i=0; i<nbinsAgl; i++) {
                ThisWorkPDAgl->SetPoint(i,AglDB.EkPerMassBinCent(i),ThisWork_PDAgl->GetBinContent(i+1));
                ThisWorkPDAgl->SetPointError(i,0,ThisWork_PDAgl->GetBinError(i+1));
        }
       
        ThisWorkPDAgl->SetMarkerStyle(3);
        ThisWorkPDAgl->SetMarkerColor(2);
        ThisWorkPDAgl->SetLineColor(2);
	ThisWorkPDAgl->SetMarkerSize(1.4);

	TLegend * legPD =new TLegend(0.4, 0.7,0.95,0.95);
        legPD->AddEntry(ThisWorkPDTOF,("This Work (" + mese + ") - TOF").c_str(), "ep");
	legPD->AddEntry(ThisWorkPDNaF,("This Work (" + mese + ") - NaF").c_str(), "ep");
	legPD->AddEntry(ThisWorkPDAgl,("This Work (" + mese + ") - Agl").c_str(), "ep");

	ThisWorkPDTOF->Draw("AP");
        ThisWorkPDNaF->Draw("Psame");
        ThisWorkPDAgl->Draw("Psame");


	for(uint n=0;n<PD_Graphs.size();n++){
                PD_Graphs[n] ->Draw("Psame");
                 legPD->AddEntry(PD_Graphs[n],PD_Graphs[n]->GetTitle(),"ep");
        }

	


	
	legPD->Draw("same");




        finalPlots.Add(c1);
	finalPlots.Add(c2);
	finalPlots.Add(c3);
        finalPlots.writeObjsInFolder("Comparison with others");

	return;

}
