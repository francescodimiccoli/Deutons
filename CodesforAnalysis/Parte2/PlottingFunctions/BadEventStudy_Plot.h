TH1F *  ExtractEdep_MC(TH3F * Histo,int i,int rev=0){

	TH2F * Slice = new TH2F("Slice","Slice",Histo->GetNbinsX(),0,Histo->GetXaxis()->GetBinUpEdge(Histo->GetNbinsX()+1),Histo->GetNbinsY(),0,Histo->GetYaxis()->GetBinUpEdge(Histo->GetNbinsY()+1));
	if(rev==0){
		for(int X=0;X<Histo->GetNbinsX();X++)
			for(int Y=0;Y<Histo->GetNbinsX();Y++)
				Slice->SetBinContent(X+1,Y+1,Histo->GetBinContent(X+1,Y+1,i+1));
	}
	else{
		for(int X=0;X<Histo->GetNbinsX();X++)
			for(int Y=0;Y<Histo->GetNbinsX();Y++)
				Slice->SetBinContent(X+1,Y+1,Histo->GetBinContent(Y+1,X+1,i+1));
	}

	TH1F * result = ProjectionXtoTH1F(Slice , "",0,Histo->GetNbinsY()+1);
	return result;	
}







void BadEventStudy_plot(
		TH2F * EdepUtofvsLtof_TOF   ,	
                TH2F * EdepUtofvsLtof_NaF   ,
                TH2F * EdepUtofvsLtof_Agl   ,
                TH2F * EdepUtofvsTrack_TOF  ,
                TH2F * EdepUtofvsTrack_NaF  ,
                TH2F * EdepUtofvsTrack_Agl  ,
                TH2F * EdepLtofvsTrack_TOF  ,
                TH2F * EdepLtofvsTrack_NaF  ,
		TH2F * EdepLtofvsTrack_Agl  ,
                TH3F * EdepUtofvsLtofD_TOF  ,
                TH3F * EdepUtofvsLtofD_NaF  ,
                TH3F * EdepUtofvsLtofD_Agl  ,
                TH3F * EdepUtofvsTrackD_TOF ,
                TH3F * EdepUtofvsTrackD_NaF ,
                TH3F * EdepUtofvsTrackD_Agl ,
                TH3F * EdepLtofvsTrackD_TOF ,
                TH3F * EdepLtofvsTrackD_NaF ,
                TH3F * EdepLtofvsTrackD_Agl ,
                TH2F * EdepUtofvsLtofHe_TOF ,
                TH2F * EdepUtofvsLtofHe_NaF ,
                TH2F * EdepUtofvsLtofHe_Agl ,
                TH2F * EdepUtofvsTrackHe_TOF,
                TH2F * EdepUtofvsTrackHe_NaF,
                TH2F * EdepUtofvsTrackHe_Agl,
                TH2F * EdepLtofvsTrackHe_TOF,
                TH2F * EdepLtofvsTrackHe_NaF,
                TH2F * EdepLtofvsTrackHe_Agl
	)	
{


		
	TH2F * EdepUtofvsLtofS_TOF = (TH2F *) EdepUtofvsLtofD_TOF ->Project3D("yx");
        TH2F * EdepUtofvsLtofS_NaF = (TH2F *) EdepUtofvsLtofD_NaF ->Project3D("yx");
        TH2F * EdepUtofvsLtofS_Agl = (TH2F *) EdepUtofvsLtofD_Agl ->Project3D("yx");
        TH2F * EdepUtofvsTrackS_TOF= (TH2F *) EdepUtofvsTrackD_TOF->Project3D("yx");
        TH2F * EdepUtofvsTrackS_NaF= (TH2F *) EdepUtofvsTrackD_NaF->Project3D("yx");
        TH2F * EdepUtofvsTrackS_Agl= (TH2F *) EdepUtofvsTrackD_Agl->Project3D("yx");
        TH2F * EdepLtofvsTrackS_TOF= (TH2F *) EdepLtofvsTrackD_TOF->Project3D("yx");
        TH2F * EdepLtofvsTrackS_NaF= (TH2F *) EdepLtofvsTrackD_NaF->Project3D("yx");
        TH2F * EdepLtofvsTrackS_Agl= (TH2F *) EdepLtofvsTrackD_Agl->Project3D("yx");



	TH1F *	EdepUtof_TOF = ProjectionXtoTH1F(EdepUtofvsLtof_TOF , "EdepUtof_TOF",0,EdepUtofvsLtof_TOF->GetNbinsY()+1);
	TH1F *  EdepUtof_NaF = ProjectionXtoTH1F(EdepUtofvsLtof_NaF , "EdepUtof_NaF",0,EdepUtofvsLtof_NaF->GetNbinsY()+1);
	TH1F *  EdepUtof_Agl = ProjectionXtoTH1F(EdepUtofvsLtof_Agl , "EdepUtof_Agl",0,EdepUtofvsLtof_Agl->GetNbinsY()+1);
	TH1F *  EdepUtofS_TOF = ProjectionXtoTH1F(EdepUtofvsLtofS_TOF , "EdepUtofS_TOF",0,EdepUtofvsLtofS_TOF->GetNbinsY()+1);
        TH1F *  EdepUtofS_NaF = ProjectionXtoTH1F(EdepUtofvsLtofS_NaF , "EdepUtofS_NaF",0,EdepUtofvsLtofS_NaF->GetNbinsY()+1);
        TH1F *  EdepUtofS_Agl = ProjectionXtoTH1F(EdepUtofvsLtofS_Agl , "EdepUtofS_Agl",0,EdepUtofvsLtofS_Agl->GetNbinsY()+1);
	TH1F *  EdepUtofHe_TOF = ProjectionXtoTH1F(EdepUtofvsLtofHe_TOF , "EdepUtofHe_TOF",0,EdepUtofvsLtofHe_TOF->GetNbinsY()+1);
        TH1F *  EdepUtofHe_NaF = ProjectionXtoTH1F(EdepUtofvsLtofHe_NaF , "EdepUtofHe_NaF",0,EdepUtofvsLtofHe_NaF->GetNbinsY()+1);
        TH1F *  EdepUtofHe_Agl = ProjectionXtoTH1F(EdepUtofvsLtofHe_Agl , "EdepUtofHe_Agl",0,EdepUtofvsLtofHe_Agl->GetNbinsY()+1);

	TH1F *	EdepLtof_TOF = ProjectionYtoTH1F(EdepUtofvsLtof_TOF , "EdepLtof_TOF",0,EdepUtofvsLtof_TOF->GetNbinsX()+1);
	TH1F *  EdepLtof_NaF = ProjectionYtoTH1F(EdepUtofvsLtof_NaF , "EdepLtof_NaF",0,EdepUtofvsLtof_NaF->GetNbinsX()+1);
	TH1F *  EdepLtof_Agl = ProjectionYtoTH1F(EdepUtofvsLtof_Agl , "EdepLtof_Agl",0,EdepUtofvsLtof_Agl->GetNbinsX()+1);
	TH1F *  EdepLtofS_TOF = ProjectionYtoTH1F(EdepUtofvsLtofS_TOF , "EdepLtofS_TOF",0,EdepUtofvsLtofS_TOF->GetNbinsX()+1);
        TH1F *  EdepLtofS_NaF = ProjectionYtoTH1F(EdepUtofvsLtofS_NaF , "EdepLtofS_NaF",0,EdepUtofvsLtofS_NaF->GetNbinsX()+1);
        TH1F *  EdepLtofS_Agl = ProjectionYtoTH1F(EdepUtofvsLtofS_Agl , "EdepLtofS_Agl",0,EdepUtofvsLtofS_Agl->GetNbinsX()+1);
	TH1F *  EdepLtofHe_TOF = ProjectionYtoTH1F(EdepUtofvsLtofHe_TOF , "EdepLtofHe_TOF",0,EdepUtofvsLtofHe_TOF->GetNbinsX()+1);
        TH1F *  EdepLtofHe_NaF = ProjectionYtoTH1F(EdepUtofvsLtofHe_NaF , "EdepLtofHe_NaF",0,EdepUtofvsLtofHe_NaF->GetNbinsX()+1);
        TH1F *  EdepLtofHe_Agl = ProjectionYtoTH1F(EdepUtofvsLtofHe_Agl , "EdepLtofHe_Agl",0,EdepUtofvsLtofHe_Agl->GetNbinsX()+1);

	TH1F *  EdepTrack_TOF = ProjectionYtoTH1F(EdepUtofvsTrack_TOF , "EdepTrack_TOF",0,EdepUtofvsTrack_TOF->GetNbinsX()+1);
        TH1F *  EdepTrack_NaF = ProjectionYtoTH1F(EdepUtofvsTrack_NaF , "EdepTrack_NaF",0,EdepUtofvsTrack_NaF->GetNbinsX()+1);
        TH1F *  EdepTrack_Agl = ProjectionYtoTH1F(EdepUtofvsTrack_Agl , "EdepTrack_Agl",0,EdepUtofvsTrack_Agl->GetNbinsX()+1);
        TH1F *  EdepTrackS_TOF = ProjectionYtoTH1F(EdepUtofvsTrackS_TOF , "EdepTrackS_TOF",0,EdepUtofvsTrackS_TOF->GetNbinsX()+1);
        TH1F *  EdepTrackS_NaF = ProjectionYtoTH1F(EdepUtofvsTrackS_NaF , "EdepTrackS_NaF",0,EdepUtofvsTrackS_NaF->GetNbinsX()+1);
        TH1F *  EdepTrackS_Agl = ProjectionYtoTH1F(EdepUtofvsTrackS_Agl , "EdepTrackS_Agl",0,EdepUtofvsTrackS_Agl->GetNbinsX()+1);
        TH1F *  EdepTrackHe_TOF = ProjectionYtoTH1F(EdepUtofvsTrackHe_TOF , "EdepTrackHe_TOF",0,EdepUtofvsTrackHe_TOF->GetNbinsX()+1);
        TH1F *  EdepTrackHe_NaF = ProjectionYtoTH1F(EdepUtofvsTrackHe_NaF , "EdepTrackHe_NaF",0,EdepUtofvsTrackHe_NaF->GetNbinsX()+1);
        TH1F *  EdepTrackHe_Agl = ProjectionYtoTH1F(EdepUtofvsTrackHe_Agl , "EdepTrackHe_Agl",0,EdepUtofvsTrackHe_Agl->GetNbinsX()+1);
	

		TCanvas * c  = new TCanvas("Upper TOF vs Lower TOF Edep");
		TCanvas * c1 = new TCanvas("Upper TOF vs Inner Tracker Edep");
		TCanvas * c2 = new TCanvas("Lower TOF vs Inner Tracker Edep");
		TCanvas * d  = new TCanvas("Deuterons E.dep. Upper TOF");
		TCanvas * d1 = new TCanvas("Deuterons E.dep. Lower TOF");
		TCanvas * d2 = new TCanvas("Deuterons E.dep. Inner Tracker");
		TCanvas * e  = new TCanvas("Upper TOF Energy dep.");
		TCanvas * e1  = new TCanvas("Lower TOF Energy dep.");
		TCanvas * e2  = new TCanvas("Inner Tracker Energy dep.");
		

		c->Divide(3,1);
		c->cd(1);
		EdepUtofvsLtof_TOF->SetMarkerColor(2);
		EdepUtofvsLtofS_TOF->SetMarkerColor(4);
		EdepUtofvsLtofS_TOF->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Upper TOF)");
		EdepUtofvsLtofS_TOF->GetYaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Lower TOF)");
		EdepUtofvsLtofS_TOF->Draw("CONT3");
		EdepUtofvsLtof_TOF->Draw("same");
		c->cd(2);
		EdepUtofvsLtof_NaF->SetMarkerColor(2);
		EdepUtofvsLtofS_NaF->SetMarkerColor(4);
		EdepUtofvsLtofS_NaF->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Upper TOF)");
                EdepUtofvsLtofS_NaF->GetYaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Lower TOF)");
		EdepUtofvsLtofS_NaF->Draw("CONT3");
		EdepUtofvsLtof_NaF->Draw("same");
		c->cd(3);
		EdepUtofvsLtof_Agl->SetMarkerColor(2);
		EdepUtofvsLtofS_Agl->SetMarkerColor(4);
		EdepUtofvsLtofS_Agl->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Upper TOF)");
                EdepUtofvsLtofS_Agl->GetYaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Lower TOF)");
		EdepUtofvsLtofS_Agl->Draw("CONT3");
		EdepUtofvsLtof_Agl->Draw("same");


		c1->Divide(3,1);
                c1->cd(1);
                EdepUtofvsTrack_TOF->SetMarkerColor(2);
                EdepUtofvsTrackS_TOF->SetMarkerColor(4);
                EdepUtofvsTrackS_TOF->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Upper TOF)");
                EdepUtofvsTrackS_TOF->GetYaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Inner Track)");
                EdepUtofvsTrackS_TOF->Draw("CONT3");
                EdepUtofvsTrack_TOF->Draw("same");
                c1->cd(2);
                gPad->SetLogy();
                gPad->SetLogx();
                EdepUtofvsTrack_NaF->SetMarkerColor(2);
                EdepUtofvsTrackS_NaF->SetMarkerColor(4);
                EdepUtofvsTrackS_NaF->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Upper TOF)");
                EdepUtofvsTrackS_NaF->GetYaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Inner Track)");
                EdepUtofvsTrackS_NaF->Draw("CONT3");
                EdepUtofvsTrack_NaF->Draw("same");
                c1->cd(3);
                gPad->SetLogy();
                gPad->SetLogx();
                EdepUtofvsTrack_Agl->SetMarkerColor(2);
                EdepUtofvsTrackS_Agl->SetMarkerColor(4);
                EdepUtofvsTrackS_NaF->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Upper TOF)");
                EdepUtofvsTrackS_NaF->GetYaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Inner Track)");
		EdepUtofvsTrackS_Agl->Draw("CONT3");
                EdepUtofvsTrack_Agl->Draw("same");

		c2->Divide(3,1);
                c2->cd(1);
                EdepLtofvsTrack_TOF->SetMarkerColor(2);
                EdepLtofvsTrackS_TOF->SetMarkerColor(4);
                EdepLtofvsTrackS_TOF->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Lower TOF)");
                EdepLtofvsTrackS_TOF->GetYaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Inner Track)");
                EdepLtofvsTrackS_TOF->Draw("CONT3");
                EdepLtofvsTrack_TOF->Draw("same");
                c2->cd(2);
                gPad->SetLogy();
                gPad->SetLogx();
                EdepLtofvsTrack_NaF->SetMarkerColor(2);
                EdepLtofvsTrackS_NaF->SetMarkerColor(4);
                EdepLtofvsTrackS_NaF->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Lower TOF)");
                EdepLtofvsTrackS_NaF->GetYaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Inner Track)");
                EdepLtofvsTrackS_NaF->Draw("CONT3");
                EdepLtofvsTrack_NaF->Draw("same");
                c2->cd(3);
                EdepLtofvsTrack_Agl->SetMarkerColor(2);
                EdepLtofvsTrackS_Agl->SetMarkerColor(4);
		EdepLtofvsTrackS_Agl->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Lower TOF)");
                EdepLtofvsTrackS_Agl->GetYaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Inner Track)");
                EdepLtofvsTrackS_Agl->Draw("CONT3");
                EdepLtofvsTrack_Agl->Draw("same");



		TH1F * EdepMCDTOF_Utof[6];
		TH1F * EdepMCDTOF_Ltof[6];
		TH1F * EdepMCDTOF_Track[6];
		TH1F * EdepMCDNaF_Utof[6];
		TH1F * EdepMCDNaF_Ltof[6];
		TH1F * EdepMCDNaF_Track[6];
		TH1F * EdepMCDAgl_Utof[6];
		TH1F * EdepMCDAgl_Ltof[6];
		TH1F * EdepMCDAgl_Track[6];
				



		for(int mc_type=0;mc_type<6;mc_type++){
			EdepMCDTOF_Utof[mc_type]=	ExtractEdep_MC(EdepUtofvsLtofD_TOF,mc_type);		
			EdepMCDNaF_Utof[mc_type]=  	ExtractEdep_MC(EdepUtofvsLtofD_NaF,mc_type);
			EdepMCDAgl_Utof[mc_type]=       ExtractEdep_MC(EdepUtofvsLtofD_Agl,mc_type);
			
			EdepMCDTOF_Ltof[mc_type]=       ExtractEdep_MC(EdepLtofvsTrackD_TOF,mc_type);    
                        EdepMCDNaF_Ltof[mc_type]=       ExtractEdep_MC(EdepLtofvsTrackD_NaF,mc_type);
                        EdepMCDAgl_Ltof[mc_type]=       ExtractEdep_MC(EdepLtofvsTrackD_Agl,mc_type);

			EdepMCDTOF_Track[mc_type]=       ExtractEdep_MC(EdepLtofvsTrackD_TOF,mc_type,1);    
                        EdepMCDNaF_Track[mc_type]=       ExtractEdep_MC(EdepLtofvsTrackD_NaF,mc_type,1);
                        EdepMCDAgl_Track[mc_type]=       ExtractEdep_MC(EdepLtofvsTrackD_Agl,mc_type,1);

		}

		d->Divide(3,1);
                d->cd(1);
                gPad->SetGridy();
                gPad->SetGridx();

		for(int mc_type=0;mc_type<6;mc_type++){	
			EdepMCDTOF_Utof[mc_type]->SetLineColor(mc_type+2);
			EdepMCDTOF_Utof[mc_type]->Scale(1/EdepMCDTOF_Utof[mc_type]->GetEntries());
			if(mc_type==0){ 
				EdepMCDTOF_Utof[mc_type]->Draw();
				EdepMCDTOF_Utof[mc_type]->SetTitle("Energy deposition (Upper TOF) TOF range");
				EdepMCDTOF_Utof[mc_type]->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Upper TOF)");
				}
			else EdepMCDTOF_Utof[mc_type]->Draw("same");
		}

		d->cd(2);
                gPad->SetGridy();
                gPad->SetGridx();

                for(int mc_type=0;mc_type<6;mc_type++){ 
                        EdepMCDNaF_Utof[mc_type]->SetLineColor(mc_type+2);
                        EdepMCDNaF_Utof[mc_type]->Scale(1/EdepMCDNaF_Utof[mc_type]->GetEntries());
                        if(mc_type==0){
                                EdepMCDNaF_Utof[mc_type]->Draw();
                                EdepMCDNaF_Utof[mc_type]->SetTitle("Energy deposition (Upper TOF) NaF range");
                                EdepMCDNaF_Utof[mc_type]->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Upper TOF)");
				}

			else EdepMCDNaF_Utof[mc_type]->Draw("same");
                }

		d->cd(3);
                gPad->SetGridy();
                gPad->SetGridx();
	
                for(int mc_type=0;mc_type<6;mc_type++){
                        EdepMCDAgl_Utof[mc_type]->SetLineColor(mc_type+2);
                        EdepMCDAgl_Utof[mc_type]->Scale(1/EdepMCDAgl_Utof[mc_type]->GetEntries());
                        if(mc_type==0){
                                EdepMCDAgl_Utof[mc_type]->Draw();
                                EdepMCDAgl_Utof[mc_type]->SetTitle("Energy deposition (Upper TOF) Agl range");
                                EdepMCDAgl_Utof[mc_type]->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Upper TOF)");
				}

			else EdepMCDAgl_Utof[mc_type]->Draw("same");
                }



		d1->Divide(3,1);
                d1->cd(1);
                gPad->SetGridy();
                gPad->SetGridx();

                for(int mc_type=0;mc_type<6;mc_type++){
                        EdepMCDTOF_Ltof[mc_type]->SetLineColor(mc_type+2);
                        EdepMCDTOF_Ltof[mc_type]->Scale(1/EdepMCDTOF_Utof[mc_type]->GetEntries());
                        if(mc_type==0){
                                EdepMCDTOF_Ltof[mc_type]->Draw();
                                EdepMCDTOF_Ltof[mc_type]->SetTitle("Energy deposition (Lower TOF) TOF range");
                                EdepMCDTOF_Ltof[mc_type]->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Lower TOF)");
				}
                        else EdepMCDTOF_Ltof[mc_type]->Draw("same");
                }

                d1->cd(2);
                gPad->SetGridy();
                gPad->SetGridx();

                for(int mc_type=0;mc_type<6;mc_type++){
                        EdepMCDNaF_Ltof[mc_type]->SetLineColor(mc_type+2);
                        EdepMCDNaF_Ltof[mc_type]->Scale(1/EdepMCDNaF_Utof[mc_type]->GetEntries());
                        if(mc_type==0){
                                EdepMCDNaF_Ltof[mc_type]->Draw();
                                EdepMCDNaF_Ltof[mc_type]->SetTitle("Energy deposition (Lower TOF) NaF range");
                                EdepMCDNaF_Ltof[mc_type]->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Lower TOF)");
				}

                        else EdepMCDNaF_Ltof[mc_type]->Draw("same");
                }

                d1->cd(3);
                gPad->SetGridy();
                gPad->SetGridx();

                for(int mc_type=0;mc_type<6;mc_type++){
                        EdepMCDAgl_Ltof[mc_type]->SetLineColor(mc_type+2);
                        EdepMCDAgl_Ltof[mc_type]->Scale(1/EdepMCDAgl_Utof[mc_type]->GetEntries());
                        if(mc_type==0){
                                EdepMCDAgl_Ltof[mc_type]->Draw();
                                EdepMCDAgl_Ltof[mc_type]->SetTitle("Energy deposition (Lower TOF) Agl range");
                                EdepMCDAgl_Ltof[mc_type]->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Lower TOF)");
				}

                        else EdepMCDAgl_Ltof[mc_type]->Draw("same");
                }



		d2->Divide(3,1);
                d2->cd(1);
                gPad->SetGridy();
                gPad->SetGridx();

                for(int mc_type=0;mc_type<6;mc_type++){
                        EdepMCDTOF_Track[mc_type]->SetLineColor(mc_type+2);
                        EdepMCDTOF_Track[mc_type]->Scale(1/EdepMCDTOF_Utof[mc_type]->GetEntries());
                        if(mc_type==0){
                                EdepMCDTOF_Track[mc_type]->Draw();
                                EdepMCDTOF_Track[mc_type]->SetTitle("Energy deposition (Inner Tracker) TOF range");
                                EdepMCDTOF_Track[mc_type]->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Inner Tracker)");
				}
                        else EdepMCDTOF_Track[mc_type]->Draw("same");
                }

                d2->cd(2);
                gPad->SetGridy();
                gPad->SetGridx();

                for(int mc_type=0;mc_type<6;mc_type++){
                        EdepMCDNaF_Track[mc_type]->SetLineColor(mc_type+2);
                        EdepMCDNaF_Track[mc_type]->Scale(1/EdepMCDNaF_Utof[mc_type]->GetEntries());
                        if(mc_type==0){
                                EdepMCDNaF_Track[mc_type]->Draw();
                                EdepMCDNaF_Track[mc_type]->SetTitle("Energy deposition (Inner Tracker) NaF range");
                                EdepMCDNaF_Track[mc_type]->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Inner Tracker)");
				}

                        else EdepMCDNaF_Track[mc_type]->Draw("same");
                }

                d2->cd(3);
                gPad->SetGridy();
                gPad->SetGridx();

                for(int mc_type=0;mc_type<6;mc_type++){
                        EdepMCDAgl_Track[mc_type]->SetLineColor(mc_type+2);
                        EdepMCDAgl_Track[mc_type]->Scale(1/EdepMCDAgl_Utof[mc_type]->GetEntries());
                        if(mc_type==0){
                                EdepMCDAgl_Track[mc_type]->Draw();
                                EdepMCDAgl_Track[mc_type]->SetTitle("Energy deposition (Inner Tracker) Agl range");
                                EdepMCDAgl_Track[mc_type]->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Inner Tracker)");
				}

                        else EdepMCDAgl_Track[mc_type]->Draw("same");
                }

		e->Divide(3,1);
                e->cd(1);
                gPad->SetGridy();
                gPad->SetGridx();

		EdepUtof_TOF->Scale(1/EdepUtof_TOF->GetEntries());
		EdepUtofS_TOF->Scale(1/EdepUtofS_TOF->GetEntries());
                EdepUtofHe_TOF->Scale(1/EdepUtofHe_TOF->GetEntries());
		EdepUtof_TOF->SetLineColor(2);
                EdepUtof_TOF->SetLineWidth(2);
                EdepUtofS_TOF->SetLineColor(4);
                EdepUtofS_TOF->SetLineWidth(2);
		EdepUtofHe_TOF->SetLineColor(3);
                EdepUtofHe_TOF->SetLineWidth(2);
		EdepUtofS_TOF->SetTitle("TOF Range");
		EdepUtofS_TOF->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Upper TOF)");
		EdepUtofS_TOF->Draw();
		EdepUtof_TOF->Draw("same");
		EdepUtofHe_TOF->Draw("same");

                e->cd(2);
                gPad->SetGridy();
                gPad->SetGridx();

		EdepUtof_NaF->Scale(1/EdepUtof_NaF->GetEntries());
                EdepUtofS_NaF->Scale(1/EdepUtofS_NaF->GetEntries());
                EdepUtofHe_NaF->Scale(1/EdepUtofHe_NaF->GetEntries());
		EdepUtof_NaF->SetLineColor(2);
                EdepUtof_NaF->SetLineWidth(2);
		EdepUtofS_NaF->SetLineColor(4);
                EdepUtofS_NaF->SetLineWidth(2);
                EdepUtofHe_NaF->SetLineColor(3);
                EdepUtofHe_NaF->SetLineWidth(2);
		EdepUtofS_NaF->SetTitle("NaF Range");
		EdepUtofS_NaF->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Upper TOF)");
		EdepUtofS_NaF->Draw();
		EdepUtof_NaF->Draw("same");
		EdepUtofHe_NaF->Draw("same");

                e->cd(3);
                gPad->SetGridy();
                gPad->SetGridx();

		EdepUtof_Agl->Scale(1/EdepUtof_Agl->GetEntries());
                EdepUtofS_Agl->Scale(1/EdepUtofS_Agl->GetEntries());
                EdepUtofHe_Agl->Scale(1/EdepUtofHe_Agl->GetEntries());
		EdepUtof_Agl->SetLineColor(2);
                EdepUtof_Agl->SetLineWidth(2);
                EdepUtofS_Agl->SetLineColor(4);
                EdepUtofS_Agl->SetLineWidth(2);
		EdepUtofHe_Agl->SetLineColor(3);
                EdepUtofHe_Agl->SetLineWidth(2);
		EdepUtofS_Agl->SetTitle("Agl Range");
		EdepUtofS_Agl->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Upper TOF)");
		EdepUtofS_Agl->Draw();
		EdepUtof_Agl->Draw("same");
		EdepUtofHe_Agl->Draw("same");



		e1->Divide(3,1);
                e1->cd(1);
                gPad->SetGridy();
                gPad->SetGridx();

		EdepLtof_TOF->Scale(1/EdepLtof_TOF->GetEntries());
		EdepLtofS_TOF->Scale(1/EdepLtofS_TOF->GetEntries());
                EdepLtofHe_TOF->Scale(1/EdepLtofHe_TOF->GetEntries());
		EdepLtof_TOF->SetLineColor(2);
                EdepLtof_TOF->SetLineWidth(2);
                EdepLtofS_TOF->SetLineColor(4);
                EdepLtofS_TOF->SetLineWidth(2);
		EdepLtofHe_TOF->SetLineColor(3);
                EdepLtofHe_TOF->SetLineWidth(2);
		EdepLtofS_TOF->SetTitle("TOF Range");
		EdepLtofS_TOF->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Lower TOF)");
		EdepLtofS_TOF->Draw();
		EdepLtof_TOF->Draw("same");
		EdepLtofHe_TOF->Draw("same");

                e1->cd(2);
                gPad->SetGridy();
                gPad->SetGridx();

		EdepLtof_NaF->Scale(1/EdepLtof_NaF->GetEntries());
                EdepLtofS_NaF->Scale(1/EdepLtofS_NaF->GetEntries());
                EdepLtofHe_NaF->Scale(1/EdepLtofHe_NaF->GetEntries());
		EdepLtof_NaF->SetLineColor(2);
                EdepLtof_NaF->SetLineWidth(2);
		EdepLtofS_NaF->SetLineColor(4);
                EdepLtofS_NaF->SetLineWidth(2);
                EdepLtofHe_NaF->SetLineColor(3);
                EdepLtofHe_NaF->SetLineWidth(2);
		EdepLtofS_NaF->SetTitle("NaF Range");
		EdepLtofS_NaF->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Lower TOF)");
		EdepLtofS_NaF->Draw();
		EdepLtof_NaF->Draw("same");
		EdepLtofHe_NaF->Draw("same");

                e1->cd(3);
                gPad->SetGridy();
                gPad->SetGridx();

		EdepLtof_Agl->Scale(1/EdepLtof_Agl->GetEntries());
                EdepLtofS_Agl->Scale(1/EdepLtofS_Agl->GetEntries());
                EdepLtofHe_Agl->Scale(1/EdepLtofHe_Agl->GetEntries());
		EdepLtof_Agl->SetLineColor(2);
                EdepLtof_Agl->SetLineWidth(2);
                EdepLtofS_Agl->SetLineColor(4);
                EdepLtofS_Agl->SetLineWidth(2);
		EdepLtofHe_Agl->SetLineColor(3);
                EdepLtofHe_Agl->SetLineWidth(2);
		EdepLtofS_Agl->SetTitle("Agl Range"); 
		EdepLtofS_Agl->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Inner Tracker)");
		EdepLtofS_Agl->Draw();
		EdepLtof_Agl->Draw("same");
		EdepLtofHe_Agl->Draw("same");
		
		e2->Divide(3,1);
                e2->cd(1);
                gPad->SetGridy();
                gPad->SetGridx();

		EdepTrack_TOF->Scale(1/EdepTrack_TOF->GetEntries());
		EdepTrackS_TOF->Scale(1/EdepTrackS_TOF->GetEntries());
                EdepTrackHe_TOF->Scale(1/EdepTrackHe_TOF->GetEntries());
		EdepTrack_TOF->SetLineColor(2);
                EdepTrack_TOF->SetLineWidth(2);
                EdepTrackS_TOF->SetLineColor(4);
                EdepTrackS_TOF->SetLineWidth(2);
		EdepTrackHe_TOF->SetLineColor(3);
                EdepTrackHe_TOF->SetLineWidth(2);
		EdepTrackS_TOF->SetTitle("TOF Range"); 
		EdepTrackS_TOF->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Inner Tracker)");
		EdepTrackS_TOF->Draw();
		EdepTrack_TOF->Draw("same");
		EdepTrackHe_TOF->Draw("same");

                e2->cd(2);
                gPad->SetGridy();
                gPad->SetGridx();

		EdepTrack_NaF->Scale(1/EdepTrack_NaF->GetEntries());
                EdepTrackS_NaF->Scale(1/EdepTrackS_NaF->GetEntries());
                EdepTrackHe_NaF->Scale(1/EdepTrackHe_NaF->GetEntries());
		EdepTrack_NaF->SetLineColor(2);
                EdepTrack_NaF->SetLineWidth(2);
		EdepTrackS_NaF->SetLineColor(4);
                EdepTrackS_NaF->SetLineWidth(2);
                EdepTrackHe_NaF->SetLineColor(3);
                EdepTrackHe_NaF->SetLineWidth(2);
		EdepTrackS_NaF->SetTitle("NaF Range"); 
		EdepTrackS_NaF->GetXaxis()->SetTitle("E. dep. (meas. - teo.) [# of sigmas] (Inner Tracker)");
		EdepTrackS_NaF->Draw();
		EdepTrack_NaF->Draw("same");
		EdepTrackHe_NaF->Draw("same");

                e2->cd(3);
                gPad->SetGridy();
                gPad->SetGridx();

		EdepTrack_Agl->Scale(1/EdepTrack_Agl->GetEntries());
                EdepTrackS_Agl->Scale(1/EdepTrackS_Agl->GetEntries());
                EdepTrackHe_Agl->Scale(1/EdepTrackHe_Agl->GetEntries());
		EdepTrack_Agl->SetLineColor(2);
                EdepTrack_Agl->SetLineWidth(2);
                EdepTrackS_Agl->SetLineColor(4);
                EdepTrackS_Agl->SetLineWidth(2);
		EdepTrackHe_Agl->SetLineColor(3);
                EdepTrackHe_Agl->SetLineWidth(2);
	 	EdepTrackS_Agl->SetTitle("Agl Range");		
		EdepTrackS_Agl->GetXaxis()->SetTitle("E.. dep (meas. - teo.) [# of sigmas] (Inner Tracker)");
		EdepTrackS_Agl->Draw();
		EdepTrack_Agl->Draw("same");
		EdepTrackHe_Agl->Draw("same");
		






		finalPlots.Add(d);		
		finalPlots.Add(d1);
                finalPlots.Add(d2);
		finalPlots.writeObjsInFolder("MC Results/Quality/D MC distrib");

		finalPlots.Add(c);
		finalPlots.Add(c1);
		finalPlots.Add(c2);
	        finalPlots.Add(e);
		finalPlots.Add(e1);
		finalPlots.Add(e2);
		finalPlots.writeObjsInFolder("MC Results/Background study");	


	

}

