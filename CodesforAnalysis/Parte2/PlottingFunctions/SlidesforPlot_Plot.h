
void SlidesforPlot (string filename)
{
    
   cout<<"************** Plots for slides *******************"<<endl;

        cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

   TH2F * RvsBetaTOF_P= (TH2F*) inputHistoFile->Get ("RvsBetaTOF_P");
   TH2F * RvsBetaNaF_P= (TH2F*) inputHistoFile->Get ("RvsBetaNaF_P");
   TH2F * RvsBetaAgl_P= (TH2F*) inputHistoFile->Get ("RvsBetaAgl_P");
   TH2F * RvsBetaTOF_D= (TH2F*) inputHistoFile->Get ("RvsBetaTOF_D");
   TH2F * RvsBetaNaF_D= (TH2F*) inputHistoFile->Get ("RvsBetaNaF_D");
   TH2F * RvsBetaAgl_D= (TH2F*) inputHistoFile->Get ("RvsBetaAgl_D");
   TH2F * RvsBetaTOF_He= (TH2F*) inputHistoFile->Get ("RvsBetaTOF_He");
   TH2F * RvsBetaNaF_He= (TH2F*) inputHistoFile->Get ("RvsBetaNaF_He");
   TH2F * RvsBetaAgl_He= (TH2F*) inputHistoFile->Get ("RvsBetaAgl_He");
   TH2F * EdepUTOFvsR_P= (TH2F*) inputHistoFile->Get ("EdepUTOFvsR_P");
   TH2F * EdepUTOFvsR_D= (TH2F*) inputHistoFile->Get ("EdepUTOFvsR_D");
   TH2F * EdepUTOFvsR_He= (TH2F*) inputHistoFile->Get ("EdepUTOFvsR_He");
   TH2F * EdepLTOFvsR_P= (TH2F*) inputHistoFile->Get ("EdepLTOFvsR_P");
   TH2F * EdepLTOFvsR_D= (TH2F*) inputHistoFile->Get ("EdepLTOFvsR_D");
   TH2F * EdepLTOFvsR_He= (TH2F*) inputHistoFile->Get ("EdepLTOFvsR_He");
   TH2F * EdepTrackvsR_P= (TH2F*) inputHistoFile->Get ("EdepTrackvsR_P");
   TH2F * EdepTrackvsR_D= (TH2F*) inputHistoFile->Get ("EdepTrackvsR_D");
   TH2F * EdepTrackvsR_He= (TH2F*) inputHistoFile->Get ("EdepTrackvsR_He");
   TH1F * MassTOF_P= (TH1F*) inputHistoFile->Get ("MassTOF_P");
   TH1F * MassTOF_D= (TH1F*) inputHistoFile->Get ("MassTOF_D");
   TH1F * MassNaF_P= (TH1F*) inputHistoFile->Get ("MassNaF_P");
   TH1F * MassNaF_D= (TH1F*) inputHistoFile->Get ("MassNaF_D");
   TH1F * MassAgl_P= (TH1F*) inputHistoFile->Get ("MassAgl_P");
   TH1F * MassAgl_D= (TH1F*) inputHistoFile->Get ("MassAgl_D");
   TH1F * MassTOF_PQ= (TH1F*) inputHistoFile->Get ("MassTOF_PQ");
   TH1F * MassTOF_DQ= (TH1F*) inputHistoFile->Get ("MassTOF_DQ");
   TH1F * MassNaF_PQ= (TH1F*) inputHistoFile->Get ("MassNaF_PQ");
   TH1F * MassNaF_DQ= (TH1F*) inputHistoFile->Get ("MassNaF_DQ");
   TH1F * MassAgl_PQ= (TH1F*) inputHistoFile->Get ("MassAgl_PQ");
   TH1F * MassAgl_DQ= (TH1F*) inputHistoFile->Get ("MassAgl_DQ");
   TH2F * RvsBetaTOF= (TH2F*) inputHistoFile->Get ("RvsBetaTOF");
   TH2F * RvsBetaNaF= (TH2F*) inputHistoFile->Get ("RvsBetaNaF");
   TH2F * RvsBetaAgl= (TH2F*) inputHistoFile->Get ("RvsBetaAgl");
   TH1F * MassTOF= (TH1F*) inputHistoFile->Get ("MassTOF");
   TH1F * MassNaF= (TH1F*) inputHistoFile->Get ("MassNaF");
   TH1F * MassAgl= (TH1F*) inputHistoFile->Get ("MassAgl");
   TH1F * MassTOFQ= (TH1F*) inputHistoFile->Get ("MassTOFQ");
   TH1F * MassNaFQ= (TH1F*) inputHistoFile->Get ("MassNaFQ");
   TH1F * MassAglQ= (TH1F*) inputHistoFile->Get ("MassAglQ");
   TH1F * DistTOF_P= (TH1F*) inputHistoFile->Get ("DistTOF_P");
   TH1F * DistNaF_P= (TH1F*) inputHistoFile->Get ("DistNaF_P");
   TH1F * DistAgl_P= (TH1F*) inputHistoFile->Get ("DistAgl_P");
   TH1F * DistTOF_D= (TH1F*) inputHistoFile->Get ("DistTOF_D");
   TH1F * DistNaF_D= (TH1F*) inputHistoFile->Get ("DistNaF_D");
   TH1F * DistAgl_D= (TH1F*) inputHistoFile->Get ("DistAgl_D");
   TH1F * DistTOF_He= (TH1F*) inputHistoFile->Get ("DistTOF_He");
   TH1F * DistNaF_He= (TH1F*) inputHistoFile->Get ("DistNaF_He");
   TH1F * DistAgl_He= (TH1F*) inputHistoFile->Get ("DistAgl_He");
   TH2F * LikvsDistTOF_P= (TH2F*) inputHistoFile->Get ("LikvsDistTOF_P");
   TH2F * LikvsDistNaF_P= (TH2F*) inputHistoFile->Get ("LikvsDistNaF_P");
   TH2F * LikvsDistAgl_P= (TH2F*) inputHistoFile->Get ("LikvsDistAgl_P");
   TH2F * LikvsDistTOF_D= (TH2F*) inputHistoFile->Get ("LikvsDistTOF_D");
   TH2F * LikvsDistNaF_D= (TH2F*) inputHistoFile->Get ("LikvsDistNaF_D");
   TH2F * LikvsDistAgl_D= (TH2F*) inputHistoFile->Get ("LikvsDistAgl_D");

   TH2F * EdepUTOFvsB_P   = (TH2F*) inputHistoFile->Get ("EdepUTOFvsB_P");			
   TH2F * EdepUTOFvsB_D	  = (TH2F*) inputHistoFile->Get ("EdepUTOFvsB_D");	
   TH2F * EdepUTOFvsB_He  = (TH2F*) inputHistoFile->Get ("EdepUTOFvsB_He");	      
   TH2F * EdepLTOFvsB_P	  = (TH2F*) inputHistoFile->Get ("EdepLTOFvsB_P");	
   TH2F * EdepLTOFvsB_D	  = (TH2F*) inputHistoFile->Get ("EdepLTOFvsB_D");	
   TH2F * EdepLTOFvsB_He  = (TH2F*) inputHistoFile->Get ("EdepLTOFvsB_He");	      
   TH2F * EdepTrackvsB_P  = (TH2F*) inputHistoFile->Get ("EdepTrackvsB_P");		
   TH2F * EdepTrackvsB_D  = (TH2F*) inputHistoFile->Get ("EdepTrackvsB_D");	      
   TH2F * EdepTrackvsB_He = (TH2F*) inputHistoFile->Get ("EdepTrackvsB_He");	      
   TH2F * EdepUTOFvsB  	  = (TH2F*) inputHistoFile->Get ("EdepUTOFvsB");	
   TH2F * EdepLTOFvsB  	  = (TH2F*) inputHistoFile->Get ("EdepLTOFvsB");	
   TH2F * EdepTrackvsB 	  = (TH2F*) inputHistoFile->Get ("EdepTrackvsB");	
                                                                         
   TH2F * MassTOFvsB_P	  = (TH2F*) inputHistoFile->Get ("MassTOFvsB_P");	      
   TH2F * MassTOFvsB_D	  = (TH2F*) inputHistoFile->Get ("MassTOFvsB_D");	      
   TH2F * MassNaFvsB_P    = (TH2F*) inputHistoFile->Get ("MassNaFvsB_P");	 		
   TH2F * MassNaFvsB_D	  = (TH2F*) inputHistoFile->Get ("MassNaFvsB_D");	       
   TH2F * MassAglvsB_P	  = (TH2F*) inputHistoFile->Get ("MassAglvsB_P");	       
   TH2F * MassAglvsB_D	  = (TH2F*) inputHistoFile->Get ("MassAglvsB_D");	       



   TCanvas *p1 =new TCanvas ("RvsBeta TOF MC");
   TCanvas *p2 =new TCanvas ("RvsBeta NaF MC");
   TCanvas *p3 =new TCanvas ("RvsBeta Agl MC");
   TCanvas *p4 =new TCanvas ("RvsBeta TOF H.L. data");
   TCanvas *p5 =new TCanvas ("RvsBeta NaF H.L. data");
   TCanvas *p6 =new TCanvas ("RvsBeta Agl H.L. data");
   TCanvas *p4_bis =new TCanvas ("EdepvsBeta TOF H.L. data");
   TCanvas *p5_bis =new TCanvas ("EdepvsBeta NaF H.L. data");
   TCanvas *p6_bis =new TCanvas ("EdepvsBeta Agl H.L. data");
   TCanvas *p7 =new TCanvas ("RvsEdep Upper TOF MC");
   TCanvas *p8 =new TCanvas ("RvsEdep Lower TOF MC");
   TCanvas *p9 =new TCanvas ("RvsEdep Tracker MC");
   TCanvas *p7_bis =new TCanvas ("BetavsEdep Upper TOF MC");
   TCanvas *p8_bis =new TCanvas ("BetavsEdep Lower TOF MC");
   TCanvas *p9_bis =new TCanvas ("BetavsEdep Tracker MC");

   TCanvas *p10=new TCanvas ("Mass TOF H.L. data");
   TCanvas *p11=new TCanvas ("Mass NaF H.L. data");
   TCanvas *p12=new TCanvas ("Mass Agl H.L. data");
   
   TCanvas *p13=new TCanvas ("Mass TOF MC");
   TCanvas *p14=new TCanvas ("Mass NaF MC");
   TCanvas *p15=new TCanvas ("Mass Agl MC");
   
   TCanvas *p13_bis=new TCanvas ("Mass TOF MC vs Beta");
   TCanvas *p14_bis=new TCanvas ("Mass NaF MC vs Beta");
   TCanvas *p15_bis=new TCanvas ("Mass Agl MC vs Beta");
   
  TCanvas *p=new TCanvas ("Beta for D_like Mass");



   TCanvas *p10Q=new TCanvas ("Qual. Mass TOF H.L. data");
   TCanvas *p11Q=new TCanvas ("Qual. Mass NaF H.L. data");
   TCanvas *p12Q=new TCanvas ("Qual. Mass Agl H.L. data");
   TCanvas *p13Q=new TCanvas ("Qual. Mass TOF MC");
   TCanvas *p14Q=new TCanvas ("Qual. Mass NaF MC");
   TCanvas *p15Q=new TCanvas ("Qual. Mass Agl MC");
   TCanvas *p16=new TCanvas ("Distance discr. TOF MC");
   TCanvas *p17=new TCanvas ("Distance discr. NaF MC");
   TCanvas *p18=new TCanvas ("Distance discr. Agl MC");
   
   TCanvas *p19=new TCanvas ("LikvsDist  TOF MC");
   TCanvas *p20=new TCanvas ("LikvsDist  NaF MC");
   TCanvas *p21=new TCanvas ("LikvsDist  Agl MC");



   cout<<"******************* R vs Beta plots ******************"<<endl;

   {
      p1->cd();
      gPad->SetGridx();
      gPad->SetGridy();
      RvsBetaTOF_P->SetMarkerColor (2);
      RvsBetaTOF_D->SetMarkerColor (4);
      RvsBetaTOF_He->SetMarkerColor (3);
      RvsBetaTOF_P->SetFillColor (2);
      RvsBetaTOF_D->SetFillColor (4);
      RvsBetaTOF_He->SetFillColor (3);

      RvsBetaTOF_D->SetTitle ("R vs Beta TOF (MC)");
      RvsBetaTOF_D->GetXaxis()->SetTitle ("R [GV]");
      RvsBetaTOF_D->GetYaxis()->SetTitle ("Beta TOF");
      RvsBetaTOF_He->GetZaxis()->SetRangeUser (1,100);
      RvsBetaTOF_D->Draw();
      RvsBetaTOF_P->Draw ("same");
     // RvsBetaTOF_He->Draw ("same");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (RvsBetaTOF_P,"Protons MC", "f");
      leg->AddEntry (RvsBetaTOF_D,"Deuterons MC", "f");
     // leg->AddEntry (RvsBetaTOF_He,"Q=1 from He fragm.", "ep");
      leg->Draw ("same");
   }

   {
      p2->cd();
      gPad->SetGridx();
      gPad->SetGridy();
      RvsBetaNaF_P->SetMarkerColor (2);
      RvsBetaNaF_D->SetMarkerColor (4);
      RvsBetaNaF_He->SetMarkerColor (3);
      RvsBetaNaF_P->SetFillColor (2);
      RvsBetaNaF_D->SetFillColor (4);
      RvsBetaNaF_He->SetFillColor (3);
      RvsBetaNaF_D->SetTitle ("R vs Beta NaF (MC)");
      RvsBetaNaF_D->GetXaxis()->SetTitle ("R [GV]");
      RvsBetaNaF_D->GetYaxis()->SetTitle ("Beta NaF");
      RvsBetaNaF_D->Draw();
      RvsBetaNaF_P->Draw ("same");
     // RvsBetaNaF_He->Draw ("same");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (RvsBetaTOF_P,"Protons MC", "f");
      leg->AddEntry (RvsBetaTOF_D,"Deuterons MC", "f");
    //  leg->AddEntry (RvsBetaTOF_He,"Q=1 from He fragm.", "ep");
      leg->Draw ("same");
   }

   {
      p3->cd();
      gPad->SetGridx();
      gPad->SetGridy();
      RvsBetaAgl_P->SetMarkerColor (2);
      RvsBetaAgl_D->SetMarkerColor (4);
      RvsBetaAgl_He->SetMarkerColor (3);
      RvsBetaAgl_P->SetFillColor (2);
      RvsBetaAgl_D->SetFillColor (4);
      RvsBetaAgl_He->SetFillColor (3);
      RvsBetaAgl_D->SetTitle ("R vs Beta Agl (MC)");
      RvsBetaAgl_D->GetXaxis()->SetTitle ("R [GV]");
      RvsBetaAgl_D->GetYaxis()->SetTitle ("Beta Agl");
      RvsBetaAgl_D->Draw();
      RvsBetaAgl_P->Draw ("same");
    //  RvsBetaAgl_He->Draw ("same");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (RvsBetaTOF_P,"Protons MC", "f");
      leg->AddEntry (RvsBetaTOF_D,"Deuterons MC", "f");
    //  leg->AddEntry (RvsBetaTOF_He,"Q=1 from He fragm.", "ep");
      leg->Draw ("same");
   }

   {
      p4->cd();
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogz();
      RvsBetaTOF->SetTitle ("R vs Beta TOF (H.L. ISS data)");
      RvsBetaTOF->GetXaxis()->SetTitle ("R [GV]");
      RvsBetaTOF->GetYaxis()->SetTitle ("Beta TOF");
      RvsBetaTOF->Draw ("col");
      protons->SetLineColor (2);
      deutons->SetLineColor (4);
      protons->SetLineWidth (3);
      deutons->SetLineWidth (3);
      protons->Draw ("same");
      deutons->Draw ("same");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (protons,"Protons Mass curve", "ep");
      leg->AddEntry (deutons,"Deuterons Mass curve", "ep");
      leg->Draw ("same");
   }

   {
      p5->cd();
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogz();
      RvsBetaNaF->SetTitle ("R vs Beta NaF (H.L. ISS data)");
      RvsBetaNaF->GetXaxis()->SetTitle ("R [GV]");
      RvsBetaNaF->GetYaxis()->SetTitle ("Beta NaF");
      RvsBetaNaF->Draw ("col");
      protons->SetLineColor (2);
      deutons->SetLineColor (4);
      protons->SetLineWidth (3);
      deutons->SetLineWidth (3);
      protons->Draw ("same");
      deutons->Draw ("same");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (protons,"Protons Mass curve", "ep");
      leg->AddEntry (deutons,"Deuterons Mass curve", "ep");
      leg->Draw ("same");
   }

   {
      p6->cd();
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogz();
      RvsBetaAgl->SetTitle ("R vs Beta Agl (H.L. ISS data)");
      RvsBetaAgl->GetXaxis()->SetTitle ("R [GV]");
      RvsBetaAgl->GetYaxis()->SetTitle ("Beta Agl");
      RvsBetaAgl->Draw ("col");
      protons->SetLineColor (2);
      deutons->SetLineColor (4);
      protons->SetLineWidth (3);
      deutons->SetLineWidth (3);
      protons->Draw ("same");
      deutons->Draw ("same");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (protons,"Protons Mass curve", "ep");
      leg->AddEntry (deutons,"Deuterons Mass curve", "ep");
      leg->Draw ("same");
   }

   cout<<"******************* Edep plots ******************"<<endl;
   
  {
      p4_bis->cd();
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogz();
      EdepUTOFvsB->SetTitle ("#beta vs E. dep. Upper TOF (H.L. ISS data)");
      EdepUTOFvsB->GetXaxis()->SetTitle ("#beta");
      EdepUTOFvsB->GetYaxis()->SetTitle ("E. dep. Upper TOF [MeV]");
      EdepUTOFvsB->Draw ("col");
   }

   {
      p5_bis->cd();
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogz();
      EdepLTOFvsB->SetTitle ("#beta vs E. dep. Lower TOF (H.L. ISS data)");
      EdepLTOFvsB->GetXaxis()->SetTitle ("#beta");
      EdepLTOFvsB->GetYaxis()->SetTitle ("E. dep. Lower TOF [MeV]");
      EdepLTOFvsB->Draw ("col");
   }

   {
      p6_bis->cd();
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogz();
      EdepTrackvsB->SetTitle ("#beta vs E. dep. Inner Tracker (H.L. ISS data)");
      EdepTrackvsB->GetXaxis()->SetTitle ("#beta");
      EdepTrackvsB->GetYaxis()->SetTitle ("E. dep. Inner Tracker [keV]");
      EdepTrackvsB->Draw ("col");
   }
 




   p7->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   EdepUTOFvsR_P->SetMarkerColor (2);
   EdepUTOFvsR_D->SetMarkerColor (4);
   EdepUTOFvsR_He->SetMarkerColor (3);
   EdepUTOFvsR_P->SetFillColor (2);
   EdepUTOFvsR_D->SetFillColor (4);  

   EdepUTOFvsR_D->SetTitle ("R vs E.dep. Upper TOF (MC)");
   EdepUTOFvsR_D->GetXaxis()->SetTitle ("R [GV]");
   EdepUTOFvsR_D->GetYaxis()->SetTitle ("E. dep. Upper TOF [MeV]");
   EdepUTOFvsR_He->GetZaxis()->SetRangeUser (100,10000);
   EdepUTOFvsR_D->GetZaxis()->SetRangeUser (10,1000);
   EdepUTOFvsR_P->GetZaxis()->SetRangeUser (150,10000);
   EdepUTOFvsR_D->GetXaxis()->SetRangeUser (0,4);
   EdepUTOFvsR_D->GetYaxis()->SetRangeUser (0,20);
   EdepUTOFvsR_D->Draw();
   EdepUTOFvsR_P->Draw ("same");
   EdepUTOFvsR_He->Draw ("same");
   {
   TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (EdepUTOFvsR_P ,"Protons MC", "f");
      leg->AddEntry (EdepUTOFvsR_D,"Deuterons MC", "f");
      leg->AddEntry (EdepUTOFvsR_He,"Helium MC", "f");	
      leg->Draw ("same");
   }

   
   p7_bis->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   EdepUTOFvsB_P->SetMarkerColor (2);
   EdepUTOFvsB_D->SetMarkerColor (4);
   EdepUTOFvsB_He->SetMarkerColor (3);
   EdepUTOFvsB_P->SetFillColor (2);
   EdepUTOFvsB_D->SetFillColor (4);  
   EdepUTOFvsB_He->SetFillColor (3);

   EdepUTOFvsB_D->SetTitle ("#beta vs E.dep. Upper TOF (MC)");
   EdepUTOFvsB_D->GetXaxis()->SetTitle ("#beta");
   EdepUTOFvsB_D->GetYaxis()->SetTitle ("E. dep. Upper TOF [MeV]");
   EdepUTOFvsB_He->GetZaxis()->SetRangeUser (100,10000);
   EdepUTOFvsB_D->GetZaxis()->SetRangeUser (10,1000);
   EdepUTOFvsB_P->GetZaxis()->SetRangeUser (150,10000);
   EdepUTOFvsB_D->GetXaxis()->SetRangeUser (0.3,1);
   EdepUTOFvsB_D->GetYaxis()->SetRangeUser (0,20);
   EdepUTOFvsB_D->Draw();
   EdepUTOFvsB_P->Draw ("same");
   EdepUTOFvsB_He->Draw ("same");
   {
   TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (EdepUTOFvsB_P ,"Protons MC", "f");
      leg->AddEntry (EdepUTOFvsB_D,"Deuterons MC", "f");
      leg->AddEntry (EdepUTOFvsB_He,"Helium MC", "f");	
      leg->Draw ("same");
   }

   
   p8->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   EdepLTOFvsR_P->SetMarkerColor (2);
   EdepLTOFvsR_D->SetMarkerColor (4);
   EdepLTOFvsR_He->SetMarkerColor (3);
   EdepLTOFvsR_P->SetFillColor (2);
   EdepLTOFvsR_D->SetFillColor (4);
   
   EdepLTOFvsR_D->SetTitle ("R vs E.dep. Lower TOF (MC)");
   EdepLTOFvsR_D->GetXaxis()->SetTitle ("R [GV]");
   EdepLTOFvsR_D->GetYaxis()->SetTitle ("E. dep. Upper TOF [MeV]");
   EdepLTOFvsR_He->GetZaxis()->SetRangeUser (100,10000);
   EdepLTOFvsR_D->GetZaxis()->SetRangeUser (10,1000);
   EdepLTOFvsR_P->GetZaxis()->SetRangeUser (150,10000);
   EdepLTOFvsR_D->GetXaxis()->SetRangeUser (0,4);
   EdepLTOFvsR_D->GetYaxis()->SetRangeUser (0,20);
   EdepLTOFvsR_D->Draw();
   EdepLTOFvsR_P->Draw ("same");
   EdepLTOFvsR_He->Draw ("same");
   {
   TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (EdepLTOFvsR_P ,"Protons MC", "f");
      leg->AddEntry (EdepLTOFvsR_D,"Deuterons MC", "f");
      leg->AddEntry (EdepLTOFvsR_He,"Helium MC", "f");           
      leg->Draw ("same");
   }
	
   p8_bis->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   EdepLTOFvsB_P->SetMarkerColor (2);
   EdepLTOFvsB_D->SetMarkerColor (4);
   EdepLTOFvsB_He->SetMarkerColor (3);
   EdepLTOFvsB_P->SetFillColor (2);
   EdepLTOFvsB_D->SetFillColor (4);
   EdepLTOFvsB_He->SetFillColor (3); 
   
   EdepLTOFvsB_D->SetTitle ("#beta vs E.dep. Lower TOF (MC)");
   EdepLTOFvsB_D->GetXaxis()->SetTitle ("#beta");
   EdepLTOFvsB_D->GetYaxis()->SetTitle ("E. dep. Upper TOF [MeV]");
   EdepLTOFvsB_He->GetZaxis()->SetRangeUser (100,10000);
   EdepLTOFvsB_D->GetZaxis()->SetRangeUser (10,1000);
   EdepLTOFvsB_P->GetZaxis()->SetRangeUser (150,10000);
   EdepLTOFvsB_D->GetXaxis()->SetRangeUser (0.3,1);
   EdepLTOFvsB_D->GetYaxis()->SetRangeUser (0,20);
   EdepLTOFvsB_D->Draw();
   EdepLTOFvsB_P->Draw ("same");
   EdepLTOFvsB_He->Draw ("same");
   {
   TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (EdepLTOFvsB_P ,"Protons MC", "f");
      leg->AddEntry (EdepLTOFvsB_D,"Deuterons MC", "f");
      leg->AddEntry (EdepLTOFvsB_He,"Helium MC", "f");           
      leg->Draw ("same");
   }
	

   p9->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   EdepTrackvsR_P->SetMarkerColor (2);
   EdepTrackvsR_D->SetMarkerColor (4);
   EdepTrackvsR_He->SetMarkerColor (3);
   EdepTrackvsR_P->SetFillColor (2);
   EdepTrackvsR_D->SetFillColor (4);
   EdepTrackvsR_He->SetFillColor (3);
   EdepTrackvsR_D->SetTitle ("R vs E.dep. Tracker (MC)");
   EdepTrackvsR_D->GetXaxis()->SetTitle ("R [GV]");
   EdepTrackvsR_D->GetYaxis()->SetTitle ("E. dep. Upper TOF [MeV]");
   EdepTrackvsR_He->GetZaxis()->SetRangeUser (100,10000);
   EdepTrackvsR_D->GetZaxis()->SetRangeUser (10,1000);
   EdepTrackvsR_P->GetZaxis()->SetRangeUser (150,10000);
   EdepTrackvsR_D->GetXaxis()->SetRangeUser (0,4);
   EdepTrackvsR_D->GetYaxis()->SetRangeUser (0,1.22);
   EdepTrackvsR_D->Draw();
   EdepTrackvsR_P->Draw ("same");
   EdepTrackvsR_He->Draw ("same");
   {
   TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (EdepTrackvsR_P ,"Protons MC", "f");
      leg->AddEntry (EdepTrackvsR_D,"Deuterons MC", "f");
      leg->AddEntry (EdepTrackvsR_He,"Helium MC", "f");           
      leg->Draw ("same");
   }

   p9_bis->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   EdepTrackvsB_P->SetMarkerColor (2);
   EdepTrackvsB_D->SetMarkerColor (4);
   EdepTrackvsB_He->SetMarkerColor (3);
   EdepTrackvsB_P->SetFillColor (2);
   EdepTrackvsB_D->SetFillColor (4);
   EdepTrackvsB_He->SetFillColor (3);
   EdepTrackvsB_D->SetTitle ("#beta vs E.dep. Tracker (MC)");
   EdepTrackvsB_D->GetXaxis()->SetTitle ("#beta");
   EdepTrackvsB_D->GetYaxis()->SetTitle ("E. dep. Upper TOF [MeV]");
   EdepTrackvsB_He->GetZaxis()->SetRangeUser (100,10000);
   EdepTrackvsB_D->GetZaxis()->SetRangeUser (10,1000);
   EdepTrackvsB_P->GetZaxis()->SetRangeUser (150,10000);
   EdepTrackvsB_D->GetXaxis()->SetRangeUser (0.3,1);
   EdepTrackvsB_D->GetYaxis()->SetRangeUser (0,1.22);
   EdepTrackvsB_D->Draw();
   EdepTrackvsB_P->Draw ("same");
   EdepTrackvsB_He->Draw ("same");
   {
   TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (EdepTrackvsB_P ,"Protons MC", "f");
      leg->AddEntry (EdepTrackvsB_D,"Deuterons MC", "f");
      leg->AddEntry (EdepTrackvsB_He,"Helium MC", "f");           
      leg->Draw ("same");
   }



   cout<<"******************* Mass plots ******************"<<endl;

   p10->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassTOF->SetLineColor (1);
   MassTOF->SetLineWidth (2);
   MassTOF->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassTOF->GetXaxis()->SetTitleSize (0.045);
   MassTOF->SetTitle ("Mass TOF H.L. data");
   MassTOF->Draw();

   p11->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassNaF->SetLineColor (1);
   MassNaF->SetLineWidth (2);
   MassNaF->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassNaF->SetTitle ("Mass NaF H.L. data");
   MassNaF->GetXaxis()->SetTitleSize (0.045);
   MassNaF->Draw();

   p12->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassAgl->SetLineColor (1);
   MassAgl->SetLineWidth (2);
   MassAgl->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassAgl->SetTitle ("Mass Agl H.L. data");
   MassAgl->GetXaxis()->SetTitleSize (0.045);
   MassAgl->Draw();

   p13->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassTOF_P->SetLineColor (2);
   MassTOF_P->SetLineWidth (2);
   MassTOF_D->SetLineColor (4);
   MassTOF_D->SetLineWidth (2);
   MassTOF_D->SetFillColor (4);
   MassTOF_P->SetFillColor (2);
   MassTOF_P->SetFillStyle (3001);
   MassTOF_D->SetFillStyle (3002);
   MassTOF_P->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassTOF_P->GetXaxis()->SetTitleSize (0.045);
   MassTOF_P->SetTitle ("Mass TOF MC");
   MassTOF_P->Draw();
   MassTOF_D->Draw ("same");

   p13_bis->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   MassTOFvsB_P->SetMarkerColor (2);
   MassTOFvsB_P->SetLineWidth (2);
   MassTOFvsB_D->SetMarkerColor (4);
   MassTOFvsB_D->SetLineWidth (2);
   MassTOFvsB_D->SetFillColor (4);
   MassTOFvsB_P->SetFillColor (2);
   MassTOFvsB_D->GetYaxis()->SetTitle ("#beta");
   MassTOFvsB_D->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassTOFvsB_D->GetXaxis()->SetTitleSize (0.045);
   MassTOFvsB_D->SetTitle ("Mass TOF MC");
   MassTOFvsB_D->Draw();
   MassTOFvsB_P->Draw ("same");
   {
   TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry ( MassTOFvsB_P,"Protons MC ", "f");
      leg->AddEntry ( MassTOFvsB_D,"Deuterons MC ", "f");
      leg->Draw ("same");
   }	

   

   p14->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassNaF_P->SetLineColor (2);
   MassNaF_P->SetLineWidth (2);
   MassNaF_D->SetLineColor (4);
   MassNaF_D->SetLineWidth (2);
   MassNaF_D->SetFillColor (4);
   MassNaF_P->SetFillColor (2);
   MassNaF_P->SetFillStyle (3001);
   MassNaF_D->SetFillStyle (3002);
   MassNaF_P->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassNaF_P->SetTitle ("Mass NaF MC");
   MassNaF_P->GetXaxis()->SetTitleSize (0.045);
   MassNaF_P->Draw();
   MassNaF_D->Draw ("same");

   p14_bis->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   MassNaFvsB_P->SetMarkerColor (2);
   MassNaFvsB_P->SetLineWidth (2);
   MassNaFvsB_D->SetMarkerColor (4);
   MassNaFvsB_D->SetLineWidth (2);
   MassNaFvsB_D->SetFillColor (4);
   MassNaFvsB_P->SetFillColor (2);
   MassNaFvsB_D->GetYaxis()->SetTitle ("#beta");
   MassNaFvsB_D->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassNaFvsB_D->GetXaxis()->SetTitleSize (0.045);
   MassNaFvsB_D->SetTitle ("Mass TOF MC");
   MassNaFvsB_D->Draw();
   MassNaFvsB_P->Draw ("same");
   {
   TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry ( MassNaFvsB_P,"Protons MC ", "f");
      leg->AddEntry ( MassNaFvsB_D,"Deuterons MC ", "f");
      leg->Draw ("same");
   }




   p15->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassAgl_P->SetLineColor (2);
   MassAgl_P->SetLineWidth (2);
   MassAgl_D->SetLineColor (4);
   MassAgl_D->SetLineWidth (2);
   MassAgl_D->SetFillColor (4);
   MassAgl_P->SetFillColor (2);
   MassAgl_P->SetFillStyle (3001);
   MassAgl_D->SetFillStyle (3002);
   MassAgl_P->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassAgl_P->SetTitle ("Mass Agl MC");
   MassAgl_P->GetXaxis()->SetTitleSize (0.045);
   MassAgl_P->Draw();
   MassAgl_D->Draw ("same");

   p15_bis->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   MassAglvsB_P->SetMarkerColor (2);
   MassAglvsB_P->SetLineWidth (2);
   MassAglvsB_D->SetMarkerColor (4);
   MassAglvsB_D->SetLineWidth (2);
   MassAglvsB_D->SetFillColor (4);
   MassAglvsB_P->SetFillColor (2);
   MassAglvsB_D->GetYaxis()->SetTitle ("#beta");
   MassAglvsB_D->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassAglvsB_D->GetXaxis()->SetTitleSize (0.045);
   MassAglvsB_D->SetTitle ("Mass TOF MC");
   MassAglvsB_D->Draw();
   MassAglvsB_P->Draw ("same");
   {
   TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry ( MassAglvsB_P,"Protons MC ", "f");
      leg->AddEntry ( MassAglvsB_D,"Deuterons MC ", "f");
      leg->Draw ("same");
   }
  



   p->Divide(3,1);
   p->cd(1);
   gPad->SetGridx();
   gPad->SetGridy();
   TH1F * BetaTOF_P = ProjectionYtoTH1F(MassTOFvsB_P    , "temp2",207,500);    	
   TH1F * BetaTOF_D = ProjectionYtoTH1F(MassTOFvsB_D    , "temp2",207,500);
   BetaTOF_P -> Scale(1/BetaTOF_P->GetEntries() );  
   BetaTOF_D -> Scale(1/BetaTOF_D->GetEntries() );
   BetaTOF_P -> SetLineColor(2);
   BetaTOF_D -> SetLineColor(4);
   BetaTOF_P -> SetLineWidth(2);
   BetaTOF_D -> SetLineWidth(2);
   BetaTOF_P -> SetFillStyle(3001);
   BetaTOF_D -> SetFillStyle(3001); 
   BetaTOF_P -> SetFillColor(2);
   BetaTOF_D -> SetFillColor(4);
   BetaTOF_P -> GetXaxis()->SetTitle ("#beta");
   BetaTOF_P -> SetTitleSize (0.045);
   BetaTOF_P ->	Draw();
   BetaTOF_D -> Draw("same");
   {
   TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (BetaTOF_P ,"Protons MC (rec. mass > 1.875)", "f");
      leg->AddEntry (BetaTOF_D ,"Deuterons MC (rec. mass > 1.875)", "f");
      leg->Draw ("same");
   }
   
    p->cd(2);
   gPad->SetGridx();
   gPad->SetGridy();
   TH1F * BetaNaF_P = ProjectionYtoTH1F(MassNaFvsB_P    , "temp2",207,500);    	
   TH1F * BetaNaF_D = ProjectionYtoTH1F(MassNaFvsB_D    , "temp2",207,500);
   BetaNaF_P -> Scale(1/BetaTOF_P->GetEntries() );  
   BetaNaF_D -> Scale(1/BetaTOF_D->GetEntries() );
   BetaNaF_P -> SetLineColor(2);
   BetaNaF_D -> SetLineColor(4);
   BetaNaF_P -> SetLineWidth(2);
   BetaNaF_D -> SetLineWidth(2);
   BetaNaF_P -> SetFillStyle(3001);
   BetaNaF_D -> SetFillStyle(3001); 
   BetaNaF_P -> SetFillColor(2);
   BetaNaF_D -> SetFillColor(4);
   BetaNaF_P -> GetXaxis()->SetTitle ("#beta");
   BetaNaF_P -> SetTitleSize (0.045);
   BetaNaF_P ->	Draw();
   BetaNaF_D -> Draw("same");
   {
   TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (BetaNaF_P ,"Protons MC (rec. mass > 1.875)", "f");
      leg->AddEntry (BetaNaF_D ,"Deuterons MC (rec. mass > 1.875)", "f");
      leg->Draw ("same");
   }

    p->cd(3);
   gPad->SetGridx();
   gPad->SetGridy();
   TH1F * BetaAgl_P = ProjectionYtoTH1F(MassAglvsB_P    , "temp2",207,500);    	
   TH1F * BetaAgl_D = ProjectionYtoTH1F(MassAglvsB_D    , "temp2",207,500);
   BetaAgl_P -> Scale(1/BetaTOF_P->GetEntries() );  
   BetaAgl_D -> Scale(1/BetaTOF_D->GetEntries() );
   BetaAgl_P -> SetLineColor(2);
   BetaAgl_D -> SetLineColor(4);
   BetaAgl_P -> SetLineWidth(2);
   BetaAgl_D -> SetLineWidth(2);
   BetaAgl_P -> SetFillStyle(3001);
   BetaAgl_D -> SetFillStyle(3001); 
   BetaAgl_P -> SetFillColor(2);
   BetaAgl_D -> SetFillColor(4);
   BetaAgl_P -> GetXaxis()->SetTitle ("#beta");
   BetaAgl_P -> SetTitleSize (0.045);
   BetaAgl_P ->	Draw();
   BetaAgl_D -> Draw("same");
   {
   TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (BetaAgl_P ,"Protons MC (rec. mass > 1.875)", "f");
      leg->AddEntry (BetaAgl_D ,"Deuterons MC (rec. mass > 1.875)", "f");
      leg->Draw ("same");
   }






   cout<<"******************* Quality Mass plots ******************"<<endl;

   p10Q->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassTOFQ->SetLineColor (1);
   MassTOFQ->SetLineWidth (2);
   MassTOFQ->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassTOFQ->SetTitle ("Mass TOF H.L. data");
   MassTOFQ->GetXaxis()->SetTitleSize (0.045);
   MassTOFQ->Draw();

   p11Q->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassNaFQ->SetLineColor (1);
   MassNaFQ->SetLineWidth (2);
   MassNaFQ->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassNaFQ->SetTitle ("Mass NaF H.L. data");
   MassNaFQ->GetXaxis()->SetTitleSize (0.045);
   MassNaFQ->Draw();

   p12Q->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassAglQ->SetLineColor (1);
   MassAglQ->SetLineWidth (2);
   MassAglQ->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassAglQ->SetTitle ("Mass Agl H.L. data");
   MassAglQ->GetXaxis()->SetTitleSize (0.045);
   MassAglQ->Draw();

   p13Q->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassTOF_PQ->SetLineColor (2);
   MassTOF_PQ->SetLineWidth (2);
   MassTOF_DQ->SetLineColor (4);
   MassTOF_DQ->SetLineWidth (2);
   MassTOF_DQ->SetFillColor (4);
   MassTOF_PQ->SetFillColor (2);
   MassTOF_PQ->SetFillStyle (3001);
   MassTOF_DQ->SetFillStyle (3002);
   MassTOF_PQ->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassTOF_PQ->GetXaxis()->SetTitleSize (0.045);
   MassTOF_PQ->SetTitle ("Mass TOF MC");
   MassTOF_PQ->Draw();
   MassTOF_DQ->Draw ("same");

   p14Q->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassNaF_PQ->SetLineColor (2);
   MassNaF_PQ->SetLineWidth (2);
   MassNaF_DQ->SetLineColor (4);
   MassNaF_DQ->SetLineWidth (2);
   MassNaF_DQ->SetFillColor (4);
   MassNaF_PQ->SetFillColor (2);
   MassNaF_PQ->SetFillStyle (3001);
   MassNaF_DQ->SetFillStyle (3002);
   MassNaF_PQ->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassNaF_PQ->SetTitle ("Mass NaF MC");
   MassNaF_PQ->GetXaxis()->SetTitleSize (0.045);
   MassNaF_PQ->Draw();
   MassNaF_DQ->Draw ("same");

   p15Q->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassAgl_PQ->SetLineColor (2);
   MassAgl_PQ->SetLineWidth (2);
   MassAgl_DQ->SetLineColor (4);
   MassAgl_DQ->SetLineWidth (2);
   MassAgl_DQ->SetFillColor (4);
   MassAgl_PQ->SetFillColor (2);
   MassAgl_PQ->SetFillStyle (3001);
   MassAgl_DQ->SetFillStyle (3002);
   MassAgl_PQ->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassAgl_PQ->SetTitle ("Mass Agl MC");
   MassAgl_PQ->GetXaxis()->SetTitleSize (0.045);
   MassAgl_PQ->Draw();
   MassAgl_DQ->Draw ("same");

   cout<<"******************* Distance discr. plots ******************"<<endl;


   p16->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   DistTOF_P->SetLineColor (2);
   DistTOF_P->SetLineWidth (2);
   DistTOF_He->SetLineWidth (2);
   DistTOF_D->SetLineColor (4);
   DistTOF_He->SetLineColor (3);
   DistTOF_D->SetLineWidth (2);
   DistTOF_D->SetFillColor (4);
   DistTOF_P->SetFillColor (2);
   DistTOF_He->SetFillColor (3);
   DistTOF_P->SetFillStyle (3001);
   DistTOF_D->SetFillStyle (3002);
   DistTOF_He->SetFillStyle (3002);
   DistTOF_He->GetXaxis()->SetTitle ("Distance discriminant");
   DistTOF_He->GetXaxis()->SetTitleSize (0.045);
   DistTOF_He->SetTitle ("Distance discr. distribution  TOF MC");
   DistTOF_He->Draw();
   DistTOF_D->Draw ("same");
   DistTOF_P->Draw ("same");

   p17->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   DistNaF_P->SetLineColor (2);
   DistNaF_P->SetLineWidth (2);
   DistNaF_He->SetLineWidth (2);
   DistNaF_D->SetLineColor (4);
   DistNaF_He->SetLineColor (3);
   DistNaF_D->SetLineWidth (2);
   DistNaF_D->SetFillColor (4);
   DistNaF_P->SetFillColor (2);
   DistNaF_He->SetFillColor (3);
   DistNaF_P->SetFillStyle (3001);
   DistNaF_D->SetFillStyle (3002);
   DistNaF_He->SetFillStyle (3002);
   DistNaF_He->GetXaxis()->SetTitle ("Distance discriminant");
   DistNaF_He->GetXaxis()->SetTitleSize (0.045);
   DistNaF_He->SetTitle ("Distance discr. distribution NaF MC");
   DistNaF_He->Draw();
   DistNaF_D->Draw ("same");
   DistNaF_P->Draw ("same");

   p18->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   DistAgl_P->SetLineColor (2);
   DistAgl_P->SetLineWidth (2);
   DistAgl_He->SetLineWidth (2);
   DistAgl_D->SetLineColor (4);
   DistAgl_He->SetLineColor (3);
   DistAgl_D->SetLineWidth (2);
   DistAgl_D->SetFillColor (4);
   DistAgl_P->SetFillColor (2);
   DistAgl_He->SetFillColor (3);
   DistAgl_P->SetFillStyle (3001);
   DistAgl_D->SetFillStyle (3002);
   DistAgl_He->SetFillStyle (3002);
   DistAgl_He->GetXaxis()->SetTitle ("Distance discriminant");
   DistAgl_He->GetXaxis()->SetTitleSize (0.045);
   DistAgl_He->SetTitle ("Distance discr. distribution Agl MC");
   DistAgl_He->Draw();
   DistAgl_D->Draw ("same");
   DistAgl_P->Draw ("same");

   cout<<"******************* R vs Dist plots ******************"<<endl;
  /* p19->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   RvsDistTOF_P->SetMarkerColor (2);
   RvsDistTOF_D->SetMarkerColor (4);
   RvsDistTOF_D->SetTitle ("R vs Distance discr. TOF (MC)");
   RvsDistTOF_D->GetXaxis()->SetTitle ("R [GV]");
   RvsDistTOF_D->GetYaxis()->SetTitle ("Distance discr. TOF");
   RvsDistTOF_P->GetZaxis()->SetRangeUser (10,4000);
   RvsDistTOF_D->GetZaxis()->SetRangeUser (1,4000);
   RvsDistTOF_D->Draw();
   RvsDistTOF_P->Draw ("same");


   p20->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   RvsDistNaF_P->SetMarkerColor (2);
   RvsDistNaF_D->SetMarkerColor (4);
   RvsDistNaF_D->SetTitle ("R vs Distance discr. NaF (MC)");
   RvsDistNaF_D->GetXaxis()->SetTitle ("R [GV]");
   RvsDistNaF_D->GetYaxis()->SetTitle ("Distance discr. NaF");
   RvsDistNaF_P->GetZaxis()->SetRangeUser (5,20);
   RvsDistNaF_D->GetZaxis()->SetRangeUser (1,20);
   RvsDistNaF_D->Draw();
   RvsDistNaF_P->Draw ("same");

   p21->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   RvsDistAgl_P->SetMarkerColor (2);
   RvsDistAgl_D->SetMarkerColor (4);
   RvsDistAgl_D->SetTitle ("R vs Distance discr. Agl (MC)");
   RvsDistAgl_D->GetXaxis()->SetTitle ("R [GV]");
   RvsDistAgl_D->GetYaxis()->SetTitle ("Distance discr. Agl");
   RvsDistAgl_P->GetZaxis()->SetRangeUser (40,400);
   RvsDistAgl_D->GetZaxis()->SetRangeUser (1,400);
   RvsDistAgl_D->Draw();
   RvsDistAgl_P->Draw ("same");
   */


   p19->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   LikvsDistTOF_P->SetMarkerColor (2);
   LikvsDistTOF_D->SetMarkerColor (4);
   LikvsDistTOF_D->GetXaxis()->SetTitle ("Likelihood Discriminant (-log(1-L))");
   LikvsDistTOF_D->GetYaxis()->SetTitle ("Distance from D TOF");
   LikvsDistTOF_D -> SetContour(15);
   LikvsDistTOF_P->SetMarkerStyle (8);
   LikvsDistTOF_P->SetMarkerSize(0.15);
   LikvsDistTOF_P->GetXaxis()->SetRangeUser(0,2.6);	
   LikvsDistTOF_D -> SetLineWidth(3);
   LikvsDistTOF_D->Smooth(4);
   LikvsDistTOF_D->SetTitle("Distance vs Likelihood (TOF range)");
   	
   LikvsDistTOF_D->Draw ("CONT6");
   LikvsDistTOF_P->Draw("same");

   p20->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   LikvsDistNaF_P->SetMarkerColor (2);
   LikvsDistNaF_D->SetMarkerColor (4);
   LikvsDistNaF_D->SetTitle ("Distance vs Likelihood (NaF range)");
   LikvsDistNaF_D->GetXaxis()->SetTitle ("Likelihood Discriminant (-log(1-L))");
   LikvsDistNaF_D->GetYaxis()->SetTitle ("Distance from D NaF");
   LikvsDistNaF_P->SetMarkerStyle (8);
   LikvsDistNaF_P->SetMarkerSize(0.15);	
   LikvsDistNaF_D -> SetContour(15);
   LikvsDistNaF_D -> SetLineWidth(3);
   LikvsDistNaF_P->GetXaxis()->SetRangeUser(0,4);	
   LikvsDistNaF_D->Smooth(4);
   LikvsDistNaF_D->Draw ("CONT6");
   LikvsDistNaF_P->Draw("same");


   p21->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   LikvsDistAgl_P->SetMarkerColor (2);
   LikvsDistAgl_D->SetMarkerColor (4);
   LikvsDistAgl_D->SetTitle ("Distance vs Likelihood (Agl range)");
   LikvsDistAgl_D->GetXaxis()->SetTitle ("Likelihood Discriminat (-log(1-L))");
   LikvsDistAgl_D->GetYaxis()->SetTitle ("Distance from D Agl");
   LikvsDistAgl_P->SetMarkerStyle (8);
   LikvsDistAgl_P->SetMarkerSize(0.15);
   LikvsDistAgl_D -> SetContour(15);
   LikvsDistAgl_D -> SetLineWidth(3);
   LikvsDistAgl_P->GetXaxis()->SetRangeUser(0,5);	
   LikvsDistAgl_D->Smooth(4);
   LikvsDistAgl_D->Draw ("CONT6");
   LikvsDistAgl_P->Draw("same");





   //fileFinalPlots->Flush();
   //fileFinalPlots->Close();
   
   
   finalPlots.Add(p  );
   finalPlots.Add(p1  );
   finalPlots.Add(p2  );
   finalPlots.Add(p3  );
   finalPlots.Add(p4  );
   finalPlots.Add(p5  );
   finalPlots.Add(p6  );
   finalPlots.Add(p4_bis);
   finalPlots.Add(p5_bis);
   finalPlots.Add(p6_bis);
   finalPlots.Add(p7  );
   finalPlots.Add(p8  );
   finalPlots.Add(p9  );
   finalPlots.Add(p7_bis);
   finalPlots.Add(p8_bis);
   finalPlots.Add(p9_bis);
   finalPlots.Add(p10 );
   finalPlots.Add(p11 );
   finalPlots.Add(p12 );
   finalPlots.Add(p13_bis );
   finalPlots.Add(p14_bis );
   finalPlots.Add(p15_bis );
   finalPlots.Add(p13 );
   finalPlots.Add(p14 );
   finalPlots.Add(p15 );
   finalPlots.Add(p10Q);
   finalPlots.Add(p11Q);
   finalPlots.Add(p12Q);
   finalPlots.Add(p13Q);
   finalPlots.Add(p14Q);
   finalPlots.Add(p15Q);
   finalPlots.Add(p16 );
   finalPlots.Add(p17 );
   finalPlots.Add(p18 );
   finalPlots.Add(p19 );
   finalPlots.Add(p20 );
   finalPlots.Add(p21 );
   finalPlots.writeObjsInFolder("Common plots for slides");
   
}
