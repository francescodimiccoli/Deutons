
void   MCControlsamplecuteff_Plot(	
    
    	   TH1*  EffCSCMCP_TH1F   	  
    	){




   TCanvas *c5	=new TCanvas("Control sample cuts Efficiency (R bins)");
   
   c5->cd(2);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph *EffMCCSC= new TGraph();
   for(int i=0; i<nbinsr; i++) EffMCCSC->SetPoint(i,PRB.RigBinCent(i), EffCSCMCP_TH1F->GetBinContent(i+1));
   EffMCCSC->SetMarkerColor(2);
   EffMCCSC->SetMarkerStyle(8);
   EffMCCSC->SetLineColor(2);
   EffMCCSC->SetLineWidth(2);
   EffMCCSC->SetTitle("Control sample cuts Efficiency (R bins)");
   EffMCCSC->GetXaxis()->SetTitle("R [GV]");
   EffMCCSC->GetYaxis()->SetTitle("Efficiency");
   EffMCCSC->GetXaxis()->SetTitleSize(0.045);
   EffMCCSC->GetYaxis()->SetTitleSize(0.045);
   EffMCCSC->Draw("APC"); 

   finalPlots.Add(c5);
   
   finalPlots.writeObjsInFolder("MC Results/Quality");

   return;	
}
