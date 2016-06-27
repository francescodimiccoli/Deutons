

void DeutonsTemplateFits_Plot(TemplateFIT *FitTOF_Dbins, 
                              TemplateFIT *   FitNaF_Dbins ,
                              TemplateFIT *   FitAgl_Dbins,
                                                             
                              TemplateFIT *   FitTOFgeo_Dbins,
                              TemplateFIT *   FitNaFgeo_Dbins,
                              TemplateFIT *   FitAglgeo_Dbins,
                                                           
                              TemplateFIT *   FitTOF_Pbins,
                              TemplateFIT *   FitNaF_Pbins ,
                              TemplateFIT *   FitAgl_Pbins,
			      string varname	
        ){
 


   TCanvas *c30_TOF[2][nbinsToF];
   TCanvas *c30_NaF[2][nbinsNaF];
   TCanvas *c30_Agl[2][nbinsAgl];

   TCanvas *c30_TOFgeo[11];
   TCanvas *c30_NaFgeo[11];
   TCanvas *c30_Aglgeo[11];

   //Primaries
   for(int bin=0; bin <nbinsToF; bin++) {
      c30_TOF[0][bin] = new TCanvas(("TOF bin:" + to_string(bin) + varname).c_str());
      c30_TOF[0][bin]->cd();
      FitTOF_Dbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin);
   }
   for(int bin=0; bin <nbinsNaF; bin++) {
      c30_NaF[0][bin] = new TCanvas(("NaF bin:" + to_string(bin) + varname).c_str());
      c30_NaF[0][bin]->cd();
      FitNaF_Dbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin);
   }
   for(int bin=0; bin <nbinsAgl; bin++) {
      c30_Agl[0][bin] = new TCanvas(("Agl bin:" + to_string(bin) + varname).c_str());
      c30_Agl[0][bin]->cd();
      FitAgl_Dbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin);
   }

   for(int bin=0; bin <nbinsToF; bin++) {
      c30_TOF[1][bin] = new TCanvas(("TOF P bin:" + to_string(bin) + varname).c_str());
      c30_TOF[1][bin]->cd();
      FitTOF_Pbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin);
   }
   for(int bin=0; bin <nbinsNaF; bin++) {
      c30_NaF[1][bin] = new TCanvas(("NaF P bin:" + to_string(bin) + varname).c_str());
      c30_NaF[1][bin]->cd();
      FitNaF_Pbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin);
   }
   for(int bin=0; bin <nbinsAgl; bin++) {
      c30_Agl[1][bin] = new TCanvas(("Agl P bin:" + to_string(bin) + varname).c_str());
      c30_Agl[1][bin]->cd();
      FitAgl_Pbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin);
   }
   //Geo. Zones
   for(int lat=1; lat<11; lat++) {


      c30_TOFgeo[lat] = new TCanvas(("TOF Geo. Zone:" + to_string(lat) + varname).c_str());
      c30_NaFgeo[lat] = new TCanvas(("NaF Geo. Zone:" + to_string(lat) + varname).c_str());
      c30_Aglgeo[lat] = new TCanvas(("Agl Geo. Zone:" + to_string(lat) + varname).c_str());

      c30_TOFgeo[lat] -> Divide(6,3);
      c30_NaFgeo[lat] -> Divide(6,3);
      c30_Aglgeo[lat] -> Divide(6,3);


      for(int bin=0; bin <nbinsToF; bin++) {
         c30_TOFgeo[lat] -> cd (bin +1);
         FitTOFgeo_Dbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin,lat);
      }
      for(int bin=0; bin <nbinsNaF; bin++) {
         c30_NaFgeo[lat] -> cd (bin +1);
         FitNaFgeo_Dbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin,lat);
      }
      for(int bin=0; bin <nbinsAgl; bin++) {
         c30_Aglgeo[lat] -> cd (bin +1);
         FitAglgeo_Dbins -> TemplateFitPlot(gPad,"Mass [GeV/C^2]",bin,lat);
      }
   }


   //Primaries


   for(int bin=0; bin <nbinsToF; bin++)
      finalPlots.Add(c30_TOF[0][bin]);
   finalPlots.writeObjsInFolder("Mass Template Fits/TOF/TOF Primaries/Dbins");

   for(int bin=0; bin <nbinsNaF; bin++)
      finalPlots.Add(c30_NaF[0][bin]);
   finalPlots.writeObjsInFolder("Mass Template Fits/NaF/NaF Primaries/Dbins");

   for(int bin=0; bin <nbinsAgl; bin++)
      finalPlots.Add(c30_Agl[0][bin]);
   finalPlots.writeObjsInFolder("Mass Template Fits/Agl/Agl Primaries/Dbins");

   for(int bin=0; bin <nbinsToF; bin++)
      finalPlots.Add(c30_TOF[1][bin]);
   finalPlots.writeObjsInFolder("Mass Template Fits/TOF/TOF Primaries/Pbins");

   for(int bin=0; bin <nbinsNaF; bin++)
      finalPlots.Add(c30_NaF[1][bin]);
   finalPlots.writeObjsInFolder("Mass Template Fits/NaF/NaF Primaries/Pbins");

   for(int bin=0; bin <nbinsAgl; bin++)
      finalPlots.Add(c30_Agl[1][bin]);
   finalPlots.writeObjsInFolder("Mass Template Fits/Agl/Agl Primaries/Pbins");
   //Geom. Zones
   for(int lat=1; lat<11; lat++)
      finalPlots.Add(c30_TOFgeo[lat]);
   finalPlots.writeObjsInFolder("Mass Template Fits/TOF/TOF Geom. Zones");

   for(int lat=1; lat<11; lat++)
      finalPlots.Add(c30_NaFgeo[lat]);
   finalPlots.writeObjsInFolder("Mass Template Fits/NaF/NaF Geom. Zones ");

   for(int lat=1; lat<11; lat++)
      finalPlots.Add(c30_Aglgeo[lat]);
   finalPlots.writeObjsInFolder("Mass Template Fits/Agl/Agl Geom. Zones  ");


return;

}
