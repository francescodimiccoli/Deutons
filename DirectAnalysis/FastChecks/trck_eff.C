const int Z = 1;
const double M[2] = {0.938,0.932*4}; 
const double K[2] = {2,4};
const double tof_charge_sigma[2] = {0.3,0.3};
const double trd_charge_sigma[2] = {0.3,0.3};
const double trk_charge_inn_left[2] = {0.3,0.3};
const double trk_charge_inn_right[2] = {0.5,0.5}; 
const double trk_charge_ul1_left[2] = {0.3,0.3};
const double trk_charge_ul1_right[2] = {0.5,0.5};

void trck_eff() {

  // gROOT->ProcessLine(".L ../build/libntp.so"); // to not to have complaints
  // ROOT::EnableImplicitMT(); 
  ROOT::RDataFrame df("Compact","/eos/ams/group/dbar/release_v7/e0_vdev_200421/neg/ISS.B1130/pass7/1520601133.root");
  // ROOT::RDataFrame df("Compact","/eos/ams/group/dbar/release_v7/e0_vdev_200421/neg/ISS.B1130/pass7/1407*.root"); 
  std::vector<std::string> colNames = df.GetColumnNames();
  for(int i=0; i<colNames.size(); i++) {
    std::cout << colNames[i] << std::endl;
  }

  // define more variables (and arrays as values ... very heavy, is there any other solution?)
  auto df_define = df
    .Define("sa_tof_q_up"    ,[](ROOT::VecOps::RVec<float>& sa_tof_q_lay){ return 0.5*(sa_tof_q_lay[0]+sa_tof_q_lay[1]); },{"sa_tof_q_lay[4]"})
    .Define("sa_tof_q_dw"    ,[](ROOT::VecOps::RVec<float>& sa_tof_q_lay){ return 0.5*(sa_tof_q_lay[2]+sa_tof_q_lay[3]); },{"sa_tof_q_lay[4]"})
    .Define("sa_tof_rig"     ,[](float sa_tof_beta) { return M[Z-1]*sa_tof_beta/sqrt(1-sa_tof_beta*sa_tof_beta)/Z; },{"sa_tof_beta"})
    .Define("log_sa_tof_rig" ,"log10(sa_tof_rig)")
    .Define("sa_tof_clsn"    ,[](ROOT::VecOps::RVec<short int>& sa_tof_clsn){ return sa_tof_clsn[0]+sa_tof_clsn[2]; },{"sa_tof_clsn[4]"})
    .Define("sa_exthit_dl1_x",[](ROOT::VecOps::RVec<float>& sa_exthit_dl1){ return sa_exthit_dl1[0];},{"sa_exthit_dl1[2]"})
    .Define("sa_exthit_dl1_y",[](ROOT::VecOps::RVec<float>& sa_exthit_dl1){ return sa_exthit_dl1[1];},{"sa_exthit_dl1[2]"})
    .Define("log_sa_ecal_rig",[](float sa_ecal_edepd) { return log10(sqrt(pow(K[Z-1]*sa_ecal_edepd,2)-pow(M[Z-1],2))); },{"sa_ecal_edepd"})
    .Define("trk_ny_inn"     ,[](short int trk_patty) { int ny=0; for (int i=0;i<7;i++) ny+=((trk_patty>>(i+1))&1); return ny; },{"trk_patty"})
    .Define("trk_nxy_inn"    ,[](short int trk_pattxy) { int nxy=0; for (int i=0;i<7;i++) nxy+=((trk_pattxy>>(i+1))&1); return nxy; },{"trk_pattxy"})
    .Define("trk_chi2y_inn"  ,[](ROOT::VecOps::RVec<float>& trk_chisqn){ return trk_chisqn[3*2+1]; },{"trk_chisqn[5][2]"})
    .Define("trk_rig_inn"    ,[](ROOT::VecOps::RVec<float>& trk_rig){ return trk_rig[3]; },{"trk_rig[5]"})
    .Define("log_trk_rig_inn","log10(trk_rig_inn)");

  // tracking efficiency denominator
  auto cut_sa_tof_beta_ncl = [](short int sa_tof_beta_ncl) { return sa_tof_beta_ncl==4; };
  auto cut_sa_tof_chisqtn = [](float sa_tof_chisqtn) { return sa_tof_chisqtn < 8; };
  auto cut_sa_tof_clsn = [](int sa_tof_clsn) { return sa_tof_clsn<=4; };
  auto cut_sa_tof_charge = [](double sa_tof_charge) { return fabs(sa_tof_charge-Z)<tof_charge_sigma[Z-1]; }; 
  auto cut_sa_trd_chi2 = [](float sa_trd_chi2) { return sa_trd_chi2<10; };
  auto cut_sa_trd_fiducial = [](short int sa_trd_fiducial) { return (sa_trd_fiducial&0xff)==0xff; }; 
  auto cut_sa_trd_charge = [](float sa_trd_q) { return fabs(sa_trd_q-Z)<trd_charge_sigma[Z-1]; };
  auto cut_sa_exthit_ql1 = [](float sa_exthit_ql1){ return (sa_exthit_ql1>Z-trk_charge_ul1_left[Z-1])&&(sa_exthit_ql1<Z+trk_charge_ul1_right[Z-1]); }; 
  auto cut_sa_exthit_dl1 = [](float x, float y){ return (fabs(x)<1)&&(fabs(y)<1); };
  auto cut_sa_exthit_sl1 = [](short int sa_exthit_status_l1){ return (sa_exthit_status_l1&0x10013D)==0; }; 

  // not checking this ... but in principle we should
  auto df_trk_base = df_define
    .Filter(cut_sa_tof_beta_ncl, {"sa_tof_beta_ncl"}, "cut_sa_tof_beta_ncl")
    .Filter(cut_sa_tof_chisqtn,  {"sa_tof_chisqtn"},  "cut_sa_tof_chisqtn")
    .Filter(cut_sa_tof_clsn,     {"sa_tof_clsn"},     "cut_sa_tof_clsn")
    .Filter(cut_sa_trd_chi2,     {"sa_trd_chi2"},     "cut_sa_trd_chi2")
    .Filter(cut_sa_trd_fiducial, {"sa_trd_fiducial"}, "cut_sa_trd_fiducial");

  // check charge cuts
  auto h_sa_tof_q_up = df_trk_base.Histo1D({"h_sa_tof_q_up","; Q_{TOF,up}",200,0,Z+3},"sa_tof_q_up");
  auto h_sa_tof_q_dw = df_trk_base.Histo1D({"h_sa_tof_q_dw","; Q_{TOF,up}",200,0,Z+3},"sa_tof_q_dw");
  auto h_sa_trd_q    = df_trk_base.Histo1D({"h_sa_trd_q"   ,"; Q_{TRD}",200,0,Z+3}   ,"sa_trd_q");

  // add charge selection 
  auto df_trk_base_more = df_trk_base
    .Filter(cut_sa_tof_charge, {"sa_tof_q_up"}, "cut_sa_tof_charge_up")
    .Filter(cut_sa_tof_charge, {"sa_tof_q_dw"}, "cut_sa_tof_charge_dw")
    .Filter(cut_sa_trd_charge, {"sa_trd_q"}, "cut_sa_trd_charge");

  // check L1 spatial and charge distributions
  auto h_sa_exthit_l1_dx_vs_dy = df_trk_base
    .Filter(cut_sa_exthit_ql1, {"sa_exthit_ql1"})
    .Filter(cut_sa_exthit_sl1, {"sa_exthit_status_l1"})
    .Histo2D({"h_sa_exthit_l1_dx_vs_dy","; #Deltax [cm]; #Deltay [cm]",200,-10,10,200,-10,10},"sa_exthit_dl1_x","sa_exthit_dl1_y");
  auto h_sa_exthit_ql1 = df_trk_base
    .Filter(cut_sa_exthit_dl1, {"sa_exthit_dl1_x","sa_exthit_dl1_y"})
    .Filter(cut_sa_exthit_sl1, {"sa_exthit_status_l1"})
    .Histo1D({"h_sa_exthit_ql1","; Unbiased L1 Charge",200,0,Z+5}, "sa_exthit_ql1"); 

  auto df_trk_den = df_trk_base_more
    .Filter(cut_sa_exthit_ql1, {"sa_exthit_ql1"}, "cut_sa_exthit_ql1")
    .Filter(cut_sa_exthit_dl1, {"sa_exthit_dl1_x","sa_exthit_dl1_y"}, "cut_sa_exthit_dl1")
    .Filter(cut_sa_exthit_sl1, {"sa_exthit_status_l1"}, "cut_sa_exthit_sl1");

  df_trk_den.Report()->Print();

  auto cut_trk_patty_inn = [](short int trk_patty) { return ((trk_patty&0x2)!=0)&&((trk_patty&0xc)!=0)&&((trk_patty&0x30)!=0)&&((trk_patty&0xc0)!=0); }; 
  auto cut_trk_nhity_inn = [](int trk_ny_inn) { return trk_ny_inn>4; };
  auto cut_trk_nhitxy_inn = [](int trk_nxy_inn) { return trk_nxy_inn>3; };
  auto cut_trk_chi2 = [](float trk_chi2) { return (trk_chi2<10)&&(trk_chi2>0); };
  auto cut_trk_q_inn = [](float trk_q_inn) { return (trk_q_inn>Z-trk_charge_inn_left[Z-1])&&(trk_q_inn<Z+trk_charge_inn_right[Z-1]); }; 
  
  auto df_trk_num = df_trk_den 
    .Filter(cut_trk_patty_inn,  {"trk_patty"},     "cut_trk_patty_inn")
    .Filter(cut_trk_nhity_inn,  {"trk_ny_inn"},    "cut_trk_nhity_inn")
    .Filter(cut_trk_nhitxy_inn, {"trk_nxy_inn"},   "cut_trk_nhity_inn")
    .Filter(cut_trk_chi2,       {"trk_chi2y_inn"}, "cut_trk_chi2_inn")
    .Filter(cut_trk_q_inn,      {"trk_q_inn"},     "cut_trk_q_inn");

  df_trk_num.Report()->Print();

  // check alternative rigidity measurements (cutoff not available!!!) 
  auto h_sa_tof_rig_vs_rig = df_trk_num
    .Filter("(trk_rig_inn>0)&&(sa_tof_rig>0)") 
    .Histo2D({"h_sa_tof_rig_vs_rig","; log_{10}(R/GV); log_{10}(R_{TOF}/GV}",200,-1,3,200,-1,3},"log_trk_rig_inn","log_sa_tof_rig");

  auto h_sa_ecal_rig_vs_rig = df_trk_num
    .Filter("(trk_rig_inn>0)&&(sa_ecal_edepd>0)")
    .Histo2D({"h_sa_ecal_rig_vs_rig","; log_{10}(R/GV); log_{10}(R_{ECAL}/GV}",200,-1,3,200,-1,3},"log_trk_rig_inn","log_sa_ecal_rig");

  auto h_den = df_trk_den.Filter("(sa_tof_rig>0)").Histo1D({"h_den","; log_{10}(R_{TOF}/GV}",200,-1,3},"log_sa_tof_rig");
  auto h_num = df_trk_num.Filter("(sa_tof_rig>0)").Histo1D({"h_num","; log_{10}(R_{TOF}/GV}",200,-1,3},"log_sa_tof_rig");

  TCanvas canvas1;
  h_sa_tof_q_up->Draw("HIST"); 
  h_sa_tof_q_dw->Draw("HIST SAME");
  h_sa_trd_q->Draw("HIST SAME");
  h_sa_exthit_ql1->Draw("HIST SAME");
  canvas1.Update();

  TCanvas canvas2;
  h_sa_tof_rig_vs_rig->Draw("COLZ");
  canvas2.Update();

  TCanvas canvas3;
  h_sa_exthit_l1_dx_vs_dy->Draw("COLZ");

  TCanvas canvas4;
  h_sa_ecal_rig_vs_rig->Draw("COLZ");

  TCanvas canvas5;
  h_num->Divide(h_den.GetPtr());  
  h_num->Draw("HIST E");

  TCanvas canvas;
  canvas.WaitPrimitive();
}
