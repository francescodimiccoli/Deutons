#include <bitset>
#include <string>

#include "Classifier.h"

// D) 0:My training; 1:Francesca's training.
// C) 0:TOF region; 1:NaF region; 2:Aerogel region.
// B) 0:TOF; 1:TRD; 2:Traker; 3:RICH.
// A) 0, ...,7 different trainings 
// classifier_type = DCBA
const int nclassifier = 31;
const int classifier_index[nclassifier] = {
    0,  1,  2,  3,  4, // tof_tof
   20, 21, 22, 23, 24, // tof_trk 
  120,121,122,123,124, // naf_trk
  130,131,132,133,134, // naf_rich
  220,221,222,223,224, // agl_trk
  230,231,232,234,235, // agl_rich
  1231 // fg 
};
const char* classifier_name[nclassifier][3] = {
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_tof_tof1_trai00_test00_Z1_N500000" ,"BDTG","Multi"}, 
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_tof_tof1_trai01_test00_Z1_N500000" ,"BDTG","Multi"}, 
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_tof_tof1_trai02_test00_Z1_N500000" ,"BDTG","Multi"},
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_tof_tof1_trai03_test00_Z1_N500000" ,"BDTG","Multi"},
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_tof_tof1_trai04_test00_Z1_N500000" ,"BDTG","Multi"},

  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_tof_trk1_trai00_test00_Z1_N500000" ,"BDTG","Multi"}, 
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_tof_trk1_trai01_test00_Z1_N500000" ,"BDTG","Multi"}, 
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_tof_trk1_trai02_test00_Z1_N500000" ,"BDTG","Multi"},
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_tof_trk1_trai03_test00_Z1_N500000" ,"BDTG","Multi"},
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_tof_trk1_trai04_test00_Z1_N500000" ,"BDTG","Multi"},
  
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_naf_trk1_trai00_test00_Z1_N500000" ,"BDTG","Multi"}, 
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_naf_trk1_trai01_test00_Z1_N500000" ,"BDTG","Multi"}, 
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_naf_trk1_trai02_test00_Z1_N500000" ,"BDTG","Multi"},
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_naf_trk1_trai03_test00_Z1_N500000" ,"BDTG","Multi"},
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_naf_trk1_trai04_test00_Z1_N500000" ,"BDTG","Multi"},

  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_naf_rich1_trai00_test00_Z1_N500000","BDTG","Multi"}, 
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_naf_rich1_trai01_test00_Z1_N500000","BDTG","Multi"}, 
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_naf_rich1_trai02_test00_Z1_N500000","BDTG","Multi"},
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_naf_rich1_trai03_test00_Z1_N500000","BDTG","Multi"},
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_naf_rich1_trai04_test00_Z1_N500000","BDTG","Multi"},

  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_agl_trk1_trai00_test00_Z1_N500000" ,"BDTG","Multi"}, 
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_agl_trk1_trai01_test00_Z1_N500000" ,"BDTG","Multi"}, 
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_agl_trk1_trai02_test00_Z1_N500000" ,"BDTG","Multi"},
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_agl_trk1_trai03_test00_Z1_N500000" ,"BDTG","Multi"},
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_agl_trk1_trai04_test00_Z1_N500000" ,"BDTG","Multi"},

  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_agl_rich1_trai00_test00_Z1_N500000","BDTG","Multi"}, 
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_agl_rich1_trai01_test00_Z1_N500000","BDTG","Multi"}, 
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_agl_rich1_trai02_test00_Z1_N500000","BDTG","Multi"},
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_agl_rich1_trai03_test00_Z1_N500000","BDTG","Multi"},
  {"/eos/ams/group/dbar/data/tmva_190415/TMVA_agl_rich1_trai04_test00_Z1_N500000","BDTG","Multi"},

  {"/afs/cern.ch/work/f/fgiovacc/public/perAlberto/weights_TMVAr80_2000EvtTrain","BDT","Single"}
};

ClassifierManager::ClassifierManager() : initialized(false) {
  data = new ClassifierData();
  pdf2db_transform = 0;
  pdf2db_mass = 0;
  initialized = false;
  current_run = -1;
  current_event = -1;
}

ClassifierManager::~ClassifierManager(){
  delete data;
  readers.clear(); // to-be-done: do it correctly 
  initialized = false;
}

void ClassifierManager::Init(){
  if (initialized) return;
  TMVA::Tools::Instance();
  TString opt = "V:Color:!Silent"; //:CreateMVAPdfs"; //TString opt = "Silent";
  for (int iclassifier=0; iclassifier<nclassifier; iclassifier++) {
    string classifier_file_name; 
    if (classifier_index[iclassifier]==1231) classifier_file_name.append("/afs/cern.ch/work/f/fgiovacc/public/perAlberto/weights_TMVAr80_2000EvtTrain/TMVAClassification_BDT.weights.xml");
    else                                     classifier_file_name.append(Form("%s/TMVAClassification_%s.weights.xml",classifier_name[iclassifier][0],classifier_name[iclassifier][1]));
    TFile* file = TFile::Open(classifier_file_name.c_str(),"read");
    if ((!file)||(file->IsZombie())) continue;
    file->Close(); 
    std::cout << "ClassifierManager::Init-read weight file " << classifier_file_name.c_str() << " for classifier with index " << classifier_index[iclassifier] << std::endl;
    readers[classifier_index[iclassifier]] = new TMVA::Reader(opt); 
    SetReaderVariables(readers[classifier_index[iclassifier]],classifier_index[iclassifier],(classifier_index[iclassifier]>1000));
    names[classifier_index[iclassifier]].append(classifier_name[iclassifier][1]); 
    multi[classifier_index[iclassifier]].append(classifier_name[iclassifier][2]);
    readers[classifier_index[iclassifier]]->BookMVA(classifier_name[iclassifier][1],classifier_file_name);
  }
  // trasformation 
  TFile* file_transform = TFile::Open("root://eosams.cern.ch//eos/ams/group/dbar/data/likelihood/release_v5_prod_10_pdf2db.root","read");
  pdf2db_transform = (PDF2DB*) file_transform->Get("pdf2db");
  file_transform->Close();
  const char* det_name[3] = {"tof","naf","agl"};
  data->_tof_logchisqtn             .SetName("tof_logchisqtn");
  data->_tof_logchisqcn             .SetName("tof_logchisqcn");
  data->_tof_logzprob               .SetName("tof_logzprob"  );
  data->_tof_logqasym               .SetName("tof_logqasym"  );
  data->_tof_logdbeta               .SetName("tof_logdbeta"  );
  data->_tof_logisolat              .SetName("tof_logisolat" );
  data->_trd_nhit                   .SetName("trd_nhit"      );
  data->_trd_vertex                 .SetName("trd_vertex"    );
  data->_trd_logphe                 .SetName("trd_logphe"    );
  data->_trd_logep                  .SetName("trd_logep"     );
  data->_trk_nyhits                 .SetName("trk_nyhits"     );
  data->_trk_nxhits                 .SetName("trk_nxhits"     );
  data->_trk_n                      .SetName("trk_n"          );
  data->_trk_logchisqkax            .SetName("trk_logchisqkax");
  data->_trk_logchisqkay            .SetName("trk_logchisqkay");
  data->_trk_logchisqchx            .SetName("trk_logchisqchx");
  data->_trk_logchisqchy            .SetName("trk_logchisqchy");
  data->_trk_logscatt34x            .SetName("trk_logscatt34x");
  data->_trk_logscatt34y            .SetName("trk_logscatt34y");
  data->_trk_logscatt56x            .SetName("trk_logscatt56x");
  data->_trk_logscatt56y            .SetName("trk_logscatt56y");
  data->_trk_logdinvru              .SetName("trk_logdinvru"  );
  data->_trk_logdinvrd              .SetName("trk_logdinvrd"  );
  data->_trk_logcc                  .SetName("trk_logcc"      );
  data->_trk_logdinvr               .SetName("trk_logdinvr"   );
  data->_trk_logdinvrms             .SetName("trk_logdinvrms" );
  data->_trk_logdinvrck             .SetName("trk_logdinvrck" );
  data->_trk_r                      .SetName("trk_r"          );
  data->_trk_q                      .SetName("trk_q"          );
  data->_trk_logqyasym              .SetName("trk_logqyasym"  );
  data->_trk_logedep_l2_offtrack    .SetName("trk_logl2off"   );
  data->_trk_logedep_inn_offtrack[0].SetName("trk_loginnoff0" );
  data->_trk_logedep_inn_offtrack[1].SetName("trk_loginnoff1" );
  data->_rich_npmt                  .SetName("rich_npmt"        );
  data->_rich_nhit_uncorr           .SetName("rich_nhit"        );
  data->_rich_nhit                  .SetName("rich_correct"     );
  data->_rich_nhit_refl             .SetName("rich_reflect"     );
  data->_rich_nhit_hyp[0]           .SetName("rich_hyp_dir"     );
  data->_rich_nhit_hyp[1]           .SetName("rich_hyp_ref"     );
  data->_rich_nhit_notused[0]       .SetName("rich_notused1"    );
  data->_rich_nhit_notused[1]       .SetName("rich_notused2"    );
  data->_rich_nclus                 .SetName("rich_nclus"       );
  data->_rich_tot_hit[0][0]         .SetName("rich_nhit_00"     );
  data->_rich_tot_hit[0][1]         .SetName("rich_nhit_01"     );
  data->_rich_tot_hit[0][2]         .SetName("rich_nhit_02"     );
  data->_rich_tot_hit[0][3]         .SetName("rich_nhit_03"     );
  data->_rich_tot_hit[0][4]         .SetName("rich_nhit_04"     );
  data->_rich_tot_hit[1][0]         .SetName("rich_nhit_10"     );
  data->_rich_tot_hit[1][1]         .SetName("rich_nhit_11"     );
  data->_rich_tot_hit[1][2]         .SetName("rich_nhit_12"     );
  data->_rich_tot_hit[1][3]         .SetName("rich_nhit_13"     );
  data->_rich_tot_hit[1][4]         .SetName("rich_nhit_14"     );
  data->_rich_np_exp                .SetName("rich_np_exp"      );
  data->_rich_exp_res               .SetName("rich_exp_res"     );
  data->_rich_exp_rms               .SetName("rich_exp_rms"     );
  data->_rich_lognp_min             .SetName("rich_lognp_min"   );
  data->_rich_logdbeta_lip          .SetName("rich_logdbeta_lip");
  data->_rich_logqpmtcons           .SetName("rich_logqpmtcons" );
  data->_rich_q                     .SetName("rich_q"           );
  data->_rich_prob                  .SetName("rich_prob"        );
  data->_rich_logdq_lip             .SetName("rich_logdq_lip"   );
  data->_rich_ratio                 .SetName("rich_ratio"       );
  data->_rich_logdbeta              .SetName("rich_logdbeta"    );
  data->_rich_logdist               .SetName("rich_logdist"     );
  for (int id=0; id<3; id++) { 
    // TOF 
    data->_tof_logchisqtn             .AddTransformation(id,pdf2db_transform->Get(Form("tof_logchisqtn_vs_logkrec_%s",det_name[id])),&(data->_logk[id]),2);
    data->_tof_logchisqcn             .AddTransformation(id,pdf2db_transform->Get(Form("tof_logchisqcn_vs_logkrec_%s",det_name[id])),&(data->_logk[id]),2);
    data->_tof_logzprob               .AddTransformation(id,pdf2db_transform->Get(Form("tof_logzprob_vs_logkrec_%s"  ,det_name[id])),&(data->_logk[id]),2);
    data->_tof_logqasym               .AddTransformation(id,pdf2db_transform->Get(Form("tof_logqasym_vs_logkrec_%s"  ,det_name[id])),&(data->_logk[id]),2);
    data->_tof_logdbeta               .AddTransformation(id,pdf2db_transform->Get(Form("tof_logdbeta_vs_logkrec_%s"  ,det_name[id])),&(data->_logk[id]),2);
    data->_tof_logisolat              .AddTransformation(id,pdf2db_transform->Get(Form("tof_logisolat_vs_logkrec_%s" ,det_name[id])),&(data->_logk[id]),0);
    // TRD
    data->_trd_nhit                   .AddTransformation(id,pdf2db_transform->Get(Form("trd_nhit_vs_logkrec_%s"      ,det_name[id])),&(data->_logk[id]),2);
    data->_trd_vertex                 .AddTransformation(id,pdf2db_transform->Get(Form("trd_vertex_vs_logkrec_%s"    ,det_name[id])),&(data->_logk[id]),0);
    data->_trd_logphe                 .AddTransformation(id,pdf2db_transform->Get(Form("trd_logphe_vs_logkrec_%s"    ,det_name[id])),&(data->_logk[id]),2);
    data->_trd_logep                  .AddTransformation(id,pdf2db_transform->Get(Form("trd_logep_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),2);
    // Tracker
    data->_trk_nyhits                 .AddTransformation(id,pdf2db_transform->Get(Form("trk_nyhits_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),0);
    data->_trk_nxhits                 .AddTransformation(id,pdf2db_transform->Get(Form("trk_nxhits_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),0);
    data->_trk_n                      .AddTransformation(id,pdf2db_transform->Get(Form("trk_n_vs_logkrec_%s"          ,det_name[id])),&(data->_logk[id]),1);
    data->_trk_logchisqkax            .AddTransformation(id,pdf2db_transform->Get(Form("trk_logchisqkax_vs_logrrec_%s",det_name[id])),&(data->_logr[id]),2);
    data->_trk_logchisqkay            .AddTransformation(id,pdf2db_transform->Get(Form("trk_logchisqkay_vs_logrrec_%s",det_name[id])),&(data->_logr[id]),2);
    data->_trk_logchisqchx            .AddTransformation(id,pdf2db_transform->Get(Form("trk_logchisqchx_vs_logrrec_%s",det_name[id])),&(data->_logr[id]),2);
    data->_trk_logchisqchy            .AddTransformation(id,pdf2db_transform->Get(Form("trk_logchisqchy_vs_logrrec_%s",det_name[id])),&(data->_logr[id]),2);
    data->_trk_logscatt34x            .AddTransformation(id,pdf2db_transform->Get(Form("trk_logscatt34x_vs_logrrec_%s",det_name[id])),&(data->_logr[id]),2);
    data->_trk_logscatt34y            .AddTransformation(id,pdf2db_transform->Get(Form("trk_logscatt34y_vs_logrrec_%s",det_name[id])),&(data->_logr[id]),2);
    data->_trk_logscatt56x            .AddTransformation(id,pdf2db_transform->Get(Form("trk_logscatt56x_vs_logrrec_%s",det_name[id])),&(data->_logr[id]),2);
    data->_trk_logscatt56y            .AddTransformation(id,pdf2db_transform->Get(Form("trk_logscatt56y_vs_logrrec_%s",det_name[id])),&(data->_logr[id]),2);
    data->_trk_logdinvru              .AddTransformation(id,pdf2db_transform->Get(Form("trk_logdinvru_vs_logrrec_%s"  ,det_name[id])),&(data->_logr[id]),2);
    data->_trk_logdinvrd              .AddTransformation(id,pdf2db_transform->Get(Form("trk_logdinvrd_vs_logrrec_%s"  ,det_name[id])),&(data->_logr[id]),2);
    data->_trk_logcc                  .AddTransformation(id,pdf2db_transform->Get(Form("trk_logcc_vs_logrrec_%s"      ,det_name[id])),&(data->_logr[id]),2);
    data->_trk_logdinvr               .AddTransformation(id,pdf2db_transform->Get(Form("trk_logdinvr_vs_logrrec_%s"   ,det_name[id])),&(data->_logr[id]),2);
    data->_trk_logdinvrms             .AddTransformation(id,pdf2db_transform->Get(Form("trk_logdinvrms_vs_logrrec_%s" ,det_name[id])),&(data->_logr[id]),2);
    data->_trk_logdinvrck             .AddTransformation(id,pdf2db_transform->Get(Form("trk_logdinvrck_vs_logrrec_%s" ,det_name[id])),&(data->_logr[id]),2);
    data->_trk_r                      .AddTransformation(id,pdf2db_transform->Get(Form("trk_r_vs_logkrec_%s"          ,det_name[id])),&(data->_logk[id]),2);
    data->_trk_q                      .AddTransformation(id,pdf2db_transform->Get(Form("trk_q_vs_logkrec_%s"          ,det_name[id])),&(data->_logk[id]),1);
    data->_trk_logqyasym              .AddTransformation(id,pdf2db_transform->Get(Form("trk_logqyasym_vs_logkrec_%s"  ,det_name[id])),&(data->_logk[id]),1);
    data->_trk_logedep_l2_offtrack    .AddTransformation(id,pdf2db_transform->Get(Form("trk_logl2off_vs_logkrec_%s"   ,det_name[id])),&(data->_logk[id]),0);
    data->_trk_logedep_inn_offtrack[0].AddTransformation(id,pdf2db_transform->Get(Form("trk_loginnoff0_vs_logkrec_%s" ,det_name[id])),&(data->_logk[id]),0);
    data->_trk_logedep_inn_offtrack[1].AddTransformation(id,pdf2db_transform->Get(Form("trk_loginnoff1_vs_logkrec_%s" ,det_name[id])),&(data->_logk[id]),0);
    // RICH
    data->_rich_npmt                  .AddTransformation(id,pdf2db_transform->Get(Form("rich_npmt_vs_logkrec_%s"        ,det_name[id])),&(data->_logk[id]),1);
    data->_rich_nhit_uncorr           .AddTransformation(id,pdf2db_transform->Get(Form("rich_nhit_vs_logkrec_%s"        ,det_name[id])),&(data->_logk[id]),1);
    data->_rich_nhit                  .AddTransformation(id,pdf2db_transform->Get(Form("rich_correct_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),1);
    data->_rich_nhit_refl             .AddTransformation(id,pdf2db_transform->Get(Form("rich_reflect_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),0);
    data->_rich_nhit_hyp[0]           .AddTransformation(id,pdf2db_transform->Get(Form("rich_hyp_dir_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),0);
    data->_rich_nhit_hyp[1]           .AddTransformation(id,pdf2db_transform->Get(Form("rich_hyp_ref_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),0);
    data->_rich_nhit_notused[0]       .AddTransformation(id,pdf2db_transform->Get(Form("rich_notused1_vs_logkrec_%s"    ,det_name[id])),&(data->_logk[id]),0);
    data->_rich_nhit_notused[1]       .AddTransformation(id,pdf2db_transform->Get(Form("rich_notused2_vs_logkrec_%s"    ,det_name[id])),&(data->_logk[id]),0);
    data->_rich_nclus                 .AddTransformation(id,pdf2db_transform->Get(Form("rich_nclus_vs_logkrec_%s"       ,det_name[id])),&(data->_logk[id]),0);
    data->_rich_tot_hit[0][0]         .AddTransformation(id,pdf2db_transform->Get(Form("rich_nhit_00_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),0);
    data->_rich_tot_hit[0][1]         .AddTransformation(id,pdf2db_transform->Get(Form("rich_nhit_01_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),0);
    data->_rich_tot_hit[0][2]         .AddTransformation(id,pdf2db_transform->Get(Form("rich_nhit_02_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),0);
  //data->_rich_tot_hit[0][3]         .AddTransformation(id,pdf2db_transform->Get(Form("rich_nhit_03_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),0);
    data->_rich_tot_hit[0][4]         .AddTransformation(id,pdf2db_transform->Get(Form("rich_nhit_04_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),1);
    data->_rich_tot_hit[1][0]         .AddTransformation(id,pdf2db_transform->Get(Form("rich_nhit_10_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),0);
    data->_rich_tot_hit[1][1]         .AddTransformation(id,pdf2db_transform->Get(Form("rich_nhit_11_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),0);
    data->_rich_tot_hit[1][2]         .AddTransformation(id,pdf2db_transform->Get(Form("rich_nhit_12_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),0);
    data->_rich_tot_hit[1][3]         .AddTransformation(id,pdf2db_transform->Get(Form("rich_nhit_13_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),0);
    data->_rich_tot_hit[1][4]         .AddTransformation(id,pdf2db_transform->Get(Form("rich_nhit_14_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),1);
    data->_rich_np_exp                .AddTransformation(id,pdf2db_transform->Get(Form("rich_np_exp_vs_logkrec_%s"      ,det_name[id])),&(data->_logk[id]),2);
    data->_rich_exp_res               .AddTransformation(id,pdf2db_transform->Get(Form("rich_exp_res_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),2);
    data->_rich_exp_rms               .AddTransformation(id,pdf2db_transform->Get(Form("rich_exp_rms_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),2);
    data->_rich_lognp_min             .AddTransformation(id,pdf2db_transform->Get(Form("rich_lognp_min_vs_logkrec_%s"   ,det_name[id])),&(data->_logk[id]),2);
    data->_rich_logdbeta_lip          .AddTransformation(id,pdf2db_transform->Get(Form("rich_logdbeta_lip_vs_logkrec_%s",det_name[id])),&(data->_logk[id]),1);
    data->_rich_logqpmtcons           .AddTransformation(id,pdf2db_transform->Get(Form("rich_logqpmtcons_vs_logkrec_%s" ,det_name[id])),&(data->_logk[id]),2);
    data->_rich_q                     .AddTransformation(id,pdf2db_transform->Get(Form("rich_q_vs_logkrec_%s"           ,det_name[id])),&(data->_logk[id]),2);
    data->_rich_prob                  .AddTransformation(id,pdf2db_transform->Get(Form("rich_prob_vs_logkrec_%s"        ,det_name[id])),&(data->_logk[id]),0);
    data->_rich_logdq_lip             .AddTransformation(id,pdf2db_transform->Get(Form("rich_logdq_lip_vs_logkrec_%s"   ,det_name[id])),&(data->_logk[id]),2);
    data->_rich_ratio                 .AddTransformation(id,pdf2db_transform->Get(Form("rich_ratio_vs_logkrec_%s"       ,det_name[id])),&(data->_logk[id]),1);
    data->_rich_logdbeta              .AddTransformation(id,pdf2db_transform->Get(Form("rich_logdbeta_vs_logkrec_%s"    ,det_name[id])),&(data->_logk[id]),2);
    data->_rich_logdist               .AddTransformation(id,pdf2db_transform->Get(Form("rich_logdist_vs_logkrec_%s"     ,det_name[id])),&(data->_logk[id]),0);
  }
  // mass likelihood
  TFile* file_mass = TFile::Open("root://eosams.cern.ch//eos/ams/group/dbar/data/likelihood/release_v5_prod_10_mass_pdf2db.root","read");
  pdf2db_mass = (PDF2DB*) file_mass->Get("mass_pdf2db");
  file_mass->Close();
  int mass_interp    = PDF2::kIDF; // PDF2::kFast
  int mass_bell_type = PDF1::kRaw; // PDF1::kAKDE
  int mass_fit_type  = PDF1::kSpline;
  const char* par_name[2] = {"prot","deut"};
  string name;
  for (int im=0; im<2; im++) {
    for (int id=0; id<3; id++) {
      printf("Load likelihood_mass_%s_%s ...\n",det_name[id],par_name[im]);
      likelihood_mass[id][im] = new Likelihood(pdf2db_mass);
      name.assign(Form("tof_qu_unc_%s_%s_logkrec"     ,det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logk[id]),&(data->_tof_qu_unc.fInput[id])     ,name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("tof_qu_unc_%s_%s_logrrec"     ,det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logr[id]),&(data->_tof_qu_unc.fInput[id])     ,name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("tof_ql_unc_%s_%s_logkrec"     ,det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logk[id]),&(data->_tof_ql_unc.fInput[id])     ,name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("tof_ql_unc_%s_%s_logrrec"     ,det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logr[id]),&(data->_tof_ql_unc.fInput[id])     ,name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trd_logphe_%s_%s_logkrec"     ,det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logk[id]),&(data->_trd_logphe.fInput[id])     ,name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trd_logphe_%s_%s_logrrec"     ,det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logr[id]),&(data->_trd_logphe.fInput[id])     ,name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trd_logep_%s_%s_logkrec"      ,det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logk[id]),&(data->_trd_logep.fInput[id])      ,name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trd_logep_%s_%s_logrrec"      ,det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logr[id]),&(data->_trd_logep.fInput[id])      ,name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trd_ampl_%s_%s_logkrec"       ,det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logk[id]),&(data->_trd_amplonpath.fInput[id]) ,name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trd_ampl_%s_%s_logrrec"       ,det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logr[id]),&(data->_trd_amplonpath.fInput[id]) ,name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logchisqkax_%s_%s_logkrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logk[id]),&(data->_trk_logchisqkax.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logchisqkax_%s_%s_logrrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logr[id]),&(data->_trk_logchisqkax.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logchisqkay_%s_%s_logkrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logk[id]),&(data->_trk_logchisqkay.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logchisqkay_%s_%s_logrrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logr[id]),&(data->_trk_logchisqkay.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logchisqchx_%s_%s_logkrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logk[id]),&(data->_trk_logchisqchx.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logchisqchx_%s_%s_logrrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logr[id]),&(data->_trk_logchisqchx.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logchisqchy_%s_%s_logkrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logk[id]),&(data->_trk_logchisqchy.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logchisqchy_%s_%s_logrrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logr[id]),&(data->_trk_logchisqchy.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logscatt34x_%s_%s_logkrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logk[id]),&(data->_trk_logscatt34x.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logscatt34x_%s_%s_logrrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logr[id]),&(data->_trk_logscatt34x.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logscatt34y_%s_%s_logkrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logk[id]),&(data->_trk_logscatt34y.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logscatt34y_%s_%s_logrrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logr[id]),&(data->_trk_logscatt34y.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logscatt56x_%s_%s_logkrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logk[id]),&(data->_trk_logscatt56x.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logscatt56x_%s_%s_logrrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logr[id]),&(data->_trk_logscatt56x.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logscatt56y_%s_%s_logkrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logk[id]),&(data->_trk_logscatt56y.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_logscatt56y_%s_%s_logrrec",det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logr[id]),&(data->_trk_logscatt56y.fInput[id]),name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_q_unc_%s_%s_logkrec"      ,det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logk[id]),&(data->_trk_qinn_unc.fInput[id])   ,name,"","",mass_interp,mass_bell_type,mass_fit_type);
      name.assign(Form("trk_q_unc_%s_%s_logrrec"      ,det_name[id],par_name[im])); likelihood_mass[id][im]->AddVariable(name,&(data->_logr[id]),&(data->_trk_qinn_unc.fInput[id])   ,name,"","",mass_interp,mass_bell_type,mass_fit_type);
    }
  }
  initialized = true;
}

void ClassifierManager::SetReaderVariables(TMVA::Reader* reader, int classifier_type, bool fgflag) {
  int var = int(classifier_type/10)%10;
  int id = int(classifier_type/100)%10;
  if (var==0) { // TOF
    reader->AddVariable("tof_logchisqtn" ,&(data->_tof_logchisqtn.fOutput[id]));
    reader->AddVariable("tof_logchisqcn" ,&(data->_tof_logchisqcn.fOutput[id]));
    reader->AddVariable("tof_logzprob"   ,&(data->_tof_logzprob  .fOutput[id]));
    reader->AddVariable("tof_logqasym"   ,&(data->_tof_logqasym  .fOutput[id]));
    reader->AddVariable("tof_logdbeta"   ,&(data->_tof_logdbeta  .fOutput[id]));
    reader->AddVariable("tof_logisolat"  ,&(data->_tof_logisolat .fOutput[id]));
  }
  else if (var==1) { // TRD
    reader->AddVariable("trd_nhit"       ,&(data->_trd_nhit  .fOutput[id]));
    reader->AddVariable("trd_vertex"     ,&(data->_trd_vertex.fOutput[id]));
    reader->AddVariable("trd_logphe"     ,&(data->_trd_logphe.fOutput[id]));
    reader->AddVariable("trd_logep"      ,&(data->_trd_logep .fOutput[id]));
  }
  else if (var==2) { // Tracker
    reader->AddVariable("trk_nyhits"     ,&(data->_trk_nyhits                 .fOutput[id]));
    reader->AddVariable("trk_nxhits"     ,&(data->_trk_nxhits                 .fOutput[id]));
    reader->AddVariable("trk_n"          ,&(data->_trk_n                      .fOutput[id]));
    reader->AddVariable("trk_logchisqkax",&(data->_trk_logchisqkax            .fOutput[id]));
    reader->AddVariable("trk_logchisqkay",&(data->_trk_logchisqkay            .fOutput[id]));
    reader->AddVariable("trk_logchisqchx",&(data->_trk_logchisqchx            .fOutput[id]));
    reader->AddVariable("trk_logchisqchy",&(data->_trk_logchisqchy            .fOutput[id]));
    reader->AddVariable("trk_logscatt34x",&(data->_trk_logscatt34x            .fOutput[id]));
    reader->AddVariable("trk_logscatt34y",&(data->_trk_logscatt34y            .fOutput[id]));
    reader->AddVariable("trk_logscatt56x",&(data->_trk_logscatt56x            .fOutput[id]));
    reader->AddVariable("trk_logscatt56y",&(data->_trk_logscatt56y            .fOutput[id]));
    reader->AddVariable("trk_logdinvru"  ,&(data->_trk_logdinvru              .fOutput[id]));
    reader->AddVariable("trk_logdinvrd"  ,&(data->_trk_logdinvrd              .fOutput[id]));
    reader->AddVariable("trk_logcc"      ,&(data->_trk_logcc                  .fOutput[id]));
    reader->AddVariable("trk_logdinvr"   ,&(data->_trk_logdinvr               .fOutput[id]));
    reader->AddVariable("trk_logdinvrms" ,&(data->_trk_logdinvrms             .fOutput[id]));
    reader->AddVariable("trk_logdinvrck" ,&(data->_trk_logdinvrck             .fOutput[id]));
    reader->AddVariable("trk_r"          ,&(data->_trk_r                      .fOutput[id]));
    reader->AddVariable("trk_q"          ,&(data->_trk_q                      .fOutput[id]));
    reader->AddVariable("trk_logqyasym"  ,&(data->_trk_logqyasym              .fOutput[id]));
    reader->AddVariable("trk_logl2off"   ,&(data->_trk_logedep_l2_offtrack    .fOutput[id]));
    reader->AddVariable("trk_loginnoff0" ,&(data->_trk_logedep_inn_offtrack[0].fOutput[id]));
    reader->AddVariable("trk_loginnoff1" ,&(data->_trk_logedep_inn_offtrack[1].fOutput[id]));
  }
  else if (var>=3) { // RICH
    if (!fgflag) {
      reader->AddVariable("rich_npmt"        ,&(data->_rich_npmt           .fOutput[id]));
      reader->AddVariable("rich_nhit"        ,&(data->_rich_nhit_uncorr    .fOutput[id]));
      reader->AddVariable("rich_correct"     ,&(data->_rich_nhit           .fOutput[id]));
      reader->AddVariable("rich_reflect"     ,&(data->_rich_nhit_refl      .fOutput[id]));
      reader->AddVariable("rich_hyp_dir"     ,&(data->_rich_nhit_hyp[0]    .fOutput[id]));
      reader->AddVariable("rich_hyp_ref"     ,&(data->_rich_nhit_hyp[1]    .fOutput[id]));
      reader->AddVariable("rich_notused1"    ,&(data->_rich_nhit_notused[0].fOutput[id]));
      reader->AddVariable("rich_notused2"    ,&(data->_rich_nhit_notused[1].fOutput[id]));
      reader->AddVariable("rich_nclus"       ,&(data->_rich_nclus          .fOutput[id]));
      reader->AddVariable("rich_nhit_00"     ,&(data->_rich_tot_hit[0][0]  .fOutput[id]));
      reader->AddVariable("rich_nhit_01"     ,&(data->_rich_tot_hit[0][1]  .fOutput[id]));
      reader->AddVariable("rich_nhit_02"     ,&(data->_rich_tot_hit[0][2]  .fOutput[id]));
    //reader->AddVariable("rich_nhit_03"     ,&(data->_rich_tot_hit[0][3]  .fOutput[id]));
      reader->AddVariable("rich_nhit_04"     ,&(data->_rich_tot_hit[0][4]  .fOutput[id]));
      reader->AddVariable("rich_nhit_10"     ,&(data->_rich_tot_hit[1][0]  .fOutput[id]));
      reader->AddVariable("rich_nhit_11"     ,&(data->_rich_tot_hit[1][1]  .fOutput[id]));
      reader->AddVariable("rich_nhit_12"     ,&(data->_rich_tot_hit[1][2]  .fOutput[id]));
      reader->AddVariable("rich_nhit_13"     ,&(data->_rich_tot_hit[1][3]  .fOutput[id]));
      reader->AddVariable("rich_nhit_14"     ,&(data->_rich_tot_hit[1][4]  .fOutput[id]));
      reader->AddVariable("rich_np_exp"      ,&(data->_rich_np_exp         .fOutput[id]));
      reader->AddVariable("rich_exp_res"     ,&(data->_rich_exp_res        .fOutput[id]));
      reader->AddVariable("rich_exp_rms"     ,&(data->_rich_exp_rms        .fOutput[id]));
      reader->AddVariable("rich_lognp_min"   ,&(data->_rich_lognp_min      .fOutput[id]));
      reader->AddVariable("rich_logdbeta_lip",&(data->_rich_logdbeta_lip   .fOutput[id]));
      reader->AddVariable("rich_logqpmtcons" ,&(data->_rich_logqpmtcons    .fOutput[id]));
      reader->AddVariable("rich_q"           ,&(data->_rich_q              .fOutput[id]));
      reader->AddVariable("rich_prob"        ,&(data->_rich_prob           .fOutput[id]));
      reader->AddVariable("rich_logdq_lip"   ,&(data->_rich_logdq_lip      .fOutput[id]));
      reader->AddVariable("rich_ratio"       ,&(data->_rich_ratio          .fOutput[id]));
      reader->AddVariable("rich_logdbeta"    ,&(data->_rich_logdbeta       .fOutput[id]));
      reader->AddVariable("rich_logdist"     ,&(data->_rich_logdist        .fOutput[id]));
    } 
    else { 
      reader->AddVariable("rich_npmt"        ,&(data->_rich_npmt           .fInput[2]));
      reader->AddVariable("rich_nhit"        ,&(data->_rich_nhit_uncorr    .fInput[2]));
      reader->AddVariable("rich_hyp_dir"     ,&(data->_rich_nhit_hyp[0]    .fInput[2]));
      reader->AddVariable("rich_hyp_ref"     ,&(data->_rich_nhit_hyp[1]    .fInput[2]));
      reader->AddVariable("rich_notused1"    ,&(data->_rich_nhit_notused[0].fInput[2]));
      reader->AddVariable("rich_notused2"    ,&(data->_rich_nhit_notused[1].fInput[2]));
      reader->AddVariable("rich_nclus"       ,&(data->_rich_nclus          .fInput[2]));
      reader->AddVariable("rich_np_exp"      ,&(data->_rich_np_exp         .fInput[2]));
      reader->AddVariable("rich_logdbeta_lip",&(data->_rich_logdbeta_lip   .fInput[2]));
      reader->AddVariable("rich_logqpmtcons" ,&(data->_rich_logqpmtcons    .fInput[2]));
      reader->AddVariable("rich_q"           ,&(data->_rich_q              .fInput[2]));
      reader->AddVariable("rich_prob"        ,&(data->_rich_prob           .fInput[2]));
      reader->AddVariable("rich_logdq_lip"   ,&(data->_rich_logdq_lip      .fInput[2]));
      reader->AddVariable("rich_ratio"       ,&(data->_rich_ratio          .fInput[2]));
      reader->AddVariable("rich_logdbeta"    ,&(data->_rich_logdbeta       .fInput[2]));
      Float_t a,b,c,d,e; 
      reader->AddSpectator("S_rig:=trk_rigidity",&a);
      reader->AddSpectator("S_Richbeta:=rich_beta",&b);
      reader->AddSpectator("S_Mass:=fabs(trk_rigidity*sqrt(1.-rich_beta*rich_beta)/rich_beta/0.93827231)",&c);
      reader->AddSpectator("S_Ekin:=((1./sqrt(1.-rich_beta*rich_beta))-1)*0.93827231",&d);
      reader->AddSpectator("S_Richqc2:=rich_prob>=0.01 && rich_npmt>=3 && rich_ratio>0.4 && rich_logqpmtcons<=1. && rich_np_exp>=2 && rich_logdbeta_lip<=-2.3 ",&e);
    }
  }
}

double ClassifierManager::GetClassifier(Event* event, int classifier_type) {
  if (!initialized) Init();
  if (readers.find(classifier_type)==readers.end()) return -100; 
  if (!((current_run==event->SHeader->run)&&(current_event==event->SHeader->event))) { 
    data->FillData(event,0);
    data->FillData(event,1);
    data->FillData(event,2);
    data->ApplyTransformations(); 
    current_run = event->SHeader->run;
    current_event = event->SHeader->event; 
  }
  float value = (multi[classifier_type].compare("Multi")==0) ? 
    readers[classifier_type]->EvaluateMulticlass(0,names[classifier_type].c_str()) :
    readers[classifier_type]->EvaluateMVA(names[classifier_type].c_str());
  // if (multi[classifier_type].compare("Multi")==0) {
  //   std::cout << readers[classifier_type]->EvaluateMulticlass(0,names[classifier_type].c_str()) << "  " 
  //        << readers[classifier_type]->EvaluateMulticlass(1,names[classifier_type].c_str()) << "  "
  //        << readers[classifier_type]->EvaluateMulticlass(2,names[classifier_type].c_str()) << "  "
  //        << readers[classifier_type]->EvaluateMulticlass(3,names[classifier_type].c_str()) << std::endl; }
  if (names[classifier_type].compare("MLP")==0) value = value*2-1; 
  return value;
}

const double LOGVARMIN = -10;
const double LOGVARMAX =  10;
static double logvar(double var) {
  if      (var>=pow(10,LOGVARMAX)) return LOGVARMAX;
  else if (var<=pow(10,LOGVARMIN)) return LOGVARMIN;
  return log10(var);
}

double get_logk(double beta) {
  double b = fabs(beta);
  b = (b>=0.999999) ? 0.999999 : b;
  double g = (b>=0.999999) ? 1./sqrt(1-0.999999*0.999999) : 1./sqrt(1-b*b);
  return (g>1) ? log10(1.87561293*(g-1)/2) : -10;
}

// beta tuning 
static double get_beta_tuned(double beta, double bgen, int id) {
  static double beta_tuning[4][3] = { // B1119
    {1.00102,    0.00118,    0.14196},
    {0.99975,   -0.00069,   -0.01122},
    {1.00001,   -0.00011,   -0.04516},
    {1.00102,    0.00118,    0.14196}
  };
  // TOF
  if ( (id==0)||(id==3) ) {
    if (bgen<=0) return 1/((1/beta)/beta_tuning[id][0]);
    double discard_defa = (1/beta-1/bgen)*bgen;
    double discard_tune = (discard_defa-beta_tuning[id][1])/(1-beta_tuning[id][2]);
    return bgen/(1+discard_tune);
  }
  // RICH/NaF
  if (bgen<=0) return beta/beta_tuning[id][0];
  double discard_defa = (beta-bgen)/bgen;
  double discard_tune = (discard_defa-beta_tuning[id][1])/(1-beta_tuning[id][2]);
  return bgen*(1+discard_tune);
}

void ClassifierData::FillData(Event* event, int id) { 
  if (!event) return; 

  // corrections to be applied to beta/rigidity should go here, if needed ... 

  NtpHeader  *Header  = event->Header; 
  NtpTrd     *Trd     = event->Trd;
  NtpTof     *Tof     = event->Tof;
  NtpTracker *Tracker = event->Tracker;
  NtpRich    *Rich    = event->Rich; 
  _logk[id] = get_logk(get_beta_tuned((id==0)?Tof->beta:Rich->beta_corrected,0,id));
  _logr[id] = (fabs(Tracker->rig[1][(id==0)?1:2])>0) ? log10(fabs(Tracker->rig[1][(id==0)?1:2])) : -10;
  float betaR = (id==0) ? fabs(Tof->beta)*Tracker->rig[1][1] : fabs(Rich->beta)*Tracker->rig[1][2];
  // TOF
  { 
    double tof_q_lay_min  =  1e+30;
    double tof_q_lay_max  = -1e+30;
    for (int ilay=0; ilay<4; ilay++) {
      if (int(Tof->flagp[ilay]/100)!=0) continue; // drop bad clusters
      if (Tof->q_lay[ilay]<=0) continue; // drop bad charges
      tof_q_lay_min = TMath::Min(tof_q_lay_min,(double)Tof->q_lay[ilay]);
      tof_q_lay_max = TMath::Max(tof_q_lay_max,(double)Tof->q_lay[ilay]);
    }
    double tof_isolation_num = 0;
    double tof_isolation_den = 0;
    for (int ilay=0; ilay<1; ilay++) {
      tof_isolation_num += Tof->edep[ilay][1]+Tof->edep[ilay][2]-Tof->edep[ilay][4]; // excluding overlap clusters
      tof_isolation_den += Tof->edep[ilay][0];
    }
    double tof_logqasym   = logvar(fabs(tof_q_lay_max-tof_q_lay_min));
    double tof_logchisqtn = logvar(Tof->chisqtn);
    double tof_logchisqcn = logvar(Tof->chisqcn);
    double tof_logzprob   = logvar(1-exp(-Tof->z_like));
    double tof_logdbeta   = logvar(fabs(fabs(Tof->evgeni_beta)-fabs(Tof->beta)));
    double tof_logisolat  = logvar((fabs(tof_isolation_den)>0)?tof_isolation_num/tof_isolation_den:0);
    // set
    _tof_logqasym  .SetValue(tof_logqasym  ,id);  
    _tof_logchisqtn.SetValue(tof_logchisqtn,id);
    _tof_logchisqcn.SetValue(tof_logchisqcn,id);
    _tof_logzprob  .SetValue(tof_logzprob  ,id);  
    _tof_logdbeta  .SetValue(tof_logdbeta  ,id);  
    _tof_logisolat .SetValue(tof_logisolat ,id); 
    int n; double rms;
    _tof_qu_unc    .SetValue(event->Tof->GetQ(n,rms,0x3,false,true),id); 
    _tof_ql_unc    .SetValue(event->Tof->GetQ(n,rms,0xc,false,true),id);
  }
  // TRD
  { 
    int trd_index = 0;
    if (Trd->trdk_like_nhit[0]<Trd->trdk_like_nhit[2]) trd_index = 2;
    int    trd_nhit   = Trd->trdk_like_nhit[trd_index];
    double trd_logphe = logvar(Trd->trdk_like_phe(trd_index));
    double trd_logep  = logvar(Trd->trdk_like_ep(trd_index));
    double trd_amplonpath = 0;
    double trd_nhit_ampl = 0;
    for (int ihit=0; ihit<20; ihit++) {
      if ( (Trd->trdk_ampl[ihit]<=0)||(Trd->trdk_path[ihit]<=0) ) continue;
      trd_amplonpath += Trd->trdk_ampl[ihit]/Trd->trdk_path[ihit];
      trd_nhit_ampl += 1;
    }
    if (trd_nhit_ampl>0) trd_amplonpath = sqrt(trd_amplonpath/trd_nhit_ampl);
    int trd_vertex = Trd->vertex_nsegx+Trd->vertex_nsegy;
    // set
    _trd_nhit      .SetValue(trd_nhit      ,id);
    _trd_vertex    .SetValue(trd_vertex    ,id);
    _trd_logphe    .SetValue(trd_logphe    ,id);
    _trd_logep     .SetValue(trd_logep     ,id);
    _trd_amplonpath.SetValue(trd_amplonpath,id);
  }
  // Tracker
  { 
    const int which_trk_charge = 0;
    double trk_logchisqkax = logvar(Tracker->chisqn[1][0]);
    double trk_logchisqkay = logvar(Tracker->chisqn[1][1]);
    double trk_logchisqchx = logvar(Tracker->chisqn[6][0]);
    double trk_logchisqchy = logvar(Tracker->chisqn[6][1]);
    double trk_logscatt34x = logvar(fabs(Tracker->scat_rigi_theta[0][0]*betaR));
    double trk_logscatt34y = logvar(fabs(Tracker->scat_rigi_theta[0][1]*betaR));
    double trk_logscatt56x = logvar(fabs(Tracker->scat_rigi_theta[1][0]*betaR));
    double trk_logscatt56y = logvar(fabs(Tracker->scat_rigi_theta[1][1]*betaR));
    double trk_logdinvrck  = logvar((fabs(Tracker->rig[1][1])>0)?fabs(1-Tracker->rig[8][1]/Tracker->rig[1][1]):-1);
    double trk_logdinvrms  = logvar((fabs(Tracker->rig[1][1])>0)?fabs(1-Tracker->rig[7][1]/Tracker->rig[1][1]):-1);
    double trk_logdinvr    = logvar((fabs(Tracker->rig[1][1])>0)?fabs(1-Tracker->rig[6][1]/Tracker->rig[1][1]):-1);
    double trk_logdinvru   = logvar((fabs(Tracker->rig[2][1])>0)?fabs(1-Tracker->rig[1][1]/Tracker->rig[2][1]):-1);
    double trk_logdinvrd   = logvar((fabs(Tracker->rig[3][1])>0)?fabs(1-Tracker->rig[1][1]/Tracker->rig[3][1]):-1);
    int    trk_nyhits      = 0;
    int    trk_nxhits      = 0;
    int    trk_n           = Tracker->q_inn_nhit[which_trk_charge];
    double trk_q           = Tracker->q_inn[which_trk_charge];
    double trk_qy_lay_max  = -1e+10;
    double trk_qy_lay_min  =  1e+10;
    double trk_qx_lay_max  = -1e+10;
    double trk_qx_lay_min  =  1e+10;
    double trk_qinn_unc    = 0;
    double trk_qyl1        = Tracker->q_lay[1][0];
    for (int ilay=1; ilay<8; ilay++) {
      if ( (Tracker->q_lay[1][ilay]>0.1)&&((Tracker->q_clu_status[1][ilay]&0x13D)==0)&&((Tracker->q_clu_status[0][ilay]&0x100)==0) ) {
        trk_qy_lay_max = TMath::Max(trk_qy_lay_max,(double)Tracker->q_lay[1][ilay]);
        trk_qy_lay_min = TMath::Min(trk_qy_lay_min,(double)Tracker->q_lay[1][ilay]);
        trk_qinn_unc += Tracker->q_lay_uncorr[1][ilay];
        trk_nyhits += 1;
      }
      if ( (Tracker->q_lay[0][ilay]>0.1)&&((Tracker->q_clu_status[0][ilay]&0x13D)==0)&&((Tracker->q_clu_status[1][ilay]&0x100)==0) ) {
        trk_qx_lay_max = TMath::Max(trk_qx_lay_max,(double)Tracker->q_lay[0][ilay]);
        trk_qx_lay_min = TMath::Min(trk_qx_lay_min,(double)Tracker->q_lay[0][ilay]);
        trk_qinn_unc += Tracker->q_lay_uncorr[0][ilay];
        trk_nxhits += 1;
      }
    }
    if (trk_nyhits+trk_nxhits>0) trk_qinn_unc /= trk_nyhits+trk_nxhits;
    double trk_logqyasym = logvar(trk_qy_lay_max-trk_qy_lay_min);
    double trk_r         = (Tracker->q_inn[which_trk_charge]>0) ? Tracker->q_inn_rms[which_trk_charge]/Tracker->q_inn[which_trk_charge] : 0;
    double trk_logcc     = -10;
    for (int ilay=0; ilay<7; ilay++) {
      if (Tracker->q_lay[1][ilay+1]<=0.1) continue;
      if (fabs(Tracker->rig[6][1])<=0) continue;
      if (fabs(Tracker->inn_rig[ilay])<=0) continue;
      trk_logcc = TMath::Max(trk_logcc,(double)logvar(fabs(Tracker->inn_rig[ilay]/Tracker->rig[6][1]-1)));
    }
//  double trk_logfeet56 = logvar(TMath::Min(Tracker->feet_dist[3],Tracker->feet_dist[4]));
    double trk_logedep_l2_offtrack = logvar(Tracker->edep_lay[1][1][6]/1000);
    double trk_logedep_inn_offtrack[2] = {0,0};
    int    trk_nlay[2]                 = {0,0};
    for (int i=0; i<2; i++) {
      trk_logedep_inn_offtrack[i] = 0;
      trk_nlay[i] = 0;
    }
    for (int ilay=1; ilay<8; ilay++) {
      int index = 0;
      if (Tracker->q_lay[1][ilay]<=0.1) index = 1;
      trk_logedep_inn_offtrack[index] += Tracker->edep_lay[1][ilay][2]/1000; // at 1 cm
      trk_nlay[index] += 1;
    }
    for (int i=0; i<2; i++) trk_logedep_inn_offtrack[i] = logvar((trk_nlay[i]>0)?trk_logedep_inn_offtrack[i]/trk_nlay[i]:0);
    // set
    _trk_logchisqkax            .SetValue(trk_logchisqkax            ,id);
    _trk_logchisqkay            .SetValue(trk_logchisqkay            ,id);
    _trk_logchisqchx            .SetValue(trk_logchisqchx            ,id);
    _trk_logchisqchy            .SetValue(trk_logchisqchy            ,id);
    _trk_logscatt34x            .SetValue(trk_logscatt34x            ,id);
    _trk_logscatt34y            .SetValue(trk_logscatt34y            ,id);
    _trk_logscatt56x            .SetValue(trk_logscatt56x            ,id);
    _trk_logscatt56y            .SetValue(trk_logscatt56y            ,id);
    _trk_logdinvru              .SetValue(trk_logdinvru              ,id);
    _trk_logdinvrd              .SetValue(trk_logdinvrd              ,id);
    _trk_logcc                  .SetValue(trk_logcc                  ,id);
    _trk_logdinvr               .SetValue(trk_logdinvr               ,id);
    _trk_logdinvrms             .SetValue(trk_logdinvrms             ,id);
    _trk_logdinvrck             .SetValue(trk_logdinvrck             ,id);
    _trk_nyhits                 .SetValue(trk_nyhits                 ,id);
    _trk_nxhits                 .SetValue(trk_nxhits                 ,id);
    _trk_n                      .SetValue(trk_n                      ,id);
    _trk_r                      .SetValue(trk_r                      ,id);
    _trk_q                      .SetValue(trk_q                      ,id);
    _trk_qyl1                   .SetValue(trk_qyl1                   ,id);
    _trk_qinn_unc               .SetValue(trk_qinn_unc               ,id);
    _trk_logqyasym              .SetValue(trk_logqyasym              ,id);
    _trk_logedep_l2_offtrack    .SetValue(trk_logedep_l2_offtrack    ,id);
    _trk_logedep_inn_offtrack[0].SetValue(trk_logedep_inn_offtrack[0],id);
    _trk_logedep_inn_offtrack[1].SetValue(trk_logedep_inn_offtrack[1],id);
  }
  // RICH
  { 
    double rich_window[2] = {4e-3,1.2e-3};
//  double rich_radx        = Rich->rad_coo[0];
//  double rich_rady        = Rich->rad_coo[1];
//  double rich_radr        = sqrt(pow(rich_radx,2)+pow(rich_rady,2));
//  double rich_radd        = TMath::Max(fabs(rich_radx),fabs(rich_rady));
    double rich_np_exp      = Rich->np_exp;
    double rich_exp_res     = Rich->beta_res;
    double rich_exp_rms     = Rich->beta_rms;
    int    rich_npmt        = Rich->npmt;
    double rich_q           = Rich->q;
    double rich_prob        = Rich->prob;
    int    rich_nhit_uncorr = Rich->nhit_uncorr;
    int    rich_nhit        = Rich->nhit;
    int    rich_nhit_refl   = Rich->nhit_refl;
    double rich_logdist     = logvar(fabs(Rich->distance_tile_border));
    vector<int> rich_pmt_crossed[2]; // primary/secondary
    int rich_nhit_crossed[2] = {0};
    for (int i=0; i<2; i++) rich_nhit_crossed[i] = 0; // primary/secondary
    for (int ipmt=0; ipmt<5; ipmt++) {
      if (Rich->pmt_np_uncorr[ipmt]<5.) continue; // is crossed 
      int index = (Rich->pmt_dist[ipmt]>3.5) ? 1 : 0;
      rich_pmt_crossed[index].push_back(Rich->pmt_pmt[ipmt]);
      rich_nhit_crossed[index] += Rich->pmt_nhit_uncorr[ipmt];
    }
    int rich_nhit_crossed_used[2] = {0};
    for (int i=0; i<2; i++) rich_nhit_crossed_used[i] = 0; // primary/secondary
    double rich_lognp_min = LOGVARMAX;
    double rich_lognp_max = LOGVARMIN;
    int rich_nhit_hyp[2] = {0,0};
    for (int i=0; i<2; i++) rich_nhit_hyp[i] = 0; // direct/reflected
    for (int ihit=0; ihit<TMath::Min(30,Header->nrichhit); ihit++) {
      if (Rich->hit_np_uncorr[ihit]<=0) continue;
      bool on_ring = ( (Rich->hit_used[ihit]>=0)&&(Rich->hit_used[ihit]<=1) );
      for (int is=0; is<2; is++) {
        bool crossed = ( (find(rich_pmt_crossed[is].begin(),rich_pmt_crossed[is].end(),int(Rich->hit_chan[ihit]/16))!=rich_pmt_crossed[is].end()) );
        bool on_ring = ( (Rich->hit_used[ihit]>=0)&&(Rich->hit_used[ihit]<=1) );
        if (on_ring&&crossed) rich_nhit_crossed_used[is] += 1;
      }
      bool crossed_carlos = ( ((Rich->hit_stat[ihit]>>30)&0x1)==0x1 );
      bool crossed_me =
        ( (find(rich_pmt_crossed[0].begin(),rich_pmt_crossed[0].end(),int(Rich->hit_chan[ihit]/16))!=rich_pmt_crossed[0].end()) )||
        ( (find(rich_pmt_crossed[1].begin(),rich_pmt_crossed[1].end(),int(Rich->hit_chan[ihit]/16))!=rich_pmt_crossed[1].end()) );
      bool crossed = (crossed_carlos||crossed_me);
      if (on_ring) {
        rich_lognp_min = TMath::Min(rich_lognp_min,logvar(Rich->hit_np_uncorr[ihit]));
        rich_lognp_max = TMath::Max(rich_lognp_max,logvar(Rich->hit_np_uncorr[ihit]));
      }
      else {
        for (int ibeta=0; ibeta<2; ibeta++) {
          if (crossed) continue;
          if (fabs(Rich->hit_beta[ihit][ibeta]-1)>3*rich_window[id-1]) continue;
          rich_nhit_hyp[ibeta] += 1;
        }
      }
    }
    int rich_nhit_notused[2] = {
      Header->nrichhit-Rich->nhit_uncorr,
      Header->nrichhit-rich_nhit_crossed[0]-rich_nhit_crossed[1]-(Rich->nhit_uncorr-rich_nhit_crossed_used[0]-rich_nhit_crossed_used[1])
    };
    double rich_logdbeta_lip = logvar(fabs(Rich->lip_beta-Rich->beta));
    double rich_logqpmtcons = logvar(fabs(Rich->q_consistency));
    double rich_logdq_lip = logvar(fabs(Rich->lip_q-Rich->q));
    double rich_logdbeta = logvar(fabs(Tof->beta-Rich->beta));
    double rich_ratio = (Rich->tot_p_uncorr>0) ? Rich->np_uncorr/Rich->tot_p_uncorr : 0;
    int rich_nclus = 0; for (int is=0; is<10; is++) if ((Rich->clus_mean[is]-Rich->beta)>0.01) rich_nclus += 1;
    // branching
    vector<int> rich_hit_on_ring[2]; // direct/reflected 
    vector<int> rich_hit_off_ring; // of the list
    for (int ihit=0; ihit<TMath::Min(30,Header->nrichhit); ihit++) {
      bool crossed_carlos = ( ((Rich->hit_stat[ihit]>>30)&0x1)==0x1 );
      bool crossed_me =
        ( (find(rich_pmt_crossed[0].begin(),rich_pmt_crossed[0].end(),int(Rich->hit_chan[ihit]/16))!=rich_pmt_crossed[0].end()) )||
        ( (find(rich_pmt_crossed[1].begin(),rich_pmt_crossed[1].end(),int(Rich->hit_chan[ihit]/16))!=rich_pmt_crossed[1].end()) );
      if (crossed_carlos||crossed_me) continue;
      bool on_ring = ( (Rich->hit_used[ihit]>=0)&&(Rich->hit_used[ihit]<=1) );
      if (on_ring) rich_hit_on_ring[Rich->hit_used[ihit]].push_back(ihit);
      else rich_hit_off_ring.push_back(ihit);
    }
    double rich_branch_beta = 0;
    double rich_branch_rms = 1e5;
    // loop on the two branches (direct/reflected)
    for (int iused=0; iused<2; iused++) {
      // exclude hit one-by-one from the branch, or exclude none        
      int size = (int) rich_hit_on_ring[iused].size();
      for (int ihyp=0; ihyp<size+1; ihyp++) {
        double mean = 0;
        double rms  = 0;
        int    n    = 0;
        double np   = 0;
        // average of opposite hypothesys of a branch (reflected/direct)
        for (int ih=0; ih<size; ih++) {
          if (ih==ihyp) continue;
          int ihit = rich_hit_on_ring[iused][ih];
          if (Rich->hit_beta[ihit][1-iused]<=0.5) continue;
          if (Rich->hit_beta[ihit][1-iused]>1+5*rich_window[id-1]) continue;
          mean += Rich->hit_beta[ihit][1-iused];
          rms += pow(Rich->hit_beta[ihit][1-iused],2);
          n++;
          np += Rich->hit_np_uncorr[ihit];
        }
        if ( (mean<0.5)||(n==0) ) continue;
        // try to add hits from the other branch 
        for (int ih=0; ih<(int)rich_hit_on_ring[1-iused].size(); ih++) {
          int ihit = rich_hit_on_ring[1-iused][ih];
          for (int iu=0; iu<2; iu++) {
            if (fabs(mean-Rich->hit_beta[ihit][iu])>3*rich_window[id-1]) continue;
            mean += Rich->hit_beta[ihit][iu];
            rms += pow(Rich->hit_beta[ihit][iu],2);
            n++;
            np += Rich->hit_np_uncorr[ihit];
          }
        }
        if ( (mean<0.5)||(n==0) ) continue;
        // try to add hits from outside 
        for (int ih=0; ih<(int)rich_hit_off_ring.size(); ih++) {
          int ihit = rich_hit_off_ring[ih];
          for (int iu=0; iu<2; iu++) {
            if (fabs(mean-Rich->hit_beta[ihit][iu])>3*rich_window[id-1]) continue;
            mean += Rich->hit_beta[ihit][iu];
            rms += pow(Rich->hit_beta[ihit][iu],2);
            n++;
            np += Rich->hit_np_uncorr[ihit];
          }
        }
        if (n<2) continue;
        mean /= n;
        rms /= n;
        rms = sqrt(rms-mean*mean);
        if (rms<rich_branch_rms) {
          rich_branch_beta = mean;
          rich_branch_rms = rms;
        }
      }
    }
    double rich_branch_dbeta = -10;
    if (rich_branch_beta>0) rich_branch_dbeta = Rich->beta-rich_branch_beta;
    // set
    _rich_np_exp      .SetValue(rich_np_exp      ,id);
    _rich_exp_res     .SetValue(rich_exp_res     ,id);
    _rich_exp_rms     .SetValue(rich_exp_rms     ,id);
    _rich_lognp_min   .SetValue(rich_lognp_min   ,id); 
    _rich_logdbeta_lip.SetValue(rich_logdbeta_lip,id);
    _rich_logqpmtcons .SetValue(rich_logqpmtcons ,id);
    _rich_q           .SetValue(rich_q           ,id);
    _rich_prob        .SetValue(rich_prob        ,id);
    _rich_npmt        .SetValue(rich_npmt        ,id);
    _rich_logdq_lip   .SetValue(rich_logdq_lip   ,id);
    _rich_ratio       .SetValue(rich_ratio       ,id);
    _rich_logdbeta    .SetValue(rich_logdbeta    ,id);
    _rich_nhit_uncorr .SetValue(rich_nhit_uncorr ,id);
    _rich_nhit        .SetValue(rich_nhit        ,id);
    _rich_nhit_refl   .SetValue(rich_nhit_refl   ,id);
    for (int i=0; i<2; i++) {
      _rich_nhit_hyp[i]    .SetValue(rich_nhit_hyp[i]    ,id);
      _rich_nhit_notused[i].SetValue(rich_nhit_notused[i],id);
    }
    _rich_logdist     .SetValue(rich_logdist,id);
    _rich_nclus       .SetValue(rich_nclus  ,id);
    for (int i=0; i<2; i++) 
      for (int j=0; j<5; j++)
        _rich_tot_hit[i][j].SetValue(Rich->tot_hit[i][j],id);
  }
}

void ClassifierData::ApplyTransformations() {
  // TOF
  _tof_logchisqtn.ApplyTransformation();
  _tof_logchisqcn.ApplyTransformation();
  _tof_logzprob  .ApplyTransformation();
  _tof_logqasym  .ApplyTransformation();
  _tof_logdbeta  .ApplyTransformation();
  _tof_logisolat .ApplyTransformation();
  _tof_qu_unc    .ApplyTransformation();
  _tof_ql_unc    .ApplyTransformation();
  // TRD
  _trd_logphe    .ApplyTransformation();
  _trd_logep     .ApplyTransformation();
  _trd_nhit      .ApplyTransformation();
  _trd_vertex    .ApplyTransformation();
  _trd_amplonpath.ApplyTransformation();
  // Tracker
  _trk_logchisqkax        .ApplyTransformation();
  _trk_logchisqkay        .ApplyTransformation();
  _trk_logchisqchx        .ApplyTransformation();
  _trk_logchisqchy        .ApplyTransformation();
  _trk_logscatt34x        .ApplyTransformation();
  _trk_logscatt34y        .ApplyTransformation();
  _trk_logscatt56x        .ApplyTransformation();
  _trk_logscatt56y        .ApplyTransformation();
  _trk_logdinvru          .ApplyTransformation();
  _trk_logdinvrd          .ApplyTransformation();
  _trk_logcc              .ApplyTransformation();
  _trk_logdinvr           .ApplyTransformation();
  _trk_logdinvrms         .ApplyTransformation();
  _trk_logdinvrck         .ApplyTransformation();
  _trk_nyhits             .ApplyTransformation();
  _trk_nxhits             .ApplyTransformation();
  _trk_n                  .ApplyTransformation();
  _trk_r                  .ApplyTransformation();
  _trk_q                  .ApplyTransformation();
  _trk_qyl1               .ApplyTransformation();
  _trk_qinn_unc           .ApplyTransformation();
  _trk_logqyasym          .ApplyTransformation();
  _trk_logedep_l2_offtrack.ApplyTransformation();
  for (int i=0; i<2; i++) _trk_logedep_inn_offtrack[i].ApplyTransformation();
  // RICH
  _rich_np_exp      .ApplyTransformation();
  _rich_exp_res     .ApplyTransformation();
  _rich_exp_rms     .ApplyTransformation();
  _rich_lognp_min   .ApplyTransformation();
  _rich_logdbeta_lip.ApplyTransformation();
  _rich_logqpmtcons .ApplyTransformation();
  _rich_q           .ApplyTransformation();
  _rich_prob        .ApplyTransformation();
  _rich_npmt        .ApplyTransformation();
  _rich_logdq_lip   .ApplyTransformation();
  _rich_ratio       .ApplyTransformation();
  _rich_logdbeta    .ApplyTransformation();
  _rich_nhit_uncorr .ApplyTransformation();
  _rich_nhit        .ApplyTransformation();
  _rich_nhit_refl   .ApplyTransformation();
  for (int i=0; i<2; i++) {
    _rich_nhit_hyp[i]    .ApplyTransformation();
    _rich_nhit_notused[i].ApplyTransformation();
    for (int j=0; j<5; j++) _rich_tot_hit[i][j].ApplyTransformation(); 
  }
  _rich_logdist.ApplyTransformation();
  _rich_nclus  .ApplyTransformation();
}

double ClassifierManager::GetMinusLogLikelihood(Event* event, int ll_type) {
  const double def_mll = 300;
  int im = (ll_type%10);
  int id = int(ll_type/10)%10; 
  if ( (im<0)||(im>1)||(id<0)||(id>2) ) return def_mll; 
  if (!likelihood_mass[id][im]) return def_mll;
  if (!((current_run==event->SHeader->run)&&(current_event==event->SHeader->event))) {
    data->FillData(event,0);
    data->FillData(event,1);
    data->FillData(event,2);
    data->ApplyTransformations();
    current_run = event->SHeader->run;
    current_event = event->SHeader->event;
  }
  double mll = -likelihood_mass[id][im]->EvalLog();
  return (mll>def_mll)?def_mll:mll;
} 
