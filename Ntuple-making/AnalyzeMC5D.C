#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "Selections5D.h"
#include "Commonglobals.cpp"

using namespace std;
void Trigg (TTree *albero,int i,TNtuple *ntupla);
int Trig_Num;
double Trig;
double TotalTrig=0;
double Totalevents=0;
double totaltrig=0;
double totaltrig2=0;
double response[44][44];
double norm[44];
float BetanS=0;
int Number=0;
int efficienzagenbeta_P[18]= {0};
int efficienzagenbeta_D[18][6]= {{0}};
int preselectedbeta_P[18]= {0};
int preselectedbeta_D[18][6]= {{0}};
TH1F * selezioni_PHLMC[10];
TH1F * selected_PHLMC[10];
TH2F * selezioni_DHLMC[10];
TH2F * selected_DHLMC[10];
TH1F * PrescaledMC = new TH1F("PrescaledMC","PrescaledMC",43,0,43);
TH1F * UnbiasMC= new TH1F("UMC","UMC",43,0,43);
TH1F * UPreselectedMC= new TH1F("PreselectedMC","PreselectedMC",43,0,43);
int MC_type=0;

int AssignMC_type(float Massa_gen);

int main(int argc, char * argv[])
{

   gROOT->ProcessLine("#include <vector>");
   ///////////////////////////// CALIBR.
   int control=0;
   string calib=argv[2];
   string nomecal=("/storage/gpfs_ams/ams/users/fdimicco/Deutons/CodesforAnalysis/CALIBRAZIONI/"+calib+".root");
   TFile *_file2 = TFile::Open(nomecal.c_str());
   cout<<"calibrazione: "<<_file2<<endl;
   if(!_file2) {
      nomecal=("/storage/gpfs_ams/ams/users/fdimicco/Deutons/CodesforAnalysis/2011_07.root");
      _file2 = TFile::Open(nomecal.c_str());
      control=1;
   }
   Rig = (TSpline3 *) _file2->Get("Fit Results/Splines/Rig");
   beta = (TSpline3 *) _file2->Get("Fit Results/Splines/beta");
   betaNaF = (TF1 *) _file2->Get("Fit Results/Splines/SigmaInvBetaNaF_spl");
   betaAgl = (TF1 *) _file2->Get("Fit Results/Splines/SigmaInvBetaAgl_spl");
   eL1 = (TSpline3 *) _file2->Get("Fit Results/Splines/eL1");
   etofu =  (TSpline3 *) _file2->Get("Fit Results/Splines/etofu");
   etrack =  (TSpline3 *) _file2->Get("Fit Results/Splines/etrack");
   etofd =  (TSpline3 *) _file2->Get("Fit Results/Splines/etofd");
   EdepL1beta =  (TSpline3 *) _file2->Get("Fit Results/Splines/EdepL1beta");
   EdepTOFbeta =  (TSpline3 *) _file2->Get("Fit Results/Splines/EdepTOFbeta");
   EdepTrackbeta =  (TSpline3 *) _file2->Get("Fit Results/Splines/EdepTrackbeta");
   EdepTOFDbeta =  (TSpline3 *) _file2->Get("Fit Results/Splines/EdepTOFDbeta");
   Corr_L1 =  (TSpline3 *) _file2->Get("Fit Results/Splines/Corr_L1");
   Corr_TOFU =  (TSpline3 *) _file2->Get("Fit Results/Splines/Corr_TOFU");
   Corr_Track =  (TSpline3 *) _file2->Get("Fit Results/Splines/Corr_Track");
   Corr_TOFD =  (TSpline3 *) _file2->Get("Fit Results/Splines/Corr_TOFD");
   cout<<Rig<<" "<<beta<<" "<<" "<<betaNaF<<" "<<betaAgl<<" "<<eL1<<" "<<etofu<<" "<<etrack<<" "<<etofd<<" "<<EdepL1beta<<" "<<EdepTOFbeta<<" "<<EdepTrackbeta<<" "<<EdepTOFDbeta<<" "<<Corr_L1<<" "<<Corr_TOFU<<" "<<Corr_Track<<" "<<Corr_TOFD<<endl;

   string tagli[10]= {"Trigger","3of4 TOF","TRD Segments","Rigidity exists","Chi^2 R","Matching TOF","Matching TRD","In TRD Accept.","1 Particle","1 Tr. Track"};
   string nome;

   //////////////////////// DICHIARAZIONE IST. //////////////////////
   for(int j=0; j<10; j++) {
      nome="SelezioniHL_MC"+tagli[j];
      selezioni_PHLMC[j]=new TH1F(nome.c_str(),nome.c_str(),43,0,43);
      nome="SelectedHL_MC"+tagli[j];
      selected_PHLMC[j]=new TH1F(nome.c_str(),nome.c_str(),43,0,43);
      nome="SelezioniHL_MCD"+tagli[j];
      selezioni_DHLMC[j]=new TH2F(nome.c_str(),nome.c_str(),18,0,18,6,0,6);
      nome="SelectedHL_MCD"+tagli[j];
      selected_DHLMC[j]=new TH2F(nome.c_str(),nome.c_str(),18,0,18,6,0,6);
   }

   /////////////////////////////////////////////////////////////////

   string Variables[9]= {"N. Anti-clusters","Unused TOF Clusters","|Rup-Rdown|:R","Unused Tracker layers","Tracker: Y Hits without X","Track Chi^2","|E.dep(lower TOF) - E.dep(upper TOF)|","|E.dep(layer 2)-E.dep(layer 1)|","|E. dep.(tot)-E.dep.(track)|"};
   for(int u2=0; u2<9; u2++) {
      nome="Splines/Spline: "+Variables[u2]+"_SGNL";
      Signal[u2]=(TSpline3 *) _file1->Get(nome.c_str());
      nome="Splines/Spline: "+Variables[u2]+"_BKGND";
      Bkgnd[u2]=(TSpline3 *) _file1->Get(nome.c_str());
   }
   cout<<_file1<<endl;
   string VariablesRICH[9]= {"N. Anti-clusters","Unused TOF Clusters","|Rup-Rdown|:R","Unused Tracker layers","Tracker: Y Hits without X","Track Chi^2","RICH Hits: tot - used","Rich Photoelectrons","|E. dep.(tot)-E.dep.(track)|"};
   for(int u2=0; u2<9; u2++) {
      nome="Splines/Spline NaF: "+VariablesRICH[u2]+"_SGNL";
      SignalNaF[u2]=(TSpline3 *) _file3->Get(nome.c_str());
      nome="Splines/Spline NaF: "+VariablesRICH[u2]+"_BKGND";
      BkgndNaF[u2]=(TSpline3 *) _file3->Get(nome.c_str());
   }
   cout<<_file3<<endl;
   for(int u2=0; u2<9; u2++) {
      nome="Splines/Spline Agl: "+VariablesRICH[u2]+"_SGNL";
      SignalAgl[u2]=(TSpline3 *) _file3b->Get(nome.c_str());
      nome="Splines/Spline Agl: "+VariablesRICH[u2]+"_BKGND";
      BkgndAgl[u2]=(TSpline3 *) _file3b->Get(nome.c_str());
   }
   cout<<_file3b<<endl;

   for(int qs=0; qs<9; qs++) cout<<Signal[qs]<<" ";
   cout<<endl;
   for(int qs=0; qs<9; qs++) cout<<BkgndNaF[qs]<<" ";
   cout<<endl;
   for(int qs=0; qs<9; qs++) cout<<SignalAgl[qs]<<" ";
   cout<<endl;

   cout<<"Vuoi Ntuple? (1-Sì;2-No)"<<endl;
   //cin>>scelta;
   scelta=1;
   cout<<"**************************** R BINS ***********************************"<<endl;
   for(int i=0; i<44; i++) {
      float temp=i+14;
      bin[i]=0.1*pow(10,temp/(9.5*2));
      //************** bin DAV
      /*        bin[i]=exp(E);
      	  E=E+a;*/
      cout<<bin[i]<<endl;

   }

   cout<<"**************************** BETA BINS ***********************************"<<endl;
   float B=0.4;
   float B1=0;
   float B2=0;
   float E=0.1;
   int binnum=1;
   float a=(log(0.9)-log(0.1))/18;
   float E2=exp(log(0.1)+1.5*a);
   while(B<0.825) {
      B=B+2*(pow(B,2)*beta->Eval(B));
      E=exp(log(0.1)+binnum*a);
      E2=exp(log(0.1)+(binnum+0.5)*a);
      B1=sqrt(1-1/(pow(E+1,2)));
      B2=sqrt(1-1/(pow(E2+1,2)));
      cout<<B<<" "<<binnum<<" "<<B1<<" "<<B2<<endl;
      binnum++;
   }




   TFile *file =TFile::Open("/storage/gpfs_ams/ams/users/fdimicco/MAIN/sommaMC/sommaMCB800.root");
   TTree *geo_stuff = (TTree *)file->Get("parametri_geo");
   string ARGV(argv[1]);
   string indirizzo_out="/storage/gpfs_ams/ams/users/fdimicco/Deutons/Risultati/"+calib+"/RisultatiMC_"+ARGV+".root";
   TFile * File = new TFile(indirizzo_out.c_str(), "RECREATE");
   TNtuple *grandezzequal = new TNtuple("grandezzequal","grandezzequal","Velocity:MC_type:R:NAnticluster:Clusterinutili:DiffR:fuoriX:layernonusati:Chisquare:Richtotused:RichPhEl:Cutmask:Momentogen:DistD:IsCharge1");
   TNtuple *grandezzesepd = new TNtuple("grandezzesepd","grandezzesepd","R:Beta:EdepL1:MC_type:Cutmask:Rmin:EdepTOF:EdepTrack:EdepTOFD:Momentogen:BetaRICH_new:LDiscriminant:BDT_response:Dist5D:Dist5D_P");
   TNtuple * pre = new TNtuple("pre","distr for giov","R:Beta:EdepL1:EdepTOFU:EdepTrack:EdepTOFD:EdepECAL:MC_type:Momentogen:Betagen:Dist5D:Dist5D_P:BetaRICH_new:Cutmask:BetanS");
   TNtuple * trig = new TNtuple("trig","trig","MC_type:Momento_gen:Ev_Num:Trig_Num:R_pre:Beta_pre:Cutmask:EdepL1:EdepTOFU:EdepTOFD:EdepTrack:BetaRICH:EdepECAL:Unbias");

   BDTreader();
   geo_stuff->SetBranchAddress("Momento_gen",&Momento_gen);
   geo_stuff->SetBranchAddress("Beta_gen",&Beta_gen);
   geo_stuff->SetBranchAddress("Massa_gen",&Massa_gen);
   geo_stuff->SetBranchAddress("Ev_Num",&Ev_Num);
   geo_stuff->SetBranchAddress("Trig_Num",&Trig_Num);
   geo_stuff->SetBranchAddress("Pres_Unbias",&Pres_Unbias);
   geo_stuff->SetBranchAddress("Preselected",&Preselected);
   geo_stuff->SetBranchAddress("Beta_pre",&Beta_pre);
   geo_stuff->SetBranchAddress("BetanS",&BetanS);
   geo_stuff->SetBranchAddress("R_pre",&R_pre);
   geo_stuff->SetBranchAddress("CUTMASK",&CUTMASK);
   geo_stuff->SetBranchAddress("trtrack_edep",&trtrack_edep);
   geo_stuff->SetBranchAddress("trtot_edep",&trtot_edep);
   geo_stuff->SetBranchAddress("Endep",&Endep);
   geo_stuff->SetBranchAddress("BetaRICH_new",&BetaRICH_new);
   geo_stuff->SetBranchAddress("RICHmask_new",&RICHmask_new);
   geo_stuff->SetBranchAddress("EdepECAL",&EdepECAL);
   geo_stuff->SetBranchAddress("Unbias",&Unbias);
   geo_stuff->SetBranchAddress("Beta_gen",&Beta_gen);
   geo_stuff->SetBranchAddress("Massa_gen",&Massa_gen);
   geo_stuff->SetBranchAddress("layernonusati",&layernonusati);
   geo_stuff->SetBranchAddress("NAnticluster",&NAnticluster);
   geo_stuff->SetBranchAddress("NTRDSegments",&NTRDSegments);
   geo_stuff->SetBranchAddress("NTofClusters",&NTofClusters);
   geo_stuff->SetBranchAddress("NTofClustersusati",&NTofClustersusati);
   geo_stuff->SetBranchAddress("Rup",&Rup);
   geo_stuff->SetBranchAddress("Rdown",&Rdown);
   geo_stuff->SetBranchAddress("R",&R);
   geo_stuff->SetBranchAddress("Chisquare",&Chisquare);
   geo_stuff->SetBranchAddress("ResiduiX",&ResiduiX);
   geo_stuff->SetBranchAddress("ResiduiY",&ResiduiY);
   geo_stuff->SetBranchAddress("chiq",&chiq);
   geo_stuff->SetBranchAddress("R_",&R_);
   geo_stuff->SetBranchAddress("Beta",&Beta);
   geo_stuff->SetBranchAddress("Betacorr",&Betacorr);
   geo_stuff->SetBranchAddress("NTrackHits",&NTrackHits);
   geo_stuff->SetBranchAddress("clustertrack",&clusterTrack);
   geo_stuff->SetBranchAddress("clustertottrack",&clustertotTrack);
   geo_stuff->SetBranchAddress("zonageo",&zonageo);
   geo_stuff->SetBranchAddress("Massa",&Massa);
   geo_stuff->SetBranchAddress("Richtotused",&Richtotused);
   geo_stuff->SetBranchAddress("RichPhEl",&RichPhEl);

   int giov=0;
   int nprotoni=0;
   int events =geo_stuff->GetEntries();
   int INDX=atoi(argv[1]);
   cout<<endl;
   cout<<INDX<<endl;
   cout<<endl;
   float entries=0;
   cout<<"Eventi: "<<events<<endl;
   //entries=entries/5;
   for(int i=(events/100)*INDX; i<(events/100)*(INDX+1); i++) {
      if(i%1300==0) cout<<i/(float)events*100<<"%"<<endl;
      geo_stuff->GetEvent(i);
      if(Massa_gen>1.8569&&Massa_gen<1.8571) Massa_gen=1.8570;
      if(Trig_Num<Trig) TotalTrig=TotalTrig+(double)Trig;
      totaltrig=TotalTrig+Trig;
      Trig=Trig_Num;
      MC_type = AssignMC_type(Massa_gen);
      R_corr=R;
      Particle_ID=0;
      Cutmask=CUTMASK;
      Cutmask=CUTMASK|(1<<10);
      Cutmask = Cutmask|(RICHmask_new<<11);
      if(!(((Cutmask&187)==187))) continue;
      entries++;
      if(Unbias==1) continue;
      if (Quality(geo_stuff,i)) {
         giov++;
         if(scelta==1) aggiungiantupla(geo_stuff,i,pre,0);
         //cout<<Beta<<" "<<BetanS<<endl;
         if(control==1) continue;
         //////////////////// CORR EDEP /////////////////////
         EdepTOFU=((EdepTOFU)*Corr_TOFU->Eval(Beta));
         EdepTOFD=((EdepTOFD)*Corr_TOFD->Eval(Beta));
         EdepTrack=((EdepTrack)*Corr_Track->Eval(Beta));
         /////////////////////////////////////////////////////
         //////////////// MATRICE DI RISPOSTA ///////////////
         if(Massa_gen<1&&Massa_gen>0.5) {
            for(int I=0; I<44; I++)
               if(fabs(Momento_gen)<bin[I+1]&&fabs(Momento_gen)>bin[I])
                  for(int J=0; J<44; J++)
                     if(fabs(R_corr)<bin[J+1]&&fabs(R_corr)>bin[J])  response[J][I]++;
            for(int I=0; I<44; I++)
               if(fabs(Momento_gen)<bin[I+1]&&fabs(Momento_gen)>bin[I]) norm[I]++;
         }
         ////////////////////////////////////////////////////
         Protoni(geo_stuff,i);
         if (Deutoni(geo_stuff,i)) if(scelta==1) Grandezzesepd(geo_stuff,i,grandezzesepd);
      }
      if(control==1) continue;
      if(scelta==1) Grandezzequal(geo_stuff,i,grandezzequal);
   }

   //ANALISI SET COMPLETO
   cout<<"Trigger tot: "<<totaltrig<<endl;
   TotalTrig=0;
   Trig=0;
   int z;
   for(z=(events/100)*INDX; z<(events/100)*(INDX+1); z++) {
      geo_stuff->GetEvent(z);
      if(Massa_gen>1.8569&&Massa_gen<1.8571) Massa_gen=1.8570;
      if(Trig_Num<Trig)TotalTrig=TotalTrig+(double)Trig;
      totaltrig2=TotalTrig+Trig;
      Trig=Trig_Num;
      MC_type = AssignMC_type(Massa_gen);
      if(z%100000==0) cout<<z/(float)events*100<<" "<<totaltrig2<<endl;

      Cutmask=CUTMASK;
      Cutmask=CUTMASK|(1<<10);
      Cutmask = Cutmask|(RICHmask_new<<11);

      EdepTrack=0;
      EdepTOFU=((*Endep)[0]+(*Endep)[1])/2;
      EdepTOFD=((*Endep)[2]+(*Endep)[3])/2;
      EdepTOFU=((EdepTOFU)*Corr_TOFU->Eval(Beta));
      EdepTOFD=((EdepTOFD)*Corr_TOFD->Eval(Beta));
      for(int layer=1; layer<8; layer++) EdepTrack+=(*trtot_edep)[layer];
      EdepTrack=EdepTrack/7;
      EdepTrack=((EdepTrack)*Corr_Track->Eval(Beta));

      if(scelta==1) Trigg(geo_stuff,z,trig);

   }




   cout<<"Eventi Tot: "<<z<<endl;
   cout<<"Preselezionate Tot: "<<entries<<endl;
   cout<<endl;
   cout<<"selezioni di qualità: "<<giov<<endl;
   cout<<"N. Protoni: "<<nprotoni<<endl;
   cout<<endl;
   cout<<"-----------------------"<<endl;
   cout<<"-----EFF. SELEZIONI (risp a PRE)----"<<endl;
   cout<<"Taglio anticluster: "<<f/(float)entries*100<<endl;
   cout<<"Taglio segmenti TRD: "<<g/(float)entries*100<<endl;
   cout<<"Taglio TOF Clusters inutilizzati: "<<h/(float)entries*100<<endl;
   cout<<"Taglio Confronto Rup/Rdown: "<<l/(float)entries*100<<endl;
   cout<<"Taglio Prob. Q: "<<m/(float)entries*100<<endl;
   cout<<"Taglio Fit multipli: "<<n/(float)entries*100<<endl;
   cout<<"Taglio HitX: "<<o/(float)entries*100<<endl;
   cout<<"Taglio Chi-Quadro: "<<r/(float)entries*100<<endl;
   cout<<"Taglio carica TOF: "<<t/(float)entries*100<<endl;
   cout<<"Taglio carica Track: "<<u/(float)entries*100<<endl;
   cout<<"Taglio carica TRD: "<<v/(float)entries*100<<endl;
   cout<<"Taglio layer non usati: "<<s/(float)entries*100<<endl;
   cout<<endl;
   cout<<"-----SELEZIONI DI REIEZIONE----"<<endl;
   cout<<"Cluster TOF: "<<b/(float)giov*100<<endl;
   cout<<endl;
   cout<<"Cluster Track: "<<d/(float)giov*100<<endl;
   cout<<endl;
   cout<<"Cluster TRD: "<<e/(float)giov*100<<endl;
   cout<<endl;

   cout<<"-----ANALISI EFFICIENZA (PROTONI)------------- "<<endl;
   cout<<endl;
   for(int i=0; i<43; i++)
      if(efficienzagen[i]!=0)
         cout<<"efficienza da"<<bin[i]<<" a " << bin[i+1]<<": "<<protoni[i][0]<<" nel bin " <<i<<" "<<( protoni[i][0]/(double)efficienzagen[i]*100)<<"%"<< " rispetto a preselezione:  "<<protoni[i][0]/(double)preselezionate[i][0]*100<<"%"<<endl;

      else cout<<"efficienza da "<<bin[i]<<" a "<< bin[i+1]<<": nessun evento generato nel bin"<<endl;
   cout<<endl;
   cout<<"-----ANALISI EFFICIENZA (DEUTONI)------------- "<<endl;
   cout<<endl;
   for(int i=0; i<43; i++)
      if(efficienzagen_D[i]!=0)
         cout<<"efficienza da"<<bin[i]<<" a " << bin[i+1]<<": "<<deutoni[i][0]<<" nel bin " <<i<<" "<<( deutoni[i][0]/(double)efficienzagen_D[i][0]*100)<<"%"<<endl;

      else cout<<"efficienza da "<<bin[i]<<" a "<< bin[i+1]<<": nessun evento generato nel bin"<<endl;
   cout<<endl;
   cout<<"----ANALISI FONDO----"<<endl;
   cout<<endl;
   for(int i=0; i<43; i++) {
      if(efficienzagen[i]!=0&&background[i]>0)
         cout<<"Fondo da"<<bin[i]<<" a " << bin[i+1]<<": "<<background[i]<<" nel bin " <<i<<" su "<< efficienzagen[i]<<" : "<<(background[i]/(double)efficienzagen[i]*100)<<"% ---> R.P. : "<< protoni[i][0]/background[i]<<endl;

      if(efficienzagen[i]!=0&&background[i]==0)
         cout<<"Fondo da"<<bin[i]<<" a " << bin[i+1]<<": "<<background[i]<<" nel bin " <<i<<" su "<<efficienzagen[i]<<" : "<<"Minore di: "<<(1/(double)efficienzagen[i]*100)<<"% ---> R.P. Maggiore di: "<<protoni[i][0]<<endl;

      if(efficienzagen[i]==0) cout<<"Fondo da "<<bin[i]<<" a "<< bin[i+1]<<": nessun evento generato nel bin"<<endl;
   }

   cout<<"------------------- MATRICE DI RISPOSTA -------------"<<endl;
   for(int j=0; j<43; j++) {
      for(int i=0; i<43; i++)
         cout<<response[j][i]<<" ";
      cout<<endl;
   }

   cout<<"------------------- NORM -------------"<<endl;
   for(int i=0; i<43; i++)
      cout<<norm[i]<<" ";
   cout<<endl;

   if(scelta==1) File->Write();
   File->Close();
   return 1;
}


void Trigg (TTree *albero,int i,TNtuple *ntupla)
{
   albero->GetEvent(i);
   ntupla->Fill(MC_type,Momento_gen,Ev_Num,Trig_Num,R_pre,Beta_pre,Cutmask,(*trtrack_edep)[0],EdepTOFU,EdepTOFD,EdepTrack,BetaRICH_new,EdepECAL,Unbias);

}

void aggiungiantupla (TTree *albero,int i,TNtuple *ntupla,int P_ID)
{
   albero->GetEvent(i);
   ntupla->Fill(R,Beta,(*trtrack_edep)[0],EdepTOFU,EdepTrack,EdepTOFD,EdepECAL,MC_type,Momento_gen,Beta_gen,Dist5D,Dist5D_P,BetaRICH_new,Cutmask,BetanS);

}

void Grandezzequal (TTree *albero,int i,TNtuple *ntupla)
{
   albero->GetEvent(i);
   ntupla->Fill(Velocity,MC_type,R,NAnticluster,NTofClusters-NTofClustersusati,fabs(Rup-Rdown)/R,fuoriX,layernonusati,Chisquare,Richtotused,RichPhEl,Cutmask,Momento_gen,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),IsCharge1);
}

void Grandezzesepd (TTree *albero,int i,TNtuple *ntupla)
{
   albero->GetEvent(i);
   ntupla->Fill(R,Beta,(*trtrack_edep)[0],MC_type,Cutmask,Rmin,EdepTOFU,EdepTrack,EdepTOFD,Momento_gen,BetaRICH_new,LDiscriminant,BDT_response,Dist5D,Dist5D_P);
}

int AssignMC_type(float Massa_gen)
{
   int MC_type=0;
   //protons MC
   int cursor=0;
   if(Massa_gen<1)	{
      if(Massa_gen<0.93805) MC_type = MC_type|(1<<(cursor+0));
      if(Massa_gen<0.939  ) MC_type = MC_type|(1<<(cursor+1));
      MC_type = MC_type|(1<<(cursor+2));
   }
   //deutons MCs
   cursor=8;
   if(Massa_gen>1&&Massa_gen<2) {
      MC_type = MC_type|(1<<(cursor+0));
      int moffset=18570;
      int c_s_type=(int)(10000*Massa_gen-moffset + 2);
      MC_type = MC_type|(1<<(cursor+c_s_type));
   }
   //Helium MC
   cursor=16;
   if(Massa_gen>3)	{
      MC_type = MC_type|(1<<(cursor+0));
      MC_type = MC_type|(1<<(cursor+2));
   }
   return	MC_type;
}
