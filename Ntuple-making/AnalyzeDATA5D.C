#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "Selections5D.h"
#include "Commonglobals.cpp"


void Controllodist (TTree *albero,int i,TNtuple *ntupla);
void Controllodistd (TTree *albero,int i,TNtuple *ntupla);
void Selez (TTree *albero,int i,TNtuple *ntupla);
int giov=0;
int nprotoni=0;
int entriestot=0;
int events=0;
int INDX;
float Encinp;
float Encind;
float esposizionep[18]= {0};
float esposizioned[18]= {0};
float esposizionepgeo[18][11]= {{0}};
float esposizionedgeo[18][11]= {{0}};
float esposizionepNaF[18]= {0};
float esposizionedNaF[18]= {0};
float esposizionepgeoNaF[18][11]= {{0}};
float esposizionedgeoNaF[18][11]= {{0}};
float esposizionepAgl[18]= {0};
float esposizionedAgl[18]= {0};
float esposizionepgeoAgl[18][11]= {{0}};
float esposizionedgeoAgl[18][11]= {{0}};
string numero[11]= {"0","1","2","3","4","5","6","7","8","9","10"};
float Rcut[11]= {18,18,16,14,12,10,8,6,4,2,1};
TH2F * unbias = new TH2F("Unbias","Unbias",43,0,43,11,0,11);
TH2F * UPreselected = new TH2F("Upreselected","Upreselected",43,0,43,11,0,11);
TH1F * UnbiasHL = new TH1F("UnbiasHL","UnbiasHL",43,0,43);
TH1F * UPreselectedHL = new TH1F("UpreselectedHL","Upreselected",43,0,43);
TH2F * selezioni_P[10];
TH2F * selected_P[10];
TH1F * selezioni_PHL[10];
TH1F * selected_PHL[10];
TH1F * selezioni_DHL[10];
TH1F * selected_DHL[10];

int main(int argc, char * argv[])
{

   string nome;
   string tagli[10]= {"Trigger","3of4 TOF","TRD Segments","Rigidity exists","Chi^2 R","Matching TOF","Matching TRD","In TRD Accept.","1 Particle","1 Tr. Track"};
   /////////// CALIBR.
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
   cout<<Rig<<" "<<beta<<" "<<betaNaF<<" "<<betaAgl<<" "<<eL1<<" "<<etofu<<" "<<etrack<<" "<<etofd<<" "<<EdepL1beta<<" "<<EdepTOFbeta<<" "<<EdepTrackbeta<<" "<<EdepTOFDbeta<<" "<<Corr_L1<<" "<<Corr_TOFU<<" "<<Corr_Track<<" "<<Corr_TOFD<<endl;

   //////////////////////// DICHIARAZIONE IST. //////////////////////
   for(int j=0; j<10; j++) {
      nome="Selezioni"+tagli[j];
      selezioni_P[j]=new TH2F(nome.c_str(),nome.c_str(),43,0,43,11,0,11);
      nome="Selected"+tagli[j];
      selected_P[j]=new TH2F(nome.c_str(),nome.c_str(),43,0,43,11,0,11);
      nome="SelezioniHL"+tagli[j];
      selezioni_PHL[j]=new TH1F(nome.c_str(),nome.c_str(),43,0,43);
      nome="SelectedHL"+tagli[j];
      selected_PHL[j]=new TH1F(nome.c_str(),nome.c_str(),43,0,43);
      nome="SelezioniHL_D"+tagli[j];
      selezioni_DHL[j]=new TH1F(nome.c_str(),nome.c_str(),43,0,43);
      nome="SelectedHL_D"+tagli[j];
      selected_DHL[j]=new TH1F(nome.c_str(),nome.c_str(),43,0,43);
      nome="SelezioniHL_MC"+tagli[j];
   }
   //////////////////////////////////////////////////////////////////
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
   float a=(log(1)-log(0.1))/18;
   float E2=exp(log(0.1)+1.5*a);
   float Betabins[18]= {0};
   float Betacent[18]= {0};
   while(B<0.825) {
      B=B+2*(pow(B,2)*beta->Eval(B));
      E=exp(log(0.1)+binnum*a);
      E2=exp(log(0.1)+(binnum+0.5)*a);
      B1=sqrt(1-1/(pow(E+1,2)));
      B2=sqrt(1-1/(pow(E2+1,2)));
      cout<<B<<" "<<binnum<<" "<<B1<<" "<<B2<<endl;
      Betabins[binnum-1]=B1;
      Betacent[binnum-1]=B2;
      binnum++;
   }

   float encinTOF[18];
   for(int i=0; i<18; i++) encinTOF[i]=1/pow(1-pow(Betabins[i],2),0.5)-1;

   cout<<"**************************** BETA BINS NaF***********************************"<<endl;
   float BetabinsNaF[18]= {0};
   float BetacentNaF[18]= {0};
   a=(log(4.025)-log(0.666))/18;
   E2=exp(log(0.666)+1.5*a);
   binnum=0;
   while(B1<0.98) {
      E=exp(log(0.666)+binnum*a);
      E2=exp(log(0.666)+(binnum+0.5)*a);
      B1=sqrt(1-1/(pow(E+1,2)));
      B2=sqrt(1-1/(pow(E2+1,2)));
      BetabinsNaF[binnum]=B1;
      BetacentNaF[binnum]=B2;
      binnum++;
   }
   float encinNaF[18];
   for(int i=0; i<18; i++) encinNaF[i]=1/pow(1-pow(BetabinsNaF[i],2),0.5)-1;

   cout<<"**************************** BETA BINS Agl***********************************"<<endl;
   float BetabinsAgl[18]= {0};
   float BetacentAgl[18]= {0};
   a=(log(9.01)-log(2.57))/18;
   E2=exp(log(2.57)+1.5*a);
   binnum=0;
   while(B1<0.995) {
      E=exp(log(2.57)+binnum*a);
      E2=exp(log(2.57)+(binnum+0.5)*a);
      B1=sqrt(1-1/(pow(E+1,2)));
      B2=sqrt(1-1/(pow(E2+1,2)));
      BetabinsAgl[binnum]=B1;
      BetacentAgl[binnum]=B2;
      binnum++;
   }
   float encinAgl[18];
   for(int i=0; i<18; i++) encinAgl[i]=1/pow(1-pow(BetabinsAgl[i],2),0.5)-1;

   string ARGV(argv[1]);
   bool check=true;
   string indirizzo_in="/storage/gpfs_ams/ams/users/fdimicco/MAIN/sommadati/sommadati"+ARGV+".root";
   TFile *file =TFile::Open(indirizzo_in.c_str());
   TTree *geo_stuff = (TTree *)file->Get("parametri_geo");
   std::cout<<geo_stuff<<endl;
   if(!(geo_stuff>0)) check=false;
   if(!(geo_stuff>0)) check=false;
   if(check==false) std::cout<<"Skipping file!"<<endl;
   string indirizzo_out="Results.root";
   TFile * File = new TFile(indirizzo_out.c_str(), "RECREATE");

   TNtuple *grandezzequal = new TNtuple("grandezzequal","grandezzequal","Velocity:Rcutoff:R:NAnticluster:Clusterinutili:DiffR:fuoriX:layernonusati:Chisquare:Richtotused:RichPhEl:Cutmask:LDiscriminant:DistD,IsCharge1");
   TNtuple *grandezzesepd = new TNtuple("grandezzesepd","grandezzesepd","R:Beta:EdepL1:Cutmask:Latitude:Rmin:EdepTOFU:EdepTrack:EdepTOFD:Rcutoff:BetaRICH_new:LDiscriminant:mcweight:Dist5D:Dist5D_P");
   TNtuple * pre = new TNtuple("Pre","distr for qual","R:Beta:EdepL1:EdepTOFU:EdepTOFD:EdepTrack:EdepECAL:Rcutoff:Latitude:Dist5D:Dist5D_P:BetaRICH_new:Cutmask:LDiscriminant");
   TNtuple * trig = new TNtuple("trig","trig","U_time:Latitude:Rcutoff:R_pre:Beta_pre:Cutmask:EdepL1:EdepTOFU:EdepTOFD:EdepTrack:BetaRICH:EdepECAL:Unbias");

   BDTreader();
   if(check) {
      geo_stuff->SetBranchAddress("U_time",&U_time);
      geo_stuff->SetBranchAddress("Latitude",&Latitude);
      geo_stuff->SetBranchAddress("zonageo",&zonageo);
      geo_stuff->SetBranchAddress("Rcutoff",&Rcutoff);
      geo_stuff->SetBranchAddress("Livetime",&Livetime);
      geo_stuff->SetBranchAddress("Unbias",&Unbias);
      geo_stuff->SetBranchAddress("Pres_Unbias",&Pres_Unbias);
      geo_stuff->SetBranchAddress("Preselected",&Preselected);
      geo_stuff->SetBranchAddress("R_pre",&R_pre);
      geo_stuff->SetBranchAddress("Beta_pre",&Beta_pre);
      geo_stuff->SetBranchAddress("CUTMASK",&CUTMASK);
      geo_stuff->SetBranchAddress("trtrack_edep",&trtrack_edep);
      geo_stuff->SetBranchAddress("trtot_edep",&trtot_edep);
      geo_stuff->SetBranchAddress("Endep",&Endep);
      geo_stuff->SetBranchAddress("BetaRICH_new",&BetaRICH_new);
      geo_stuff->SetBranchAddress("RICHmask_new",&RICHmask_new);
      geo_stuff->SetBranchAddress("EdepECAL",&EdepECAL);
      geo_stuff->SetBranchAddress("layernonusati",&layernonusati);
      geo_stuff->SetBranchAddress("NAnticluster",&NAnticluster);
      geo_stuff->SetBranchAddress("NTofClusters",&NTofClusters);
      geo_stuff->SetBranchAddress("NTofClustersusati",&NTofClustersusati);
      geo_stuff->SetBranchAddress("Rup",&Rup);
      geo_stuff->SetBranchAddress("Rdown",&Rdown);
      geo_stuff->SetBranchAddress("R",&R);
      geo_stuff->SetBranchAddress("chiq",&chiq);
      geo_stuff->SetBranchAddress("R_",&R_);
      geo_stuff->SetBranchAddress("Chisquare",&Chisquare);
      geo_stuff->SetBranchAddress("ResiduiX",&ResiduiX);
      geo_stuff->SetBranchAddress("ResiduiY",&ResiduiY);
      geo_stuff->SetBranchAddress("Beta",&Beta);
      geo_stuff->SetBranchAddress("Betacorr",&Betacorr);
      geo_stuff->SetBranchAddress("NTrackHits",&NTrackHits);
      geo_stuff->SetBranchAddress("zonageo",&zonageo);
      geo_stuff->SetBranchAddress("Massa",&Massa);
      geo_stuff->SetBranchAddress("Richtotused",&Richtotused);
      geo_stuff->SetBranchAddress("RichPhEl",&RichPhEl);
   }

   giov=0;
   nprotoni=0;
   if(check) events =geo_stuff->GetEntries();
   INDX=atoi(argv[1]);
   cout<<endl;
   cout<<INDX<<endl;
   cout<<endl;
   cout<<"Eventi: "<<events<<endl;
   int entries=0;
   if(entries>entriestot) entries=entriestot;

   for(int i=0; i<events; i++) {
      if(!check) break;
      if(i%1300==0) cout<<i/(float)events*100<<"%"<<endl;
      int k = geo_stuff->GetEvent(i);
      R_corr=R;
      if(i==0) { tbeg=U_time; cout <<"Tempo Iniziale: "<<tbeg<<endl;}
      if(i==events-1) {tend=U_time; cout <<"Tempo Finale: "<<tend<<endl;}
      for(int I=0; I<43; I++) if(R<bin[I+1]&&R>bin[I]) preselezionate[I][zonageo]++;
      Particle_ID=0;

      Cutmask=CUTMASK;
      Cutmask=CUTMASK|(1<<10);
      //Temp. rich bug fixing
      if(((Cutmask>>11)==0||(Cutmask>>11)==512)&&BetaRICH_new==-1) RICHmask_new=1;
      Cutmask = Cutmask|(RICHmask_new<<11);
      if(!(((Cutmask&187)==187))) continue;
      entries++;
      if(Unbias==1) continue;
      if (Quality(geo_stuff,i)) {
         giov++;
         if(scelta==1) aggiungiantupla(geo_stuff,i,pre,0);
         if(control==1) continue;
         Protoni(geo_stuff,i);
         if (Deutoni(geo_stuff,i)) if(scelta==1) Grandezzesepd(geo_stuff,i,grandezzesepd);

      }
      if(control==1) continue;
      if(scelta==1) Grandezzequal(geo_stuff,i,grandezzequal);
   }

   cout<<"------Calcolo Live Time-----"<<endl;
   zona=0;
   for(int j=0; j<43; j++) for(int i=0; i<11; i++) tempobingeo[j][i]=0;
   int contriniz=0;
   int zona1=0;
   int zona2=0;
   int contaeventi=0;
   int z=0;
   int Tempo=0;
   for(z=0; z<events; z++) {
      if(!check) break;
      int k = geo_stuff->GetEvent(z);
      contaeventi++;
      R=R_pre;

      if(z%1000000==0) cout<<z/(float)events*100<<"%"<<endl;
      //////// SCARTO EVENTI FINCHE' NON ENTRO IN UNA NUOVA ZONA GEOM
      zona2=zona1;
      for(int x=0; x<11; x++) {
         double geo= geomag[x]  ;
         double geo2=geomag[x+1];
         if(Latitude>geo && Latitude<=geo2)
            zona1=x;
      }

      if(zona1!=zona2&&z!=0) contriniz=1;
      if(contriniz!=1) {continue;}

      /////////////////////////////////////////////////////
      for(int i=0; i<11; i++) {
         double geo= geomag[i]  ;
         double geo2=geomag[i+1];
         if(Latitude>geo && Latitude<=geo2)
            zona=i;

         else tempozona[i]=U_time;
      }
      if(U_time!=tempozona[zona]) {
         tempozona[zona]=U_time;
         Time[zona]=Time[zona]+Livetime;
         contasecondi[zona]++;
      }

      ////////// CALCOLO TEMPO DI ESPOSIZIONE /////////////
      if(U_time!=Tempo) {
         for(int i=0; i<43; i++)
            if(fabs(1.2*Rcutoff)<bin[i+1]&&Rcutoff!=0) {
               tempobin[i]=tempobin[i]+Livetime;
               tempobingeo[i][zonageo]=tempobingeo[i][zonageo]+Livetime;
            }
         if(Rcutoff!=0) Tempo=U_time;
         for(int m=0; m<18; m++) {
            Encind=pow(((1+pow((fabs(1.2*Rcutoff)/1.875),2))),0.5)-1;
            Encinp=pow(((1+pow((fabs(1.2*Rcutoff)/0.938),2))),0.5)-1;
            if(Encinp<encinTOF[m]&&Rcutoff!=0) {
               esposizionep[m]=esposizionep[m]+Livetime;
               esposizionepgeo[m][zonageo]=esposizionepgeo[m][zonageo]+Livetime;
            }
            if(Encinp<encinNaF[m]&&Rcutoff!=0) {
               esposizionepNaF[m]=esposizionepNaF[m]+Livetime;
               esposizionepgeoNaF[m][zonageo]=esposizionepgeoNaF[m][zonageo]+Livetime;
            }
            if(Encinp<encinAgl[m]&&Rcutoff!=0) {
               esposizionepAgl[m]=esposizionepAgl[m]+Livetime;
               esposizionepgeoAgl[m][zonageo]=esposizionepgeoAgl[m][zonageo]+Livetime;
            }

            if(Encind<encinTOF[m]&&Rcutoff!=0) {
               esposizioned[m]=esposizioned[m]+Livetime;
               esposizionedgeo[m][zonageo]=esposizionedgeo[m][zonageo]+Livetime;
            }
            if(Encind<encinNaF[m]&&Rcutoff!=0) {
               esposizionedNaF[m]=esposizionedNaF[m]+Livetime;
               esposizionedgeoNaF[m][zonageo]=esposizionedgeoNaF[m][zonageo]+Livetime;
            }
            if(Encind<encinAgl[m]&&Rcutoff!=0) {
               esposizionedAgl[m]=esposizionedAgl[m]+Livetime;
               esposizionedgeoAgl[m][zonageo]=esposizionedgeoAgl[m][zonageo]+Livetime;
            }

         }
      }
      /////////////////////////////////////////////////////

      Cutmask=CUTMASK;
      Cutmask=CUTMASK|(1<<10);
      Cutmask = Cutmask|(RICHmask_new<<11);

      EdepTrack=0;
      EdepTOFU=((*Endep)[0]+(*Endep)[1])/2;
      EdepTOFD=((*Endep)[2]+(*Endep)[3])/2;

      for(int layer=1; layer<8; layer++) EdepTrack+=(*trtot_edep)[layer];
      EdepTrack=EdepTrack/7;

      if(scelta==1) Trigg(geo_stuff,z,trig);
   }
   /////////////////////////////////////////////////////

   for(int i=0; i<11; i++) {
      for(int j=0; j<43; j++) std::cout<<tempobingeo[j][i]<<" ";
      std::cout<<endl;
   }

   for(int i=1; i<11; i++)   { for(int m=0; m<18; m++) cout<<esposizionedgeo[m][i]<<" "; cout<<endl;}
   cout<<endl;
   for(int i=1; i<11; i++)   { for(int m=0; m<18; m++) cout<<esposizionepgeo[m][i]<<" "; cout<<endl;}
   cout<<"Eventi Tot: "<<events<<endl;
   cout<<"Preselezionate Tot: "<<entriestot<<endl;
   cout<<endl;
   cout<<"Eventi Analizzati: "<<contaeventi<<endl;
   cout<<"Preselezionate Analizzate: "<<entries<<endl;
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
   cout<<"Taglio Chi-Quadro traccia: "<<r/(float)entries*100<<endl;
   cout<<"Taglio carica TOF: "<<t/(float)entries*100<<endl;
   cout<<"Taglio carica Track: "<<u/(float)entries*100<<endl;
   cout<<"Taglio carica TRD: "<<v/(float)entries*100<<endl;
   cout<<"Taglio layer non usati: "<<s/(float)entries*100<<endl;
   cout<<endl;
   cout<<"-----SELEZIONI DI REIEZIONE----"<<endl;
   cout<<"Edep TOF (R): "<<a/(float)giov*100<<endl;
   cout<<"Edep TOF (Beta): "<<alpha/(float)giov*100<<endl;
   cout<<"Cluster TOF: "<<b/(float)giov*100<<endl;
   cout<<endl;
   cout<<"Edep Track (R): "<<c/(float)giov*100<<endl;
   cout<<"Edep Track (Beta): "<<delta/(float)giov*100<<endl;
   cout<<"Cluster Track: "<<d/(float)giov*100<<endl;
   cout<<endl;
   cout<<"Edep TRD (R): "<<e/(float)giov*100<<endl;
   cout<<"Edep TRD (Beta): "<<Gamma/(float)giov*100<<endl;
   cout<<endl;

   cout<<"-----ANALISI TEMPI------------- "<<endl;
   cout<<endl;
   for(int i=0; i<11; i++) {
      cout<<"Tempo trascorso in zona "<<i<<": "<<contasecondi[i];
      if(contasecondi[i]!=0)
         cout<<" sec: Percentuale di tempo vivo: "<<Time[i]/(double)contasecondi[i]*100<<"% "<<endl;
      else cout<<" sec"<<endl;
   }

   cout<<"Tempo di Volo: "<<tend-tbeg<<" sec."<<endl;
   int Contasecondi = 0;
   for (int y=0; y<11; y++)
      Contasecondi=Contasecondi+contasecondi[y];
   cout<<"Tempo di Misura: "<<Contasecondi<<" sec."<<endl;

   for(int i=0; i<43; i++)
      cout<<"Tempo di esposizione bin "<<i<<": "<<tempobin[i]<<endl;


   TH1F * Tempi=new TH1F("Tempi","Tempi",11,0,11);
   for(int i=1; i<11; i++) Tempi->SetBinContent(i,Time[i]);

   TH2F * Esposizionegeo=new TH2F("esposizionegeo","esposizionegeo",43,0,43,11,0,11);
   for(int i=1; i<11; i++) for(int j=0; j<43; j++) Esposizionegeo->SetBinContent(j+1,i,tempobingeo[j][i]);

   TH2F * Esposizionepgeo=new TH2F("esposizionepgeo","esposizionepgeo",18,0,18,11,0,11);
   for(int i=1; i<11; i++) for(int j=0; j<18; j++) Esposizionepgeo->SetBinContent(j+1,i,esposizionepgeo[j][i]);
   TH2F * EsposizionepgeoNaF=new TH2F("esposizionepgeoNaF","esposizionepgeoNaF",18,0,18,11,0,11);
   for(int i=1; i<11; i++) for(int j=0; j<18; j++) EsposizionepgeoNaF->SetBinContent(j+1,i,esposizionepgeoNaF[j][i]);
   TH2F * EsposizionepgeoAgl=new TH2F("esposizionepgeoAgl","esposizionepgeoAgl",18,0,18,11,0,11);
   for(int i=1; i<11; i++) for(int j=0; j<18; j++) EsposizionepgeoAgl->SetBinContent(j+1,i,esposizionepgeoAgl[j][i]);

   TH2F * Esposizionedgeo=new TH2F("esposizionedgeo","esposizionedgeo",18,0,18,11,0,11);
   for(int i=1; i<11; i++) for(int j=0; j<18; j++) Esposizionedgeo->SetBinContent(j+1,i,esposizionedgeo[j][i]);
   TH2F * EsposizionedgeoNaF=new TH2F("esposizionedgeoNaF","esposizionedgeoNaF",18,0,18,11,0,11);
   for(int i=1; i<11; i++) for(int j=0; j<18; j++) EsposizionedgeoNaF->SetBinContent(j+1,i,esposizionedgeoNaF[j][i]);
   TH2F * EsposizionedgeoAgl=new TH2F("esposizionedgeoAgl","esposizionedgeoAgl",18,0,18,11,0,11);
   for(int i=1; i<11; i++) for(int j=0; j<18; j++) EsposizionedgeoAgl->SetBinContent(j+1,i,esposizionedgeoAgl[j][i]);
   Tempi->Write();
   Esposizionegeo->Write();
   Esposizionepgeo->Write();
   EsposizionepgeoNaF->Write();
   EsposizionepgeoAgl->Write();
   Esposizionedgeo->Write();
   EsposizionedgeoNaF->Write();
   EsposizionedgeoAgl->Write();

   if(scelta==1) File->Write();
   if(scelta==1) File->Close();

   return 1;
}

void Trigg (TTree *albero,int i,TNtuple *ntupla)
{
   albero->GetEvent(i);
   ntupla->Fill(U_time,Latitude,Rcutoff,R_pre,Beta_pre,Cutmask,(*trtrack_edep)[0],EdepTOFU,EdepTOFD,EdepTrack,BetaRICH_new,EdepECAL,Unbias);

}


void aggiungiantupla (TTree *albero,int i,TNtuple *ntupla,int P_ID)
{
   albero->GetEvent(i);
   ntupla->Fill(R,Beta,(*trtrack_edep)[0],EdepTOFU,EdepTOFD,EdepTrack,EdepECAL,Rcutoff,Latitude,Dist5D,Dist5D_P,BetaRICH_new,Cutmask,LDiscriminant);
}

void Grandezzequal (TTree *albero,int i,TNtuple *ntupla)
{
   albero->GetEvent(i);
   ntupla->Fill(Velocity,Rcutoff,R,NAnticluster,NTofClusters-NTofClustersusati,fabs(Rup-Rdown)/R,fuoriX,layernonusati,Chisquare,Richtotused,RichPhEl,Cutmask,LDiscriminant,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),IsCharge1);
}

void Grandezzesepd (TTree *albero,int i,TNtuple *ntupla)
{
   albero->GetEvent(i);
   ntupla->Fill(R,Beta,(*trtrack_edep)[0],Cutmask,Latitude,Rmin,EdepTOFU,EdepTrack,EdepTOFD,Rcutoff,BetaRICH_new,LDiscriminant,BDT_response,Dist5D,Dist5D_P);
}


