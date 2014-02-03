#include "TROOT.h"
#include "TFile.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include <TVector3.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cstring>
#include <vector>
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "Functions.h"
#include "Selections.h"


using namespace std;

int Ev_Num;
int Number=0;

int main(){

    gROOT->ProcessLine("#include <vector>");

for(int i=0;i<24;i++)
{
        float temp=i+7;
        bin[i]=0.1*pow(10,temp/9.5);
        //************** bin DAV
/*        bin[i]=exp(E);
        E=E+a;*/
        cout<<bin[i]<<endl;

}

//string Nome_file= "/afs/ams.cern.ch/user/fdimicco/Dimiccoli/istogrammiMC/" + file_name; 
TFile *file =TFile::Open("sommaMC.root");
TTree *cut_stuff = (TTree *)file->Get("grandezze_tagli");
TTree *geo_stuff = (TTree *)file->Get("parametri_geo");
//string Nome_Output="/afs/ams.cern.ch/user/fdimicco/Dimiccoli/Analisi/ntupleMC/Risultati_"+file_name;
TFile * File = new TFile("./risultati/RisultatiMC.root", "RECREATE");

TNtuple *grandezzequal = new TNtuple("grandezzequal","grandezzequal","R:Beta:Massa_gen:Massa:NAnticluster:Clusterinutili:DiffR:ProbQ:residuoY:residuoX:fuoriX:fuoriY:layernonusati:Chisquare");
TNtuple *grandezzesep = new TNtuple("grandezzesep","grandezzesep","R:Beta:Massa_gen:Massa:clusterTOFfuori:clusterTrackfuori:EdepTOF:EdepTOFteo:EdepTrack:EdepTrackteo:EdepTRD:EdepTRDteo");
TNtuple *grandezzesepd = new TNtuple("grandezzesepd","grandezzesepd","R:Beta:Massa_gen:Massa:clusterTOFfuori:clusterTrackfuori:EdepTOF:EdepTOFteo:EdepTrack:EdepTrackteo:EdepTRD:EdepTRDteo");

TNtuple * protons = new TNtuple("protons","distr for pre","R:Beta:Betacorr:BetaRICH:EdepTOF:CaricaTOF:EdepTRD:CaricaTRD:EdepTrack:CaricaTrack:Massa:Massagen:Momentogen:Betagen");
TNtuple * deutons = new TNtuple("deutons","distr for pre","R:Beta:Betacorr:BetaRICH:EdepTOF:CaricaTOF:EdepTRD:CaricaTRD:EdepTrack:CaricaTrack:Massa:Massagen:Momentogen:Betagen");
TNtuple * pre = new TNtuple("pre","distr for pre","R:Beta:Betacorr:BetaRICH:EdepTOF:CaricaTOF:EdepTRD:CaricaTRD:EdepTrack:CaricaTrack:Massa:Massagen:Momentogen:Betagen:Qbest");
TNtuple * qual = new TNtuple("qual","distr for giov","R:Beta:Betacorr:BetaRICH:EdepTOF:CaricaTOF:EdepTRD:CaricaTRD:EdepTrack:CaricaTrack:Massa:Massagen:Momentogen:Betagen:Qbest");
TNtuple * Background= new TNtuple("Background","distr for pre","R:Beta:Betacorr:BetaRICH:EdepTOF:CaricaTOF:EdepTRD:CaricaTRD:EdepTrack:CaricaTrack:Massa:Massagen:Momentogen:Betagen:Qbest");

geo_stuff->SetBranchAddress("Momento_gen",&Momento_gen);
geo_stuff->SetBranchAddress("Beta_gen",&Beta_gen);
geo_stuff->SetBranchAddress("Massa_gen",&Massa_gen);
geo_stuff->SetBranchAddress("Ev_Num",&Ev_Num);
cut_stuff->SetBranchAddress("Ev_Num",&Ev_Num);
cut_stuff->SetBranchAddress("Momento_gen",&Momento_gen);
cut_stuff->SetBranchAddress("Beta_gen",&Beta_gen);
cut_stuff->SetBranchAddress("Massa_gen",&Massa_gen);
cut_stuff->SetBranchAddress("CaricaTOF",&CaricaTOF);
cut_stuff->SetBranchAddress("CaricaTRD",&CaricaTRD);
cut_stuff->SetBranchAddress("CaricaTrack",&CaricaTrack);
cut_stuff->SetBranchAddress("ProbQ",&ProbQ);
cut_stuff->SetBranchAddress("Qbest",&Qbest);
cut_stuff->SetBranchAddress("Endep",&Endep);
cut_stuff->SetBranchAddress("layernonusati",&layernonusati);
cut_stuff->SetBranchAddress("NAnticluster",&NAnticluster);
cut_stuff->SetBranchAddress("NTRDSegments",&NTRDSegments);
cut_stuff->SetBranchAddress("NTofClusters",&NTofClusters);
cut_stuff->SetBranchAddress("NTofClustersusati",&NTofClustersusati);
cut_stuff->SetBranchAddress("Rup",&Rup);
cut_stuff->SetBranchAddress("Rdown",&Rdown);
cut_stuff->SetBranchAddress("R",&R);
cut_stuff->SetBranchAddress("Chisquare",&Chisquare);
cut_stuff->SetBranchAddress("ResiduiX",&ResiduiX);
cut_stuff->SetBranchAddress("ResiduiY",&ResiduiY);
cut_stuff->SetBranchAddress("chiq",&chiq);
cut_stuff->SetBranchAddress("R_",&R_);
cut_stuff->SetBranchAddress("Beta",&Beta);
cut_stuff->SetBranchAddress("Betacorr",&Betacorr);
cut_stuff->SetBranchAddress("BetaRICH",&BetaRICH);
cut_stuff->SetBranchAddress("TRDclusters",&TRDclusters);
cut_stuff->SetBranchAddress("NTRDclusters",&NTRDclusters);
cut_stuff->SetBranchAddress("endepostatrack",&endepostatrack);
cut_stuff->SetBranchAddress("NTrackHits",&NTrackHits);
cut_stuff->SetBranchAddress("clusterTrack",&clusterTrack);
cut_stuff->SetBranchAddress("EdepTRD",&EdepTRD);
cut_stuff->SetBranchAddress("zonageo",&zonageo);
cut_stuff->SetBranchAddress("Massa",&Massa);
int giov=0;
int nprotoni=0;
int entries=cut_stuff->GetEntries();
int events =geo_stuff->GetEntries();


cout<<"Eventi: "<<events<<endl;
cout<<"Preselezionate: "<<entries<<endl;
//cut_stuff->Draw("Betacorr:R>>rigvsbeta(1000,0,15,100,0,1)");
entries=11500000;
//if(entries>entriestot) entries=entriestot; 
//entries=entriestot;
for(int i=0;i<entries;i++)
{
	if(i%100000==0) cout<<i/(float)entries*100<<"%"<<endl;
	int k = cut_stuff->GetEvent(i);	
        Number=Ev_Num;
	//aggiungiantupla(cut_stuff,i,pre);
	for(int I=0;I<23;I++) if(fabs(R)<bin[I+1]&&fabs(R)>bin[I]) preselezionate[I][0]++;
	//Grandezzequal(cut_stuff,i,grandezzequal);
	if (Quality(cut_stuff,i)){
 		giov++;
		//aggiungiantupla(cut_stuff,i,qual);
		for(int J=0;J<23;J++) if(fabs(R)<bin[J+1]&&fabs(R)>bin[J]) quality[J][0]++;
                //Grandezzesep(cut_stuff,i,grandezzesep);                  			
                //Grandezzesepd(cut_stuff,i,grandezzesepd);             

		if (Protoni(cut_stuff,i)) {
					   nprotoni++;
					   //aggiungiantupla(cut_stuff,i,protons);	
					   for(int K=0;K<23;K++) if(fabs(R)<bin[K+1]&&fabs(R)>bin[K]) protoni[K][0]++; 
						}				 
		if (Deutoni(cut_stuff,i)) {
                                           	if(Massa_gen>1.5){
							//aggiungiantupla(cut_stuff,i,deutons);
					   		for(int L=0;L<23;L++) if(fabs(R)<bin[L+1]&&fabs(R)>bin[L]) deutoni[L][0]++;  }
						if(Massa_gen<1){
							//aggiungiantupla(cut_stuff,i,Background);
                                                        for(int L=0;L<23;L++) if(fabs(R)<bin[L+1]&&fabs(R)>bin[L]) background[L]++;  }

					  }	 


		
		}
}

//ANALISI EFFICIENZA DI GENERAZIONE
int z;
for(z=0;z<events;z++)
{
	int k = geo_stuff->GetEvent(z);
	if(z%100000==0) cout<<z/(float)events*100<<endl;
	for(int M=0;M<23;M++) if(fabs(Momento_gen)<bin[M+1]&&fabs(Momento_gen)>bin[M]) efficienzagen[M]++;
	if(Ev_Num==Number) break;
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
cout<<"Taglio Chi-Quadro: "<<r/(float)entries*100<<endl;
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
cout<<"Edep Track (Beta): "<<Gamma/(float)giov*100<<endl;
cout<<"Cluster Track: "<<d/(float)giov*100<<endl;
cout<<endl;
cout<<"Edep TRD (R): "<<e/(float)giov*100<<endl;
cout<<"Edep TRD (Beta): "<<delta/(float)giov*100<<endl;
cout<<endl;

cout<<"-----ANALISI EFFICIENZA------------- "<<endl;
cout<<endl;
 for(int i=0; i<23;i++)
        if(efficienzagen[i]!=0)
		cout<<"efficienza da"<<bin[i]<<" a " << bin[i+1]<<": "<<protoni[i][0]<<" nel bin " <<i<<" "<<( protoni[i][0]/(double)efficienzagen[i]*100)<<"%"<< " rispetto a preselezione:  "<<protoni[i][0]/(double)preselezionate[i][0]*100<<"%"<<endl;

        else cout<<"efficienza da "<<bin[i]<<" a "<< bin[i+1]<<": nessun evento generato nel bin"<<endl;
cout<<endl;        
cout<<"----ANALISI FONDO----"<<endl;
cout<<endl;
for(int i=0; i<23;i++) {
	 if(efficienzagen[i]!=0&&background[i]>0)
         	   cout<<"Fondo da"<<bin[i]<<" a " << bin[i+1]<<": "<<background[i]<<" nel bin " <<i<<" su "<< efficienzagen[i]<<" : "<<(background[i]/(double)efficienzagen[i]*100)<<"%"<<endl;

	  if(efficienzagen[i]!=0&&background[i]==0)
	cout<<"Fondo da"<<bin[i]<<" a " << bin[i+1]<<": "<<background[i]<<" nel bin " <<i<<" su "<<efficienzagen[i]<<" : "<<"Minore di: "<<(1/(double)efficienzagen[i]*100)<<"%"<<endl;

        if(efficienzagen[i]==0) cout<<"Fondo da "<<bin[i]<<" a "<< bin[i+1]<<": nessun evento generato nel bin"<<endl;
}

//File->Write();
File->Close();

string numero[11]={"0","1","2","3","4","5","6","7","8","9","10"};
{
FILE *fp;

        string stringa="./risultati/protoni0.dat";
        fp=fopen(stringa.c_str(),"w");
        for(int j=0;j<23;j++)
        fprintf(fp,"%7d",protoni[j][0]);
        fprintf(fp,"\n");
        fclose(fp);
}
{
FILE *fp;


        string stringa="./risultati/deutoni0.dat";
        fp=fopen(stringa.c_str(),"w");
        for(int j=0;j<23;j++)
        fprintf(fp,"%7d",deutoni[j][0]);
        fprintf(fp,"\n");
        fclose(fp);
}
{
FILE *fp;


        string stringa="./risultati/preselezionate0.dat";
        fp=fopen(stringa.c_str(),"w");
        for(int j=0;j<23;j++)
        fprintf(fp,"%7d",preselezionate[j][0]);
        fprintf(fp,"\n");
        fclose(fp);

}
{
FILE *fp;


        string stringa="./risultati/quality0.dat";
        fp=fopen(stringa.c_str(),"w");
        for(int j=0;j<23;j++)
        fprintf(fp,"%7d",quality[j][0]);
        fprintf(fp,"\n");
        fclose(fp);

}
{
FILE *fp;


        string stringa="./risultati/tempi0.dat";
        fp=fopen(stringa.c_str(),"w");
        for(int j=0;j<23;j++)
        fprintf(fp,"%7d",efficienzagen[j]);
        fprintf(fp,"\n");
        fclose(fp);

}
{
FILE *fp;


        string stringa="./risultati/background.dat";
        fp=fopen(stringa.c_str(),"w");
        for(int j=0;j<23;j++)
        fprintf(fp,"%7d",background[j]);
        fprintf(fp,"\n");
        fclose(fp);

}

return 1;
}


void aggiungiantupla (TTree *albero,int i,TNtuple *ntupla)
{
	 int k = albero->GetEvent(i);
	 ntupla->Fill(R,Beta,Betacorr,BetaRICH,((*Endep)[0]+(*Endep)[1]+(*Endep)[2]+(*Endep)[3])/4,CaricaTOF,EdepTRD/NTRDclusters,CaricaTRD,endepostatrack/NTrackHits,CaricaTrack,Massa,Massa_gen,Momento_gen,Beta_gen,Qbest);
}


void Grandezzequal (TTree *albero,int i,TNtuple *ntupla)
{
         int k = albero->GetEvent(i);

ntupla->Fill(R,Beta,Massa_gen,Massa,NAnticluster,NTofClusters-NTofClustersusati,fabs(Rup-Rdown)/R,ProbQ,residuoY,residuoX,fuoriX,fuoriY,layernonusati,Chisquare);
}


void Grandezzesep (TTree *albero,int i,TNtuple *ntupla)
{
        clusterTOFfuori=0;
        clusterTrackfuori=0;
        int k = albero->GetEvent(i);
        for (int j=0;j<4;j++)
        if (fabs((*Endep)[j]-fverap->Eval(R))>2) clusterTOFfuori++;
        float EndepTOF=((*Endep)[0]+(*Endep)[1]+(*Endep)[2]+(*Endep)[3])/4;
        for (int j=0;j<NTrackHits;j++)
        if(fabs((*clusterTrack)[j]-fveratrp->Eval(R))>40) clusterTrackfuori++;

ntupla->Fill(R,Beta,Massa_gen,Massa,clusterTOFfuori,clusterTrackfuori,EndepTOF,fverap->Eval(R),endepostatrack/NTrackHits,fveratrp->Eval(R),EdepTRD/NTRDclusters,fveraTRDp->Eval(R));
}

void Grandezzesepd (TTree *albero,int i,TNtuple *ntupla)
{
        clusterTOFfuori=0;
        clusterTrackfuori=0;
        int k = albero->GetEvent(i);
        for (int j=0;j<4;j++)
        if (fabs((*Endep)[j]-fvera->Eval(R))>2) clusterTOFfuori++;
        float EndepTOF=((*Endep)[0]+(*Endep)[1]+(*Endep)[2]+(*Endep)[3])/4;
        for (int j=0;j<NTrackHits;j++)
        if(fabs((*clusterTrack)[j]-fveratr->Eval(R))>40) clusterTrackfuori++;

ntupla->Fill(R,Beta,Massa_gen,Massa,clusterTOFfuori,clusterTrackfuori,EndepTOF,fvera->Eval(R),endepostatrack/NTrackHits,fveratr->Eval(R),EdepTRD/NTRDclusters,fveraTRD->Eval(R));
}

