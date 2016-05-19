using namespace std;

TH2F *RvsBetaTOF_P=new TH2F("RvsBetaTOF_P","RvsBetaTOF_P",500,0,6,500,0.4,1);
TH2F *RvsBetaNaF_P=new TH2F("RvsBetaNaF_P","RvsBetaNaF_P",500,1,10,500,0.75,1);
TH2F *RvsBetaAgl_P=new TH2F("RvsBetaAgl_P","RvsBetaAgl_P",500,3,15,500,0.95,1);
TH2F *RvsBetaTOF_D=new TH2F("RvsBetaTOF_D","RvsBetaTOF_D",500,0,6,500,0.4,1);
TH2F *RvsBetaNaF_D=new TH2F("RvsBetaNaF_D","RvsBetaNaF_D",500,1,10,500,0.75,1);
TH2F *RvsBetaAgl_D=new TH2F("RvsBetaAgl_D","RvsBetaAgl_D",500,3,15,500,0.95,1);
TH2F *RvsBetaTOF_He=new TH2F("RvsBetaTOF_He","RvsBetaTOF_He",500,0,6,500,0.4,1);
TH2F *RvsBetaNaF_He=new TH2F("RvsBetaNaF_He","RvsBetaNaF_He",500,1,10,500,0.75,1);
TH2F *RvsBetaAgl_He=new TH2F("RvsBetaAgl_He","RvsBetaAgl_He",500,3,15,500,0.95,1);

TH1F *MassTOF_P=new TH1F("MassTOF_P","MassTOF_P",500,0,4.5);
TH1F *MassTOF_D=new TH1F("MassTOF_D","MassTOF_D",500,0,4.5);
TH1F *MassNaF_P=new TH1F("MassNaF_P","MassNaF_P",500,0,4.5);
TH1F *MassNaF_D=new TH1F("MassNaF_D","MassNaF_D",500,0,4.5);
TH1F *MassAgl_P=new TH1F("MassAgl_P","MassAgl_P",500,0,4.5);
TH1F *MassAgl_D=new TH1F("MassAgl_D","MassAgl_D",500,0,4.5);

TH1F *MassTOF_PQ=new TH1F("MassTOF_PQ","MassTOF_PQ",500,0,4.5);
TH1F *MassTOF_DQ=new TH1F("MassTOF_DQ","MassTOF_DQ",500,0,4.5);
TH1F *MassNaF_PQ=new TH1F("MassNaF_PQ","MassNaF_PQ",500,0,4.5);
TH1F *MassNaF_DQ=new TH1F("MassNaF_DQ","MassNaF_DQ",500,0,4.5);
TH1F *MassAgl_PQ=new TH1F("MassAgl_PQ","MassAgl_PQ",500,0,4.5);
TH1F *MassAgl_DQ=new TH1F("MassAgl_DQ","MassAgl_DQ",500,0,4.5);


TH2F *EdepUTOFvsR_P=new TH2F("EdepUTOFvsR_P","EdepUTOFvsR_P",500,0,10,500,0,40);
TH2F *EdepUTOFvsR_D=new TH2F("EdepUTOFvsR_D","EdepUTOFvsR_D",500,0,10,500,0,40);
TH2F *EdepUTOFvsR_He=new TH2F("EdepUTOFvsR_He","EdepUTOFvsR_He",500,0,10,500,0,40);
TH2F *EdepLTOFvsR_P=new TH2F("EdepLTOFvsR_P","EdepLTOFvsR_P",500,0,10,500,0,40);
TH2F *EdepLTOFvsR_D=new TH2F("EdepLTOFvsR_D","EdepLTOFvsR_D",500,0,10,500,0,40);
TH2F *EdepLTOFvsR_He=new TH2F("EdepLTOFvsR_He","EdepLTOFvsR_He",500,0,10,500,0,40);
TH2F *EdepTrackvsR_P=new TH2F("EdepTrackvsR_P","EdepTrackvsR_P",500,0,10,500,0,4);
TH2F *EdepTrackvsR_D=new TH2F("EdepTrackvsR_D","EdepTrackvsR_D",500,0,10,500,0,4);
TH2F *EdepTrackvsR_He=new TH2F("EdepTrackvsR_He","EdepTrackvsR_He",500,0,10,500,0,4);


TH2F *RvsBetaTOF=new TH2F("RvsBetaTOF","RvsBetaTOF",500,0,6,500,0.4,1);
TH2F *RvsBetaNaF=new TH2F("RvsBetaNaF","RvsBetaNaF",500,1,10,500,0.75,1);
TH2F *RvsBetaAgl=new TH2F("RvsBetaAgl","RvsBetaAgl",500,3,15,500,0.95,1);


TH1F *MassTOF=new TH1F("MassTOF","MassTOF",500,0,4.5);
TH1F *MassNaF=new TH1F("MassNaF","MassNaF",500,0,4.5);
TH1F *MassAgl=new TH1F("MassAgl","MassAgl",500,0,4.5);

TH1F *MassTOFQ=new TH1F("MassTOFQ","MassTOFQ",500,0,4.5);
TH1F *MassNaFQ=new TH1F("MassNaFQ","MassNaFQ",500,0,4.5);
TH1F *MassAglQ=new TH1F("MassAglQ","MassAglQ",500,0,4.5);


TH2F *LikvsDistTOF_P=new TH2F("LikvsDistTOF_P","LikvsDistTOF_P",500,0,6,500,-1,1);
TH2F *LikvsDistNaF_P=new TH2F("LikvsDistNaF_P","LikvsDistNaF_P",500,0,6,500,-1,1);
TH2F *LikvsDistAgl_P=new TH2F("LikvsDistAgl_P","LikvsDistAgl_P",500,0,6,500,-1,1);
TH2F *LikvsDistTOF_D=new TH2F("LikvsDistTOF_D","LikvsDistTOF_D",500,0,6,500,-1,1);
TH2F *LikvsDistNaF_D=new TH2F("LikvsDistNaF_D","LikvsDistNaF_D",500,0,6,500,-1,1);
TH2F *LikvsDistAgl_D=new TH2F("LikvsDistAgl_D","LikvsDistAgl_D",500,0,6,500,-1,1);

TH2F *RvsDistTOF_P=new TH2F("RvsDistTOF_P","RvsDistTOF_P",500,0,6,500,-1,1);
TH2F *RvsDistNaF_P=new TH2F("RvsDistNaF_P","RvsDistNaF_P",500,1,10,500,-1,1);
TH2F *RvsDistAgl_P=new TH2F("RvsDistAgl_P","RvsDistAgl_P",500,2,19,500,-1,1);
TH2F *RvsDistTOF_D=new TH2F("RvsDistTOF_D","RvsDistTOF_D",500,0,6,500,-1,1);
TH2F *RvsDistNaF_D=new TH2F("RvsDistNaF_D","RvsDistNaF_D",500,1,10,500,-1,1);
TH2F *RvsDistAgl_D=new TH2F("RvsDistAgl_D","RvsDistAgl_D",500,2,19,500,-1,1);
TH2F *RvsDistTOF_He=new TH2F("RvsDistTOF_He","RvsDistTOF_He",500,0,6,500,-1,1);
TH2F *RvsDistNaF_He=new TH2F("RvsDistNaF_He","RvsDistNaF_He",500,1,10,500,-1,1);
TH2F *RvsDistAgl_He=new TH2F("RvsDistAgl_He","RvsDistAgl_He",500,2,19,500,-1,1);



TH1F *DistTOF_P=new TH1F("DistTOF_P","DistTOF_P",500,-1,1);
TH1F *DistNaF_P=new TH1F("DistNaF_P","DistNaF_P",500,-1,1);
TH1F *DistAgl_P=new TH1F("DistAgl_P","DistAgl_P",500,-1,1);
TH1F *DistTOF_D=new TH1F("DistTOF_D","DistTOF_D",500,-1,1);
TH1F *DistNaF_D=new TH1F("DistNaF_D","DistNaF_D",500,-1,1);
TH1F *DistAgl_D=new TH1F("DistAgl_D","DistAgl_D",500,-1,1);
TH1F *DistTOF_He=new TH1F("DistTOF_He","DistTOF_He",500,-1,1);
TH1F *DistNaF_He=new TH1F("DistNaF_He","DistNaF_He",500,-1,1);
TH1F *DistAgl_He=new TH1F("DistAgl_He","DistAgl_He",500,-1,1);


TH2F *sigmagen_bad =new TH2F("sigmagen_bad","sigmagen_bad",500,0,30,500,0,30);



void SlidesforPlot_Fill(TNtuple *ntupla, int l){
         ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	float Betagen= pow (pow (Momento_gen/Massa_gen,2) / (1+pow (Momento_gen/Massa_gen,2) ),0.5);
           if(Herejcut) {   	 
		if(Massa_gen<1&&Massa_gen>0.5){
                                RvsBetaTOF_P->Fill(R,Beta);
				if((((int)Cutmask)>>11)==512) RvsBetaNaF_P->Fill(R,BetaRICH);
				if((((int)Cutmask)>>11)==0) RvsBetaAgl_P->Fill(R,BetaRICH);
				if(Betastrongcut&&BetaRICH<0) MassTOF_P->Fill((R/Beta)*pow(1-pow(Beta,2),0.5));
				if(Betastrongcut&&(((int)Cutmask)>>11)==512) MassNaF_P->Fill((R/BetaRICH)*pow(1-pow(BetaRICH,2),0.5));
				if(Betastrongcut&&(((int)Cutmask)>>11)==0) MassAgl_P->Fill((R/BetaRICH)*pow(1-pow(BetaRICH,2),0.5));
				if(Likcut&&Distcut){
					if(Betastrongcut&&BetaRICH<0) MassTOF_PQ->Fill((R/Beta)*pow(1-pow(Beta,2),0.5));
                                	if(Betastrongcut&&(((int)Cutmask)>>11)==512) MassNaF_PQ->Fill((R/BetaRICH)*pow(1-pow(BetaRICH,2),0.5));
                                	if(Betastrongcut&&(((int)Cutmask)>>11)==0) MassAgl_PQ->Fill((R/BetaRICH)*pow(1-pow(BetaRICH,2),0.5));
					}
				if(Betastrongcut&&BetaRICH<0&&(R/Beta)*pow(1-pow(Beta,2),0.5)>2)
					sigmagen_bad->Fill(fabs(R-Momento_gen)/(pow(Momento_gen,2)*Rig->Eval(Momento_gen)),fabs(Beta-Betagen)/(pow(Beta,2)*beta->Eval(Beta)));
				}
                if(Massa_gen<2&&Massa_gen>1.5){
                                RvsBetaTOF_D->Fill(R,Beta);
				if((((int)Cutmask)>>11)==512) RvsBetaNaF_D->Fill(R,BetaRICH);
                                if((((int)Cutmask)>>11)==0) RvsBetaAgl_D->Fill(R,BetaRICH);
                 		if(Betastrongcut&&BetaRICH<0) MassTOF_D->Fill((R/Beta)*pow(1-pow(Beta,2),0.5));
                                if(Betastrongcut&&(((int)Cutmask)>>11)==512) MassNaF_D->Fill((R/BetaRICH)*pow(1-pow(BetaRICH,2),0.5));
                                if(Betastrongcut&&(((int)Cutmask)>>11)==0) MassAgl_D->Fill((R/BetaRICH)*pow(1-pow(BetaRICH,2),0.5));       		
				if(Likcut&&Distcut){
                                        if(Betastrongcut&&BetaRICH<0) MassTOF_DQ->Fill((R/Beta)*pow(1-pow(Beta,2),0.5));
                                        if(Betastrongcut&&(((int)Cutmask)>>11)==512) MassNaF_DQ->Fill((R/BetaRICH)*pow(1-pow(BetaRICH,2),0.5));
                                        if(Betastrongcut&&(((int)Cutmask)>>11)==0) MassAgl_DQ->Fill((R/BetaRICH)*pow(1-pow(BetaRICH,2),0.5));
                                        }
				}
                if(Massa_gen<4.5&&Massa_gen>2.5){
                                if(BetaRICH<0) RvsBetaTOF_He->Fill(R,Beta);
				if((((int)Cutmask)>>11)==512) RvsBetaNaF_He->Fill(R,BetaRICH);
                                if((((int)Cutmask)>>11)==0) RvsBetaAgl_He->Fill(R,BetaRICH);
                       		
				}
                }
	if(Massa_gen<1&&Massa_gen>0.5) {
				EdepUTOFvsR_P->Fill(R,EdepTOFU);
				EdepLTOFvsR_P->Fill(R,EdepTOFD);
				EdepTrackvsR_P->Fill(R,EdepTrack);	
				if(Betastrongcut&&BetaRICH<0) DistTOF_P->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
                                if(Betastrongcut&&(((int)Cutmask)>>11)==512) DistNaF_P->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
                                if(Betastrongcut&&(((int)Cutmask)>>11)==0) DistAgl_P->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));

                                if(BetaRICH<0) RvsDistTOF_P->Fill(R,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
                                if((((int)Cutmask)>>11)==512) RvsDistNaF_P->Fill(R,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
                                if((((int)Cutmask)>>11)==0) RvsDistAgl_P->Fill(R,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
	
		}
	if(Massa_gen<2&&Massa_gen>1.5) {
                                EdepUTOFvsR_D->Fill(R,EdepTOFU);
                                EdepLTOFvsR_D->Fill(R,EdepTOFD);
                                EdepTrackvsR_D->Fill(R,EdepTrack);
				if(Betastrongcut&&BetaRICH<0) DistTOF_D->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
                                if(Betastrongcut&&(((int)Cutmask)>>11)==512) DistNaF_D->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
                                if(Betastrongcut&&(((int)Cutmask)>>11)==0) DistAgl_D->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));

                                if(BetaRICH<0) RvsDistTOF_D->Fill(R,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
                                if((((int)Cutmask)>>11)==512) RvsDistNaF_D->Fill(R,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
                                if((((int)Cutmask)>>11)==0) RvsDistAgl_D->Fill(R,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));

                }
	if(Massa_gen<4.5&&Massa_gen>2.5) {
                                EdepUTOFvsR_He->Fill(R,EdepTOFU);
                                EdepLTOFvsR_He->Fill(R,EdepTOFD);
                                EdepTrackvsR_He->Fill(R,EdepTrack);
                		if(Betastrongcut&&BetaRICH<0&&R>1) DistTOF_He->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
                                if(Betastrongcut&&(((int)Cutmask)>>11)==512&&R>1) DistNaF_He->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
                                if(Betastrongcut&&(((int)Cutmask)>>11)==0&&R>1) DistAgl_He->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));

                                if(BetaRICH<0) RvsDistTOF_He->Fill(R,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
                                if((((int)Cutmask)>>11)==512) RvsDistNaF_He->Fill(R,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
                                if((((int)Cutmask)>>11)==0) RvsDistAgl_He->Fill(R,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D));
		}
	
		
}


void SlidesforPlot_D_Fill(TNtuple *ntupla, int l){
         ntupla->GetEvent(l);
        if(Beta<=0||R<=0) return;

           if(Herejcut&&Latitude>0.8) {
                                RvsBetaTOF->Fill(R,Beta);
				if((((int)Cutmask)>>11)==512) RvsBetaNaF->Fill(R,BetaRICH);
                                if((((int)Cutmask)>>11)==0) RvsBetaAgl->Fill(R,BetaRICH);
                		if(Betastrongcut&&BetaRICH<0) MassTOF->Fill((R/Beta)*pow(1-pow(Beta,2),0.5));
                                if(Betastrongcut&&(((int)Cutmask)>>11)==512) MassNaF->Fill((R/BetaRICH)*pow(1-pow(BetaRICH,2),0.5));
                                if(Betastrongcut&&(((int)Cutmask)>>11)==0) MassAgl->Fill((R/BetaRICH)*pow(1-pow(BetaRICH,2),0.5));
				if(Likcut&&Distcut){
                                        if(Betastrongcut&&BetaRICH<0) MassTOFQ->Fill((R/Beta)*pow(1-pow(Beta,2),0.5));
                                        if(Betastrongcut&&(((int)Cutmask)>>11)==512) MassNaFQ->Fill((R/BetaRICH)*pow(1-pow(BetaRICH,2),0.5));
                                        if(Betastrongcut&&(((int)Cutmask)>>11)==0) MassAglQ->Fill((R/BetaRICH)*pow(1-pow(BetaRICH,2),0.5));
                                        }

		}

}




void SlidesforPlot_Write(){
RvsBetaTOF_P->Write();
RvsBetaNaF_P->Write();
RvsBetaAgl_P->Write();
RvsBetaTOF_D->Write();
RvsBetaNaF_D->Write();
RvsBetaAgl_D->Write();
RvsBetaTOF_He->Write();
RvsBetaNaF_He->Write();
RvsBetaAgl_He->Write();
MassTOF_P->Write();
MassTOF_D->Write();
MassNaF_P->Write();
MassNaF_D->Write();
MassAgl_P->Write();
MassAgl_D->Write();
MassTOF_PQ->Write();
MassTOF_DQ->Write();
MassNaF_PQ->Write();
MassNaF_DQ->Write();
MassAgl_PQ->Write();
MassAgl_DQ->Write();
EdepUTOFvsR_P->Write();
EdepUTOFvsR_D->Write();
EdepUTOFvsR_He->Write();
EdepLTOFvsR_P->Write();
EdepLTOFvsR_D->Write();
EdepLTOFvsR_He->Write();
EdepTrackvsR_P->Write();
EdepTrackvsR_D->Write();
EdepTrackvsR_He->Write();
RvsBetaTOF->Write();
RvsBetaNaF->Write();
RvsBetaAgl->Write();
MassTOF->Write();
MassNaF->Write();
MassAgl->Write();
MassTOFQ->Write();
MassNaFQ->Write();
MassAglQ->Write();
DistTOF_P->Write();
DistNaF_P->Write();
DistAgl_P->Write();
DistTOF_D->Write();
DistNaF_D->Write();
DistAgl_D->Write();
DistTOF_He->Write();
DistNaF_He->Write();
DistAgl_He->Write();
RvsDistTOF_P->Write();
RvsDistNaF_P->Write();
RvsDistAgl_P->Write();
RvsDistTOF_D->Write();
RvsDistNaF_D->Write();
RvsDistAgl_D->Write();
RvsDistTOF_He->Write();
RvsDistNaF_He->Write();
RvsDistAgl_He->Write();
sigmagen_bad->Write();
}


void SlidesforPlot(TFile * file1){
TH2F * RvsBetaTOF_P=(TH2F*)file1->Get("RvsBetaTOF_P");
TH2F * RvsBetaNaF_P=(TH2F*)file1->Get("RvsBetaNaF_P");
TH2F * RvsBetaAgl_P=(TH2F*)file1->Get("RvsBetaAgl_P");
TH2F * RvsBetaTOF_D=(TH2F*)file1->Get("RvsBetaTOF_D");
TH2F * RvsBetaNaF_D=(TH2F*)file1->Get("RvsBetaNaF_D");
TH2F * RvsBetaAgl_D=(TH2F*)file1->Get("RvsBetaAgl_D");
TH2F * RvsBetaTOF_He=(TH2F*)file1->Get("RvsBetaTOF_He");
TH2F * RvsBetaNaF_He=(TH2F*)file1->Get("RvsBetaNaF_He");
TH2F * RvsBetaAgl_He=(TH2F*)file1->Get("RvsBetaAgl_He");
TH2F * EdepUTOFvsR_P=(TH2F*)file1->Get("EdepUTOFvsR_P");
TH2F * EdepUTOFvsR_D=(TH2F*)file1->Get("EdepUTOFvsR_D");
TH2F * EdepUTOFvsR_He=(TH2F*)file1->Get("EdepUTOFvsR_He");
TH2F * EdepLTOFvsR_P=(TH2F*)file1->Get("EdepLTOFvsR_P");
TH2F * EdepLTOFvsR_D=(TH2F*)file1->Get("EdepLTOFvsR_D");
TH2F * EdepLTOFvsR_He=(TH2F*)file1->Get("EdepLTOFvsR_He");
TH2F * EdepTrackvsR_P=(TH2F*)file1->Get("EdepTrackvsR_P");
TH2F * EdepTrackvsR_D=(TH2F*)file1->Get("EdepTrackvsR_D");
TH2F * EdepTrackvsR_He=(TH2F*)file1->Get("EdepTrackvsR_He");
TH1F * MassTOF_P=(TH1F*)file1->Get("MassTOF_P");
TH1F * MassTOF_D=(TH1F*)file1->Get("MassTOF_D");
TH1F * MassNaF_P=(TH1F*)file1->Get("MassNaF_P");
TH1F * MassNaF_D=(TH1F*)file1->Get("MassNaF_D");
TH1F * MassAgl_P=(TH1F*)file1->Get("MassAgl_P");
TH1F * MassAgl_D=(TH1F*)file1->Get("MassAgl_D");
TH1F * MassTOF_PQ=(TH1F*)file1->Get("MassTOF_PQ");
TH1F * MassTOF_DQ=(TH1F*)file1->Get("MassTOF_DQ");
TH1F * MassNaF_PQ=(TH1F*)file1->Get("MassNaF_PQ");
TH1F * MassNaF_DQ=(TH1F*)file1->Get("MassNaF_DQ");
TH1F * MassAgl_PQ=(TH1F*)file1->Get("MassAgl_PQ");
TH1F * MassAgl_DQ=(TH1F*)file1->Get("MassAgl_DQ");
TH2F * RvsBetaTOF=(TH2F*)file1->Get("RvsBetaTOF");
TH2F * RvsBetaNaF=(TH2F*)file1->Get("RvsBetaNaF");
TH2F * RvsBetaAgl=(TH2F*)file1->Get("RvsBetaAgl");
TH1F * MassTOF=(TH1F*)file1->Get("MassTOF");
TH1F * MassNaF=(TH1F*)file1->Get("MassNaF");
TH1F * MassAgl=(TH1F*)file1->Get("MassAgl");
TH1F * MassTOFQ=(TH1F*)file1->Get("MassTOFQ");
TH1F * MassNaFQ=(TH1F*)file1->Get("MassNaFQ");
TH1F * MassAglQ=(TH1F*)file1->Get("MassAglQ");
TH1F * DistTOF_P=(TH1F*)file1->Get("DistTOF_P");
TH1F * DistNaF_P=(TH1F*)file1->Get("DistNaF_P");
TH1F * DistAgl_P=(TH1F*)file1->Get("DistAgl_P");
TH1F * DistTOF_D=(TH1F*)file1->Get("DistTOF_D");
TH1F * DistNaF_D=(TH1F*)file1->Get("DistNaF_D");
TH1F * DistAgl_D=(TH1F*)file1->Get("DistAgl_D");
TH1F * DistTOF_He=(TH1F*)file1->Get("DistTOF_He");
TH1F * DistNaF_He=(TH1F*)file1->Get("DistNaF_He");
TH1F * DistAgl_He=(TH1F*)file1->Get("DistAgl_He");
TH2F * RvsDistTOF_P=(TH2F*)file1->Get("RvsDistTOF_P");
TH2F * RvsDistNaF_P=(TH2F*)file1->Get("RvsDistNaF_P");
TH2F * RvsDistAgl_P=(TH2F*)file1->Get("RvsDistAgl_P");
TH2F * RvsDistTOF_D=(TH2F*)file1->Get("RvsDistTOF_D");
TH2F * RvsDistNaF_D=(TH2F*)file1->Get("RvsDistNaF_D");
TH2F * RvsDistAgl_D=(TH2F*)file1->Get("RvsDistAgl_D");

TCanvas *p1 =new TCanvas("RvsBeta TOF MC");
TCanvas *p2 =new TCanvas("RvsBeta NaF MC");
TCanvas *p3 =new TCanvas("RvsBeta Agl MC");
TCanvas *p4 =new TCanvas("RvsBeta TOF H.L. data");
TCanvas *p5 =new TCanvas("RvsBeta NaF H.L. data");
TCanvas *p6 =new TCanvas("RvsBeta Agl H.L. data");
TCanvas *p7 =new TCanvas("RvsEdep Upper TOF MC");
TCanvas *p8 =new TCanvas("RvsEdep Lower TOF MC");
TCanvas *p9 =new TCanvas("RvsEdep Tracker MC");
TCanvas *p10=new TCanvas("Mass TOF H.L. data");
TCanvas *p11=new TCanvas("Mass NaF H.L. data");
TCanvas *p12=new TCanvas("Mass Agl H.L. data");
TCanvas *p13=new TCanvas("Mass TOF MC");
TCanvas *p14=new TCanvas("Mass NaF MC");
TCanvas *p15=new TCanvas("Mass Agl MC");
TCanvas *p10Q=new TCanvas("Qual. Mass TOF H.L. data");
TCanvas *p11Q=new TCanvas("Qual. Mass NaF H.L. data");
TCanvas *p12Q=new TCanvas("Qual. Mass Agl H.L. data");
TCanvas *p13Q=new TCanvas("Qual. Mass TOF MC");
TCanvas *p14Q=new TCanvas("Qual. Mass NaF MC");
TCanvas *p15Q=new TCanvas("Qual. Mass Agl MC");
TCanvas *p16=new TCanvas("Distance discr. TOF MC");
TCanvas *p17=new TCanvas("Distance discr. NaF MC");
TCanvas *p18=new TCanvas("Distance discr. Agl MC");
TCanvas *p19=new TCanvas("DistvsR  TOF MC");
TCanvas *p20=new TCanvas("DistvsR  NaF MC");
TCanvas *p21=new TCanvas("DistvsR  Agl MC");



cout<<"******************* R vs Beta plots ******************"<<endl;

{
p1->cd();
gPad->SetGridx();
gPad->SetGridy();
RvsBetaTOF_P->SetMarkerColor(2);
RvsBetaTOF_D->SetMarkerColor(4);
RvsBetaTOF_He->SetMarkerColor(3);
RvsBetaTOF_P->SetFillColor(2);
RvsBetaTOF_D->SetFillColor(4);
RvsBetaTOF_He->SetFillColor(3);

RvsBetaTOF_D->SetTitle("R vs Beta TOF (MC)");
RvsBetaTOF_D->GetXaxis()->SetTitle("R [GV]");
RvsBetaTOF_D->GetYaxis()->SetTitle("Beta TOF");	
RvsBetaTOF_He->GetZaxis()->SetRangeUser(1,100);
RvsBetaTOF_D->Draw();
RvsBetaTOF_P->Draw("same");
RvsBetaTOF_He->Draw("same");
TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
leg->AddEntry(RvsBetaTOF_P,"Protons MC", "ep");
leg->AddEntry(RvsBetaTOF_D,"Deutons MC", "ep");
leg->AddEntry(RvsBetaTOF_He,"Q=1 from He fragm.", "ep");
leg->Draw("same");
}

{
p2->cd();
gPad->SetGridx();
gPad->SetGridy();
RvsBetaNaF_P->SetMarkerColor(2);
RvsBetaNaF_D->SetMarkerColor(4);
RvsBetaNaF_He->SetMarkerColor(3);
RvsBetaNaF_P->SetFillColor(2);
RvsBetaNaF_D->SetFillColor(4);
RvsBetaNaF_He->SetFillColor(3);
RvsBetaNaF_D->SetTitle("R vs Beta NaF (MC)");
RvsBetaNaF_D->GetXaxis()->SetTitle("R [GV]");
RvsBetaNaF_D->GetYaxis()->SetTitle("Beta NaF"); 
RvsBetaNaF_D->Draw();
RvsBetaNaF_P->Draw("same");
RvsBetaNaF_He->Draw("same");
TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);                
leg->AddEntry(RvsBetaTOF_P,"Protons MC", "ep");
leg->AddEntry(RvsBetaTOF_D,"Deutons MC", "ep");
leg->AddEntry(RvsBetaTOF_He,"Q=1 from He fragm.", "ep");
leg->Draw("same");
}

{
p3->cd();
gPad->SetGridx();
gPad->SetGridy();
RvsBetaAgl_P->SetMarkerColor(2);
RvsBetaAgl_D->SetMarkerColor(4);
RvsBetaAgl_He->SetMarkerColor(3);
RvsBetaAgl_P->SetFillColor(2);
RvsBetaAgl_D->SetFillColor(4);
RvsBetaAgl_He->SetFillColor(3);
RvsBetaAgl_D->SetTitle("R vs Beta Agl (MC)");
RvsBetaAgl_D->GetXaxis()->SetTitle("R [GV]");
RvsBetaAgl_D->GetYaxis()->SetTitle("Beta Agl"); 
RvsBetaAgl_D->Draw();
RvsBetaAgl_P->Draw("same");
RvsBetaAgl_He->Draw("same");
TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);                
leg->AddEntry(RvsBetaTOF_P,"Protons MC", "ep");
leg->AddEntry(RvsBetaTOF_D,"Deutons MC", "ep");
leg->AddEntry(RvsBetaTOF_He,"Q=1 from He fragm.", "ep");
leg->Draw("same");
}

{p4->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogz();
RvsBetaTOF->SetTitle("R vs Beta TOF (H.L. ISS data)");
RvsBetaTOF->GetXaxis()->SetTitle("R [GV]");
RvsBetaTOF->GetYaxis()->SetTitle("Beta TOF");
RvsBetaTOF->Draw("col");
protons->SetLineColor(2);
deutons->SetLineColor(4);
protons->SetLineWidth(3);
deutons->SetLineWidth(3);
protons->Draw("same");
deutons->Draw("same");
TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);                
leg->AddEntry(protons,"Protons Mass curve", "ep");
leg->AddEntry(deutons,"Deutons Mass curve", "ep");
leg->Draw("same");
}

{
p5->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogz();
RvsBetaNaF->SetTitle("R vs Beta NaF (H.L. ISS data)");
RvsBetaNaF->GetXaxis()->SetTitle("R [GV]");
RvsBetaNaF->GetYaxis()->SetTitle("Beta NaF");
RvsBetaNaF->Draw("col");
protons->SetLineColor(2);
deutons->SetLineColor(4);
protons->SetLineWidth(3);
deutons->SetLineWidth(3);
protons->Draw("same");
deutons->Draw("same");
TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
leg->AddEntry(protons,"Protons Mass curve", "ep");
leg->AddEntry(deutons,"Deutons Mass curve", "ep");
leg->Draw("same");
}

{
p6->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogz();
RvsBetaAgl->SetTitle("R vs Beta Agl (H.L. ISS data)");
RvsBetaAgl->GetXaxis()->SetTitle("R [GV]");
RvsBetaAgl->GetYaxis()->SetTitle("Beta Agl");
RvsBetaAgl->Draw("col");
protons->SetLineColor(2);
deutons->SetLineColor(4);
protons->SetLineWidth(3);
deutons->SetLineWidth(3);
protons->Draw("same");
deutons->Draw("same");
TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
leg->AddEntry(protons,"Protons Mass curve", "ep");
leg->AddEntry(deutons,"Deutons Mass curve", "ep");
leg->Draw("same");
}

cout<<"******************* R vs Edep plots ******************"<<endl;
p7->cd();
gPad->SetGridx();
gPad->SetGridy();
EdepUTOFvsR_P->SetMarkerColor(2);
EdepUTOFvsR_D->SetMarkerColor(4);
EdepUTOFvsR_He->SetMarkerColor(3);
EdepUTOFvsR_D->SetTitle("R vs E.dep. Upper TOF (MC)");
EdepUTOFvsR_D->GetXaxis()->SetTitle("R [GV]");
EdepUTOFvsR_D->GetYaxis()->SetTitle("E. dep. Upper TOF [MeV]");
EdepUTOFvsR_He->GetZaxis()->SetRangeUser(100,10000);
EdepUTOFvsR_D->GetZaxis()->SetRangeUser(10,1000);	
EdepUTOFvsR_P->GetZaxis()->SetRangeUser(150,10000);
EdepUTOFvsR_D->GetXaxis()->SetRangeUser(0,4);
EdepUTOFvsR_D->GetYaxis()->SetRangeUser(0,20);
EdepUTOFvsR_D->Draw();
EdepUTOFvsR_P->Draw("same");
EdepUTOFvsR_He->Draw("same");

p8->cd();
gPad->SetGridx();
gPad->SetGridy();
EdepLTOFvsR_P->SetMarkerColor(2);
EdepLTOFvsR_D->SetMarkerColor(4);
EdepLTOFvsR_He->SetMarkerColor(3);
EdepLTOFvsR_D->SetTitle("R vs E.dep. Lower TOF (MC)");
EdepLTOFvsR_D->GetXaxis()->SetTitle("R [GV]");
EdepLTOFvsR_D->GetYaxis()->SetTitle("E. dep. Upper TOF [MeV]");
EdepLTOFvsR_He->GetZaxis()->SetRangeUser(100,10000);
EdepLTOFvsR_D->GetZaxis()->SetRangeUser(10,1000);
EdepLTOFvsR_P->GetZaxis()->SetRangeUser(150,10000);
EdepLTOFvsR_D->GetXaxis()->SetRangeUser(0,4);
EdepLTOFvsR_D->GetYaxis()->SetRangeUser(0,20);
EdepLTOFvsR_D->Draw();
EdepLTOFvsR_P->Draw("same");
EdepLTOFvsR_He->Draw("same");

p9->cd();
gPad->SetGridx();
gPad->SetGridy();
EdepTrackvsR_P->SetMarkerColor(2);
EdepTrackvsR_D->SetMarkerColor(4);
EdepTrackvsR_He->SetMarkerColor(3);
EdepTrackvsR_D->SetTitle("R vs E.dep. Tracker (MC)");
EdepTrackvsR_D->GetXaxis()->SetTitle("R [GV]");
EdepTrackvsR_D->GetYaxis()->SetTitle("E. dep. Upper TOF [MeV]");
EdepTrackvsR_He->GetZaxis()->SetRangeUser(100,10000);
EdepTrackvsR_D->GetZaxis()->SetRangeUser(10,1000);
EdepTrackvsR_P->GetZaxis()->SetRangeUser(150,10000);
EdepTrackvsR_D->GetXaxis()->SetRangeUser(0,4);
EdepTrackvsR_D->GetYaxis()->SetRangeUser(0,1.22);
EdepTrackvsR_D->Draw();
EdepTrackvsR_P->Draw("same");
EdepTrackvsR_He->Draw("same");


cout<<"******************* Mass plots ******************"<<endl;

p10->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogy();
MassTOF->SetLineColor(1);
MassTOF->SetLineWidth(2);
MassTOF->GetXaxis()->SetTitle("Mass [GeV/c^2]");
MassTOF->GetXaxis()->SetTitleSize(0.045);
MassTOF->SetTitle("Mass TOF H.L. data");
MassTOF->Draw();

p11->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogy();
MassNaF->SetLineColor(1);
MassNaF->SetLineWidth(2);
MassNaF->GetXaxis()->SetTitle("Mass [GeV/c^2]");
MassNaF->SetTitle("Mass NaF H.L. data");
MassNaF->GetXaxis()->SetTitleSize(0.045);
MassNaF->Draw();

p12->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogy();
MassAgl->SetLineColor(1);
MassAgl->SetLineWidth(2);
MassAgl->GetXaxis()->SetTitle("Mass [GeV/c^2]");
MassAgl->SetTitle("Mass Agl H.L. data");
MassAgl->GetXaxis()->SetTitleSize(0.045);
MassAgl->Draw();

p13->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogy();
MassTOF_P->SetLineColor(2);
MassTOF_P->SetLineWidth(2);
MassTOF_D->SetLineColor(4);
MassTOF_D->SetLineWidth(2);
MassTOF_D->SetFillColor(4);
MassTOF_P->SetFillColor(2);
MassTOF_P->SetFillStyle(3001);
MassTOF_D->SetFillStyle(3002);
MassTOF_P->GetXaxis()->SetTitle("Mass [GeV/c^2]");
MassTOF_P->GetXaxis()->SetTitleSize(0.045);
MassTOF_P->SetTitle("Mass TOF MC");
MassTOF_P->Draw();
MassTOF_D->Draw("same");

p14->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogy();
MassNaF_P->SetLineColor(2);
MassNaF_P->SetLineWidth(2);
MassNaF_D->SetLineColor(4);
MassNaF_D->SetLineWidth(2);
MassNaF_D->SetFillColor(4);
MassNaF_P->SetFillColor(2);
MassNaF_P->SetFillStyle(3001);
MassNaF_D->SetFillStyle(3002);
MassNaF_P->GetXaxis()->SetTitle("Mass [GeV/c^2]");
MassNaF_P->SetTitle("Mass NaF MC");
MassNaF_P->GetXaxis()->SetTitleSize(0.045);
MassNaF_P->Draw();
MassNaF_D->Draw("same");

p15->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogy();
MassAgl_P->SetLineColor(2);
MassAgl_P->SetLineWidth(2);
MassAgl_D->SetLineColor(4);
MassAgl_D->SetLineWidth(2);
MassAgl_D->SetFillColor(4);
MassAgl_P->SetFillColor(2);
MassAgl_P->SetFillStyle(3001);
MassAgl_D->SetFillStyle(3002);
MassAgl_P->GetXaxis()->SetTitle("Mass [GeV/c^2]");
MassAgl_P->SetTitle("Mass Agl MC");
MassAgl_P->GetXaxis()->SetTitleSize(0.045);
MassAgl_P->Draw();
MassAgl_D->Draw("same");


cout<<"******************* Quality Mass plots ******************"<<endl;

p10Q->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogy();
MassTOFQ->SetLineColor(1);
MassTOFQ->SetLineWidth(2);
MassTOFQ->GetXaxis()->SetTitle("Mass [GeV/c^2]");
MassTOFQ->SetTitle("Mass TOF H.L. data");
MassTOFQ->GetXaxis()->SetTitleSize(0.045);
MassTOFQ->Draw();

p11Q->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogy();
MassNaFQ->SetLineColor(1);
MassNaFQ->SetLineWidth(2);
MassNaFQ->GetXaxis()->SetTitle("Mass [GeV/c^2]");
MassNaFQ->SetTitle("Mass NaF H.L. data");
MassNaFQ->GetXaxis()->SetTitleSize(0.045);
MassNaFQ->Draw();

p12Q->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogy();
MassAglQ->SetLineColor(1);
MassAglQ->SetLineWidth(2);
MassAglQ->GetXaxis()->SetTitle("Mass [GeV/c^2]");
MassAglQ->SetTitle("Mass Agl H.L. data");
MassAglQ->GetXaxis()->SetTitleSize(0.045);
MassAglQ->Draw();

p13Q->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogy();
MassTOF_PQ->SetLineColor(2);
MassTOF_PQ->SetLineWidth(2);
MassTOF_DQ->SetLineColor(4);
MassTOF_DQ->SetLineWidth(2);
MassTOF_DQ->SetFillColor(4);
MassTOF_PQ->SetFillColor(2);
MassTOF_PQ->SetFillStyle(3001);
MassTOF_DQ->SetFillStyle(3002);
MassTOF_PQ->GetXaxis()->SetTitle("Mass [GeV/c^2]");
MassTOF_PQ->GetXaxis()->SetTitleSize(0.045);
MassTOF_PQ->SetTitle("Mass TOF MC");
MassTOF_PQ->Draw();
MassTOF_DQ->Draw("same");

p14Q->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogy();
MassNaF_PQ->SetLineColor(2);
MassNaF_PQ->SetLineWidth(2);
MassNaF_DQ->SetLineColor(4);
MassNaF_DQ->SetLineWidth(2);
MassNaF_DQ->SetFillColor(4);
MassNaF_PQ->SetFillColor(2);
MassNaF_PQ->SetFillStyle(3001);
MassNaF_DQ->SetFillStyle(3002);
MassNaF_PQ->GetXaxis()->SetTitle("Mass [GeV/c^2]");
MassNaF_PQ->SetTitle("Mass NaF MC");
MassNaF_PQ->GetXaxis()->SetTitleSize(0.045);
MassNaF_PQ->Draw();
MassNaF_DQ->Draw("same");

p15Q->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogy();
MassAgl_PQ->SetLineColor(2);
MassAgl_PQ->SetLineWidth(2);
MassAgl_DQ->SetLineColor(4);
MassAgl_DQ->SetLineWidth(2);
MassAgl_DQ->SetFillColor(4);
MassAgl_PQ->SetFillColor(2);
MassAgl_PQ->SetFillStyle(3001);
MassAgl_DQ->SetFillStyle(3002);
MassAgl_PQ->GetXaxis()->SetTitle("Mass [GeV/c^2]");
MassAgl_PQ->SetTitle("Mass Agl MC");
MassAgl_PQ->GetXaxis()->SetTitleSize(0.045);
MassAgl_PQ->Draw();
MassAgl_DQ->Draw("same");

cout<<"******************* Distance discr. plots ******************"<<endl;


p16->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogy();
DistTOF_P->SetLineColor(2);
DistTOF_P->SetLineWidth(2);
DistTOF_He->SetLineWidth(2);
DistTOF_D->SetLineColor(4);
DistTOF_He->SetLineColor(3);
DistTOF_D->SetLineWidth(2);
DistTOF_D->SetFillColor(4);
DistTOF_P->SetFillColor(2);
DistTOF_He->SetFillColor(3);
DistTOF_P->SetFillStyle(3001);
DistTOF_D->SetFillStyle(3002);
DistTOF_He->SetFillStyle(3002);
DistTOF_He->GetXaxis()->SetTitle("Distance discriminant");
DistTOF_He->GetXaxis()->SetTitleSize(0.045);
DistTOF_He->SetTitle("Distance discr. distribution  TOF MC");
DistTOF_He->Draw();
DistTOF_D->Draw("same");
DistTOF_P->Draw("same");

p17->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogy();
DistNaF_P->SetLineColor(2);
DistNaF_P->SetLineWidth(2);
DistNaF_He->SetLineWidth(2);
DistNaF_D->SetLineColor(4);
DistNaF_He->SetLineColor(3);
DistNaF_D->SetLineWidth(2);
DistNaF_D->SetFillColor(4);
DistNaF_P->SetFillColor(2);
DistNaF_He->SetFillColor(3);   
DistNaF_P->SetFillStyle(3001);
DistNaF_D->SetFillStyle(3002);
DistNaF_He->SetFillStyle(3002);
DistNaF_He->GetXaxis()->SetTitle("Distance discriminant");
DistNaF_He->GetXaxis()->SetTitleSize(0.045);
DistNaF_He->SetTitle("Distance discr. distribution NaF MC");
DistNaF_He->Draw();
DistNaF_D->Draw("same");
DistNaF_P->Draw("same");

p18->cd();
gPad->SetGridx();
gPad->SetGridy();
gPad->SetLogy();
DistAgl_P->SetLineColor(2);
DistAgl_P->SetLineWidth(2);
DistAgl_He->SetLineWidth(2);
DistAgl_D->SetLineColor(4);
DistAgl_He->SetLineColor(3);
DistAgl_D->SetLineWidth(2);
DistAgl_D->SetFillColor(4);
DistAgl_P->SetFillColor(2);
DistAgl_He->SetFillColor(3);
DistAgl_P->SetFillStyle(3001);
DistAgl_D->SetFillStyle(3002);
DistAgl_He->SetFillStyle(3002);
DistAgl_He->GetXaxis()->SetTitle("Distance discriminant");
DistAgl_He->GetXaxis()->SetTitleSize(0.045);
DistAgl_He->SetTitle("Distance discr. distribution Agl MC");
DistAgl_He->Draw();
DistAgl_D->Draw("same");
DistAgl_P->Draw("same");

cout<<"******************* R vs Dist plots ******************"<<endl;
p19->cd();
gPad->SetGridx();
gPad->SetGridy();
RvsDistTOF_P->SetMarkerColor(2);
RvsDistTOF_D->SetMarkerColor(4);
RvsDistTOF_D->SetTitle("R vs Distance discr. TOF (MC)");
RvsDistTOF_D->GetXaxis()->SetTitle("R [GV]");
RvsDistTOF_D->GetYaxis()->SetTitle("Distance discr. TOF");
RvsDistTOF_P->GetZaxis()->SetRangeUser(10,4000);
RvsDistTOF_D->GetZaxis()->SetRangeUser(1,4000);
RvsDistTOF_D->Draw();
RvsDistTOF_P->Draw("same");


p20->cd();
gPad->SetGridx();
gPad->SetGridy();
RvsDistNaF_P->SetMarkerColor(2);
RvsDistNaF_D->SetMarkerColor(4);
RvsDistNaF_D->SetTitle("R vs Distance discr. NaF (MC)");
RvsDistNaF_D->GetXaxis()->SetTitle("R [GV]");
RvsDistNaF_D->GetYaxis()->SetTitle("Distance discr. NaF");
RvsDistNaF_P->GetZaxis()->SetRangeUser(5,20);
RvsDistNaF_D->GetZaxis()->SetRangeUser(1,20);
RvsDistNaF_D->Draw(); 
RvsDistNaF_P->Draw("same");

p21->cd();
gPad->SetGridx();
gPad->SetGridy();
RvsDistAgl_P->SetMarkerColor(2);
RvsDistAgl_D->SetMarkerColor(4);
RvsDistAgl_D->SetTitle("R vs Distance discr. Agl (MC)");
RvsDistAgl_D->GetXaxis()->SetTitle("R [GV]");
RvsDistAgl_D->GetYaxis()->SetTitle("Distance discr. Agl");
RvsDistAgl_P->GetZaxis()->SetRangeUser(40,400);
RvsDistAgl_D->GetZaxis()->SetRangeUser(1,400);
RvsDistAgl_D->Draw();
RvsDistAgl_P->Draw("same");



        cout<<"*** Updating Results file ***"<<endl;
        string nomefile="./Final_plots/"+mese+".root";
        TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
        f_out->mkdir("Common plots for slides");
        f_out->cd("Common plots for slides");
         p1->Write();
         p2->Write();
         p3 ->Write();
	 p4 ->Write();
         p5 ->Write();
	 p6 ->Write();
         p7 ->Write();
         p8 ->Write();
         p9 ->Write();
         p10->Write();
         p11->Write();
         p12->Write();
         p13->Write();
         p14->Write();
         p15->Write();
         p10Q->Write();
         p11Q->Write();
         p12Q->Write(); 
	 p13Q->Write();
         p14Q->Write();
         p15Q->Write();
         p16->Write();
         p17->Write();
         p18->Write();
         p19->Write();
         p20->Write();
         p21->Write();
	 f_out->Write();
	 f_out->Close();
}
