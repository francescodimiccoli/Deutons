//Author P. Zuccon -- MIT
#include "VPreScaler.h"


void VPreScaler::node::PrintNode(int rec) {
  printf("Node %s Prescaling: %ld ", name.c_str(), ps);
  if (pso >= 0) printf(" Ovveridden to: %ld", pso);
  printf("\n");
  if (dr) {
    printf("   %s -> %s;\n", name.c_str(), dr->name.c_str());
    if (rec)dr->PrintNode(rec);
  }
  if (dl) {
    printf("   %s -> %s;\n", name.c_str(), dl->name.c_str());
    if (rec)dl->PrintNode(rec);
  }
}



void VPreScaler::node::PrintfNode(FILE* ff, int flag) {
  std::string color = "black";
  if (!dr && !dl) color = "blue";
  if (flag == 0) {
    fprintf(ff, " \"%s\" [  label =<<table border=\"0\" cellborder=\"0\" cellpadding=\"3\" bgcolor=\"white\"><tr><td align=\"center\" port=\"r0\"><font color=\"%s\"> %s </font></td></tr><tr><td align=\"left\" port=\"r1\"><font color=\"%s\">  PS: %ld </font> </td></tr>", name.c_str(), color.c_str(), basename.c_str(), color.c_str(), ps);
    if (pso >= 0) fprintf(ff, "<tr><td align=\"center\" port=\"r2\"><font color=\"red\">Overridden: %ld </font></td></tr>", pso);
    fprintf(ff, "</table>>];\n");
  }
  if (dr) {
    if (flag == 1)
      fprintf(ff, "   %s -> %s [fontcolor=\"blue\" label=\"%ld\"];\n", name.c_str(), dr->name.c_str(), psdr);
    dr->PrintfNode(ff, flag);
  }
  if (dl) {
    if (flag == 1)
      fprintf(ff, "   %s -> %s[fontcolor=\"blue\" label=\"%ld\"] ;\n", name.c_str(), dl->name.c_str(), psdl);
    dl->PrintfNode(ff, flag);
  }
  return;
}


VPSCategory* VPreScaler::FindCat(long int evpatt) const {
  auto it = categ.begin();
  for (; it != categ.end(); it++) {
    vpbitset test = (evpatt & it->second->catmask.to_ulong());
    if (it->first == test.to_ulong())  break;
  }
  if (it == categ.end()) return 0;
  else return it->second;
}



std::string  VPreScaler::GetName(vpbitset _cid, vpbitset _cmas) const {
  std::string nameout;
  int first = 1;
  for (int ii = 0; ii < cond.size(); ii++) {
    std::string namecut = "(";
    namecut += cond[ii].name;
    namecut += ")";
    if (!first) namecut += "&&";
    if (_cmas[ii]) {
      first = 0;
      if (_cid[ii]) {
        nameout.insert(0, " " + namecut);
      } else {
        nameout.insert(0, "!" + namecut);
      }
    }
  }
  return nameout;
}

long int   VPreScaler::GetPrs(vpbitset _cid, vpbitset _cmas) const {
  long int prsout = 1;
  for (int ii = 0; ii < cond.size(); ii++) {
    if (_cmas[ii]) {
      if (_cid[ii]) {
        prsout *= cond[ii].prsFactor1;
      } else {
        prsout *= cond[ii].prsFactor0;
      }
    }
  }
  return prsout;
}


void VPreScaler::AddCategory(vpbitset cid, vpbitset cmas, long int override) {
  std::string name = GetName( cid, cmas);
  vpbitset cat2 = cid ;
  if (categ.find(cat2.to_ulong()) != categ.end()) {
    printf("VPreScaler::AddCategory-E- the condition that you are entering %s %s  overlap with a previous one!\n", cid.to_string().c_str(), cmas.to_string().c_str());
    categ[cat2.to_ulong()]->Print();
  }


  categ[cat2.to_ulong()] = new   VPSCategory(name, cid, cmas, GetPrs(cid, cmas));

  if (override >= 0)
    categ[cat2.to_ulong()]->OverridePrf(override);
}



bool  VPreScaler::CheckConsistency(int verbose) {
  TestConsistency = true;
  int ncond = cond.size();
  long int ii, max = pow(2, ncond);
  int nlost = 0;
  for (ii = 0; ii < max; ii++) {
    vpbitset aa = ii;
    int found = 0;
    for (auto it = categ.begin(); it != categ.end(); it++) {
      vpbitset test =  (aa & it->second->catmask);
      if (it->first == test.to_ulong()) found++;
    }
    if (verbose || found>1) printf("VPreScaler::CheckConsistency-I-  Combination  %s found %d times\n", aa.to_string().c_str(), found);
    if(found==0) {printf("VPreScaler::CheckConsistency-W- Combination  %s  NOT FOUND!\n", aa.to_string().c_str()); nlost++;}
  }
  if (nlost == 0) Consistency = true;
  return (nlost == 0);

}


long int  VPreScaler::PrescaleEvent(const VPSEV& ev)  {
  long int  evpatt = 0;
  long int keep = -1;
  if (!TestConsistency) CheckConsistency();
  for (int ii = 0; ii < cond.size(); ii++)
    if (cond[ii].Eval(ev))
      evpatt |= 1 << ii;

  VPSCategory * aa = FindCat(evpatt);
  if (!aa) {
    printf("\n\nPrescaler::PrescaleEvent-E- Event has unknown category: %ld!  This should never happen !!\n\n ",evpatt);
    keep=evpatt*-1;
    return keep;
  }
  aa->Scaler++;
  if( aa->GetPrf()==0){
    aa->ScalerLost++;
    keep=evpatt*-1;
  }else
  if (aa->Scaler >= aa->GetPrf()) {
    aa->Scaler = 0;
    aa->ScalerGot++;
    keep = evpatt;
  } else
    aa->ScalerLost++;

  return keep;
}

void VPreScaler::BuildTree() {
  int ctop = cond.size() - 1;
  node* top = new node("TOP", "TOP", cond[ctop].prsFactor0, cond[ctop].prsFactor1);
  for (auto it = categ.begin(); it != categ.end(); it++) {
    const VPSCategory* cc = it->second;
    node* ref = top;
    std::string nlev;
    char nnam[100];
    char nnam2[100];

    //printf("\n%s %s %s\n",cc.name.c_str(),cc.catid.to_string().c_str(),cc.catmask.to_string().c_str());
    for (int kk = cond.size() - 1; kk >= 0; kk--) {
      if (cc->catmask[kk]) {
        if (cc->catid[kk]) {
          nlev += "1";
          if (!ref->dl) {
            sprintf(nnam, "%s_%s", cond[kk].name.c_str(), nlev.c_str());
            sprintf(nnam2, "(%s)", cond[kk].name.c_str());
            ref->AddNodeLeft( new node(nnam, nnam2, cond[kk - 1].prsFactor0, cond[kk - 1].prsFactor1));
            //printf("addning left %s\n",nnam);
          }
          ref = ref->dl;

        } else {
          nlev += "0";
          if (!ref->dr) {
            sprintf(nnam, "N_%s_%s", cond[kk].name.c_str(), nlev.c_str());
            sprintf(nnam2, "!(%s)", cond[kk].name.c_str());
            ref->AddNodeRight( new node(nnam, nnam2, cond[kk - 1].prsFactor0, cond[kk - 1].prsFactor1));
            // printf("addning Right %s\n",nnam);
          }
          ref = ref->dr;
        }
      }
      else{
	if(ref){
	  ref->psdr= cond[kk - 1].prsFactor0;
	  ref->psdl= cond[kk - 1].prsFactor1;
	}
      }
    }
    long int b1 = ref->ps;
    long int b2 = cc->GetPrf();
    if (b1 != b2) ref->pso = b2;
  }

  PrintTree(top);
  delete top;
  return;
}


void VPreScaler::PrintTree(node* top) const {

  FILE* ff = fopen("pp.dot", "w");
  fprintf(ff, "digraph GVGraph {\n");
  fprintf(ff, "   size=\"80,100\"\n");
  fprintf(ff, "   node [label=\"\\N\"]\n");

  top->PrintfNode(ff, 0);
  top->PrintfNode(ff, 1);
  fprintf(ff, "}\n");
  fclose(ff);
  return;


}
