#include "rundb.h"
#include <cstdlib>
#include <stdio.h>
#include <cstring>



int  rundb::readdb(const char * fname){

  

  FILE * fin;
  fin=fopen(fname,"r");
  if(!fin) {printf("cant'open %s\n",fname);return -1;}

  

  char* fp=0;
  size_t ssize=0;

  while(1){
    getline(&fp,&ssize,fin);
    if(feof(fin)) break;
    if(ferror(fin)) return -2;
    if(fp[0]=='#') continue;
    int st[8]  ={ 0,  11, 17,  25,  35,  45, 55, 66 };
    int len[8] ={10,   6,  7,   10,  10,  10, 11,  10 };
    char dest[8][20];
    memset(dest,0,8*20*sizeof(char));

    for (int ii=0;ii<8;ii++){
      strncpy(dest[ii],&(fp[st[ii]]),len[ii]);
      dest[ii][len[ii]]='\0';
    }

   

    rundb_el el;
    el.run=atoi(dest[0]);
    el.Rpmin=atof(dest[1]);
    el.Rpmax=atof(dest[2]);
    el.RTrig=atoi(dest[3]);
    el.pmin=atof(dest[4]);
    el.Trig=atoi(dest[6]);
    el.Events=atoi(dest[7]);
    Add(el);

  }  

  if(fp) delete fp;
  fclose(fin);
  return 1;

}
