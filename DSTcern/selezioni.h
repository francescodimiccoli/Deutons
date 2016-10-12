#include <root.h>


void selezioni();


AMSEventR* minimumbias(AMSEventR * ev,int fit);
void minimumbiasTOF(AMSEventR *ev,bool n[]);
bool minimumbiasTRIGG(AMSEventR *ev);
bool minimumbiasTRD(AMSEventR* ev);
bool minimumbiasTRACKER(AMSEventR *pev,int fit);



AMSEventR* golden(AMSEventR * ev, int s, int fit); 
bool goldenTRACKER(AMSEventR *ev,int s,int fit);
void goldenTOF(AMSEventR *ev,int s,int fit,bool m[]);
bool goldenTRD(AMSEventR* ev, int s,int fit);



bool InTrdAcceptance(AMSEventR* ev, int s,int fit);

AMSEventR* preselect ( AMSEventR* ev); 
AMSEventR* giovacchiniRICH (AMSEventR* ev);
		
int RichQual_new(RichRingR *prich);
int RichQual (AMSEventR* ev) ;


