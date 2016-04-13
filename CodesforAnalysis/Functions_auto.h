#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cstring>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <cstdlib>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF2.h"
#include "TVector3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSpline.h"
#include "TFractionFitter.h"
#include "THStack.h"
#include "TNtuple.h"
#include "TObject.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"


extern const int nbinsr=43;
extern const int nbinsbeta=18;
extern const int nbinsToF=18;
extern const int nbinsNaF=18;
extern const int nbinsAgl=18;

////////////// VALORI CENTRALI BINS //////////////////
double Beta_cent[30]={0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,0.71,0.73,0.75,0.77,0.79,0.81,0.83,0.85,0.87,0.89,0.91,0.93,0.95,0.97,0.99,};
double PiccoBeta[30]={0,0.01,0.02,0.03,0.440917,0.366032,0.415973,0.455373,0.491159,0.521921,0.551876,0.580046,0.605327,0.632911,0.656168,0.681199,0.704225,0.728863,0.750751,0.769231,0.793651,0.814332,0.83612,0.859106,0.877193,0.902527,0.922509,0.943396,0.957854,0.980392,};
double betacent[30]={0,0,0,0,0.483333,0.503333,0.523333,0.543333,0.563333,0.583333,0.603333,0.623333,0.643333,0.663333,0.683333,0.703333,0.723333,0.743333,0.763333,0.783333,0.803333,0.823333,0.843333,0.863333,0.883333,0.903333,0.923333,0.943333,0.963333,0.983333,};
double valorecent[24]={0.620376,0.79053,1.28364,1.63572,2.08435,2.65604,3.38452,4.31281,5.49571,7.00304,8.9238,11.3714,14.4903,18.4646,23.5289,29.9823,38.2058,48.6846,62.0376,79.053,100.735,128.364,162.378,0,};
///////////////////////////////////////////////////////

////////////////// CORREZIONE R Mis:Gen //////////////////////
double R_gen[34]={0.5633,       0.6490, 0.8449, 1.0408, 1.1755, 1.3347, 1.5306, 1.7510, 1.9469, 2.1551, 2.3510, 2.5347, 2.7551, 2.9878, 3.2204, 3.5143, 3.7347, 3.9429, 4.1755, 4.4204, 4.6653, 4.9469, 5.1918, 5.5224, 5.7918, 6.0122, 7.0000, 8.0000, 9.0000, 10.0000,        11.0000,        20.0000,        50.0000,        100.0000};
double R_mis[34]={0.7306,       1.0000, 1.3367, 1.6902, 1.9764, 2.1279, 2.2626, 2.3805, 2.4310, 2.5320, 2.7172, 2.8182, 3.0034, 3.1886, 3.3232, 3.5926, 3.8114, 4.0303, 4.2660, 4.4680, 4.6532, 4.9226, 5.1414, 5.4613, 5.7643, 6.0337, 7.0000, 8.0000, 9.0000, 10.0000,        11.0000,        20.0000,        50.0000,        100.0000};
//////////////////////////////////////////////////////////////
////////////////// CORREZIONE RICH //////////////////////
double R_rich[25]={3.18055,     3.57491,        4.01817,        4.51638,        5.07637,        5.70579,        6.41325,        7.20844,        8.10221,        9.10681,        10.236,         11.5051,        12.9317,        14.5351,        16.3373,        18.3629,        20.6398,        23.1989,        26.0753,        29.3084,        32.9424,        37.0269,        41.6179,        46.7781,        52.5781};
double Corr_rich[25]={1.00213,  1.00045,        1.00006,        0.99960,        0.99930,        0.99921,        0.99903,        0.99895,        0.99894,        0.99888,        0.99890,        0.99887,        0.99885,        0.99882,        0.99887,        0.99883,        0.99881,        0.99883,        0.99886,        0.99885,        0.99880,        0.99888,        0.99884,        0.99883,        0.99888};
//////////////////////////////////////////////////////////////
////////////// SIGMA INVERSE //////////////////
double sigmaEL1Uinv[30]={0.967736,1.03912,1.10202,1.14133,1.23905,1.26,1.30519,1.34853,1.44947,1.46982,1.56834,1.73692,1.66074,1.80898,1.86769,2.03014,2.19262,2.2068,2.31205,2.53633,2.57257,2.74351,2.86411,3.11201,3.10839,3.29759,3.316,3.33046,3.37105,3.43659,};
double sigmaEtofUinv[30]={0.0100976,0.0110663,0.0119915,0.0130737,0.0144239,0.0158222,0.016933,0.0186202,0.0204328,0.0220403,0.0235844,0.0255532,0.0279154,0.0300694,0.0326443,0.0355004,0.0386063,0.0416717,0.0447282,0.0482057,0.0516991,0.0553607,0.0594391,0.0638552,0.0702149,0.0751617,0.0722792,0.0690407,0.0686127,0.0684623,};
double sigmaEtrackinv[30]={0.395739,0.421681,0.446726,0.472621,0.499512,0.526704,0.555722,0.577731,0.610716,0.646882,0.677966,0.710215,0.747245,0.779576,0.808805,0.841817,0.883018,0.922074,0.948715,1.00304,1.03319,1.06437,1.12453,1.15809,1.18681,1.21149,1.23365,1.24104,1.25902,1.24932,};
double sigmaETofDinv[30]={0.00840365,0.00891004,0.0100159,0.0114633,0.0126798,0.0137551,0.0157324,0.0172159,0.0189864,0.0214713,0.0235472,0.025618,0.0283304,0.0309723,0.0333629,0.0367026,0.0393653,0.0433426,0.0462506,0.0503075,0.0525681,0.0579311,0.0616576,0.0652398,0.0746341,0.0803164,0.0816279,0.0764055,0.0760838,0.0762046,};
double sigmabetainv[30]={0,0,0,0,0,0.0929367,0.056543,0.0437849,0.0379535,0.0341409,0.0316607,0.0297379,0.0288931,0.0279226,0.0275195,0.0276231,0.0277531,0.0282105,0.0286759,0.0292834,0.029913,0.0308043,0.031344,0.0321792,0.0330781,0.0337626,0.0343831,0.0352896,0.0358694,0.0360589,};
double sigmaRinv[24]={0.348281,0.199019,0.127102,0.0903468,0.0657225,0.0492223,0.0376547,0.0291273,0.0229646,0.0181335,0.0145401,0.0118079,0.00979204,0.00828369,0.00714761,0.00634645,0.00567329,0.00529938,0.00503379,0.00483305,0.00457456,0.00448186,0.00443432,0,};
///////////////////////////////////////////////////////

////////////// CURVE TEORICHE //////////////////
double EL1[30]={0.245637,0.237801,0.22797,0.219124,0.211125,0.202325,0.193408,0.186255,0.179249,0.170916,0.162978,0.156784,0.150011,0.143709,0.136972,0.130807,0.126058,0.120904,0.115025,0.110849,0.106122,0.101712,0.0966255,0.0935516,0.0900914,0.0872063,0.0861243,0.0857394,0.0856475,0.0857541,};
double ETOFU[30]={6.798,6.45514,6.11129,5.79702,5.49895,5.2031,4.93988,4.68787,4.46327,4.24062,4.04368,3.84512,3.67525,3.51634,3.36878,3.216,3.08505,2.98135,2.86505,2.75327,2.64612,2.54435,2.46038,2.35904,2.25442,2.15184,2.10065,2.0454,2.03427,2.02858,};
double ETrack[30]={0.455533,0.423873,0.395205,0.369995,0.347445,0.32806,0.310335,0.294561,0.280993,0.269133,0.257275,0.245638,0.234709,0.225833,0.216852,0.209108,0.202228,0.195626,0.189467,0.183724,0.179296,0.174042,0.168874,0.163991,0.160435,0.156167,0.154189,0.15319,0.152808,0.152689,};
double ETOFD[30]={9.42654,8.4053,7.58894,6.93,6.38088,5.91862,5.50694,5.14038,4.81146,4.52634,4.27067,4.05544,3.84172,3.66187,3.47696,3.32249,3.17592,3.04209,2.91392,2.79447,2.68009,2.59572,2.49149,2.38706,2.28299,2.21518,2.1336,2.10057,2.04695,2.02921,};
TF1 *protons = new TF1("f1","pow((pow(x,2)/pow(0.938,2)/(1 + pow(x,2)/pow(0.938,2))),0.5)",0.1,100);
TF1 *deutons = new TF1("f1","pow((pow(x,2)/pow(1.875,2)/(1 + pow(x,2)/pow(1.875,2))),0.5)",0.1,100);
///////////////////////////////////////////////////////

////////////// CORREZIONE E.DEP. MC //////////////////
double CorrL1[30]={1.02913,1.01709,1.01083,0.996921,0.987542,0.982316,0.981171,0.973789,0.963865,0.964013,0.953305,0.941808,0.946517,0.943413,0.939304,0.932319,0.927261,0.928456,0.921824,0.918713,0.917363,0.914923,0.910774,0.902295,0.907529,0.90471,0.906166,0.904332,0.899394,0.895222,};
double CorrTOFU[30]={0.92277,0.92649,0.92627,0.923377,0.927777,0.933831,0.935906,0.936918,0.943672,0.94863,0.950819,0.952972,0.958021,0.962392,0.966266,0.970695,0.975242,0.981886,0.985758,0.989501,0.994808,0.998525,1.00233,1.00368,1.00721,1.01387,1.02325,1.01933,1.01011,1.00364,};
double CorrTrack[30]={1.24315,1.20935,1.1797,1.15242,1.12857,1.11066,1.09591,1.08397,1.07584,1.07164,1.0672,1.0639,1.05909,1.05632,1.05426,1.05126,1.04966,1.0485,1.04613,1.04542,1.04371,1.04116,1.03774,1.03417,1.02951,1.02518,1.02073,1.01271,1.00428,0.998837,};
double CorrTOFD[30]={0.999484,1.04346,0.908516,0.913861,0.924209,0.926476,0.930559,0.938475,0.942968,0.951335,0.962027,0.967232,0.972022,0.976821,0.981199,0.984842,0.987938,0.992611,0.995943,0.999435,1.00331,1.00616,1.00969,1.01104,1.0145,1.02572,1.03999,1.03676,1.02362,1.00867,};
///////////////////////////////////////////////////////

////////////// DEFINIZIONE SPLINES //////////////////
TSpline3 *Rig;
TSpline3 *beta;
TF1 *betaNaF;
TF1 *betaAgl;
TSpline3 *eL1;
TSpline3 *etofu;
TSpline3 *etrack;
TSpline3 *etofd;
TSpline3 *EdepL1beta;
TSpline3 *EdepTOFbeta;
TSpline3 *EdepTrackbeta;
TSpline3 *EdepTOFDbeta;
TSpline3 *Corr_L1;
TSpline3 *Corr_TOFU ;
TSpline3 *Corr_Track;
TSpline3 *Corr_TOFD;

///////////////////////////////////////////////////////

void Functions_auto(){}
