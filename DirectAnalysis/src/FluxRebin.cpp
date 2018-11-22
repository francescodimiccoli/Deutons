#include <vector>
#include <algorithm>

#include "FluxRebin.h"
#include <TSpline.h>

#include <iostream>

std::vector<Double_t> integrate(Double_t * values, const Double_t * bins, Int_t n) {
	std::vector<Double_t> cumsum(n+1);
	cumsum[0] = 0;
	for( Int_t i = 0 ; i < n; i++) 
		cumsum[i+1] = cumsum[i] + values[i] * (bins[i+1] - bins[i]);
	return cumsum;
}


// Rebins a TGrahAsymmError, assuming that X errorbars denote bin widths.
// Y-axis errorbar values are treated properly
//
TGraphAsymmErrors * flux_rebin(TGraphAsymmErrors * g, Double_t bins[], Int_t n) 
{
	Int_t gnum = g->GetN();

	// All bin edges
	std::vector<Double_t> x;
	x.push_back( g->GetX()[0] - g->GetEXlow()[0]) ;
	for(Int_t i = 0; i < gnum; i++ ) 
		x.push_back(g->GetX()[i] + g->GetEXhigh()[i]);

	// Integrals of graph values/error
	std::vector<Double_t> y    = integrate(g->GetY     () , x.data() , gnum); 
	std::vector<Double_t> high = integrate(g->GetEYhigh() , x.data() , gnum); 
	std::vector<Double_t> low  = integrate(g->GetEYlow () , x.data() , gnum); 

	// Creating splines
	TSpline3 spline_y   ( "spline_y"    , x.data() ,    y.data() , x.size() );
	TSpline3 spline_high( "spline_high" , x.data() , high.data() , x.size() );
	TSpline3 spline_low ( "spline_low"  , x.data() ,  low.data() , x.size() );

	//Creating graph data
	std::vector<Double_t> nx, ny, exl, exh, eyl, eyh;
	for(Int_t i = 0; i < n - 1; i++) {
		Double_t a = bins[i] , b = bins[i+1];
		Double_t bcenter = (a + b) / 2;
		Double_t bwidth  = b - a;
		nx. push_back( bcenter     );
		exl.push_back( bcenter - a );
		exh.push_back( b - bcenter );        

		ny. push_back( ( spline_y.   Eval(b) - spline_y.   Eval(a) ) / bwidth );
		eyl.push_back( ( spline_high.Eval(b) - spline_high.Eval(a) ) / bwidth );
		eyh.push_back( ( spline_low. Eval(b) - spline_low. Eval(a) ) / bwidth );
	}

	return new TGraphAsymmErrors(n-1, nx.data(), ny.data(), exl.data(), exh.data(), eyl.data(), eyh.data());
}

