#include <iostream>
#include <vector>
#include <random>
#include <string>

#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TGraph.h"
#include "TF1.h"

int main()
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    int N = 1000000;
    int nDim = 15;
    std::map<int,double> accepted;
    double r;
    for( int d=2; d<=nDim; d++ )
	accepted[d] = 0;


    for( auto& it: accepted )
	{

	    int dim = it.first;

	    for( int i=0; i<N; i++ )
		{
		    r = 0;
		    for( int d=0; d<dim; d++ )
			r += pow( dis(gen), 2. );
		    r = sqrt(r);

		    if( r <= 1. )
			accepted[dim] = accepted[dim] + 1.;
		}
	    accepted[dim] *= 100./(double)N;
	    
	    std::cout << "Acceptance for Dim-" << dim << ": " << accepted[dim] << std::endl;
	}

    TGraph* gr = new TGraph();
    gr->SetName("Acceptance-vs-Dim");
    gr->SetTitle("Acceptance-vs-Dim");
    gr->GetXaxis()->SetTitle("Dimension");
    gr->GetYaxis()->SetTitle("Acceptance [%]");
    gr->SetMarkerStyle(22);
    
    for( auto& it: accepted )
	gr->AddPoint( it.first, it.second );


    TGraph* gr_th = new TGraph();
    gr_th->SetName("TheoreticalAcceptance-vs-Dim");
    gr_th->SetTitle("TheoreticalAcceptance-vs-Dim");
    for( auto& it: accepted )
	{
	    double d = it.first;
	    // Volume of a n-dim sphere, divided by volume of n-dimensional cube
	    double y = 100. * pow( M_PI, 0.5 * d ) / TMath::Gamma( 0.5 * d + 1. ) / pow( 2., d );
	    gr_th->AddPoint( d, y );
	}
    gr_th->SetMarkerStyle(23);
    gr_th->SetMarkerColor(2);
    
    gStyle->SetLabelFont(132,"XY");
    gStyle->SetTitleFont(132,"XY");// Set title font for x-y axes
    gStyle->SetTitleFont(132,"");// Set title font for all TPAds
    gStyle->SetTitleFontSize(0.1);
    gStyle->SetFrameLineWidth(0);
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->cd();
    gr->Draw("AP");
    gr_th->Draw("P SAME");

    
    app->Run(kTRUE);
    
    return 0;
}
