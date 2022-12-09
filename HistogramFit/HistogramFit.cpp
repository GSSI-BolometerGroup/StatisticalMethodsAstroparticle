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
#include "TF1.h"

double GaussianPlusFlat( double *x, double* par )
{
    // par[0] = integral
    // par[1] = mean
    // par[2] = sigma
    // par[3] = flat background
    // par[4] = bin width
    // par[5] = fit range
    return par[0] / sqrt( 2. * M_PI ) / par[2] * exp( - pow( x[0] - par[1], 2. ) / 2. / pow( par[2], 2. ) ) + par[3] / par[5];
}

int main()
{

    double mu    = 2039.;// keV, Q-value of 76Ge (for neutrinoless double beta decay)
    double sigma = 1.5;// keV, typical resolution of a germanium detector
    double min   = 2000.;
    double max   = 2080.;
    int nBins    = 10. * ( max-min );// 0.1 keV bins, much thinner than energy resolution
    TH1D* histo  = new TH1D("histo","histo",nBins,min,max);

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(min,max);
    std::normal_distribution<> gaus(mu,sigma);

    int B = 10000;// Number of background events (uniformly distributed)
    int S = 10000;// Number of signal events (Gaussian distributed around mu)

    // Fill histogram with background and signal events
    for( int b=0; b<B; b++ )
	histo->Fill( dis(gen) );
    for( int s=0; s<S; s++ )
	histo->Fill( gaus(gen) );

    // Define fitting function
    TF1* func = new TF1( "func", GaussianPlusFlat, min, max, 6 );
    func->SetParameter(0,10);
    func->SetParameter(1,2040);
    func->SetParameter(2,1);
    func->SetParameter(3,10);
    func->FixParameter(4,histo->GetBinWidth(1));
    func->FixParameter(5,max-min);
    
    
    // Draw
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    histo->Draw();
    histo->Fit(func,"R");
    histo->Fit(func,"RLL");
    histo->Fit(func,"R");
    app->Run(kTRUE);
    
    return 0;
}
