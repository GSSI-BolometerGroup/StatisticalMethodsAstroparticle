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
#include "TGaxis.h"
#include "TPad.h"
#include "TFrame.h"

// Smarter implementation of hte Binomial,
// which minimizes the number of calculations and includes the log-based
// calculation of the binomial coefficient reported above.
// It has a relative precision of ~1e-15 with respect to the straightforward implementation.
double SmartBinomial( double k, double n, double p )
{
    // Binomial coefficient
    double sum_1 = 0;
    double sum_2 = 0;
    for( double i=(int)(n-k+1); i<=n; i+=1. )
	sum_1 += log(i);
    for( double i=1.; i<=k; i+=1. )
	sum_2 += log(i);
    // Exponential terms
    double exp_1 = k * log(p);
    double exp_2 = (n-k) * log(1.-p);

    return exp( sum_1 - sum_2 + exp_1 + exp_2 );
}

int main()
{
    double n = 100;
    double k = 10;
    double p = k/n;

    // Create Binomial distribution with n=100 and p=0.1
    double mink  = -0.5;
    double maxk  = 0.5+n;
    int    nBinsK = n+1;
    TH1D* binomialHisto = new TH1D( "Binomial", "Binomial", nBinsK, mink, maxk );
    binomialHisto->GetXaxis()->SetTitle("k");
    binomialHisto->GetYaxis()->SetTitle("P(k|n,p)");
    for( int b=1; b<=nBinsK; b++ )
	{
	    double x = binomialHisto->GetBinCenter(b);
	    binomialHisto->SetBinContent( b, SmartBinomial( x, n, p ) );
	}

    // Now create posterior of p, given a measured k=10
    double binWidthP = 0.01 / n;
    double minp = -0.5*binWidthP;
    double maxp = 1. + 0.5*binWidthP;
    int nBinsP = (maxp-minp)/binWidthP;
    TH1D* posterior = new TH1D( "Posterior", "Posterior", nBinsP, minp, maxp );
    posterior->GetXaxis()->SetTitle("p");
    posterior->GetYaxis()->SetTitle("P(p|n,k)");
    posterior->SetLineColor(2);
    for( int b=1; b<nBinsP; b++ )
	{
	    double x = posterior->GetBinCenter(b);
	    if( x > 0. && x < 1. )
		posterior->SetBinContent( b, SmartBinomial( k, n, x ) );
	}

    std::cout << "Posterior mean: " << posterior->GetMean() << std::endl;
    std::cout << "Posterior variance:  " << pow(posterior->GetStdDev(),2.) << std::endl;

    
    // Draw
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->Divide(1,2);
    can->cd(1);
    binomialHisto->Draw();
    can->cd(2);
    posterior->Draw();

    app->Run(kTRUE);
    
    return 0;
}
