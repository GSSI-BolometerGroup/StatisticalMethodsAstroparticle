#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <iomanip>

#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TPaveStats.h"

// In this example, we compute the binomial distributions with n=150 and p={0.1,0.3,0.5,0.7,0.9},
// and plot them all together.
// We do not use any predefined method for the binomial distribution, or for the factorial,
// with the goal of showing the numerical issues connected with their calculation.
// To cross check the validity of your implementation, repeat with n=500.

// Works only up to n=170, after which it hits the max value allowed for a double and returns nan.
double Factorial( double n ){
    if( n <= 1. )
	return 1.;
    else
	return n * Factorial(n-1.);
}

// This is the most straightforward implementation of the binomial coefficient.
// It works only up to n=170, after which Factorial(n) returns nan.
double BinomialCoefficient( double n, double k )
{
    return Factorial(n) / Factorial(k) / Factorial(n-k);
}

// Smarter implementation of the Binomial coefficient,
// which minimizes the number of calculations, and performs the sum over logs
// to avoid incurring in the limited precision of doubles.
// The use of log and exp induces a relative precision of ~1e-13
// with respect to the straightforward implementation.
double SmartBinomialCoefficient( double n, double k )
{

    double sum_1 = 0;
    double sum_2 = 0;
    for( double i=(int)(n-k+1); i<=n; i+=1. )
	sum_1 += log(i);
    for( double i=1.; i<=k; i+=1. )
	sum_2 += log(i);

    return exp( sum_1 - sum_2 );
}

// Straightforward implementation of the binomial.
// Works up to n=170.
double Binomial( double k,
		 double n,
		 double p )
{
    return BinomialCoefficient( n, k ) * pow( p, k ) * pow( 1.-p, n-k );
}

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

    // Define n
    double n=500;
    // Prepare and fill histograms (one for each value of p)
    std::vector<TH1D*> histo;
    int nBins = n+1;
    double min = -0.5;
    double max = 0.5 + n;
    // Loop over possible values of p
    for( double p=0.1; p<=0.9; p+=0.2 )
	{
	    std::string name( "p=" + std::to_string(p) );
	    histo.push_back( new TH1D( name.c_str(), name.c_str(), nBins, min, max ) );
	    histo.back()->GetXaxis()->SetTitle("k");
	    histo.back()->GetYaxis()->SetTitle("P(k)");
	    for( int b=1; b<=nBins; b++ )
		{
		    double k = histo.back()->GetBinCenter(b);
		    // Set the content of each bin with the corresponding value of the Binomial distribution
		    histo.back()->SetBinContent( b, SmartBinomial( k, n, p ) );
		}
	}

    // Some graphical options
    for( size_t i=0; i<histo.size(); i++ )
	{
	    histo[i]->SetLineColor(kGreen-i);
	    histo[i]->SetLineWidth(2);
	    histo[i]->SetLineStyle(1+i);
	}

    // Draw all histograms together
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->cd();
    histo[0]->Draw();
    for( long unsigned int i=1; i<histo.size(); i++ )
	histo[i]->Draw("same");
    can->SaveAs("Binomial.jpg");
    can->SaveAs("Binomial.pdf");
    app->Run(kTRUE);
    
    return 0;
}
