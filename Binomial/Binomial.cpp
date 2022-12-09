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
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    //std::cout << Binomial( 0, 2, 0.5 ) << std::endl;

    double n=150;
    std::vector<TH1D*> histo;
    int nBins = n+1;
    double min = -0.5;
    double max = 0.5 + n;
    for( double p=0.1; p<=0.9; p+=0.2 )
	{
	    std::string name( "n=" + std::to_string(p) );
	    histo.push_back( new TH1D( name.c_str(), name.c_str(), nBins, min, max ) );
	    for( int b=1; b<=nBins; b++ )
		{
		    double k = histo.back()->GetBinCenter(b);
		    histo.back()->SetBinContent( b, SmartBinomial( k, n, p ) );
		}
	}

    for( size_t i=0; i<histo.size(); i++ )
	{
	    histo[i]->SetLineColor(kGreen-i);
	    histo[i]->SetLineWidth(2);
	    histo[i]->SetLineStyle(1+i);
	}
    
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->cd();
    histo[0]->Draw();
    for( long unsigned int i=1; i<histo.size(); i++ )
	histo[i]->Draw("same");
    app->Run(kTRUE);
    
    return 0;
}
