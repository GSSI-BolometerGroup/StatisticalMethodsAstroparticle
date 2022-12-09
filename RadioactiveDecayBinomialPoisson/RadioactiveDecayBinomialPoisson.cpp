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

// Smarter implementation of hte Binomial,
// which minimizes the number of calculations and includes the log-based
// calculation of the binomial coefficient reported above.
// It has a relative precision of ~1e-15 with respect to the straightforward implementation.
double SmartBinomial( double k, double n, double p, double A=1 )
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

    return A * exp( sum_1 - sum_2 + exp_1 + exp_2 );
}

// Computes the logarithm of the factorial.
// Useful if the factorial is divided by a similarly large number, as in the Poisson case.
double LogFactorial( double n )
{
    if( n <= 1. )
	return 0.;
    return log(n) + LogFactorial(n-1.);
}

// Smarter way to implement the Poisson distribution, which works for n>170.
double SmartPoisson( double n,
		     double lambda,
		     double A=1 )
{
    return A * exp( -lambda + n * log(lambda) - LogFactorial(n) );
}

int main()
{
    // Prepare random number machinery
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    double halflife = 138.4;// Halflife of 210Po (in days)
    std::exponential_distribution dis(log(2.)/halflife);


    double n = 100.;// Initial number of nuclei
    std::vector<double> dt{0.1,0.2,0.5,1,2,3,4,5,6};
    //for( double x=1.; x<=6.; x+=1. )
    //	dt.push_back( x * halflife );
    size_t nCases = dt.size();
    int M=10000;// Number of toy-MC experiments

    std::vector<TH1D*> dataHisto(nCases);
    std::vector<TH1D*> binomialHisto(nCases);
    std::vector<TH1D*> poissonHisto(nCases);

    TF1* expfunc = new TF1( "expfunc", "[0]*exp(-[0]*x)", 0, 10. * halflife );
    expfunc->SetParameter( 0, log(2.)/halflife);
    expfunc->SetNpx(100000);
    
    for( size_t c=0; c<nCases; c++ )
	{

	    // Compute expected number of decays in given dt period
	    double p = expfunc->Integral(0.,dt[c]*halflife);
	    double lambda = p * n;
	    // Prepare histograms
	    double min = std::max( -0.5+(int)(lambda-7.*sqrt(lambda)), -0.5 );
	    double max = 0.5+(int)(lambda+7.*sqrt(lambda));
	    int nBins = max-min;
	    // Histograms for toy-MC data
	    std::string name( "Data N=" + std::to_string((int)n) + " dt=" + std::to_string(dt[c]) + " T1/2" );
		dataHisto[c] = new TH1D( name.c_str(), name.c_str(), nBins, min, max );
		// Histogram for binomial distro
		name = "Binomial N=" + std::to_string((int)n) + " dt=" + std::to_string(dt[c]) + "T1/2";
		binomialHisto[c] = new TH1D( name.c_str(), name.c_str(), nBins, min, max );
		// Histogram for Poisson distro
		name = "Poisson N=" + std::to_string((int)n) + " dt=" + std::to_string(dt[c]) + "T1/2";
		poissonHisto[c] = new TH1D( name.c_str(), name.c_str(), nBins, min, max );
		// Graphic settings
		dataHisto[c]->GetXaxis()->SetTitle("k");
		dataHisto[c]->GetXaxis()->SetTitle("Entries");
		binomialHisto[c]->SetLineColor(2);
		poissonHisto[c]->SetLineColor(4);

		// Loop over toy-MC experiments
		for( int i=0; i<M; i++ )
		    {
			// Generate n time values
			std::vector<double> time(n);
			for( auto& t: time )
			    t = dis(gen);
			// count cases with t<dT
			int k=0;
			for( auto& t: time )
			    if( t < dt[c]*halflife ) k++;
			// Fill data histogram
			dataHisto[c]->Fill(k);
		    }

		// Fill binomial and Poisson distributions
		for( int b=1; b<nBins; b++ )// The bin number of ROOT TH1D starts from 1.
		    {
			double k = binomialHisto[c]->GetBinCenter(b);
			double tmp = SmartBinomial(k,n,p,M);
			if( !std::isnan(tmp) )// Just to avoid some rare failure of SmartBinomial()
			    binomialHisto[c]->SetBinContent( b, tmp );
			poissonHisto[c]->SetBinContent( b, SmartPoisson(k,lambda,M) );
		    }
	}

    // Draw
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->Divide(3,3);
    for( size_t c=0; c<nCases; c++ )
	{
	    can->cd(c+1);
	    dataHisto[c]->Draw();
	    binomialHisto[c]->Draw("SAME");
	    poissonHisto[c]->Draw("SAME");
	}
    app->Run(kTRUE);
    
    return 0;
}
