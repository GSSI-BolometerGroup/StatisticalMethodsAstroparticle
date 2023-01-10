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

// In this example, we compute the distribution of the squared y of a variable x, y=x^2,
// where x has two possible distributions:
// a) Gaussian with mu=2 and sigma=1
// b) Poisson with lambda=3.
// We proceed in two ways:
// 1) Generate 10^6 random values of x, and histogram the corresponding value of y
// 2) Compute the theoretical distribution of y.
// In the end, we plot the distributions obtained with methods 1) and 2) and see if they overlap.

// These are just functions necessary to compute the Gaussian and Poisson distributions
double Factorial( double n ){
    if( n <= 1. )
	return 1.;
    else
	return n * Factorial(n-1.);
}

double Gaussian( double *x, double* par )
{
    // par[0] = integral
    // par[1] = mean
    // par[2] = sigma
    return par[0] / sqrt( 2. * M_PI ) / par[2] * exp( - pow( x[0] - par[1], 2. ) / 2. / pow( par[2], 2. ) );
}

double GaussianSquared( double *y, double* par )
{
    // par[0] = integral
    // par[1] = mean
    // par[2] = sigma
    double plus_x = sqrt(y[0]);
    double minus_x = -sqrt(y[0]);
    return 0.5 * ( Gaussian( &plus_x, par ) + Gaussian( &minus_x, par ) ) / plus_x;
}

double Poisson( double *x, double* par )
{
    // par[1] = lambda
    double integral = par[0];
    double lambda = par[1];
    return  integral * TMath::Poisson( x[0], lambda );//integral * exp( -lambda ) * pow( lambda, x[0] ) / Factorial( x[0] );
}

double PoissonSquared( double *y, double* par )
{
    // par[0] = integral
    // par[1] = lambda
    double plus_x = sqrt(y[0]);
    double minus_x = -sqrt(y[0]);
    return 0.5 * ( Poisson( &plus_x, par ) + Poisson( &minus_x, par ) ) / plus_x;
}

int main()
{
    // Number of events to be generated
    int N = 1000000;
	
    // Prepare histograms for Gaussian case
    double mu    = 2.;
    double sigma = 1.;
    double min = mu - 3.*sigma;
    double max = mu + 3.*sigma;
    double minSquared = std::min( pow(min,2.), pow(max,2.) );
    minSquared = std::min( 0., minSquared );
    double maxSquared = std::max( pow(min,2.), pow(max,2.) );
    int nBins = 10000;
    TH1D* gaussian_histo = new TH1D( "Gaussian", "Gaussian", nBins, min, max );
    TH1D* gaussianSquared_histo = new TH1D( "GaussianSquared", "GaussianSquared", nBins, minSquared, maxSquared );
    double binWidth = gaussian_histo->GetBinWidth(1);
    double binWidthSquared = gaussianSquared_histo->GetBinWidth(1);

    // Gaussian distribution
    TF1* gaussian = new TF1( "Gaussian", Gaussian, min, max, 3 );
    gaussian->SetParameter( 0, (double)N*binWidth );
    gaussian->SetParameter( 1, mu );
    gaussian->SetParameter( 2, sigma );

    // Theoretical distribution of y, with a Gaussian-distributed x
    TF1* gaussianSquared = new TF1( "GaussianSquared", GaussianSquared, 0, pow(max,2.), 3 );
    gaussianSquared->SetParameter( 0, (double)N*binWidthSquared );
    gaussianSquared->SetParameter( 1, mu );
    gaussianSquared->SetParameter( 2, sigma );

    // Prepare histogram for Poisson case
    double lambda = 3.;
    min = 0.;
    max = lambda + 5 * sqrt(lambda);
    TH1D* poisson_histo = new TH1D( "Poisson", "Poisson", nBins, min, max );

    // Poisson distribution
    TF1* poisson = new TF1( "Poisson", Poisson, min+0.5*binWidth, max-0.5*binWidth, 2 );
    poisson->SetParameter( 0, 1 );
    poisson->SetParameter( 1, lambda );

    minSquared = 0.;
    maxSquared = pow( max, 2. );
    TH1D* poissonSquared_histo = new TH1D( "PoissonSquared", "PoissonSquared", nBins, minSquared, maxSquared );
    binWidth = poisson_histo->GetBinWidth(1);
    binWidthSquared = poissonSquared_histo->GetBinWidth(1);

    // Theoretical distribution of y, with a Poisson-distributed y
    TF1* poissonSquared = new TF1( "PoissonSquared", PoissonSquared,
				   minSquared+0.5*binWidthSquared,
				   maxSquared-0.5*binWidthSquared, 2 );
    poissonSquared->SetParameter( 0, (double)N*binWidthSquared );
    poissonSquared->SetParameter( 1, lambda );
    
    // Populate the histograms with random values of x,
    // and the corresponding value of y.
    for( int i=0; i<N; i++ )
	{
	    // Generate random Gaussian x
	    double x = gaussian->GetRandom();
	    // Compute corresponding y
	    double y = pow(x,2.);
	    // Populate histograms
	    gaussian_histo->Fill(x);
	    gaussianSquared_histo->Fill(y);

	    // Repeat for the Poisson case
	    x = poisson->GetRandom();
	    y = pow(x,2.);
	    poisson_histo->Fill(x);
	    poissonSquared_histo->Fill(y);
	}

    // Draw the 4 distributions (x and y, Gaussian and Poisson case)
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->Divide(2,2);
    can->cd(1);
    gaussian_histo->Draw();
    gaussian->SetNpx(10000);
    gaussian->Draw("same");
    can->cd(3);
    gaussianSquared_histo->Draw();
    gaussianSquared->SetNpx(10000);
    gaussianSquared->Draw("same");
    can->cd(2);
    poisson_histo->Draw();
    poisson->SetParameter( 0, (double)N*binWidth );
    poisson->SetNpx(10000);
    poisson->Draw("same");
    can->cd(4);
    poissonSquared_histo->Draw();
    poissonSquared->SetNpx(10000);
    poissonSquared->Draw("same");
    
    
    app->Run(kTRUE);
    
    return 0;
}
