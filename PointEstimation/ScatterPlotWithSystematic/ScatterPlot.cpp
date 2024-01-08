#include <iostream>
#include <vector>
#include <random>
#include <string>

#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLegend.h"

double ChiSquare( double x2, double n )
{
    double x = sqrt(x2);
    double p = pow( 2., -0.5 * n );
    p *= pow( x, n-2. );
    p *= exp( -0.5 * x2 );
    p /= TMath::Gamma(0.5*n);
    return p;
}

int main()
{
    // We will generate N points distributed along a quadratic curve f(x)=a+bx+cx^2,
    // shifting them randomly according to a Gaussian distribution centered at 0 and with a sigma
    // also randomly generated with a flat distribution between 10 and 100.
    // We will fit the scatter plot with a straight line g(x)=a+bx, systematically ignoring the 2nd order term.
    // In other words:
    // x_i     = equally spaced points
    // sigma_i = flat random number in [10,100]
    // dy_i    = Gaussian random number with mean 0 and sigma=sigma_i
    // y_i     = f(x_i) + dy_i
    // To each point y_i, we assign the uncertainty sigma_i.
    // At this point, we fit the (x,y) scatterplot with g(x) and extract the minimum Chi2 from the minimizer.
    // We repeat this M times, populate a histogram with the min-Chi2 values,
    // and compare it with the theoretical Chi2 distribution for N-2 degrees of freedom (because g(x) has 2 parameters).
    // We repeat this a few times for different values of c, so for different amplitudes of the quadratic term in f(x).
    
    // Set the min and max of x, and the step
    double min = 0.;
    double max = 1000.;
    double step = 100.;
    int N = (max-min)/step + 1;
    // Set the parameters of f(x)
    double offset = 3.;// Just a random value
    double slope = 1.5;// Just a random value
    // Set the number of toy experiments
    int M = 10000;

    // We'll store X and Y on vectors. Xerr will be filled with zeroes.
    std::vector<double>* X    = new std::vector<double>(N);
    std::vector<double>* Xerr = new std::vector<double>(N);
    std::vector<double>* Y    = new std::vector<double>(N);
    std::vector<double>* Yerr = new std::vector<double>(N);
    // Copy this to a TGraphErrors object. At the moment all values are zero, but we'll update later on for each toy experiment.
    TGraphErrors* scatterplot = new TGraphErrors( N, &(X->at(0)), &(Y->at(0)), &(Xerr->at(0)), &(Yerr->at(0)) );
    scatterplot->GetXaxis()->SetTitle("Amplitude [mV]");
    scatterplot->GetYaxis()->SetTitle("Energy [keV]");

    // Prepare a random number generator with flat distribution to generate sigma_i
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    double minsigma = 10;
    double maxsigma = 100;
    std::uniform_real_distribution<> dis(minsigma,maxsigma);

    std::vector<double> C;
    for( double c=0.; c<=5.e-4; c+=1.e-4 )
	C.push_back(c);
    size_t nC = C.size();
    
    // Define the fitting function
    TF1* func = new TF1( "func", "[0]+[1]*x", min, max );
    // Prepare the histo for the min-Chi2 distribution
    double maxChi2=N+10*sqrt(N);
    int nBinsChi2 = 100;
    std::vector<TH1D*> h_Chi2(nC);
    for( size_t c=0; c<nC; c++ )
	{
	    std::string name( "c=" + std::to_string(C[c]) );
	    h_Chi2[c] = new TH1D(name.c_str(),name.c_str(),nBinsChi2,0,maxChi2);
	    h_Chi2[c]->GetXaxis()->SetTitle("#Chi^{2}");
	}
    h_Chi2[0]->SetLineColor(1);
    h_Chi2[1]->SetLineColor(kBlue+1);
    h_Chi2[2]->SetLineColor(kGreen+2);
    h_Chi2[3]->SetLineColor(kCyan+1);
    h_Chi2[4]->SetLineColor(kOrange);
    h_Chi2[5]->SetLineColor(kMagenta+2);
    
    // This is just ROOT stuff needed to retrieve the min-Chi2 of each fit
    TFitResultPtr fitresultPtr;
    TFitResult* fitresult;

    // Loop over C
    for( size_t c=0; c<nC; c++ )
	// Loop over the toy experiments
	for( int m=0; m<M; m++ )
	    {
		// Loop over the points
		for( int i=0; i<N; i++ )
		    {
			// Generate random sigma
			double sigma = dis(gen);
			// Create Gaussian random generator
			std::normal_distribution<> dissigma(0,sigma);
			// Generate the point
			X->at(i)    = min + step * i;
			Xerr->at(i) = 0.;
			Y->at(i)    = offset
			    + X->at(i) * slope
			    + dissigma(gen)
			    + C[c] * std::pow(X->at(i),2.);
			Yerr->at(i) = sigma;
			// Update the TGraphErrors
			scatterplot->SetPoint( i, X->at(i), Y->at(i) );
			scatterplot->SetPointError( i, Xerr->at(i), Yerr->at(i) );
		    }
		// Fit and get the min-Chi2 for the current fit.
		// Fit options:
		// N --> does not draw
		// Q --> does not print fit results on the terminal
		// R --> Use range specified in function declaration
		// S --> Save fit result to TFitResult object (necessary to retrieve min-Chi2)
		fitresultPtr = scatterplot->Fit( func, "NQRS" );
		fitresult = fitresultPtr.Get();
		// Fill the histogram of min-Chi2 values
		h_Chi2[c]->Fill( fitresult->Chi2() );
	    }

    // Compute the theoretical Chi2 for (N-2) degrees of freedom
    int nBinsTheorChi2 = 1000;
    TH1D* theoreticalChi2 = new TH1D("theoreticalChi2", "theoreticalChi2", nBinsTheorChi2,0,maxChi2 );
    for( int b=1; b<=1000; b++ )
	theoreticalChi2->SetBinContent( b, M * nBinsTheorChi2 / nBinsChi2 * theoreticalChi2->GetBinWidth(1)*ChiSquare( theoreticalChi2->GetBinCenter(b), N-2 ) );
    theoreticalChi2->SetLineColor(2);
    theoreticalChi2->SetLineWidth(2);

    // Draw
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->Divide(1,2);
    can->cd(1);
    // Draw the last scatter plot and the corresponding function
    scatterplot->Draw("AP");
    func->Draw("same");
    can->cd(2);
    // Draw the min-Chi2 histo and the corrisponding theoretical Chi2 distribution
    
    h_Chi2[0]->Draw();
    for( size_t c=1; c<nC; c++ )
	h_Chi2[c]->Draw("SAME");
    theoreticalChi2->Draw("SAME");

    TLegend* leg = new TLegend(0.75,0.5,0.91,0.91);
    for( auto& h: h_Chi2 )
	leg->AddEntry( h, h->GetName(), "L" );
    leg->Draw();
	
    can->SaveAs("ScatterPlotWithSystematic.jpg");    
    app->Run();
    
    return 0;
}
