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

double Gaussian( double x,
		 double mean,
		 double sigma )
{
    return exp( -0.5 * pow( (x-mean)/sigma, 2.) ) / sqrt(2.*M_PI) / sigma;
}

int main()
{

    // Prepare random number machinery
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

    // Number of points to generate in both cases
    int M = 100001;// Add one event because we'll skip the last one

    // --------------------------------
    // Gaussian random number generator
    double mean = 3;
    double sigma = 2;
    std::normal_distribution<> gaus(mean,sigma);

    double x;
    double prevx = 0;

    // Prepare histograms
    double nBins = 1000;
    double minx = mean - 7. * sigma;
    double maxx = mean + 7. * sigma;
    double mindiff = mean - 14. * sigma;
    double maxdiff = mean + 14. * sigma;
    TH1D* random_x = new TH1D( "random_x", "random_x", nBins, minx, maxx );
    TH1D* random_diff = new TH1D( "random_diff", "random_diff", nBins, mindiff, maxdiff );
    // Populate histograms with Gaussian distro, and with the difference
    // of number x[i] and x[i-1]
    for( int i=0; i<M; i++ )
	{
	    x = gaus(gen);
	    if( i!= 0 )
		{
		    random_x ->Fill( x );
		    random_diff->Fill( x - prevx );
		}
	    
	    prevx = x;
	}

    // -------------------
    // Metropolis-Hastings

    
    // Define the proposal function
    // Caveat: this proposal function is just guessed, and was not optimized in any way.
    // Real-life algorithms have some non-trivial internal method to optimize the parameters of the proposal function.
    double proposal_mean = 0;
    double proposal_sigma = 0.3;
    std::normal_distribution<> proposal(proposal_mean,proposal_sigma);

    // Prepare histograms
    TH1D* mcmc_x_accepted    = new TH1D( "mcmc_x_accepted",    "mcmc_x_accepted",    nBins, minx, maxx );
    mindiff = proposal_mean - 7. * proposal_sigma;
    maxdiff = proposal_mean + 7. * proposal_sigma;
    TH1D* mcmc_diff_accepted = new TH1D( "mcmc_diff_accepted", "mcmc_diff_accepted", nBins, mindiff, maxdiff );

    // Define test variable r
    std::uniform_real_distribution<> test(0,1);
    // Start with a random point in the parameter space
    std::uniform_real_distribution<> mcmc_first(-1,1);//(minx,maxx);
    prevx = mcmc_first(gen);
    double prevy = Gaussian( prevx, mean, sigma );
    double y, r;
    // Generate following points using proposal function
    bool accepted;
    int nAccepted = 0;
    int nBurnIn = 1000;
    //for( int i=1; i<M; i++ )
    // Generate points until we accept M of them.
    // Skip the first ones (nBurnIn) to guarantee some sort of indepencence from the starting point
    while( nAccepted < M + nBurnIn )
	{
	    accepted = false;

	    // Generate new point
	    x = prevx + proposal(gen);
	    y = Gaussian( x, mean, sigma );

	    // Check if we should accept the point or not
	    r = test(gen);
	    if( y/prevy > r )
		accepted = true;
	    
	    if( accepted )
		{
		    // Fill histograms of accepted points (provided we passed the burn-in phase)
		    if( nAccepted > nBurnIn )
			{
			    mcmc_x_accepted->Fill(x);
			    mcmc_diff_accepted->Fill(x-prevx);
			}
		    // Update previous point to the current one
		    nAccepted ++;
		    prevx = x;
		    prevy = y;

		}
		    
	}


    random_x->GetXaxis()->SetTitle("x_{random} [bananas]");
    random_x->GetYaxis()->SetTitle("Entries");
    random_diff->GetXaxis()->SetTitle("#delta x_{random} [bananas]");
    random_diff->GetYaxis()->SetTitle("Entries");

    mcmc_x_accepted->GetXaxis()->SetTitle("x_{mcmc} [bananas]");
    mcmc_x_accepted->GetYaxis()->SetTitle("Entries");
    mcmc_diff_accepted->GetXaxis()->SetTitle("#delta x_{mcmc} [bananas]");
    mcmc_diff_accepted->GetYaxis()->SetTitle("Entries");
    
    // Draw
    gStyle->SetOptTitle(0);
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->Divide(2,2);
    can->cd(1);
    random_x->Draw();
    can->cd(3);
    random_diff->Draw();
    can->cd(2);
    mcmc_x_accepted->Draw();
    can->cd(4);
    mcmc_diff_accepted->Draw();

    can->SaveAs("RandomNumberCorrelation.jpg");
    
    app->Run(kTRUE);
    return 0;
}
