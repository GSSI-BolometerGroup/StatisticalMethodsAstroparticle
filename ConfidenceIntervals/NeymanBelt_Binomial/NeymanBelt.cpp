#include <iostream>
#include <vector>
#include <random>
#include <string>

#include "TH2D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TF1.h"

// In this example, we construct the Neyman belt for a Binomial distribution with N generated events.
// This applies to the determination of the trigger efficiency using pulser-injected events:
// suppose we inject N events into our detector, measure k events, and want to evaluate the trigger efficiency p.
// We proceed as follows:
// 1- create a 2-dim histogram (TH2D)
// 2- for each value of p, populate k with the corresponding Binomial distribution
// 3- for each value of p, extract the 68% shortest interval for k
// 4- invert the confidence belt
// 5- compute the coverage for each value of p

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
    int N = 10;
    double minK   = -0.5;// This is just so that the first bin is centered at 0
    double maxK   = +0.5+N;
    int    nBinsK = N+1;// k goes from 0 to N
    double dP     = 0.001;// Let's use a fine binning for p, which goes from 0 to 1
    double minP   = 0.;
    double maxP   = 1.;
    int    nBinsP = (maxP-minP)/dP+0.5;

    // Prepare the 2-dim histogram
    TH2D* belt = new TH2D( "Neyman belt", "Neyman belt", nBinsK, minK, maxK, nBinsP, minP, maxP );
    belt->GetXaxis()->SetTitle("k [counts]");
    belt->GetYaxis()->SetTitle("p");

    // Populate the Neyman belt with the theoretical Binomial distribution.
    // Alternatively we could use random numbers, but this is faster.
    for( int bp=1; bp<=nBinsP; bp++ )
	{
	    double p = 0.5 * dP + dP * (bp-1);
	    for( int bk=1; bk<=nBinsK; bk++ )
		{
		    double k = bk-1;
		    belt->SetBinContent( bk, bp, SmartBinomial( k, N, p ) );
		}
	}

    // Prepare a 2-dim histogram to draw the 68% Neyman belt
    TH2D* belt68 = new TH2D( "Neyman belt 68%", "Neyman belt 68%", nBinsK, minK, maxK, nBinsP, minP, maxP );
    belt68->GetXaxis()->SetTitle("k [counts]");
    belt68->GetYaxis()->SetTitle("p");

    // Populate the 68% Neyman belt
    double quantile = 0.6826895;
    double integral;
    // Loop over p (bp is the bin along p, bk is the bin along k)
    for( int bp=1; bp<=nBinsP; bp++ )
	{

	    // Let's start by finding the mode for each fixed value of p,
	    integral = 0.;
	    int maxBK = -1;
	    for( int bk=1; bk<=nBinsK; bk++ )
		if( belt->GetBinContent(bk,bp) > integral )
		    {
			integral = belt->GetBinContent(bk,bp);
			maxBK = bk;
		    }

	    // Compute the shortest interval
	    int binLeft = maxBK;
	    int binRight = maxBK;
	    while( integral < quantile)
		{
		    double intLeft  = belt->GetBinContent( binLeft-1, bp );
		    double intRight = belt->GetBinContent( binRight+1, bp);
		    if( intLeft > intRight )
			{
			    integral += intLeft;
			    binLeft --;
			}
		    else if( intLeft < intRight )
			{
			    integral += intRight;
			    binRight ++;
			}
		    else// It will probably never happen that intLeft == intRight
			{
			    integral += intLeft + intRight;
			    binLeft --;
			    binRight ++;
			}
		}

	    // Populate the 68% confidence belt
	    for( int bk=binLeft; bk<=binRight; bk++ )
		belt68->SetBinContent( bk, bp, 1 );
	}

    // Compute the coverage
    TH1D* coverage = new TH1D( "Coverage(p)", "Coverage(p)", nBinsP, minP, maxP );
    coverage->GetXaxis()->SetTitle("p");
    coverage->GetYaxis()->SetTitle("Coverage(p)");
    coverage->GetYaxis()->SetRangeUser(0,1);
    for( int bp=1; bp<=nBinsP; bp++ )
	for( int bk=1; bk<=nBinsK; bk++ )
	    if( belt68->GetBinContent( bk, bp ) == 1 )
		coverage->AddBinContent( bp, belt->GetBinContent(bk,bp) );


    // Draw
    gStyle->SetOptStat(0);
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 2400, 900 );
    can->Divide(3,1);
    can->cd(1);
    belt->Draw("colz");
    can->cd(2);
    belt68->Draw("colz");
    can->cd(3);
    coverage->Draw();
    app->Run();
    
    return 0;
}
