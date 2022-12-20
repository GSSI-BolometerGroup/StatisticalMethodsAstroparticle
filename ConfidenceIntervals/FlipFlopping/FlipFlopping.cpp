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

// In this example, we construct the 90% Neyman belt for a variable x
// which is Gaussian distributed with mean mu and sigma=1.
// The observable x can take any value, even negative ones,
// but we assume that the underlying parameter mu can only be non-negative.
// If the measured value of x is negative or close to zero, say if x<3*sigma,
// we decide to quote an upper limit on mu, instead than a central value with symmetric uncertainties,
// because this would go into the unphysical region.
// This produces the so-called flip-flopping problem.
// Infact, if x<3*sigma, instead of building the 90% central Neyman belt
// and inverting it to obtain an interval for mu,
// we directly place a limit on mu by computing (vertically!) the 90% quantile of PDF(mu).
// As a result, the belt has a discontinuity for x=3.
// This leads to a situation of undercoverage for 1.3<x<4.7.

double Gaussian( double x,
		 double mu,
		 double sigma )
{
    return exp( -0.5 * std::pow( (x-mu)/sigma, 2. ) ) / sqrt(2.*M_PI) / sigma;
}

int main()
{
    double dX      = 0.01;
    double minX    = -3.-0.5*dX;
    double maxX    = 11.+0.5*dX;
    int    nBinsX  = (maxX-minX)/dX+0.5;
    double dMu     = 0.01;
    double minMu   = 0.-0.5*dMu;
    double maxMu   = 7.+0.5*dMu;
    int    nBinsMu = (maxMu-minMu)/dMu+0.5;
    double sigma   = 1.;
    
    // Prepare the 2-dim histogram
    TH2D* belt = new TH2D( "Neyman belt", "Neyman belt", nBinsX, minX, maxX, nBinsMu, minMu, maxMu );
    belt->GetXaxis()->SetTitle("x [bananas]");
    belt->GetYaxis()->SetTitle("#mu [bananas]");


    // Populate the Neyman belt
    for( int bmu=1; bmu<=nBinsMu; bmu++ )
	{
	    double mu = minMu + dMu * ( 0.5 + bmu - 1 );

	    for( int bx=1; bx<=nBinsX; bx++ )
		{
		    double x = minX + dX * ( 0.5 + bx - 1 );
		    belt->SetBinContent( bx, bmu, Gaussian( x, mu, sigma ) );
		}
	}

    // Construct the 90% central confidence belt
    TH2D* beltcentral = new TH2D( "Neyman belt central", "Neyman belt central", nBinsX, minX, maxX, nBinsMu, minMu, maxMu );
    beltcentral->GetXaxis()->SetTitle("x [bananas]");
    beltcentral->GetYaxis()->SetTitle("#mu [bananas]");
    double quantile = 0.9;
    
    for( int bmu=1; bmu<=nBinsMu; bmu++ )
	{
	    double integral = 0;
	    for( int bx=1; bx<=nBinsX; bx++ )
		integral += belt->GetBinContent(bx,bmu);

	    double partial = 0;
	    int    modeBin = -1;
	    for( int bx=1; bx<=nBinsX; bx++ )
		if( belt->GetBinContent(bx,bmu) > partial )
		    {
			partial = belt->GetBinContent(bx,bmu);
			modeBin = bx;
		    }

	    int binLeft  = modeBin;
	    int binRight = modeBin;
	    while( partial < quantile * integral )
		{
		    double intLeft  = belt->GetBinContent( binLeft-1, bmu );
		    double intRight = belt->GetBinContent( binRight+1, bmu);
		    if( intLeft > intRight )
			{
			    partial += intLeft;
			    binLeft --;
			}
		    else if( intLeft < intRight )
			{
			    partial += intRight;
			    binRight ++;
			}
		    else// It will probably never happen that intLeft == intRight
			{
			    partial += intLeft + intRight;
			    binLeft --;
			    binRight ++;
			}
		}

	    for( int bx=binLeft; bx<=binRight; bx++ )
		beltcentral->SetBinContent( bx, bmu, 1 );
	}
    
    // Construct the 90% combined confidence belt
    TH2D* beltcombined = new TH2D( "Neyman belt combined", "Neyman belt combined", nBinsX, minX, maxX, nBinsMu, minMu, maxMu );
    beltcombined->GetXaxis()->SetTitle("x [bananas]");
    beltcombined->GetYaxis()->SetTitle("#mu [bananas]");    

    bool upper;
    double limit;
    for( int bx=1; bx<=nBinsX; bx++ )// Notice that we start looping over bx now!
	{
	    double x = minX + dX * ( 0.5 + bx - 1 );

	    if( x / ( 3. * sigma ) < 1 )
		{
		    upper = true;
		    limit = std::max( 0., x );
		    limit += 1.282 * sigma;//corresponds to the 90% upper limit on a Gaussian
		}
	    else
		upper = false;

	    
	    for( int bmu=1; bmu<=nBinsMu; bmu++ )
		{
		    double mu = minMu + dMu * ( 0.5 + bmu - 1 );
		    if( upper )
			{
			    if( mu <= limit )
				beltcombined->SetBinContent( bx, bmu, 1 );
			}
		    else
			beltcombined->SetBinContent( bx, bmu, beltcentral->GetBinContent(bx,bmu) );
		}
	}

    // Compute the coverage of x. Yes, the coverage of x, not of mu!
    TH1D* coverage = new TH1D("Coverage(#mu)","Coverage(#mu)",nBinsMu,minMu,maxMu);
    coverage->GetXaxis()->SetTitle("#mu [bananas]");
    coverage->GetYaxis()->SetTitle("Coverage(#mu)");
    coverage->GetYaxis()->SetRangeUser(0,1);
    for( int bmu=1; bmu<=nBinsMu; bmu++ )
	for( int bx=1; bx<=nBinsX; bx++ )
	    if( beltcombined->GetBinContent( bx, bmu ) == 1 )
		coverage->AddBinContent( bmu, belt->GetBinContent(bx,bmu)*dX );
    
    // Draw
    gStyle->SetOptStat(0);
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->Divide(2,2);
    can->cd(1);
    belt->SetStats(0);
    belt->Draw("colz");
    can->cd(2);
    beltcentral->Draw("colz");
    can->cd(3);
    beltcombined->Draw("colz");
    can->cd(4);
    coverage->Draw();
    app->Run();
        
    return 0;
}
