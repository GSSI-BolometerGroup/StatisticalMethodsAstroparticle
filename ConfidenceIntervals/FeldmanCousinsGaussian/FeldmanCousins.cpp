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

// In this example, we construct the 90% Feldman-Cousins belt for a variable x
// which is Gaussian distributed with mean mu and sigma=1.
// The observable x can take any value, even negative ones,
// but we assume that the underlying parameter mu can only be non-negative.
// The Feldman-Cousins method, introduced in Phys. Rev. D 57 (1998) 3873-3889,
// solves the flip-flopping problem by using an ordering principle
// based on the likelihood-ratio instead than on the likelihood itself.

double Gaussian( double x,
		 double mu,
		 double sigma )
{
    return exp( -0.5 * std::pow( (x-mu)/sigma, 2. ) ) / sqrt(2.*M_PI) / sigma;
}

double Fmax( double* x,
	     double* p )
{
    double xx = x[0];
    if( xx < 0. )
	return exp( - 0.5 * std::pow( xx, 2. ) ) / sqrt( 2. * M_PI );
    else
	return 1. / sqrt( 2. * M_PI );
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

    // Extract the maximum of f(mu) for fixed values of x from the belt itself.
    // This works well until we reach the upper cutoff of the Neyman belt for x=7.
    TH1D* fmax = new TH1D( "max{f(x|#mu)}", "max{f(x|#mu)}", nBinsX, minX, maxX );
    fmax->GetXaxis()->SetTitle("x");
    fmax->GetYaxis()->SetTitle("max{f(x|#mu)}");
    fmax->SetFillColor(kBlue-5);
    for( int bx=1; bx<=nBinsX; bx++ )
	{
	    double max = 0.;
	    for( int bmu=1; bmu<=nBinsMu; bmu++ )
		if( belt->GetBinContent(bx,bmu) > max )
		    max = belt->GetBinContent(bx,bmu);

	    fmax->SetBinContent( bx, max );
	}
    // Theoretical calculation of the max(f(x)).
    // We will use this to compute the likelihood ratio distribution
    TF1* theoretical_fmax = new TF1( "theoretical_fmax", Fmax, minX, maxX, 0 );

    // Compute the distribution of the likelihood ratio
    TH2D* Lratio = new TH2D( "Likelihood ratio", "Likelihood ratio", nBinsX, minX, maxX, nBinsMu, minMu, maxMu );
    Lratio->GetXaxis()->SetTitle("x [bananas]");
    Lratio->GetYaxis()->SetTitle("#mu [bananas]");
    
    for( int bx=1; bx<=nBinsX; bx++ )
	for( int bmu=1; bmu<=nBinsMu; bmu++ )
	    {
		double x = minX + dX * ( 0.5 + bx - 1 );
		Lratio->SetBinContent( bx, bmu, belt->GetBinContent(bx,bmu) / theoretical_fmax->Eval(x) );
	    }

    // Construct the Feldman-Cousins confidence belt
    TH2D* FCbelt = new TH2D( "Feldman-Cousins belt", "Feldman-Cousins belt", nBinsX, minX, maxX, nBinsMu, minMu, maxMu );
    FCbelt->GetXaxis()->SetTitle("x [bananas]");
    FCbelt->GetYaxis()->SetTitle("#mu [bananas]");

    double quantile = 0.9;
    // Loop over mu
    for( int bmu=1; bmu<=nBinsMu; bmu++ )
	{
	    // Compute the integral of f(x) for a fixed value of mu.
	    // This is not one because of the non-unitary bin-width of the histogram used for the Neyman belt.
	    double integral = 0.;
	    for( int bx=1; bx<=nBinsX; bx++ )
		integral += belt->GetBinContent(bx,bmu);

	    // Find the mode of the likelihood-ratio for a fixed value of mu.
	    // At the same time, start computing the partial integral.
	    double modeBin = -1;
	    double partial = 0.;
	    double maxLratio = 0;
	    for( int bx=1; bx<=nBinsX; bx++ )
		if( Lratio->GetBinContent(bx,bmu) > maxLratio )
		    {
			maxLratio = Lratio->GetBinContent(bx,bmu);
			partial   = belt->GetBinContent(bx,bmu);
			modeBin   = bx;
		    }
	    double binLeft  = modeBin;
	    double binRight = modeBin;

	    // Compute the partial integral using the likelihood-ratio as ordering principle.
	    while( partial < quantile * integral )
		{
		    double lambdaLeft  = Lratio->GetBinContent( binLeft-1, bmu );
		    double lambdaRight = Lratio->GetBinContent( binRight+1, bmu );
		    double intLeft  = belt->GetBinContent( binLeft-1, bmu );
		    double intRight = belt->GetBinContent( binRight+1, bmu );

		    bool moveLeft  = false;
		    bool moveRight = false;
		    if( lambdaLeft  >= lambdaRight && binLeft > 1 ) moveLeft = true;
		    if( lambdaRight >= lambdaLeft  && binRight < nBinsX ) moveRight = true;

		    if( moveLeft )
			{
			    partial += intLeft;
			    binLeft --;
			}

		    if( moveRight )
			{
			    partial += intRight;
			    binRight ++;
			}
		    
		    if( !moveLeft && !moveRight )
			{
			    std::cout << "Cazzo: " << bmu << "\t" << partial/integral << std::endl;
			    break;
			}
		}

	    // Fill the Feldman-Cousins confidence belt
	    for( int bx=binLeft; bx<=binRight; bx++ )
		FCbelt->SetBinContent( bx, bmu, 1 );
	}

    // Compute the coverage of mu
    TH1D* coverage = new TH1D("Coverage(#mu)","Coverage(#mu)",nBinsMu,minMu,maxMu);
    coverage->GetXaxis()->SetTitle("#mu [bananas]");
    coverage->GetYaxis()->SetTitle("Coverage(#mu)");
    coverage->GetYaxis()->SetRangeUser(0,1);
    for( int bmu=1; bmu<=nBinsMu; bmu++ )
	for( int bx=1; bx<=nBinsX; bx++ )
	    if( FCbelt->GetBinContent( bx, bmu ) == 1 )
		coverage->AddBinContent( bmu, belt->GetBinContent(bx,bmu)*dX );

    // Compute the coverage for mu, this time using toy-MC experiments
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<> gaus(0,sigma);

    // Copy the lower and upper limit of the FC belt to a vector
    // to speed up the comparison later on.
    std::vector<double> lowMu(nBinsX);
    std::vector<double> highMu(nBinsX);
    for( int bx=0; bx<nBinsX; bx++ )
	{
	    bool foundLow = false;
	    for( int bmu=0; bmu<nBinsMu; bmu++ )
		if( FCbelt->GetBinContent( bx+1, bmu+1 ) == 1 )
		    {
			if( !foundLow )
			    {
				lowMu[bx] = minMu + dMu * ( 0.5 + bmu );
				foundLow = true;
			    }
			else
			    highMu[bx] = minMu + dMu * ( 0.5 + bmu );
		    }
	}

    TH1D* toycoverage = new TH1D("Toy-MC coverage(#mu)","Toy-MC coverage(#mu)",nBinsMu,minMu,maxMu);
    toycoverage->GetXaxis()->SetTitle("#mu [bananas]");
    toycoverage->GetYaxis()->SetTitle("Coverage(#mu)");
    toycoverage->GetYaxis()->SetRangeUser(0,1);
    
    int M = 10000;
    for( int bmu=0; bmu<nBinsMu; bmu++ )
	{
	    double mu = minMu + dMu * ( 0.5 + bmu );
	    double nAccepted = 0.;
	    for( int m=0; m<M; m++ )
		{
		    double x = mu + gaus(gen);
		    int bx = ( x - minX ) / dX;
		
		    if( mu >= lowMu[bx] &&
			mu <= highMu[bx] )
			nAccepted += 1.;
		}
	    toycoverage->SetBinContent( bmu, nAccepted/M );

	}
    

    // Draw
    gStyle->SetOptStat(0);
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->Divide(3,2);
    can->cd(1);
    belt->Draw("colz");
    can->cd(2);
    fmax->Draw();
    theoretical_fmax->Draw("same");
    can->cd(3);
    Lratio->Draw("colz");
    can->cd(4);
    FCbelt->Draw("colz");
    can->cd(5);
    coverage->Draw();
    can->cd(6);
    toycoverage->Draw();

    can->SaveAs("FeldmanCousinsGaussian.jpg");
    can->SaveAs("FeldmanCousinsGaussian.pdf");
	
    app->Run();
        
    return 0;
}
