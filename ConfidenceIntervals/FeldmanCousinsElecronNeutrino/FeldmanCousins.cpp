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

// In this example, we reproduce the KATRIN limit on the electron neutrino mass
// reported in Phys. Rev. Lett. 123 (2019) 221802.

double Gaussian( double x,
		 double mu,
		 double sigma )
{
    return exp( -0.5 * std::pow( (x-mu*mu)/sigma, 2. ) ) / sqrt(2.*M_PI) / sigma;
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
    double dMu     = 0.01;
    double minMu   = 0.-0.5*dMu;
    double maxMu   = 4.+0.5*dMu;
    int    nBinsMu = (maxMu-minMu)/dMu+0.5;
    double dX      = 0.01;
    double minX    = -5.-0.5*dX;
    double maxX    = 20.+0.5*dX;
    int    nBinsX  = (maxX-minX)/dX+0.5;

    double sigma   = 1.;
    
    // Prepare the 2-dim histogram
    TH2D* belt = new TH2D( "Neyman belt", "Neyman belt", nBinsX, minX, maxX, nBinsMu, minMu, maxMu );
    belt->GetXaxis()->SetTitle("m^{2} [eV^{2}]");
    belt->GetYaxis()->SetTitle("m [eV]");


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
    TH1D* fmax = new TH1D( "max{f(m^{2}|m)}", "max{f(m^{2}|m)}", nBinsX, minX, maxX );
    fmax->GetXaxis()->SetTitle("m^{2}");
    fmax->GetYaxis()->SetTitle("max{f(m^{2}|m)}");
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
    Lratio->GetXaxis()->SetTitle("m^{2} [eV^{2}]");
    Lratio->GetYaxis()->SetTitle("m [eV]");
    
    for( int bx=1; bx<=nBinsX; bx++ )
	for( int bmu=1; bmu<=nBinsMu; bmu++ )
	    {
		double x = minX + dX * ( 0.5 + bx - 1 );
		Lratio->SetBinContent( bx, bmu, belt->GetBinContent(bx,bmu) / theoretical_fmax->Eval(x) );
	    }


    // Construct the Feldman-Cousins confidence belt
    TH2D* FCbelt = new TH2D( "Feldman-Cousins belt", "Feldman-Cousins belt", nBinsX, minX, maxX, nBinsMu, minMu, maxMu );
    FCbelt->GetXaxis()->SetTitle("m^{2} [eV^{2}]");
    FCbelt->GetYaxis()->SetTitle("m [eV]");

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

    // Get the actual limit
    double bestFitM2     = -1.0;// eV
    int    bestFitBinX  = ( bestFitM2 - minX ) / dX + 0.5;
    int    bestFitBinMu = 1;
    while( FCbelt->GetBinContent( bestFitBinX, bestFitBinMu ) == 1 )
	bestFitBinMu ++;
    double limit = FCbelt->GetYaxis()->GetBinLowEdge( bestFitBinMu );
    std::cout << "90% C.L. limit on m: " << limit << " eV" << std::endl;
    
    // Compute the coverage of x. Yes, the coverage of x, not of mu!
    TH1D* coverage = new TH1D("Coverage(m)","Coverage(m)",nBinsMu,minMu,maxMu);
    coverage->GetXaxis()->SetTitle("m [eV]");
    coverage->GetYaxis()->SetTitle("Coverage(m)");
    coverage->GetYaxis()->SetRangeUser(0,1);
    for( int bmu=1; bmu<=nBinsMu; bmu++ )
	for( int bx=1; bx<=nBinsX; bx++ )
	    if( FCbelt->GetBinContent( bx, bmu ) == 1 )
		coverage->AddBinContent( bmu, belt->GetBinContent(bx,bmu)*dX );

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
    can->SaveAs("FeldmanCousinsKATRIN.jpg");
     
    app->Run();
        
    return 0;
}
