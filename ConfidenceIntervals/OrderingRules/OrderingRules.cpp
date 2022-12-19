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

// In this example, we suppose to have a posterior-PDF for a parameter lambda
// describing the number of signal counts measured by some experiment (with zero background).
// Assuming n=3, the posterior for lambda will be asymmetric, with mode at lambda=3.
// We will compute the 68% interval in three different ways:
// 1- The central values with symmetric error bars containing a 68% quantile;
// 2- Integrate from the left and from the right, until we reach a 16% quantile on each side;
// 3- Compute the shortest interval containing the 68% quantile.
// As a cross check, we will print the coverage of all three methods.

double LogFactorial( double n )
{
    if( n <= 1. )
	return 0.;
    return log(n) + LogFactorial(n-1.);
}

double SmartPoisson( double n,
		     double lambda,
		     double A=1 )
{
    return A * exp( -lambda + n * log(lambda) - LogFactorial(n) );
}

int main()
{
    // Create posterior for n=3
    double n=3;
    double min = 0.;// Lambda defined >= 0
    double max = 15.;
    int nBins = 10000;// Reduce bin size to visualize continuous behaviour
    double binWidth = (max-min) / nBins;
    TH1D* posterior = new TH1D("Posterior","Posterior", nBins, min, max );
    posterior->GetXaxis()->SetTitle("#lambda [events]");
    posterior->GetYaxis()->SetTitle("P(#lambda|n=3)");
    // Fill the posterior
    for( int b=1; b<=nBins; b++ )
	{
	    double x = posterior->GetBinCenter(b);
	    posterior->SetBinContent( b, SmartPoisson( n, x ) );
	}
    // Normalize the posterior to 1, so that we can ignore the bin width.
    posterior->Scale( 1./posterior->Integral(1,nBins) );
    // Set all bin error to zero just for drawing reasons (error bars get messed up by the rescaling)
    for( int b=1; b<=nBins; b++ )
	posterior->SetBinError(b,0);
    
    // Find the mode
    double modeBin = posterior->GetMaximumBin();
    double mode    = posterior->GetBinCenter( modeBin );
    std::cout << "Mode: " << mode << std::endl << std::endl;
    // This is just for drawing purposes
    TH1D* modehisto = new TH1D("mode","mode",nBins,min,max );
    modehisto->SetLineColor(kRed);
    modehisto->SetBinContent( modeBin, posterior->GetBinContent(modeBin) );
    
    // Define the quantile required for the interval
    double quantile = 0.6826895;// Corresponds to +-1 sigma for a Gaussian distribution
    
    // Central interval
    TH1D* central = new TH1D( "Central interval", "Central interval", nBins, min, max );
    central->SetFillColor(kBlue-5);
    double integral = posterior->GetBinContent( modeBin );
    int binShift = 1;
    while( integral < quantile )
	{
	    integral += posterior->GetBinContent( modeBin + binShift );
	    integral += posterior->GetBinContent( modeBin - binShift );
	    binShift ++;
	}
    double errorCI = binWidth * ( 2. * binShift + 1 );
    std::cout << "Symmetric error with central interval: " << errorCI << std::endl;
    std::cout << "Coverage for central interval: "
	      << posterior->Integral( modeBin-binShift, modeBin+binShift )
	      << std::endl << std::endl;
    
    for( int b=modeBin-binShift; b<=modeBin+binShift; b++ )
	central->SetBinContent( b, posterior->GetBinContent(b) );

    // Equal areas
    TH1D* equalAreas = new TH1D( "Equal areas", "Equal areas", nBins, min, max );
    equalAreas->SetFillColor(kBlue-5);
    integral = 0.;
    double integralThr = 0.5 * ( 1. - quantile );
    int b = 1;
    while( integral < integralThr )
	{
	    integral += posterior->GetBinContent(b);
	    b++;
	}
    int lowEAErrorBin = b-1;
    integral = 0.;
    b = nBins;
    while( integral < integralThr )
	{
	    integral += posterior->GetBinContent(b);
	    b--;
	}
    int highEAErrorBin = b+1;
    for( b=lowEAErrorBin; b<=highEAErrorBin; b++ )
	equalAreas->SetBinContent( b, posterior->GetBinContent(b) );
    
    double errEALow  = binWidth * ( modeBin - lowEAErrorBin );
    double errEAHigh = binWidth * ( highEAErrorBin - modeBin);
    std::cout << "Low error for equal-areas method: " << errEALow << std::endl;
    std::cout << "High error for equal-areas method: " << errEAHigh << std::endl;
    std::cout << "Coverage for equal-areas method: "
	      << posterior->Integral( lowEAErrorBin, highEAErrorBin )
	      << std::endl << std::endl;
    
    // Shortest interval
    TH1D* shortest = new TH1D( "Shortest", "Shortest", nBins, min, max );
    shortest->SetFillColor(kBlue-5);
    int binLeft = modeBin;
    int binRight = modeBin;
    integral = posterior->GetBinContent(modeBin);
    while( integral < quantile )
	{
	    double intLeft  = posterior->GetBinContent( binLeft - 1 );
	    double intRight = posterior->GetBinContent( binRight + 1 );
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
    for( b=binLeft; b<=binRight; b++ )
	shortest->SetBinContent( b, posterior->GetBinContent(b) );
    
    double errShortLeft  = binWidth * ( modeBin - binLeft );
    double errShortRight = binWidth * ( binRight - modeBin );
    std::cout << "Low error for shortest interval method: " << errShortLeft << std::endl;
    std::cout << "High error for shortest interval method: " << errShortRight << std::endl;
    std::cout << "Coverage for equal-areas method: "
	      << posterior->Integral( binLeft, binRight )
	      << std::endl << std::endl;
    
    // Draw
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 2400, 900 );
    can->Divide(3,1);
    can->cd(1);
    posterior->Draw();
    central->Draw("same");
    modehisto->Draw("same");
    posterior->Draw("AXIS SAME");
    can->cd(2);
    posterior->Draw();
    equalAreas->Draw("same");
    modehisto->Draw("same");
    posterior->Draw("AXIS SAME");
    can->cd(3);
    posterior->Draw();
    shortest->Draw("same");
    modehisto->Draw("same");
    posterior->Draw("AXIS SAME");
    app->Run();
    
    return 0;
}
