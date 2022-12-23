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

// In this example, we generate two variables x and y with standard normal distributions
// and compute the distribution of their ratio z=x/y with two methods:
// 1) generating random numbers for x and y, and filling up a histogram with z=x/y;
// 2) generating the theoretical Cauchy distribution for z.

int main()
{
    // Random number machinery
    std::random_device rd;
    std::mt19937 gen(rd());

    // Number of events to be generated
    int N = 100000;
    int nBins = 1000;

    // Prepare Gaussian distribution for x
    double muX    = 0.;
    double sigmaX = 1.;
    double minX  = muX - 5.*sigmaX;
    double maxX  = muX + 5.*sigmaX;
    std::normal_distribution<> gausX( muX, sigmaX );
    TH1D* hX = new TH1D( "x", "x", nBins, minX, maxX );
    hX->GetXaxis()->SetTitle( "x [bananas]" );
    hX->GetYaxis()->SetTitle( "P(x)" );
    // This is the theoretical distribution for x, which is used exclusively for plotting.
    TF1* fX = new TF1( "fX", "[0]*[3]/sqrt(2.*TMath::Pi()) /[2] * exp( -0.5 * pow( ( x-[1] )/[2], 2. ) )", minX, maxX );
    fX->SetParameter( 0, N );
    fX->SetParameter( 1, muX );
    fX->SetParameter( 2, sigmaX );
    fX->SetParameter( 3, hX->GetBinWidth(1) );

    // Prepare distribution for y
    double muY    = 0.;
    double sigmaY = 1.;
    double minY  = muY - 5.*sigmaY;
    double maxY  = muY + 5.*sigmaY;
    std::normal_distribution<> gausY( muY, sigmaY );
    TH1D* hY = new TH1D( "y", "y", nBins, minY, maxY );
    hY->GetXaxis()->SetTitle( "y [mangos]" );
    hY->GetYaxis()->SetTitle( "P(y)" );
    TF1* fY = new TF1( "fY", "[0]*[3]/sqrt(2.*TMath::Pi()) /[2] * exp( -0.5 * pow( ( x-[1] )/[2], 2. ) )", minY, maxY );
    fY->SetParameter( 0, N );
    fY->SetParameter( 1, muY );
    fY->SetParameter( 2, sigmaY );
    fY->SetParameter( 3, hY->GetBinWidth(1) );

    // Prepare distribution for z=x/y
    double minZ = -10.;
    double maxZ = +10.;
    TH1D* hZ = new TH1D( "z=x/y", "z=x/y", nBins, minZ, maxZ );
    hZ->GetXaxis()->SetTitle( "z [bananas/mangos]" );
    hZ->GetYaxis()->SetTitle( "P(z)" );
    // The theoretical distribution for z is a Cauchy
    TF1* fZ = new TF1( "fZ", "[0]*[1]/TMath::Pi()/( 1. + pow( x, 2. ) )", minZ, maxZ );
    fZ->SetParameter( 0, N );
    fZ->SetParameter( 1, hZ->GetBinWidth(1) );
    fZ->SetNpx(nBins);

    // Generate random values for x and y, and compute z.
    // Fill all three histograms.
    for( int n=0; n<N; n++ )
	{
	    double x = gausX(gen);
	    double y = gausY(gen);
	    double z = x / y;

	    hX->Fill(x);
	    hY->Fill(y);
	    hZ->Fill(z);
	}

    // Draw
    gStyle->SetOptStat(0);
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->Divide(3,1);
    can->cd(1);
    hX->Draw();
    fX->Draw("same");
    can->cd(2);
    hY->Draw();
    fY->Draw("same");
    can->cd(3);
    hZ->Draw();
    fZ->Draw("same");
    app->Run();
    
    return 0;
}
