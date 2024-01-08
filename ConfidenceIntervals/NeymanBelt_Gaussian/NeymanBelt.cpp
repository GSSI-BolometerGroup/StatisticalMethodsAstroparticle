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

// In this example, we construct the Neyman belt for a variable x which depends
// on a parameter theta, with relation x=theta/2,
// and distributed according to a Gaussian PDF with sigma=1.
// We proceed as follows:
// 1- create a 2-dim histogram (TH2D)
// 2- for each value of theta, populate x with N Gaussian-distributed random numbers
// 3- invert the confidence belt for a fixed value of x
// 4- fit the confidence belt of theta.

int main()
{
    double dx    = 0.01;// bin width along x
    double dTheta = 2.*dx;// bin width along theta
    double minX  = 0.-0.5*dx;
    double maxX  = 10.+0.5*dx;
    double minTheta = 2. * minX;
    double maxTheta = 2. * maxX;
    int    nBins = (maxX-minX)/dx+0.5;
    double sigma = 1.;
    int    N     = 100000;// number of x values generated for each value of theta

    // Prepare the 2-dim histogram
    TH2D* belt = new TH2D( "Neyman belt", "Neyman belt", nBins, minX, maxX, nBins, minTheta, maxTheta );
    belt->GetXaxis()->SetTitle("x [bananas]");
    belt->GetYaxis()->SetTitle("#theta [bananas]");

    // Prepare random number machinery
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<> dis(0.,sigma);

    // Populate the Neyman belt
    double shift;
    for( int b=1; b<=nBins; b++ )
	{
	    double mu = dx * (b-1);
	    for( int i=0; i<N; i++ )
		{
		    shift = dis(gen);
		    belt->Fill( mu+shift, 2.*mu );
		}
	}

    // We'll invert the Neyman confidence belt for the middle bin.
    int projectionBin = nBins/2 + 1;
    // Invert the confidence belt using the "ProjectionY" method of the TH1D class
    TH1D* theta = belt->ProjectionY( "theta", projectionBin, projectionBin );
    // Prepare Gaussian function to fit the inverted confidence belt
    TF1* func = new TF1("func","[0]*[3]/sqrt(2.*TMath::Pi())/[2]*exp( -0.5*pow( (x-[1])/[2], 2. ))",minTheta,maxTheta);
    func->SetParameter(0,0.9*N);
    func->SetParameter(1,0.55*(minTheta+maxTheta));
    func->SetParameter(2,1.3);
    func->FixParameter(3,dTheta);
    func->SetParName(0,"Integral");
    func->SetParName(1,"Mean");
    func->SetParName(2,"Sigma");
    func->SetParName(3,"Bin width");
    
    // Draw
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->Divide(1,2);
    can->cd(1);
    belt->SetStats(0);
    belt->Draw("colz");
    can->cd(2);
    theta->Draw();
    theta->Fit(func,"RL");

    can->SaveAs("NeymanBeltGaussian.jpg");
    
    app->Run();
    
    return 0;
}
