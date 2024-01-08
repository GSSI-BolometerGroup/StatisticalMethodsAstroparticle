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

int main()
{
    // Prepare random number machinery
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    double minT = 0.;
    double maxT = 1.e3;// This is just a "random" choice, and has no influence on the result
    std::uniform_real_distribution<> dis(minT,maxT);

    // Generate N time values uniformly distributed between min and max
    int N=1000000;
    std::vector<double> time(N);
    for( auto& t: time )
	t = dis(gen);

    // Sort the list of time values
    std::sort( time.begin(), time.end() );

    // Compute the delta-t between subsequent events
    std::vector<double> deltaT(N);
    double prevT = 0;
    for( int i=0; i<N; i++ )
	{
	    deltaT[i] = time[i] - prevT;
	    prevT = time[i];
	}

    // Create and fill up histogram with delta-t values
    int nBins = 1000;
    double minDeltaT = 0;
    double maxDeltaT = 1.01 * ( *std::max_element( deltaT.begin(), deltaT.end() ) );// Set the histogram range to accomodate all delta-t entries
    TH1D* histo = new TH1D( "DeltaT", "DeltaT", nBins, minDeltaT, maxDeltaT );
    histo->GetXaxis()->SetTitle("#delta t [bananas]");
    histo->GetYaxis()->SetTitle("Entries");
    for( auto& t:deltaT )
	histo->Fill(t);

    // Define an exponential function to fit the histogram: r_1*N*w*exp(-t/r_2)
    // Parameter 0 is the rate (the one outside the exponential)
    // Parameter 1 is fixed to the number of generated events N
    // Parameter 2 is fixed to the bin width
    // Parameter 3 is again the rate, but this time inside the exponential.
    // In practice, we implement a function with two parameters (0 and 3) corresponding to the rate,
    // and we want to check whether the fit provides the correct number of values for them.
    TF1* func = new TF1( "func", "[0]*[1]*[2]*exp(-[3]*x)", minDeltaT, maxDeltaT );
    // Let's initialize the function parameters to some reasonable value
    func->SetParameter(0,1);
    func->FixParameter(1,(double)N);
    func->FixParameter(2,histo->GetBinWidth(1));
    func->SetParameter(3,100);
    // Define parameters name (just for plotting)
    func->SetParName( 0, "r_{1}" );
    func->SetParName( 1, "N" );
    func->SetParName( 2, "bin width" );
    func->SetParName( 3, "r_{2}" );
    // Set number of points associated to the function (just for plotting)
    func->SetNpx(10000);

    // Compute theoretical value for event rate
    double DeltaT = maxT - minT;
    double r = (double)N/DeltaT;
    std::cout << "Theoretical value for rate: " << r << std::endl;
    
    // Draw and fit
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->cd();
    histo->Draw();
    histo->Fit(func,"LLR" );
    can->SaveAs("ExponentialFromUniform.jpg");
    app->Run(kTRUE);

    
    return 0;
}
