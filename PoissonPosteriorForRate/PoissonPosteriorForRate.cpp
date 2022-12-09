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
#include "TGaxis.h"

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


    // Create Poisson distribution with lambda=5
    double lambda = 5.;
    double min = -0.5;
    double max = 25.5;
    int nBins = max-min;
    TH1D* poissonHisto = new TH1D("Poisson", "Poisson", nBins, min, max );
    poissonHisto->GetXaxis()->SetTitle("n");
    poissonHisto->GetYaxis()->SetTitle("P(n|#lambda=5)");
    poissonHisto->GetYaxis()->SetRangeUser(0,0.2);
    for( int b=1; b<=nBins; b++ )
	poissonHisto->SetBinContent( b, SmartPoisson( poissonHisto->GetBinCenter(b), lambda ) );

    // Now create posterior for n=5
    double n=5;
    min = 0.;// Lambda defined >= 0
    max = 25.;
    nBins = 10000;// Reduce bin size to visualize continuous behaviour
    TH1D* posterior = new TH1D("Posterior","Posterior", nBins, min, max );
    posterior->SetLineColor(2);
    for( int b=1; b<=nBins; b++ )
	{
	    double x = posterior->GetBinCenter(b);
	    posterior->SetBinContent( b, SmartPoisson( n, x ) );
	}
    std::cout << "Posterior mean: " << posterior->GetMean() << std::endl;
    std::cout << "Posterior variance:  " << pow(posterior->GetStdDev(),2.) << std::endl;

	
    // Draw
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->cd();
    poissonHisto->Draw();
    can->Update();

    // Stuff for drawing second set of axes
    double xmax = gPad->GetUxmax();///poissonHisto->GetXaxis()->GetXmax();
    double ymax = 0.2;//1.1*posterior->GetMaximum();
    //posterior->Scale(scale);
    //poissonHisto->GetYaxis()->SetRangeUser(ymin,ymax);

    
    posterior->Draw("SAME");

    TGaxis *xaxis = new TGaxis( gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax(), 0, xmax, 510, "-" );
    xaxis->SetLineColor(2);
    xaxis->SetLabelColor(2);
    xaxis->SetTitleColor(2);
    xaxis->SetTitle("#lambda");
    xaxis->Draw();
    
    TGaxis *yaxis = new TGaxis( gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), 0, ymax, 510, "+L" );
    yaxis->SetLineColor(2);
    yaxis->SetLabelColor(2);
    yaxis->SetTitleColor(2);
    yaxis->SetTitle("P(#lambda|n=5)");
    yaxis->Draw();

    app->Run(kTRUE);
    
    return 0;
}
