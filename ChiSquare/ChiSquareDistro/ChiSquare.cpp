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

double ChiSquare( double x2, double n )
{
    double x = sqrt(x2);
    double p = pow( 2., -0.5 * n );
    p *= pow( x, n-2. );
    p *= exp( -0.5 * x2 );
    p /= TMath::Gamma(0.5*n);
    return p;
}

int main()
{
    // Prepare random number machinery
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<> dis(0,1);

    std::vector<double> N{1,2,3,5,10,20};
    size_t nCases = N.size();
    int M=100000;
    std::vector<TH1D*> theoreticalChi2(nCases);
    std::vector<TH1D*> dataChi2(nCases);
    for( int i=0; i<nCases; i++ )
	{
	    std::string name = "Toy-MC n=" + std::to_string((int)N[i]);
	    double min = 0;
	    double max = N[i] + 5.*sqrt(2.*N[i]);// The variance of the Chi2 distro is 2N
	    double nBins = 30. * (max-min);
	    dataChi2[i] = new TH1D( name.c_str(), name.c_str(), nBins, min, max );
	    dataChi2[i]->GetXaxis()->SetTitle("#Chi^{2}");
	    std::string axistitle = "P(#Chi^{2}|n=" + std::to_string((int)N[i]) + ")";
	    dataChi2[i]->GetYaxis()->SetTitle( axistitle.c_str() );

	    for( int m=0; m<M; m++ )
		{
		    double chi2 = 0;
		    for( int n=0; n<N[i]; n++ )
			chi2 += pow( dis(gen), 2. );
		    double weight = 1. / M / dataChi2[i]->GetBinWidth(1);
		    dataChi2[i]->Fill(chi2,weight);
		}

	    name = "Theoretical n=" + std::to_string((int)N[i]);
	    theoreticalChi2[i] = new TH1D( name.c_str(), name.c_str(), nBins, min, max );
	    theoreticalChi2[i]->SetLineColor(2);
	    for( int b=1; b<=nBins; b++ )
		theoreticalChi2[i]->SetBinContent( b, ChiSquare( theoreticalChi2[i]->GetBinCenter(b), N[i] ) );

		    
	}

    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->Divide(3,2);
    for( int i=0; i<nCases; i++ )
	{
	    can->cd(i+1);
	    dataChi2[i]->Draw();
	    theoreticalChi2[i]->Draw("SAME");
	}
    can->SaveAs("Chi2.jpg");
    app->Run(kTRUE);
    return 0;
}
