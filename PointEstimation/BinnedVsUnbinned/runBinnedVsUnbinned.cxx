#include <chrono>
#include <string>

#include "TApplication.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TStyle.h"

#include <BAT/BCLog.h>

#include "BinnedVsUnbinned.h"

// In this example, we perform a maximum-likelihood fit of data composed
// by S Gaussian-distributed signal events and B flat-distributed background events in the 2000-2080 keV range.
// We perform the fit twice:
// 1) as an unbinned fit using the extended-likelihood method;
// 2) as a binned fit, using a Poisson likelihood for each bin.

int main()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::error, BCLog::error);

    std::vector<int>* S = new std::vector<int>{5,20,50,100};
    std::vector<int>* B = new std::vector<int>{10,30,100,200};

    size_t nCases = S->size() * B->size();
    int M = 1000;
    
    // create new BinnedVsUnbinned object
    
    BCLog::OutSummary("Test model created");

    std::vector<double>* N    = new std::vector<double>(nCases);
    std::vector<double>* Nerr = new std::vector<double>(nCases);
    std::vector<double>* R    = new std::vector<double>(nCases);
    std::vector<double>* Rerr = new std::vector<double>(nCases);
    
    int n = 0;
    double nBins;
    for( auto& s: *S )
	for( auto& b: *B )
	    {

		std::vector<double> ratio(M);
		for( int m=0; m<M; m++ )
		    {
			BinnedVsUnbinned model("BinnedVsUnbinned");
			// populate the data
			//model.SetData( S[s], B[b] );
			model.SetData( s, b );
			nBins = model.GetNBins();
			
			// Extended-likelihood fit
			auto startUnbinned = std::chrono::high_resolution_clock::now();
			model.SetType( BinnedVsUnbinned::Type::kUnbinned );
			model.FindMode();
			auto endUnbinned = std::chrono::high_resolution_clock::now();
			double durationUnbinned = std::chrono::duration_cast<std::chrono::microseconds>(endUnbinned - startUnbinned).count();

			auto startBinned = std::chrono::high_resolution_clock::now();
			model.ResetResults();
			model.SetType( BinnedVsUnbinned::Type::kBinned );
			model.FindMode();
			auto endBinned = std::chrono::high_resolution_clock::now();
			double durationBinned = std::chrono::duration_cast<std::chrono::microseconds>(endBinned - startBinned).count();

			ratio[m] = durationUnbinned/durationBinned;

		    }

		double averageRatio = std::reduce( ratio.begin(), ratio.end() ) / M;

		double variance = 0.;
		for( auto& r: ratio )
		    variance += std::pow( r - averageRatio, 2. );
		double errorRatio = sqrt( variance / M );
		
		N->at(n)    = (double)(s+b)/nBins;//S[s]+B[b];
		Nerr->at(n) = 0;
		R->at(n)    = averageRatio;
		Rerr->at(n) = errorRatio;
	  		
		n++;
	    }

    TGraphErrors* gr_ratio = new TGraphErrors( nCases, &(N->at(0)), &(R->at(0)), &(Nerr->at(0)), &(Rerr->at(0)) );
    gr_ratio->GetXaxis()->SetTitle( "#frac{N_{events}}{N_{bins}}" );
    gr_ratio->GetYaxis()->SetTitle( "#frac{t_{unbinned}}{t_{binned}}" );
    gr_ratio->GetYaxis()->SetRangeUser( 0, 3 );
    
    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    // Draw
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->cd();
    can->SetGridy();
    gr_ratio->Draw("AP");
    app->Run();
    
    return 0;
}
