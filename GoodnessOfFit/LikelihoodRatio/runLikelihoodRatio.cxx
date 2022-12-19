#include <BAT/BCLog.h>

#include "TH1D.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMath.h"
#include "TStyle.h"
#include "TF1.h"

#include "LikelihoodRatio.h"

// In this example, we will generate M toy-histograms populated with a Gaussian signal over a flat background.
// We will fit them with a 2*log[(likelihood)/(max-likelihood)] following the method proposed by
// Bakers and Cousins, Nucl. Instrum. Methods 221 (1984) 437-442.
// We will approximate the max-likelihood by imposing that the expectation value for each bin
// is equal to the number of data events in that bin.
// This works if each bin has at least a few events, but fails if just one bin has zero entries.
// Maximinzing 2*log[(likelihood)/(max-likelihood)] is equivalent to maximizing the likelihood itself,
// given that max-likelihood depends only on the data and not on the parameters.
// According to Wilk's theorem, the quantity -2*log[(likelihood)/(max-likelihood)]
// follows a Chi2 distribution with a number of degrees of freedom
// equal to the number of bins minus the number of fit parameters.
// For each toy experiment, we will extract such Chi2 and compare it
// with the theoretical distribution for that number of DOF.



// This will be used to produce the theoretical Chi2 distribution
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
    // open log file, printing only errors
    BCLog::OpenLog("log.txt", BCLog::error, BCLog::error);

    // Number of toys
    int M = 100000;

    // Prepare Chi2 histogram
    int nBins = 1000;
    double maxChi2 = 100.;
    TH1D* histo_Chi2 = new TH1D( "Chi2", "Chi2", nBins, 0, maxChi2 );
    histo_Chi2->GetXaxis()->SetTitle("#Chi^{2}");
    
    int ndf;// Number of degrees of freedom

    // Histogram for one data spectrum
    TH1D* histo_data;
    // Best fitting function for one data spectrum
    TF1* bestfit;
    
    // Loop over toys
    for( int m=0; m<M; m++ )
	{
	    // create ne1w LikelihoodRatio object
	    LikelihoodRatio model("LikelihoodRatio");
	    // Create the data
	    model.SetData();
	    // Fit with Minuit using custom function defined in LikelihoodRatio class
	    model.FindMode();

	    // Fill Chi2 histogram
	    histo_Chi2->Fill(  -model.LogLikelihood(model.GetBestFitParameters()) );
	    // Get number of degrees of freedom (just for the first toy, given that it's the same for all)
	    if( m == 0 )
		{
		    ndf = model.GetNDF();
		    histo_data = model.GetDataHisto();
		    bestfit = model.GetFittingFunction();
		}
	}

    // Compute the theoretical Chi2 for (N-2) degrees of freedom
    TH1D* theoreticalChi2 = new TH1D("theoreticalChi2", "theoreticalChi2", nBins,0,maxChi2 );
    for( int b=1; b<=1000; b++ )
	theoreticalChi2->SetBinContent( b, M*theoreticalChi2->GetBinWidth(1)*ChiSquare( theoreticalChi2->GetBinCenter(b), ndf ) );
    theoreticalChi2->SetLineColor(2);
    theoreticalChi2->SetLineWidth(2);
    
    // Draw
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->Divide(1,2);
    can->cd(1);
    // Draw data and fitting function for first fit
    histo_data->Draw();
    bestfit->Draw("same");
    can->cd(2);
    // Draw Chi distributions
    histo_Chi2->Draw();
    theoreticalChi2->Draw("SAME");
    
    app->Run();

    
    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
