// ***************************************************************
// This file was created using the bat-project script
// for project Efficiency.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************
#include <vector>
#include <string>

#include <BAT/BCLog.h>

#include "TH1D.h"
#include "TCanvas.h"
#include "TApplication.h"

#include "Efficiency.h"

// In this example, we simulate a pulse generator injecting artificial pulses into some detector
// and we run a Bayesian fit with different likelihood functions to evaluate the trigger efficiency.
// We inject pulses between 1 and 20 keV, with 1 keV steps.
// For each energy, we inject 1e3 events.
// Specifically, we generate 1e3 values with a uniform distribution between 0 and 1,
// compare them with the theoretical efficiency curve, implemented as an error function with mu=3keV and sigma=1keV,
// and count the number of accepted events k.
// Then we fit k as a function of energy using three types of fits:
// 1) Likelihood with a binomial term for each point
// 2) Likelihood with a Poisson term for each point
// 3) Chi2
// Finally, we extract the posterior for the efficiency obtained from each of the 3 fits,
// and plot them all together.

int main()
{
    // open log file
    //BCLog::OpenLog("log.txt", BCLog::error, BCLog::error);
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // Number of injected pulser events
    int N = 1000;
    // Input trigger efficiency (asymptotic value)
    double P = 0.99;
    
    
    // create new Efficiency object
    Efficiency model( "Efficiency", N, P );
			
    // set precision of Markov Chain
    model.SetPrecision(BCEngineMCMC::kMedium);
    BCLog::OutSummary("Test model created");
    
    // First fit: binomial likelihood
    model.SetFitMethod( Efficiency::FitMethod::kBinomial );
    // Run Metropolis-Hastings
    model.MarginalizeAll(BCIntegrate::kMargMetropolis);
    // Find global mode with Minuit, initializing it to the approximate global mode
    // obtained from Metropolis-Hastings
    model.FindMode(model.GetBestFitParameters());
    // Print some information on the fit
    model.PrintSummary();
    // Retrieve the marginalized distribution of the efficiency
    TH1D* h_Eff_binomial = (TH1D*)model.GetMarginalizedHistogram("Efficiency")->Clone();
    // Reset the model
    model.ResetResults();

    // Repeat with Poisson likelihood
    model.SetFitMethod( Efficiency::FitMethod::kPoisson );
    model.MarginalizeAll(BCIntegrate::kMargMetropolis);
    model.FindMode(model.GetBestFitParameters());
    model.PrintSummary();
    TH1D* h_Eff_poisson = (TH1D*)model.GetMarginalizedHistogram("Efficiency")->Clone();
    model.ResetResults();

    // Repeat with Chi2 fit
    model.SetFitMethod( Efficiency::FitMethod::kChiSquare );
    model.MarginalizeAll(BCIntegrate::kMargMetropolis);
    model.FindMode(model.GetBestFitParameters());
    model.PrintSummary();
    TH1D* h_Eff_ChiSquare = (TH1D*)model.GetMarginalizedHistogram("Efficiency")->Clone();
    model.ResetResults();

    // Draw
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->cd();
    h_Eff_binomial->SetLineColor(1);
    h_Eff_poisson->SetLineColor(2);
    h_Eff_ChiSquare->SetLineColor(3);
    h_Eff_binomial->Draw();
    h_Eff_poisson->Draw("same");
    h_Eff_ChiSquare->Draw("same");
    app->Run(kTRUE);

    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
