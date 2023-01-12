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

int main()
{
    // open log file
    //BCLog::OpenLog("log.txt", BCLog::error, BCLog::error);
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    int N = 1000;
    double P = 0.99;
    
    
    // create new Efficiency object
    Efficiency model( "Efficiency", N, P );
			
    // set precision of Markov Chain
    model.SetPrecision(BCEngineMCMC::kMedium);
			
    BCLog::OutSummary("Test model created");
    model.SetFitMethod( Efficiency::FitMethod::kBinomial );
    model.MarginalizeAll(BCIntegrate::kMargMetropolis);
    // Find global mode with Minuit, initializing it to the approximate global mode
    // obtained from Metropolis-Hastings
    model.FindMode(model.GetBestFitParameters());
    model.PrintSummary();
    TH1D* h_Eff_binomial = (TH1D*)model.GetMarginalizedHistogram("Efficiency")->Clone();
    model.ResetResults();
    
    model.SetFitMethod( Efficiency::FitMethod::kPoisson );
    model.MarginalizeAll(BCIntegrate::kMargMetropolis);
    model.FindMode(model.GetBestFitParameters());
    model.PrintSummary();
    TH1D* h_Eff_poisson = (TH1D*)model.GetMarginalizedHistogram("Efficiency")->Clone();
    model.ResetResults();
    
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
    //////////////////////////////
    // perform your analysis here

    // Normalize the posterior by integrating it over the full parameter space
    // model.Normalize();

    // Write Markov Chain to a ROOT file as a TTree
    // model.WriteMarkovChain(model.GetSafeName() + "_mcmc.root", "RECREATE");

    // run MCMC, marginalizing posterior
    //model.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // run mode finding; by default using Minuit
    //model.FindMode(model.GetBestFitParameters());

    // draw all marginalized distributions into a PDF file
    //model.PrintAllMarginalized(model.GetSafeName() + "_plots.pdf");

    // print summary plots
    // model.PrintParameterPlot(model.GetSafeName() + "_parameters.pdf");
    // model.PrintCorrelationPlot(model.GetSafeName() + "_correlation.pdf");
    // model.PrintCorrelationMatrix(model.GetSafeName() + "_correlationMatrix.pdf");
    // model.PrintKnowledgeUpdatePlots(model.GetSafeName() + "_update.pdf");

    // print results of the analysis into a text file
    //model.PrintSummary();

    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
