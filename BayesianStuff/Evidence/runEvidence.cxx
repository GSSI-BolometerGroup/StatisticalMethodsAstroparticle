// ***************************************************************
// This file was created using the bat-project script
// for project Evidence.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"

#include "Evidence.h"

int main()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    int S = 5;
    int B = 25;
    
    // create new Evidence object
    Evidence m( "Evidence", S, B );
    
    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);

    BCLog::OutSummary("Test model created");

    // Fit with bkg only
    m.SetModel( Evidence::Model::kBkgOnly );
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);
    m.FindMode(m.GetBestFitParameters());
    m.PrintAllMarginalized(m.GetSafeName() + "_BkgOnly_plots.pdf");
    m.PrintSummary();
    m.ComputeIntegral();
    
    m.ResetResults();
    
    // Fit with signal + background
    m.SetModel( Evidence::Model::kSgnPlusBkg );
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);
    m.FindMode(m.GetBestFitParameters());
    m.PrintAllMarginalized(m.GetSafeName() + "_SgnPlusBkg_plots.pdf");
    m.PrintSummary();
    m.ComputeIntegral();

    m.ComputeBayesFactor();
    
    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
