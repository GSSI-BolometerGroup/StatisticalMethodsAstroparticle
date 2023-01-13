#include <BAT/BCLog.h>

#include "SimultaneousFit.h"

// In this example, we perform a combined fit between two data sets:
// 1- Pulser injected data, produced in the same way as in the "Efficiency" example,
//    i.e. by simulating a pulse generator that injects 1e3 events at 20 different energies
//    between 1 and 20 keV, with 1 keV steps, and simulating a trigger efficiency
//    described by a Gaussian CDF with mu=3keV, sigma=1keV and asymptotic efficiency=0.9.
// 2- Physics data, recorded between 10 and 20 keV, and featuring:
//    a) a Gaussian signal peak with mu=15keV and sigma=1.2keV
//    b) a flat background.
//    We generate 20 signal and 20 background events, and apply the efficiency curve to them.
// We then perform a simultaneous fit on the two datasets, with the efficiency as a common parameter.

int main()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    int    nInjected = 1000;
    double eff = 0.9;
    int    S = 20;
    int    B = 20;
    
    // create new SimultaneousFit object
    SimultaneousFit m( "SimultaneousFit", nInjected, eff, S, B );

    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);

    BCLog::OutSummary("Test model created");

    // run MCMC, marginalizing posterior
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // run mode finding; by default using Minuit
    m.FindMode(m.GetBestFitParameters());

    // draw all marginalized distributions into a PDF file
    m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

    // print results of the analysis into a text file
    m.PrintSummary();

    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
