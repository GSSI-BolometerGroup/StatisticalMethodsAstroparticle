#include "TFile.h"
#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"

#include <BAT/BCLog.h>

#include "RadioactiveDecayFit.h"

void Normalize( TH1D* h )
{
    double integral = h->Integral( 1, h->GetNbinsX() );
    h->Scale( 1./integral );
    for( int b=1; b<h->GetNbinsX(); b++ )
	h->SetBinError(b,0);
    return;
}

int main()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // create new RadioactiveDecayFit object
    RadioactiveDecayFit m("RadioactiveDecayFit");
    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);
    BCLog::OutSummary("Test model created");


    // -----------------------
    // Run Binomial likelihood
    m.SetFitMethod(RadioactiveDecayFit::FitMethod::kBinomial);
    BCLog::OutSummary("Running Binomial likelihood");
    // Run the Metropolis-Hastings algorithm
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);
    // Run Minuit to find the global mode, initializing the parameters to the "best fit"
    // from Metropolis-Hastings. This helps the convergence in some cases.
    m.FindMode(m.GetBestFitParameters());
    // Produce a pdf with the marginalized distributions
    m.PrintAllMarginalized(m.GetSafeName() + "_Binomial_plots.pdf");
    m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_Binomial_update.pdf");
    // Get posterior distribution for parameter N
    TH1D* h_N_binomial = (TH1D*)m.GetMarginalizedHistogram("N")->Clone();
    // Normalize the posterior of N to 1
    Normalize( h_N_binomial );
    // Get the best fit curve
    TF1* bestfit_binomial = m.GetBestFit();
    bestfit_binomial->SetLineColor(1);
    bestfit_binomial->SetLineWidth(2);
    // Print some infos on the fit outome, interval estimation, etcetera
    m.PrintSummary();
    m.ResetResults();
    // -----------------------
    // Run Multinomial likelihood
    m.SetFitMethod(RadioactiveDecayFit::FitMethod::kMultinomial);
    BCLog::OutSummary("Running Multinomial likelihood");
    // Run the Metropolis-Hastings algorithm
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);
    // Run Minuit to find the global mode, initializing the parameters to the "best fit"
    // from Metropolis-Hastings. This helps the convergence in some cases.
    m.FindMode(m.GetBestFitParameters());
    // Produce a pdf with the marginalized distributions
    m.PrintAllMarginalized(m.GetSafeName() + "_Multinomial_plots.pdf");
    m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_Multinomial_update.pdf");
    // Get posterior distribution for parameter N
    TH1D* h_N_multinomial = (TH1D*)m.GetMarginalizedHistogram("N")->Clone();
    // Normalize the posterior of N to 1
    Normalize( h_N_multinomial );
    h_N_multinomial->SetLineColor(4);
    // Get the best fit curve
    TF1* bestfit_multinomial = m.GetBestFit();
    bestfit_multinomial->SetLineColor(4);
    bestfit_multinomial->SetLineWidth(1);
    // Print some infos on the fit outome, interval estimation, etcetera
    m.PrintSummary();
    m.ResetResults();
    // ----------------------
    // Run Poisson likelihood
    //
    // Do the same as above, but using a Poisson extended-likelihood term
    m.SetFitMethod(RadioactiveDecayFit::FitMethod::kPoisson);
    BCLog::OutSummary("Running Poisson likelihood");
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);
    m.FindMode(m.GetBestFitParameters());
    m.PrintAllMarginalized(m.GetSafeName() + "_Poisson_plots.pdf");
    m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_Poisson_update.pdf");
    TH1D* h_N_poisson = (TH1D*)m.GetMarginalizedHistogram("N")->Clone();
    Normalize( h_N_poisson );
    h_N_poisson->SetLineColor(2);
    TF1* bestfit_poisson = m.GetBestFit();
    bestfit_poisson->SetLineStyle(7);
    m.PrintSummary();

    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    // ----
    // Draw

    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->Divide(2,1);
    can->cd(1);
    h_N_binomial->Draw();
    h_N_multinomial->Draw("same");
    h_N_poisson->Draw("same");
    can->cd(2);
    m.GetData()->Draw();
    TF1* func = new TF1("func","[0]*[2]*log(2.)/[1]*exp(-log(2.)*x/[1])",0,m.GetData()->GetXaxis()->GetXmax());
    func->SetParameter(0,1000*m.GetData()->GetBinWidth(1));
    func->SetParameter(1,138.376);
    func->FixParameter(2,m.GetData()->GetBinWidth(1));
    func->SetLineColor(3);
    m.GetData()->Fit(func,"RL");	
    bestfit_binomial->Draw("same");
    bestfit_multinomial->Draw("same");
    bestfit_poisson->Draw("same");
    app->Run();

    return 0;
}
