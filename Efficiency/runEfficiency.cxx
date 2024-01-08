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
    BCLog::OpenLog("log.txt", BCLog::error, BCLog::error);
    //BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    std::vector<int>    N{100,1000,10000};
    std::vector<double> P{0.9,0.95,0.99};
    int M = 1000;// number of toy-MC experiments
    
    size_t nCases = N.size() * P.size();
    std::vector<TH1D*> h_eff_binomial( nCases );
    std::vector<TH1D*> h_eff_poisson ( nCases );
    std::vector<TH1D*> h_eff_chi2    ( nCases );
    std::vector<TH1D*> h_threshold_binomial( nCases );
    std::vector<TH1D*> h_threshold_poisson ( nCases );
    std::vector<TH1D*> h_threshold_chi2    ( nCases );
    std::vector<TH1D*> h_sigma_binomial( nCases );
    std::vector<TH1D*> h_sigma_poisson ( nCases );
    std::vector<TH1D*> h_sigma_chi2    ( nCases );
    
    int i=0;
    std::string h_name;
    int nBins = 1000;
    for( auto& n: N )
	for( auto& p: P )
	    {
		std::cout << n << "\t" << p << std::endl;
 
		h_name = "Eff-binomial n=" + std::to_string(n) + " p=" + std::to_string(p);
		h_eff_binomial[i] = new TH1D( h_name.c_str(), h_name.c_str(), nBins, 0.8, 1 );
		h_name = "Eff-Poisson n=" + std::to_string(n) + " p=" + std::to_string(p);
		h_eff_poisson[i]  = new TH1D( h_name.c_str(), h_name.c_str(), nBins, 0.8, 1 );
		h_name = "Eff-Chi2 n=" + std::to_string(n) + " p=" + std::to_string(p);
		h_eff_chi2[i]     = new TH1D( h_name.c_str(), h_name.c_str(), nBins, 0.8, 1 );
		h_name = "Threshold-binomial n=" + std::to_string(n) + " p=" + std::to_string(p);
		h_threshold_binomial[i] = new TH1D( h_name.c_str(), h_name.c_str(), nBins, 0.8, 1 );
		h_name = "Threshold-Poisson n=" + std::to_string(n) + " p=" + std::to_string(p);
		h_threshold_poisson[i]  = new TH1D( h_name.c_str(), h_name.c_str(), nBins, 0.8, 1 );
		h_name = "Threshold-Chi2 n=" + std::to_string(n) + " p=" + std::to_string(p);
		h_threshold_chi2[i]     = new TH1D( h_name.c_str(), h_name.c_str(), nBins, 0.8, 1 );
		h_name = "Sigma-binomial n=" + std::to_string(n) + " p=" + std::to_string(p);
		h_sigma_binomial[i] = new TH1D( h_name.c_str(), h_name.c_str(), nBins, 0.8, 1 );
		h_name = "Sigma-Poisson n=" + std::to_string(n) + " p=" + std::to_string(p);
		h_sigma_poisson[i]  = new TH1D( h_name.c_str(), h_name.c_str(), nBins, 0.8, 1 );
		h_name = "Sigma-Chi2 n=" + std::to_string(n) + " p=" + std::to_string(p);
		h_sigma_chi2[i]     = new TH1D( h_name.c_str(), h_name.c_str(), nBins, 0.8, 1 );

		h_eff_binomial[i]      ->SetLineColor(1);
		h_threshold_binomial[i]->SetLineColor(1);
		h_sigma_binomial[i]    ->SetLineColor(1);
		h_eff_poisson[i]      ->SetLineColor(2);
		h_eff_poisson[i]->SetFillColor(2);
		h_threshold_poisson[i]->SetLineColor(2);
		h_sigma_poisson[i]    ->SetLineColor(2);
		h_eff_chi2[i]      ->SetLineColor(4);
		h_threshold_chi2[i]->SetLineColor(4);
		h_sigma_chi2[i]    ->SetLineColor(4);

		double eff_binomial_weight = 0;
		double eff_poisson_weight  = 0;
		double eff_chi2_weight     = 0;
		
		for( int m=0; m<M; m++ )	    
		    {
			
			// create new Efficiency object
			Efficiency model( "Efficiency", n, p );
			
			// set precision
			//model.SetPrecision(BCEngineMCMC::kHigh);
			
			//BCLog::OutSummary("Test model created");
			model.SetFitMethod( Efficiency::FitMethod::kBinomial );
			//model.MarginalizeAll(BCIntegrate::kMargMetropolis);
			model.FindMode();
			//model.FindMode(model.GetBestFitParameters());
			//model.PrintSummary();
			h_eff_binomial[i]      ->Fill( model.GetBestFitParameters()[0],
						       1./pow(model.GetBestFitParameterErrors()[0],2.) );
			h_threshold_binomial[i]->Fill( model.GetBestFitParameters()[1],
						       1./pow(model.GetBestFitParameterErrors()[1],2.) );
			h_sigma_binomial[i]    ->Fill( model.GetBestFitParameters()[2],
						       1./pow(model.GetBestFitParameterErrors()[2],2.) );
			eff_binomial_weight += 1./pow(model.GetBestFitParameterErrors()[0],2.);
		    
			model.SetFitMethod( Efficiency::FitMethod::kPoisson );
			model.FindMode();
			//model.FindMode(model.GetBestFitParameters());
			//model.PrintSummary();
			h_eff_poisson[i]      ->Fill( model.GetBestFitParameters()[0],
						      1./pow(model.GetBestFitParameterErrors()[0],2.) );
			h_threshold_poisson[i]->Fill( model.GetBestFitParameters()[1],
						      1./pow(model.GetBestFitParameterErrors()[1],2.) );
			h_sigma_poisson[i]    ->Fill( model.GetBestFitParameters()[2],
						      1./pow(model.GetBestFitParameterErrors()[2],2.) );
			eff_poisson_weight += 1./pow(model.GetBestFitParameterErrors()[0],2.);
			
			model.SetFitMethod( Efficiency::FitMethod::kChiSquare );
			model.FindMode();
			//model.FindMode(model.GetBestFitParameters());
			//model.PrintSummary();
			h_eff_chi2[i]      ->Fill( model.GetBestFitParameters()[0],
						   1./pow(model.GetBestFitParameterErrors()[0],2.) );
			h_threshold_chi2[i]->Fill( model.GetBestFitParameters()[1],
						   1./pow(model.GetBestFitParameterErrors()[1],2.) );
			h_sigma_chi2[i]    ->Fill( model.GetBestFitParameters()[2],
						   1./pow(model.GetBestFitParameterErrors()[2],2.) );
			eff_chi2_weight += 1./pow(model.GetBestFitParameterErrors()[0],2.);
		    }

		h_eff_binomial[i]->Scale(1./eff_binomial_weight);
		h_eff_poisson[i]->Scale(1./eff_poisson_weight);
		h_eff_chi2[i]->Scale(1./eff_chi2_weight);
		i++;
	    }

    // Draw
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->Divide(P.size(),N.size());
    for( size_t c=0; c<nCases; c++ )
	{
	    can->cd(c+1);
	    for( int b=1; b<=nBins; b++ )
		{
		    h_eff_poisson[c]->SetBinError(b,0);
		    h_eff_binomial[c]->SetBinError(b,0);
		    h_eff_chi2[c]->SetBinError(b,0);
		}
	    h_eff_poisson[c]->Draw("");
	    h_eff_binomial[c]->Draw("same");
	    h_eff_chi2[c]->Draw("same");
	}

    can->SaveAs("Efficiency.jpg");
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
