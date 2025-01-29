#include <iostream>
#include <vector>
#include <random>
#include <string>

#include "TH1D.h"
#include "TH2D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TF1.h"

// In this example, we compare the performance of the likelihood, Neyman-Chi2 and Pearson-Chi2 fits
// on a binned energy spectrum (histogram) containing S Gaussian-distributed signal events
// over B uniformly-distributed background events.
// We consider several combinations of S and B, going from 10 to 10000.
// For each combination, we generate M toy-MC spectra, and fit them with the 3 above-mentioned methods.
// We then plot the distributions of S an B and S vs B for all methods.

double GaussianPlusFlat( double *x, double* par )
{
    // par[0] = integral
    // par[1] = mean
    // par[2] = sigma
    // par[3] = flat background
    // par[4] = bin width
    // par[5] = fit range
    return par[0] * par[4] / sqrt( 2. * M_PI ) / par[2] * exp( - pow( x[0] - par[1], 2. ) / 2. / pow( par[2], 2. ) ) + par[3] / par[5] * par[4];
}

int main()
{

    double mu    = 2039.;// keV, Q-value of 76Ge (for neutrinoless double beta decay)
    double sigma = 1.5;// keV, typical resolution of a germanium detector
    double min   = 2000.;
    double max   = 2080.;
    int nBins    = 10. * ( max-min );// 0.1 keV bins, much thinner than energy resolution
    
    // Prepare histogram, which will be reset and filled for each toy-MC experiment
    TH1D* histo  = new TH1D("histo","histo",nBins,min,max);

    // Prepare fitting function
    TF1* func = new TF1( "func", GaussianPlusFlat, min, max, 6 );
    func->SetParName( 0, "Integral" );
    func->SetParName( 1, "Mean" );
    func->SetParName( 2, "Sigma" );
    func->SetParName( 3, "Background" );
    func->SetParName( 4, "Bin width" );
    func->SetParName( 5, "Fit range" );
    func->SetNpx(10000);

    // Random number machinery
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(min,max);// This is for the background distribution
    std::uniform_real_distribution<> smear(0,1);// This is used to initialize the parameters of the fitting function so that they're wrong by some amount
    std::normal_distribution<> gaus(mu,sigma);// This is for the signal distribution

    // List of signal and background counts to be simulated
    std::vector<int> B{10,100,1000,10000};
    std::vector<int> S{10,100,1000,10000};

    // Lists of 1-dim and 2-dim histograms for S, B and S-vs-B
    size_t nCases = B.size() * S.size();
    std::vector<TH1D*> h_s_Likelihood(nCases);
    std::vector<TH1D*> h_s_NeymanChi2(nCases);
    std::vector<TH1D*> h_s_PearsonChi2(nCases);
    std::vector<TH1D*> h_b_Likelihood(nCases);
    std::vector<TH1D*> h_b_NeymanChi2(nCases);
    std::vector<TH1D*> h_b_PearsonChi2(nCases);
    std::vector<TH2D*> h_sb_Likelihood(nCases);
    std::vector<TH2D*> h_sb_NeymanChi2(nCases);
    std::vector<TH2D*> h_sb_PearsonChi2(nCases);
    int M = 300;

    int c=0;
    std::cout << "s\tb"<< std::endl;
    std::cout << "-------------" << std::endl;
    std::string h_name;
    
    // Loop over S and B
    for( auto& b: B )
	for( auto& s: S )
	    {

		std::cout << s << "\t" << b << std::endl;

		// Estimate the combined uncertainty to define histogram ranges for S and B
		double err = sqrt( s + b );
		double nErr = 7.;
		double mins = std::max( 0., -nErr*err+s );
		double maxs = nErr*err+s;
		double minb = std::max( 0., -nErr*err+b );
		double maxb = nErr*err+b;

		// Create histograms for this specific combination of S and B
		int nBins = 1000;
		h_name = "S Likelihood S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_s_Likelihood[c] = new TH1D( h_name.c_str(), h_name.c_str(), nBins, mins, maxs );
		h_name = "S NeymanChi2 S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_s_NeymanChi2[c] = new TH1D( h_name.c_str(), h_name.c_str(), nBins, mins, maxs );
		h_name = "S PearsonChi2 S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_s_PearsonChi2[c] = new TH1D( h_name.c_str(), h_name.c_str(), nBins, mins, maxs );
		h_s_NeymanChi2[c]->SetLineColor(2);
		h_s_PearsonChi2[c]->SetLineColor(4);
		
		h_name = "B Likelihood S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_b_Likelihood[c] = new TH1D( h_name.c_str(), h_name.c_str(), nBins, minb, maxb );
		h_name = "B NeymanChi2 S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_b_NeymanChi2[c] = new TH1D( h_name.c_str(), h_name.c_str(), nBins, minb, maxb );
		h_name = "B PearsonChi2 S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_b_PearsonChi2[c] = new TH1D( h_name.c_str(), h_name.c_str(), nBins, minb, maxb );
		h_b_NeymanChi2[c]->SetLineColor(2);
		h_b_PearsonChi2[c]->SetLineColor(4);

		h_name = "S vs B Likelihood S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_sb_Likelihood[c] = new TH2D( h_name.c_str(), h_name.c_str(), nBins, mins, maxs, nBins, minb, maxb );
		h_name = "S vs B NeymanChi2 S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_sb_NeymanChi2[c] = new TH2D( h_name.c_str(), h_name.c_str(), nBins, mins, maxs, nBins, minb, maxb );
		h_name = "S vs B PearsonChi2 S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_sb_PearsonChi2[c] = new TH2D( h_name.c_str(), h_name.c_str(), nBins, mins, maxs, nBins, minb, maxb );
		h_sb_Likelihood[c]->GetXaxis()->SetTitle( "S" );
		h_sb_Likelihood[c]->GetYaxis()->SetTitle( "B" );
		h_sb_NeymanChi2[c]->GetXaxis()->SetTitle( "S" );
		h_sb_NeymanChi2[c]->GetYaxis()->SetTitle( "B" );
		h_sb_PearsonChi2[c]->GetXaxis()->SetTitle( "S" );
		h_sb_PearsonChi2[c]->GetYaxis()->SetTitle( "B" );
		
		// Loop over toy-MC experiments
		for( int m=0; m<M; m++ )
		    {
			histo->Reset();
			
			// Fill histogram with background and signal events
			for( int bb=0; bb<b; bb++ )
			    histo->Fill( dis(gen) );
			for( int ss=0; ss<s; ss++ )
			    histo->Fill( gaus(gen) );

			// Initialize the fitting function.
			// We fit only the parameters S and B, while all the other parameters are fixed to their known input values.
			func->SetParameter( 0, mins + (maxs-mins) * smear(gen) );// +- 10% smearing
			func->FixParameter( 1, mu );
			func->FixParameter( 2, sigma );
			func->SetParameter( 3, minb + (maxb-minb) * smear(gen) );//b + b * 0.1 * (maxb-minb) * smear(gen) );
			func->FixParameter( 4, histo->GetBinWidth(1) );
			func->FixParameter( 5, max-min );
			func->SetParLimits( 0, mins, maxs );
			func->SetParLimits( 3, minb, maxb );

			// At this point we can fit multiple times, passing various fit options.
			// Meaning of the fit options:
			// R --> Use the range specified in the TF1 definition
			// N --> Do now plot the fitting function (otherwise it will keep on popping-up new windows
			// Q --> Do not print the fit result on the terminal (otherwise it will be pretty slow)
			// L --> Likelihood fit
			// P --> Pearson Chi2 fit
			// Notice that the default fit mode is the Neyman-Chi2.
			
			// Likelihood fit
			histo->Fit( func, "RNQL" );
			h_s_Likelihood[c]->Fill( func->GetParameter(0) );
			h_b_Likelihood[c]->Fill( func->GetParameter(3) );
			h_sb_Likelihood[c]->Fill( func->GetParameter(0), func->GetParameter(3) );
			// Neyman Chi2 fit
			histo->Fit( func, "RNQ" );
			h_s_NeymanChi2[c]->Fill( func->GetParameter(0) );
			h_b_NeymanChi2[c]->Fill( func->GetParameter(3) );
			h_sb_NeymanChi2[c]->Fill( func->GetParameter(0), func->GetParameter(3) );
			// Pearson Chi2 fit
			histo->Fit( func, "RNQP" );
			h_s_PearsonChi2[c]->Fill( func->GetParameter(0) );
			h_b_PearsonChi2[c]->Fill( func->GetParameter(3) );
			h_sb_PearsonChi2[c]->Fill( func->GetParameter(0), func->GetParameter(3) );
		    }		
		c++;
	    }
    
    // Draw
    gStyle->SetOptStat(0);
    TApplication* app = new TApplication( "app", NULL, 0 );
    // Distribution of fitted S for all combinations of S and B, and all 3 fit methods
    TCanvas* can1 = new TCanvas( "S", "S", 1600, 900 );
    can1->Divide( B.size(), S.size() );
    for( size_t c=0; c<nCases; c++ )
	{
	    can1->cd(c+1);
	    h_s_Likelihood[c]->Draw();
	    h_s_NeymanChi2[c]->Draw("SAME");
	    h_s_PearsonChi2[c]->Draw("SAME");
	}
    // Distribution of fitted B for all combinations of S and B, and all 3 fit methods
    TCanvas* can2 = new TCanvas( "B", "B", 1600, 900 );
    can2->Divide( B.size(), S.size() );
    for( size_t c=0; c<nCases; c++ )
	{
	    can2->cd(c+1);
	    h_b_Likelihood[c]->Draw();
	    h_b_NeymanChi2[c]->Draw("SAME");
	    h_b_PearsonChi2[c]->Draw("SAME");
	}
    
    TCanvas* can3 = new TCanvas( "S vs B Likelihood", "S vs B Likelihood", 1600, 900 );
    can3->Divide( B.size(), S.size() );
    for( size_t c=0; c<nCases; c++ )
	{
	    can3->cd(c+1);
	    h_sb_Likelihood[c]->Draw("colz");
	}
    TCanvas* can4 = new TCanvas( "S vs B Neyman Chi2", "S vs B Neyman Chi2", 1600, 900 );
    can4->Divide( B.size(), S.size() );
    for( size_t c=0; c<nCases; c++ )
	{
	    can4->cd(c+1);
	    h_sb_NeymanChi2[c]->Draw("colz");
	}
    TCanvas* can5 = new TCanvas( "S vs B Pearson Chi2", "S vs B Pearson Chi2", 1600, 900 );
    can5->Divide( B.size(), S.size() );
    for( size_t c=0; c<nCases; c++ )
	{
	    can5->cd(c+1);
	    h_sb_PearsonChi2[c]->Draw("colz");
	}

    can1->SaveAs("LikelihoodVsChi2_1.jpg");
    can2->SaveAs("LikelihoodVsChi2_2.jpg");
    can3->SaveAs("LikelihoodVsChi2_3.jpg");
    can4->SaveAs("LikelihoodVsChi2_4.jpg");
    can5->SaveAs("LikelihoodVsChi2_5.jpg");

    can1->SaveAs("LikelihoodVsChi2_1.pdf");
    can2->SaveAs("LikelihoodVsChi2_2.pdf");
    can3->SaveAs("LikelihoodVsChi2_3.pdf");
    can4->SaveAs("LikelihoodVsChi2_4.pdf");
    can5->SaveAs("LikelihoodVsChi2_5.pdf");
    
    app->Run(kTRUE);
    
    return 0;
}
