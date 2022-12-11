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
    TH1D* histo  = new TH1D("histo","histo",nBins,min,max);
    TF1* func = new TF1( "func", GaussianPlusFlat, min, max, 6 );
    func->SetParName( 0, "Integral" );
    func->SetParName( 1, "Mean" );
    func->SetParName( 2, "Sigma" );
    func->SetParName( 3, "Background" );
    func->SetParName( 4, "Bin width" );
    func->SetParName( 5, "Fit range" );
    func->SetNpx(10000);
    
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(min,max);
    std::uniform_real_distribution<> smear(0,1);
    std::normal_distribution<> gaus(mu,sigma);

    std::vector<int> B{10,100,1000,10000};
    std::vector<int> S{10,100,1000,10000};

    size_t nCases = B.size() * S.size();
    std::vector<TH1D*> h_s_Poisson(nCases);
    std::vector<TH1D*> h_s_NeymanChi2(nCases);
    std::vector<TH1D*> h_s_PearsonChi2(nCases);
    std::vector<TH1D*> h_b_Poisson(nCases);
    std::vector<TH1D*> h_b_NeymanChi2(nCases);
    std::vector<TH1D*> h_b_PearsonChi2(nCases);
    std::vector<TH2D*> h_sb_Poisson(nCases);
    std::vector<TH2D*> h_sb_NeymanChi2(nCases);
    std::vector<TH2D*> h_sb_PearsonChi2(nCases);
    int M = 300;
    //int B = 1000;// Number of background events (uniformly distributed)
    //int S = 1000;// Number of signal events (Gaussian distributed around mu)

    int c=0;
    std::string h_name;
    for( auto& b: B )
	for( auto& s: S )
	    {

		std::cout << s << "\t" << b << std::endl;

		double err = sqrt( s + b );
		double mins = std::max( 0., -5.*err+s );
		double maxs = 5.*err+s;
		double minb = std::max( 0., -5.*err+b );
		double maxb = 5.*err+b;
		
		int nBins = 1000;
		h_name = "S Poisson S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_s_Poisson[c] = new TH1D( h_name.c_str(), h_name.c_str(), nBins, mins, maxs );
		h_name = "S NeymanChi2 S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_s_NeymanChi2[c] = new TH1D( h_name.c_str(), h_name.c_str(), nBins, mins, maxs );
		h_name = "S PearsonChi2 S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_s_PearsonChi2[c] = new TH1D( h_name.c_str(), h_name.c_str(), nBins, mins, maxs );
		h_s_NeymanChi2[c]->SetLineColor(2);
		h_s_PearsonChi2[c]->SetLineColor(4);
		
		h_name = "B Poisson S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_b_Poisson[c] = new TH1D( h_name.c_str(), h_name.c_str(), nBins, minb, maxb );
		h_name = "B NeymanChi2 S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_b_NeymanChi2[c] = new TH1D( h_name.c_str(), h_name.c_str(), nBins, minb, maxb );
		h_name = "B PearsonChi2 S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_b_PearsonChi2[c] = new TH1D( h_name.c_str(), h_name.c_str(), nBins, minb, maxb );
		h_b_NeymanChi2[c]->SetLineColor(2);
		h_b_PearsonChi2[c]->SetLineColor(4);

		h_name = "S vs B Poisson S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_sb_Poisson[c] = new TH2D( h_name.c_str(), h_name.c_str(), nBins, mins, maxs, nBins, minb, maxb );
		h_name = "S vs B NeymanChi2 S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_sb_NeymanChi2[c] = new TH2D( h_name.c_str(), h_name.c_str(), nBins, mins, maxs, nBins, minb, maxb );
		h_name = "S vs B PearsonChi2 S=" + std::to_string(s) + " B=" + std::to_string(b);
		h_sb_PearsonChi2[c] = new TH2D( h_name.c_str(), h_name.c_str(), nBins, mins, maxs, nBins, minb, maxb );
		h_sb_Poisson[c]->GetXaxis()->SetTitle( "S" );
		h_sb_Poisson[c]->GetYaxis()->SetTitle( "B" );
		h_sb_NeymanChi2[c]->GetXaxis()->SetTitle( "S" );
		h_sb_NeymanChi2[c]->GetYaxis()->SetTitle( "B" );
		h_sb_PearsonChi2[c]->GetXaxis()->SetTitle( "S" );
		h_sb_PearsonChi2[c]->GetYaxis()->SetTitle( "B" );
		
		
		for( int m=0; m<M; m++ )
		    {
			histo->Reset();
			
			// Fill histogram with background and signal events
			for( int bb=0; bb<b; bb++ )
			    histo->Fill( dis(gen) );
			for( int ss=0; ss<s; ss++ )
			    histo->Fill( gaus(gen) );

			// Define fitting function

			func->SetParameter( 0, mins + (maxs-mins) * smear(gen) );// +- 10% smearing
			//func->SetParameter( 1, mu + mu * 0.0005 * smear(gen) );// +- 2 keV smearing
			func->FixParameter( 1, mu );
			//func->SetParameter( 2, sigma + sigma * 0.1 * smear(gen));
			func->FixParameter( 2, sigma );
			func->SetParameter( 3, minb + (maxb-minb) * smear(gen) );//b + b * 0.1 * (maxb-minb) * smear(gen) );
			func->FixParameter( 4, histo->GetBinWidth(1) );
			func->FixParameter( 5, max-min );
			func->SetParLimits( 0, mins, maxs );
			//func->SetParLimits( 1, mu-3*sigma, mu+3*sigma );
			//func->SetParLimits( 2, 0.5*sigma, 1.5*sigma );
			func->SetParLimits( 3, minb, maxb );

			histo->Fit( func, "RNQL" );
			h_s_Poisson[c]->Fill( func->GetParameter(0) );
			h_b_Poisson[c]->Fill( func->GetParameter(3) );
			h_sb_Poisson[c]->Fill( func->GetParameter(0), func->GetParameter(3) );
			histo->Fit( func, "RNQ" );
			h_s_NeymanChi2[c]->Fill( func->GetParameter(0) );
			h_b_NeymanChi2[c]->Fill( func->GetParameter(3) );
			h_sb_NeymanChi2[c]->Fill( func->GetParameter(0), func->GetParameter(3) );
			histo->Fit( func, "RNQP" );
			h_s_PearsonChi2[c]->Fill( func->GetParameter(0) );
			h_b_PearsonChi2[c]->Fill( func->GetParameter(3) );
			h_sb_PearsonChi2[c]->Fill( func->GetParameter(0), func->GetParameter(3) );
		    }		
		c++;
	    }
    
    // Draw
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can1 = new TCanvas( "can1", "can1", 1600, 900 );
    can1->Divide( B.size(), S.size() );
    for( size_t c=0; c<nCases; c++ )
	{
	    can1->cd(c+1);
	    h_s_Poisson[c]->Draw();
	    h_s_NeymanChi2[c]->Draw("SAME");
	    h_s_PearsonChi2[c]->Draw("SAME");
	}
    TCanvas* can2 = new TCanvas( "can2", "can2", 1600, 900 );
    can2->Divide( B.size(), S.size() );
    for( size_t c=0; c<nCases; c++ )
	{
	    can2->cd(c+1);
	    h_b_Poisson[c]->Draw();
	    h_b_NeymanChi2[c]->Draw("SAME");
	    h_b_PearsonChi2[c]->Draw("SAME");
	}
    TCanvas* can3 = new TCanvas( "can3", "can3", 1600, 900 );
    can3->Divide( B.size(), S.size() );
    for( size_t c=0; c<nCases; c++ )
	{
	    can3->cd(c+1);
	    h_sb_Poisson[c]->Draw("colz");
	}
    TCanvas* can4 = new TCanvas( "can4", "can4", 1600, 900 );
    can4->Divide( B.size(), S.size() );
    for( size_t c=0; c<nCases; c++ )
	{
	    can4->cd(c+1);
	    h_sb_NeymanChi2[c]->Draw("colz");
	}
    TCanvas* can5 = new TCanvas( "can5", "can5", 1600, 900 );
    can5->Divide( B.size(), S.size() );
    for( size_t c=0; c<nCases; c++ )
	{
	    can5->cd(c+1);
	    h_sb_PearsonChi2[c]->Draw("colz");
	}
    //histo->Draw();
    //histo->Fit(func,"R");


    //histo->Fit(func,"RLL");
    //histo->Fit(func,"R");
    app->Run(kTRUE);
    
    return 0;
}
