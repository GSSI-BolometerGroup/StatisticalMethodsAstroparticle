#include <iostream>
#include <vector>
#include <random>
#include <string>

#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TPaveStats.h"

int main()
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    int N=12;
    int K=100000000;
    int nBins = 1000;
    double min = -10.;
    double max = 10.;
    std::vector<TH1D*> histo(N);
    std::vector<TH1D*> histo_range(N);
    for (int n=1; n<=N; n++)
	{

	    std::string name( "N = " + std::to_string(n) );
	    histo[n-1] = new TH1D( name.c_str(), name.c_str(), nBins, min, max );

	    
	    for( int i=0; i<K; i++ )
		{

		    double g = 0.;
		    for( int k=0; k<n; k++ )
			g += dis(gen);
		    g -= 0.5 * n;
		    g *= sqrt( 12. / n );
		    histo[n-1]->Fill(g);
		}
	}

    for (int n=1; n<=N; n++)
	{
	    std::string name( "Range for N=" + std::to_string(n) );
	    histo_range[n-1] = new TH1D( name.c_str(), name.c_str(), nBins, min, max );
	    histo_range[n-1]->SetFillColor(18);
	    histo_range[n-1]->SetLineWidth(0);
	    
	    double left  = - sqrt( 3. * n );
	    double right = + sqrt( 3. * n );
	    double max = histo[n-1]->GetMaximum();
	    for( int b=histo_range[n-1]->FindBin(left); b<=histo_range[n-1]->FindBin(right); b++ )
		histo_range[n-1]->SetBinContent( b, 2.*max );
	}
		 

    gStyle->SetLabelFont(132,"XY");
    gStyle->SetTitleFont(132,"XY");// Set title font for x-y axes
    gStyle->SetTitleFont(132,"");// Set title font for all TPAds
    gStyle->SetTitleFontSize(0.1);
    gStyle->SetFrameLineWidth(0);

    std::vector<TPaveStats*> stat(N);
    TApplication* app = new TApplication( "app", NULL, 0 );
    TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    can->Divide(3,4);
    for( int i=0; i<N; i++ )
	{
	    can->cd(i+1);
	    can->GetPad(i+1)->SetTopMargin(0.02);
	    can->GetPad(i+1)->SetRightMargin(0.02);
	    can->GetPad(i+1)->SetLeftMargin(0.07);
	    can->GetPad(i+1)->SetBottomMargin(0.18);
	    can->GetPad(i+1)->SetLogy();
	    histo[i]->GetYaxis()->SetRangeUser(1.e-1,1.1*histo[i]->GetMaximum());
	    histo[i]->GetXaxis()->SetTitle("g");
	    histo[i]->GetXaxis()->SetLabelSize(0.1);
	    histo[i]->GetYaxis()->SetLabelSize(0.1);
	    histo[i]->GetXaxis()->SetTitleSize(0.1);
	    histo[i]->GetXaxis()->SetTitleOffset(0.9);
	    histo[i]->GetXaxis()->SetLabelFont(132);
	    histo[i]->GetYaxis()->SetLabelFont(132);
	    histo[i]->GetXaxis()->SetTitleFont(132);
	    histo[i]->Draw("AXIS");
	    gPad->Update();
	    stat[i] = (TPaveStats*)histo[i]->FindObject("stats");
	    stat[i]->SetLineWidth(0);
	    stat[i]->SetFillStyle(0);
	    stat[i]->SetTextFont(132);
	    stat[i]->SetTextSize(0.1);
	    stat[i]->SetOptStat(1101);
	    stat[i]->SetX1NDC(0.7);
	    stat[i]->SetX2NDC(0.99);
	    stat[i]->SetY1NDC(0.7);
	    stat[i]->SetY2NDC(0.99);
	    histo_range[i]->Draw("SAME");
	    histo[i]->Draw("SAME");
	}
    can->SaveAs("CentralLimit.pdf");
    can->SaveAs("CentralLimitLog.jpg");
    app->Run(kTRUE);
    
    return 0;
}
