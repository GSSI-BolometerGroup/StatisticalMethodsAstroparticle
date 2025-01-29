#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>

#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"

double NoiseSpike( double A,
		   double t,
		   double t0,
		   double f,
		   double tau )
{
    return A * exp( - ( t - t0 ) / tau ) * sin( f * ( t - t0 ) );
}

double Pulse( double E,
	      double t,
	      double baselineT,
	      double risetime,
	      double decaytime,
	      double rise_dt,
	      double MSEfraction,
	      bool   noise,
	      double A,
	      double t0,
	      double f,
	      double tau )
{
    double value;
    if( rise_dt == 0. )
	value = 0.5 * E * ( 1. + erf( ( t - baselineT ) / risetime ) ) * exp( - ( t - baselineT ) / decaytime );
    else
	value = 0.5 * E * MSEfraction * ( 1. + erf( ( t - baselineT ) / risetime ) ) * exp( - ( t - baselineT ) / decaytime )
	    + 0.5 * E * ( 1. - MSEfraction ) * ( 1. + erf( ( t - baselineT + rise_dt ) / risetime ) ) * exp( - ( t - baselineT + rise_dt ) / decaytime );

    if( noise && t > t0 )
	value += NoiseSpike( A, t, t0, f, tau );
    
    return value;
}

struct Peak{
    std::string Name;
    double Energy;
    double SingleSideFraction;
    int    NEvents;
};

int main()
{

    double ns = 1.;
    double us = 1.e3;
    double dt = 10. * ns;// Sampling time
    double Dt = 160. * us;// Time window
    double baselineT    = 40. * us;
    double baseline     = 0.;// Average baseline
    double baselineRMS  = 10.;
    double risetime     = 0.2 * us;
    double risetimeRMS  = 0.02 * us;
    double decaytime    = 150. * us;
    double decaytimeRMS = 5. * us;
    double E_RMS = 1.5;
    double risetime_DtMin = 0. * us;
    double risetime_DtMax = 0.7 * us;
    double noiseFraction = 0.1;
    double noiseA = 1000.;
    double noiseARMS = 100.;
    double noiseF = 100. / ns;
    double noiseFRMS = 10. / ns;
    double noiseTau = 200. * ns;
    double noiseTauRMS = 20. * ns;
    
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(12345); // Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<> baseline_rdm( baseline, baselineRMS );
    std::normal_distribution<> E_rdm( 0., E_RMS );
    std::normal_distribution<> decaytime_rdm( decaytime, decaytimeRMS );
    std::uniform_real_distribution<> MSEFraction_rdm(0.,1.);
    std::normal_distribution<> SSEFraction_rdm(0.,0.02);
    std::uniform_real_distribution<> PSD_rdm(0.,1.);
    std::normal_distribution<> risetime_rdm( risetime, risetimeRMS );
    std::uniform_real_distribution<> risetime_dt_rdm( risetime_DtMin, risetime_DtMax );
    std::uniform_real_distribution<> noise_rdm(0.,1.);
    std::uniform_real_distribution<> noiseT0_rdm(0.,Dt);
    std::normal_distribution<> noiseA_rdm( noiseA, noiseARMS );
    std::normal_distribution<> noiseF_rdm( noiseF, noiseFRMS );
    std::normal_distribution<> noiseTau_rdm( noiseTau, noiseTauRMS );

    std::vector<Peak> peaklistCal;
    peaklistCal.push_back( { std::string("DEP"),     1592., 1.,   3000  } );
    peaklistCal.push_back( { std::string("1620keV"), 1620., 0.2,  3000  } );
    peaklistCal.push_back( { std::string("SEP"),     2103., 0.05, 3000  } );
    peaklistCal.push_back( { std::string("FEP"),     2615., 0.1,  3000 } );
    
    int nBins = Dt / dt;
    //std::vector<double>* wf = new std::vector<double>(nBins);
    //TH1D* h_wf = new TH1D( "wf", "wf", nBins, 0., Dt );
    //gStyle->SetOptStat(0);
    //gStyle->SetOptTitle(0);
    //TApplication* app = new TApplication( "app", NULL, 0 );
    //TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
    //can->cd();

    //TFile* outCal = new TFile( "Calibration.root", "RECREATE" );
    //outCal->cd();
    //TTree* treeCal = new TTree( "Events", "Events" );
    double E;
    int    peakType;
    double decay;
    double rise;
    double MSEFraction;
    //treeCal->Branch( "E", &E, "E/D" );
    //treeCal->Branch( "Label", &peakType, "Label/I" );
    //treeCal->Branch( "DecayTime", &decay, "DecayTime/D" );
    //treeCal->Branch( "RiseTime", &rise, "RiseTime/D" );
    //treeCal->Branch( "MSEFraction", &MSEFraction, "MSEFraction/D" );

    //treeCal->Branch( "WF", wf );

    std::ofstream outCal("Calibration.txt");    
    outCal << "Energy\tLabel\tDecayTime\tRiseTime\tMSEFraction" << std::endl;   
    for( auto& peak : peaklistCal )
	for( int k=0; k<peak.NEvents; k++ )
	    {
		if( peak.Name.compare( "DEP" ) == 0 )
		    peakType = 0;
		else if( peak.Name.compare( "1620keV" ) == 0 )
		    peakType = 1;
		else if( peak.Name.compare( "SEP" ) == 0 )
		    peakType = 2;
		else if( peak.Name.compare( "FEP" ) == 0 )
		    peakType = 3;
		
		E = peak.Energy + E_rdm(gen);
		decay = decaytime_rdm(gen);
		double psd = PSD_rdm(gen);
		double rise_dt = 0.;
		MSEFraction = SSEFraction_rdm(gen);
		
		if( psd > peak.SingleSideFraction )
		    {
			rise_dt = risetime_dt_rdm(gen);
			MSEFraction += MSEFraction_rdm(gen);
		    }
		rise = risetime_rdm(gen) + rise_dt;
    

		//for( int i=0; i<nBins; i++ )
		//  {
		//	wf->at(i) = baseline_rdm(gen);
		//	wf->at(i) += Pulse( E, dt*i, baselineT, risetime, decay, rise_dt, MSEFraction, false, 0, 0, 0, 0 );
		//	h_wf->SetBinContent( i+1, wf->at(i) );
		//  }
		outCal << E << "\t"
		       << peakType << "\t"
		       << decay << "\t"
		       << rise << "\t"
		       << MSEFraction << std::endl;
		//treeCal->Fill();
		//h_wf->Draw();
		//app->Run(kTRUE);
	    }
    //outCal->Write();
    //outCal->Close();


    double minE = 2000.;
    double maxE = 2080.;
    std::uniform_real_distribution<> bkg_rdm(minE,maxE);
    
    std::vector<Peak> peaklistBkg;
    peaklistBkg.push_back( { std::string("Sgn"),  2039., 1., 3 } );
    peaklistBkg.push_back( { std::string("Bkg"),  minE, 0.2, 20 } );

    //TFile* outBkg = new TFile( "Physics.root", "RECREATE" );
    //outBkg->cd();
    //TTree* treeBkg = new TTree( "Events", "Events" );
    //treeBkg->Branch( "E", &E, "E/D" );
    //treeBkg->Branch( "DecayTime", &decay, "DecayTime/D" );
    //treeBkg->Branch( "RiseTime", &rise, "RiseTime/D" );
    //treeBkg->Branch( "MSEFraction", &MSEFraction, "MSEFraction/D" );

    std::ofstream outBkg("Physics.txt");
    outBkg << "Energy\tDecayTime\tRiseTime\tMSEFraction" << std::endl;
    for( auto& peak : peaklistBkg )
	for( int k=0; k<peak.NEvents; k++ )
	    {
		if( peak.Name.compare( "Sgn" ) == 0 )
		    E = peak.Energy + E_rdm(gen);
		else
		    E = bkg_rdm(gen);

		decay = decaytime_rdm(gen);
		double psd = PSD_rdm(gen);
		double rise_dt = 0.;
		MSEFraction = SSEFraction_rdm(gen);
		
		if( psd > peak.SingleSideFraction )
		    {
			rise_dt = risetime_dt_rdm(gen);
			MSEFraction += MSEFraction_rdm(gen);
		    }
		rise = risetime_rdm(gen) + rise_dt;

		outBkg << E << "\t"
		       << decay << "\t"
		       << rise << "\t"
		       << MSEFraction << std::endl;
		//treeBkg->Fill();
		std::cout << peak.Name << "\t" << E << std::endl;
	    }
    //outBkg->Write();
    //outBkg->Close();

    return 0;
}
