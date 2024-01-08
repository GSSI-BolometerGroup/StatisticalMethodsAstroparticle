// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__EVIDENCE__H
#define __BAT__EVIDENCE__H

#include <BAT/BCModel.h>

#include <string>
#include <vector>
#include <random>

// This is a Evidence header file.
// Model source code is located in file Evidence/Evidence.cxx

// ---------------------------------------------------------
class Evidence : public BCModel
{
public:
    enum Model{ kBkgOnly, kSgnPlusBkg };
    
private:
    int fS;
    int fB;
    int fN;
    Model fModel;
    double fMinE;
    double fMaxE;
    double fDeltaE;

    std::vector<double> fE;
    std::vector<double> fGaussianPDF;
    double fBkgPDF;
    double fMu;
    double fSigma;
    double fLogFactN;
    
    double fIntegralBkgOnly;
    double fIntegralSgnPlusBkg;
    int    fNHitMiss;
    int    fNAcceptedHitMiss;

    std::random_device fRD;
    std::mt19937 fGen;
    std::uniform_real_distribution<> fBkg;
    std::normal_distribution<>       fSgn;
    std::vector< std::uniform_real_distribution<> > fSpaceSampler;
    std::uniform_real_distribution<> fR;

    double LogFactorial( double n );
    
public:

    // Constructor
    Evidence( const std::string& name,
	      int s,
	      int b );

    // Destructor
    ~Evidence();

    void SetModel( Model model );
    
    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double>& pars);

    double LogLikelihoodBkgOnly( const std::vector<double>& pars );
    double LogLikelihoodSgnPlusBkg( const std::vector<double>& pars );

    // Overload LogAprioriProbability if not using built-in 1D priors
    double LogAPrioriProbability(const std::vector<double> & pars);

    // Overload CalculateObservables if using observables
    // void CalculateObservables(const std::vector<double> & pars);

    void ComputeIntegral();

    void ComputeBayesFactor();
};
// ---------------------------------------------------------

#endif
