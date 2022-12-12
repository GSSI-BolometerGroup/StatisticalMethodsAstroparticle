// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__RADIOACTIVEDECAYFIT__H
#define __BAT__RADIOACTIVEDECAYFIT__H

#include <BAT/BCModel.h>

#include <string>
#include <vector>

// This is a RadioactiveDecayFit header file.
// Model source code is located in file RadioactiveDecayFit/RadioactiveDecayFit.cxx

// ---------------------------------------------------------
class RadioactiveDecayFit : public BCModel
{
public:
    enum FitMethod{ kBinomial, kPoisson };
    
private:
    int    fN;// Original number of 210Po isotopes
    double fHalflife;
    double fDecayRate;
    double fDeltaT;// Measurement time
    std::vector<double> fT;// Time of each event
    size_t fK;

    double fP;// Binomial probability of events in fDeltaT
    double fLambda;// Poisson expectation for number of events in fDeltaT (=n*fP)

    FitMethod fMethod;
    
    double LogBinomial( double k,
			double n,
			double p );
    double LogFactorial( double n );
    double LogPoisson( double n,
		       double lambda );
    
public:

    // Constructor
    RadioactiveDecayFit(const std::string& name);

    // Destructor
    ~RadioactiveDecayFit();

    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double>& pars);

    // Overload LogAprioriProbability if not using built-in 1D priors
    // double LogAPrioriProbability(const std::vector<double> & pars);

    // Overload CalculateObservables if using observables
    // void CalculateObservables(const std::vector<double> & pars);

    void SetFitMethod( FitMethod method );
    
};
// ---------------------------------------------------------

#endif
