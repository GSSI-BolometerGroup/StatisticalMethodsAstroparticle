// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__EFFICIENCY__H
#define __BAT__EFFICIENCY__H

#include <BAT/BCModel.h>

#include <string>
#include <vector>
#include <map>

#include "TF1.h"

// This is a Efficiency header file.
// Model source code is located in file Efficiency/Efficiency.cxx

// ---------------------------------------------------------
class Efficiency : public BCModel
{
public:
    enum FitMethod{ kBinomial, kPoisson, kChiSquare };
    
private:
    int fN;
    double fMinE;
    double fMaxE;
    double fThreshold;
    double fSigma;
    double fEfficiency;// Asympthotic value
    TF1 fEfficiencyCurve;
    std::map<int,double> fK;// map of (energy,k)
    FitMethod fMethod;
    
    double EvaluateEfficiency( double E,
			       double asympEff,
			       double threshold,
			       double sigma );
    double LogBinomial( double k,
			double n,
			double p );
    double LogFactorial( double n );
    double LogPoisson( double n,
		       double lambda );
    double LogChi2( double n, double mu );
    
public:

    // Constructor
    Efficiency( const std::string& name,
		const int    n=100,
		const double eff=0.9 );

    // Destructor
    ~Efficiency();

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
