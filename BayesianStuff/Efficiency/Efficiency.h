#ifndef __BAT__EFFICIENCY__H
#define __BAT__EFFICIENCY__H

#include <BAT/BCModel.h>

#include <string>
#include <vector>
#include <map>

#include "TF1.h"

class Efficiency : public BCModel
{
public:
    // With this enumerator, we switch between the 3 types of fit
    enum FitMethod{ kBinomial, kPoisson, kChiSquare };
    
private:
    int fN;// Number of injected events for each energy
    double fMinE;// Minimum energy
    double fMaxE;// Maximum energy
    double fThreshold;// Trigger threshold (mean of the error function)
    double fSigma;// Sigma of the efficiency curve (sigma of the error function)
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

    void SetFitMethod( FitMethod method );
};
// ---------------------------------------------------------

#endif
