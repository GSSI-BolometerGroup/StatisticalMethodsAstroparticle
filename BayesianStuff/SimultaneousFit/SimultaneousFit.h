#ifndef __BAT__SIMULTANEOUSFIT__H
#define __BAT__SIMULTANEOUSFIT__H

#include <BAT/BCModel.h>

#include <string>
#include <vector>
#include <map>

#include "TF1.h"

class SimultaneousFit : public BCModel
{
private:
    int fNInjected;// Number of injected events for each energy
    int fS;// Number of generated signal events
    int fB;// Number of generated background events
    double fMinEEff;// Minimum energy of pulser-injected events
    double fMaxEEff;// Maximum energy of pulser-injected events
    double fThreshold;// Trigger threshold (mean of the error function)
    double fSigmaEff;// Sigma of the efficiency curve (sigma of the error function)
    double fEfficiency;// Asympthotic value
    TF1 fEfficiencyCurve;
    std::map<int,double> fK;// map of (energy,k)

    double fMinE;// Minimum of energy spectrum
    double fMaxE;// Maximum of energy spectrum
    double fDeltaE;// Width of fit range
    double fMu;// Position of signal peak
    double fSigma;// Resolution of signal peak


    std::vector<double> fE;// List of event energies
    std::vector<double> fGaussianPDF;// Value of the signal PDF for each event
    double fBkgPDF;// Value of the background PDF for all events (the bkg is flat, so the PDF is always the same)
    
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
    
public:

    // Constructor
    SimultaneousFit( const std::string& name,
		     int nInjected,
		     double eff,
		     int s,
		     int b );

    // Destructor
    ~SimultaneousFit();

    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double>& pars);

    // Overload LogAprioriProbability if not using built-in 1D priors
    // double LogAPrioriProbability(const std::vector<double> & pars);

    // Overload CalculateObservables if using observables
    // void CalculateObservables(const std::vector<double> & pars);

};
// ---------------------------------------------------------

#endif
