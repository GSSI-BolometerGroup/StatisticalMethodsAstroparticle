#ifndef __BAT__RADIOACTIVEDECAYFIT__H
#define __BAT__RADIOACTIVEDECAYFIT__H

#include <BAT/BCModel.h>

#include <string>
#include <vector>

#include "TH1D.h"
#include "TF1.h"

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

    TH1D* fDataHisto;
    
    double LogBinomial( double k,
			double n,
			double p );
    double LogFactorial( double n );
    double LogPoisson( double n,
		       double lambda );
    
public:

    RadioactiveDecayFit(const std::string& name);

    ~RadioactiveDecayFit();

    double LogLikelihood(const std::vector<double>& pars);

    void SetFitMethod( FitMethod method );

    TH1D* GetData();
    TF1* GetBestFit();
};

#endif
