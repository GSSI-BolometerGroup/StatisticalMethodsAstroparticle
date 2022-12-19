#ifndef __BAT__LIKELIHOODRATIO__H
#define __BAT__LIKELIHOODRATIO__H

#include <BAT/BCModel.h>

#include <string>
#include <vector>
#include <random>

#include "TH1D.h"
#include "TF1.h"

class LikelihoodRatio : public BCModel
{

private:
    double fMinE;
    double fMaxE;
    double fDeltaE;// Energy range
    double fDE;// Bin width
    int fNBins;// Number of bins
    std::vector<double> fK;// Vector of bin entries. We use a vector instead of a TH1D to minimize the CPU time
    std::vector<double> fGaussianPDF;// Gaussian PDF at bin k
    double fS;// Injected number of signal events
    double fB;// Injected number of background events
    double fMu;// Mean value of the signal
    double fSigma;// Sigma of the signal
    TH1D* fDataHisto;// This is used just for drawing

    // Random number machinery
    std::random_device fRD;
    std::mt19937 fGen;
    std::uniform_real_distribution<> fBkg;;
    std::normal_distribution<> fSgn;

    // Set all entries of fK to zero
    void ResetData();
    
public:

    LikelihoodRatio(const std::string& name);

    ~LikelihoodRatio();

    double LogLikelihood(const std::vector<double>& pars);

    // Create new set of toy data
    void SetData();

    // Getters
    TH1D* GetDataHisto();
    TF1* GetFittingFunction();
    int GetNDF();
};

#endif
