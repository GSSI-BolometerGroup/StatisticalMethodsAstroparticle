// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__BINNEDVSUNBINNED__H
#define __BAT__BINNEDVSUNBINNED__H

#include <BAT/BCModel.h>

#include <string>
#include <vector>
#include <random>

#include "TH1D.h"

// This is a BinnedVsUnbinned header file.
// Model source code is located in file BinnedVsUnbinned/BinnedVsUnbinned.cxx

// ---------------------------------------------------------
class BinnedVsUnbinned : public BCModel
{
public:
    enum Type{ kBinned, kUnbinned };
    
private:
    Type fType;
    
    double fMinE;
    double fMaxE;
    double fDeltaE;
    double fDE;
    int    fNBins;
    std::vector<double> fE;
    std::vector<double> fK;
    std::vector<double> fGaussianPDF;
    std::vector<double> fGaussianBinCDF;
    double fBkgPDF;
    double fBkgBinCDF;
    int    fS;
    int    fB;
    int    fN;
    double fMu;
    double fSigma;
    TH1D*  fDataHisto;

    std::random_device fRD;
    std::mt19937 fGen;
    std::uniform_real_distribution<> fBkg;
    std::normal_distribution<>       fSgn;
    
public:

    // Constructor
    BinnedVsUnbinned(const std::string& name);

    // Destructor
    ~BinnedVsUnbinned();

    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double>& pars);

    // Overload LogAprioriProbability if not using built-in 1D priors
    // double LogAPrioriProbability(const std::vector<double> & pars);

    // Overload CalculateObservables if using observables
    // void CalculateObservables(const std::vector<double> & pars);

    void ResetData();
    void SetData( int s, int b );
    void SetType( Type type );

    int GetN(){ return fN; };
    int GetNBins(){ return fNBins; };
    
};
// ---------------------------------------------------------

#endif
