// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "BinnedVsUnbinned.h"

#include "BAT/BCConstantPrior.h"
// #include <BAT/BCMath.h>

// ---------------------------------------------------------
BinnedVsUnbinned::BinnedVsUnbinned(const std::string& name)
    : BCModel(name)
{
    fType = Type::kUnbinned;
    
    fMinE   = 2000.;// keV
    fMaxE   = 2080.;// keV
    fDeltaE = fMaxE - fMinE;
    fDE     = fDeltaE / 80;
    fNBins  = fDeltaE / fDE;
    fK      = std::vector<double>(fNBins);
    for( auto&k : fK )
	k = 0.;

    fMu    = 2039.;// keV
    fSigma = 1.5;// keV

    fGaussianBinCDF = std::vector<double>( fNBins );
    for( int b=0; b<fNBins; b++ )
	{
	    double Elow = fMinE + fDE * b;
	    double Eup  = fMinE + fDE * ( b + 1 );
	    fGaussianBinCDF[b] = 0.5 * ( std::erf( ( Eup  - fMu ) / std::sqrt(2.) / fSigma ) -
					 std::erf( ( Elow - fMu ) / std::sqrt(2.) / fSigma ) );
	}
    
    fGen = std::mt19937(fRD());
    fBkg = std::uniform_real_distribution<>(fMinE,fMaxE);
    fSgn = std::normal_distribution<>(fMu,fSigma);

    fBkgPDF    = 1. / fDeltaE;
    fBkgBinCDF = fDE / fDeltaE;

    double errSB = sqrt( fS + fB );
    double nErr = 7.;
    double minS = std::max( 0., fS - nErr * errSB );
    double maxS = fS + 7. * errSB;
    AddParameter( "S", minS, maxS, "S", "[counts]" );
    GetParameters().Back().SetPrior( new BCConstantPrior() );
    GetParameters().Back().SetNbins(300);

    double minB = std::max( 0., fB - nErr * errSB );
    double maxB = fB + 7. * errSB;
    AddParameter( "B", minB, maxB, "B", "[counts]" );
    GetParameters().Back().SetPrior( new BCConstantPrior() );
    GetParameters().Back().SetNbins(300);
    
}

void BinnedVsUnbinned::ResetData()
{
    for( auto& k : fK )
	k = 0;

    for( auto& fs : fGaussianPDF )
	fs = 0;
    
    return;
}

void BinnedVsUnbinned::SetData( int s, int b )
{
    ResetData();

    fS = s;
    fB = b;
    fN = fS + fB;
    fE     = std::vector<double>(fS+fB);

    for( int s=0; s<fS; s++ )
	{
	    double e = fSgn(fGen);
	    int    i = ( e - fMinE ) / fDE;
	    fE[s]    = e;// set vector of event energies
	    fK[i]    += 1.;// fill histogram
	}

    for( int b=0; b<fB; b++ )
	{
	    double e = fBkg(fGen);
	    int    i = ( e - fMinE ) / fDE;
	    fE[fS+b] = e;
	    fK[i]    += 1.;
	}

    fGaussianPDF = std::vector<double>(fS+fB);
    for( int i=0; i<fS+fB; i++ )
	{
	    double e = fE[i];
	    fGaussianPDF[i] = exp( - pow( ( e - fMu ) / fSigma, 2. ) / 2. ) / sqrt(2.*M_PI) / fSigma;
	}
    
    return;
}

// ---------------------------------------------------------
BinnedVsUnbinned::~BinnedVsUnbinned()
{
    ;
}

// ---------------------------------------------------------
double BinnedVsUnbinned::LogLikelihood(const std::vector<double>& pars)
{

    double logL = 0.;

    double s = pars[0];
    double b = pars[1];
    
    if( fType == Type::kUnbinned )
	{
	    double lambda = s + b;

	    // Extended likelihood:
	    // logL = -lambda + n*log(lambda) - log(n!) + Sum_i { log( fs*s/lambda + fb*b/lambda ) }
	    //      = -lambda + n*log(lambda) - log(n!) + Sum_i { log( fs*s + fb*b ) - log(lambda) }
	    //      = -lambda + n*log(lambda) - log(n!) + Sum_i { log( fs*s + fb*b ) } - n*log*lambda)
	    //      = - lambda + Sum_i { log( fs*s + fb*g ) } --> log(n!) is constant so we can neglect it
	    
	    logL = -lambda;
	    for( auto& fs: fGaussianPDF )// Loop over events
		logL += log( s * fs + b * fBkgPDF );
	}
    else if( fType == Type::kBinned )
	{

	    // Binned likelihood:
	    // logL = Sum_b { -lambda_b + k_b * log(lambda_b) + log(k_b!) }
	    //      = Sum_b { -lambda_b + k_b * log(lambda_b) } --> log(k_b!) is constant, so we can neglect it
	    
	    for( int i=0; i<fNBins; i++ )// Loop over bins
		{
		    double lambda = s * fGaussianBinCDF[i] + b * fBkgBinCDF;
		    logL += -lambda + fK[i] * log(lambda);
		}
	}
    
    return logL;
}

// ---------------------------------------------------------
// double BinnedVsUnbinned::LogAPrioriProbability(const std::vector<double>& pars)
// {
//     // return the log of the prior probability p(pars)
//     // If you use built-in priors, leave this function commented out.
// }

// ---------------------------------------------------------
// void BinnedVsUnbinned::CalculateObservables(const std::vector<double>& pars)
// {
//     // Calculate and store obvserables. For example:
//     GetObservable(0) = pow(pars[0], 2);
// }

void BinnedVsUnbinned::SetType( Type type )
{
    fType = type;
    
    return;
}
