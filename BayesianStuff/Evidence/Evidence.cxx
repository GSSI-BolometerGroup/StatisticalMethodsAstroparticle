// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "Evidence.h"

// #include <BAT/BCMath.h>
#include <BAT/BCConstantPrior.h>

// ---------------------------------------------------------
Evidence::Evidence( const std::string& name,
		    int s,
		    int b )
    : BCModel(name)
{
    fS = s;
    fB = b;
    fN = s+b;
    fModel = Model::kBkgOnly;

    fMinE   = 2000.;// keV
    fMaxE   = 2080.;// keV
    fDeltaE = fMaxE - fMinE;

    fMu    = 2039.;// keV
    fSigma = 1.5;// keV
    fLogFactN = LogFactorial( fS + fB );
    
    // Generate data
    fGen = std::mt19937( fRD() );
    fBkg = std::uniform_real_distribution<>(fMinE,fMaxE);
    fSgn = std::normal_distribution<>(fMu,fSigma);

    fE     = std::vector<double>(fS+fB);
    for( int s=0; s<fS; s++ )
	fE[s] = fSgn(fGen);// set vector of event energies
    for( int b=0; b<fB; b++ )
	fE[s+b] = fBkg(fGen);

    fGaussianPDF = std::vector<double>(fS+fB);
    for( int i=0; i<fS+fB; i++ )
	{
	    double e = fE[i];
	    fGaussianPDF[i] = exp( - pow( ( e - fMu ) / fSigma, 2. ) / 2. ) / sqrt(2.*M_PI) / fSigma;
	}
    fBkgPDF    = 1. / fDeltaE;
    
    double errSB = sqrt( fS + fB );
    double nErr = 7.;

    double minB = std::max( 0., fB - nErr * errSB );
    double maxB = fB + 7. * errSB;
    AddParameter( "B", minB, maxB, "B", "[counts]" );
    GetParameters().Back().SetNbins(300);

    double minS = std::max( 0., fS - nErr * errSB );
    double maxS = fS + 7. * errSB;
    AddParameter( "S", minS, maxS, "S", "[counts]" );
    GetParameters().Back().SetNbins(300);

    // Random number machinery for Hit-And-Miss
    for( unsigned int i=0; i<GetNParameters(); i++ )
	{
	    double min = GetParameter(i).GetLowerLimit();
	    double max = GetParameter(i).GetUpperLimit();

	    fSpaceSampler.push_back( std::uniform_real_distribution<>(min,max) );
	}
    fNHitMiss = 1000000;
    fNAcceptedHitMiss = 0;

}

// ---------------------------------------------------------
Evidence::~Evidence()
{
    // destructor
}

double Evidence::LogFactorial( double n )
{
    if( n <= 1. )
	return 0.;
    else
	return log(n) + LogFactorial( n-1. );
}

void Evidence::SetModel( Model model )
{
    fModel = model;
    
    if( fModel == Model::kBkgOnly )
	GetParameter(1).Fix(0.);
    else
	GetParameter(1).Unfix();
    
    return;
}

// ---------------------------------------------------------
double Evidence::LogLikelihood(const std::vector<double>& pars)
{
    double logL = 0;
    if( fModel == Model::kBkgOnly )
	logL = LogLikelihoodBkgOnly( pars );
    else
	logL = LogLikelihoodSgnPlusBkg( pars );
    
    return logL;
}

double Evidence::LogLikelihoodBkgOnly(const std::vector<double>& pars)
{
    double logL = 0.;

    double b = pars[0];
    logL = -b + fN * log(b) - fLogFactN;
    logL += fN * log( fBkgPDF );
	    
    return logL;
}

double Evidence::LogLikelihoodSgnPlusBkg(const std::vector<double>& pars)
{
    double logL = 0.;

    double b = pars[0];
    double s = pars[1];
    double lambda = s + b;
    logL = -lambda - fLogFactN;
    for( auto& fs: fGaussianPDF )// Loop over events
	logL += log( s * fs + b * fBkgPDF );
	    
    return logL;
}

// ---------------------------------------------------------
double Evidence::LogAPrioriProbability(const std::vector<double>& pars)
{
    // return the log of the prior probability p(pars)
    // If you use built-in priors, leave this function commented out.
    double logP = 0;
    for( unsigned int i=0; i<GetNParameters(); i++ )
	if( !GetParameter(i).Fixed() )
	    logP -= log(GetParameter(i).GetRangeWidth());
    
    return logP;     
}

// ---------------------------------------------------------
// void Evidence::CalculateObservables(const std::vector<double>& pars)
// {
//     // Calculate and store obvserables. For example:
//     GetObservable(0) = pow(pars[0], 2);
// }

void Evidence::ComputeIntegral()
{

    
    double maxP = exp( LogLikelihood( GetBestFitParameters() )  + LogAPrioriProbability( GetBestFitParameters() ) );
    fR = std::uniform_real_distribution<>( 0., maxP );

    std::vector<double> pars = std::vector<double>(GetNParameters());
    double p,r;
    fNAcceptedHitMiss = 0;
    for( int i=0; i<fNHitMiss; i++ )
	{
	    for( unsigned int p=0; p<GetNParameters(); p++ )
		if( !GetParameter(p).Fixed() )
		    pars[p] = fSpaceSampler[p](fGen);

	    p = exp( LogLikelihood( pars ) + LogAPrioriProbability( pars ) );
	    r = fR(fGen);

	    if( p >= r )
		fNAcceptedHitMiss ++;
	}


    double volume = maxP;
    for( unsigned int p=0; p<GetNParameters(); p++ )
	if( !GetParameter(p).Fixed() )
	    volume *= GetParameter(p).GetRangeWidth();

    if( fModel == Model::kBkgOnly )
	{
	    fIntegralBkgOnly = (double)fNAcceptedHitMiss / fNHitMiss * volume;

	    BCLog::OutSummary( "Computed integral with Acceptance-Rejection algorithm for Bkg-only model" );
	    BCLog::OutSummary( "Max-P:      " + std::to_string(maxP) );
	    BCLog::OutSummary( "Acceptance: " + std::to_string((double)fNAcceptedHitMiss / fNHitMiss) );
	    //BCLog::OutSummary( "Integral:   " + std::to_string(fIntegralBkgOnly) );
	    std::cout << "Integral: " << fIntegralBkgOnly << std::endl;
	}
    else
	{
	    fIntegralSgnPlusBkg = (double)fNAcceptedHitMiss / fNHitMiss * volume;

	    BCLog::OutSummary( "Computed integral with Acceptance-Rejection algorithm for Sgn+Bkg model" );
	    BCLog::OutSummary( "Max-P:      " + std::to_string(maxP) );
	    BCLog::OutSummary( "Acceptance: " + std::to_string((double)fNAcceptedHitMiss / fNHitMiss) );
	    //BCLog::OutSummary( "Integral:   " + std::to_string(fIntegralSgnPlusBkg) );
	    std::cout << "Integral: " << fIntegralSgnPlusBkg << std::endl;
	}

    return;
}

void Evidence::ComputeBayesFactor()
{
    double P_BkgOnly    = 0.5;
    double P_SgnPlusBkg = 0.5;

    double BayesFactor = fIntegralSgnPlusBkg * P_SgnPlusBkg / fIntegralBkgOnly / P_BkgOnly;
    BCLog::OutSummary( "Bayes factor: " + std::to_string( BayesFactor ) );
    
    return;
}
