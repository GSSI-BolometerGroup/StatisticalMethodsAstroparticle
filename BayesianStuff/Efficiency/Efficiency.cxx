// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************
#include <random>

#include "Efficiency.h"

#include <BAT/BCConstantPrior.h>
// #include <BAT/BCMath.h>

// ---------------------------------------------------------
Efficiency::Efficiency( const std::string& name,
			const int n,
			const double eff )
    : BCModel(name)
{
    fMethod = FitMethod::kPoisson;
    
    fN          = n;
    fMinE       = 1.;// keV
    fMaxE       = 20.;// keV
    fThreshold  = 3.;// keV
    fSigma      = 1.;// keV
    fEfficiency = eff;
    // Define the input (true) efficiency curve
    fEfficiencyCurve = TF1( "EfficiencyCurve", "0.5*[0]*(1.+TMath::Erf((x-[1])/[2]))", -1.+fMinE, 1.+fMaxE );
    fEfficiencyCurve.SetParameter(0,fEfficiency);
    fEfficiencyCurve.SetParameter(1,fThreshold);
    fEfficiencyCurve.SetParameter(2,fSigma);

    // Random number machinery
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0,1);

    // Populate points
    for( int e=fMinE; e<=fMaxE; e++ )
	{
	    double k = 0;
	    for( int i=0; i<fN; i++ )
		if( dis(gen) < EvaluateEfficiency( e,
						   fEfficiency,
						   fThreshold,
						   fSigma ) )
		    k += 1.;
	    
	    fK[e] = k;
	}


    // Add parameters
    double minEff = std::max( 0., 0.95*fEfficiency );
    double maxEff = std::min( 1., 1.05*fEfficiency );
    AddParameter( "Efficiency", minEff, maxEff, "#epsilon", "" );
    GetParameters().Back().SetPrior( new BCConstantPrior() );// Flat between 0 and 1 for the efficiency
    GetParameters().Back().SetNbins(300);

    AddParameter( "Threshold", 0.9*fThreshold, 1.1*fThreshold, "E_{thr}", "[keV]" );
    GetParameters().Back().SetPrior( new BCConstantPrior() );
    GetParameters().Back().SetNbins(300);

    AddParameter( "Sigma", 0.9*fSigma, 1.1*fSigma, "#sigma", "[keV]" );
    GetParameters().Back().SetPrior( new BCConstantPrior() );
    GetParameters().Back().SetNbins(300);

    
}

// ---------------------------------------------------------
Efficiency::~Efficiency()
{
    // destructor
}

// ---------------------------------------------------------
double Efficiency::LogLikelihood(const std::vector<double>& pars)
{
    // return the log of the conditional probability p(data|pars).

    double logL = 0.;
    if( fMethod == FitMethod::kBinomial )
	{
	    double E, k, p;
	    for( auto& it: fK )
		{
		    E = it.first;
		    k = it.second;
		    p = EvaluateEfficiency( E, pars[0], pars[1], pars[2] );
		    logL += LogBinomial( k, fN, p );
		}
	}
    else if( fMethod == FitMethod::kPoisson )
	{
	    double E, k, lambda;
	    for( auto& it: fK )
		{
		    E = it.first;
		    k = it.second;
		    lambda = EvaluateEfficiency( E, pars[0], pars[1], pars[2] ) * fN;
		    logL += LogPoisson( k, lambda );
		}
	}
    else if( fMethod == FitMethod::kChiSquare )
	{
	    double E, k, mu;
	    for( auto& it: fK )
		{
		    E = it.first;
		    k = it.second;
		    mu = EvaluateEfficiency( E, pars[0], pars[1], pars[2] ) * fN;
		    logL += LogChi2( k, mu );
		}
	}
    
    return logL;
}

// ---------------------------------------------------------
// double Efficiency::LogAPrioriProbability(const std::vector<double>& pars)
// {
//     // return the log of the prior probability p(pars)
//     // If you use built-in priors, leave this function commented out.
// }

// ---------------------------------------------------------
// void Efficiency::CalculateObservables(const std::vector<double>& pars)
// {
//     // Calculate and store obvserables. For example:
//     GetObservable(0) = pow(pars[0], 2);
// }

double Efficiency::EvaluateEfficiency( double E,
				       double asympEff,
				       double threshold,
				       double sigma )
{
    return 0.5 * asympEff * ( 1. + erf( ( E-threshold ) / sigma ) );
}

double Efficiency::LogBinomial( double k, double n, double p )
{
    // To save time, we drop the Binomial coefficient, which does not depend on p.
    // This reduces the CPU time by a factor ~20
    double exp_1 = k * log(p);
    double exp_2 = (n-k) * log(1.-p);

    return exp_1 + exp_2;

    // Full formula, including the Binomial coefficient (sum1-sum2)
    //double sum_1 = 0;
    //double sum_2 = 0;
    //for( double i=(int)(n-k+1); i<=n; i+=1. )
    //    sum_1 += log(i);
    //for( double i=1.; i<=k; i+=1. )
    //    sum_2 += log(i);
    //double exp_1 = k * log(p);
    //double exp_2 = (n-k) * log(1.-p);
    //return sum_1 - sum_2 + exp_1 + exp_2;
}

// Computes the logarithm of the factorial.
double Efficiency::LogFactorial( double n )
{
    if( n <= 1. )
        return 0.;
    return log(n) + LogFactorial(n-1.);
}

double Efficiency::LogPoisson( double n,
			       double lambda )
{
    // Faster formula, neglecting the factorial term which does not depend on lambda
    return -lambda + n * log(lambda);
    
    // Full formula, including factorial term
    //return exp( -lambda + n * log(lambda) - LogFactorial(n) );
}

double Efficiency::LogChi2( double n, double mu )
{
    return -0.5 * pow( n-mu, 2. ) / mu - log( sqrt(2.*M_PI) * mu );
}

void Efficiency::SetFitMethod( FitMethod method )
{
    fMethod = method;
    
    return;
}
