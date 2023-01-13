#include <random>

#include <BAT/BCConstantPrior.h>

#include "SimultaneousFit.h"

SimultaneousFit::SimultaneousFit( const std::string& name,
				  int nInjected,
				  double eff,
				  int s,
				  int b )
    : BCModel(name)
{


    fNInjected  = nInjected;
    fEfficiency = eff;// True value for the asymptotic efficiency
    fS          = s;// True value of signal events
    fB          = b;// True value of background events
    fMinEEff    = 1.;// Minimum energy of pulser-injected events, keV
    fMaxEEff    = 20.;// Maximum energy of pulser-injected events, keV
    fThreshold  = 3.;// Trigger theshold keV
    fSigmaEff   = 1.;// Sigma of the efficiency curve, keV

    // Define the input (true) efficiency curve
    fEfficiencyCurve = TF1( "EfficiencyCurve", "0.5*[0]*(1.+TMath::Erf((x-[1])/[2]))", -1.+fMinE, 1.+fMaxE );
    fEfficiencyCurve.SetParameter(0,fEfficiency);
    fEfficiencyCurve.SetParameter(1,fThreshold);
    fEfficiencyCurve.SetParameter(2,fSigmaEff);

    // Random number machinery
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0,1);

    // Populate pulser-injected events
    for( int e=fMinEEff; e<=fMaxEEff; e++ )
	{
	    double k = 0;
	    for( int i=0; i<fNInjected; i++ )
		if( dis(gen) < EvaluateEfficiency( e,
						   fEfficiency,
						   fThreshold,
						   fSigmaEff ) )
		    k += 1.;
	    
	    fK[e] = k;
	}

    // Define the fit range
    fMinE   = 10.;// keV
    fMaxE   = 20.;// keV
    fDeltaE = fMaxE - fMinE;// Width of the energy range, keV
    // Define the signal shape (which is fixed, so we will not fit these two parameters).
    fMu     = 15.;// Mean of the Gaussian signal, keV
    fSigma  = 1.2;// Sigma of the Gaussian signal, keV

    // Random generator for signal and background distributions
    std::uniform_real_distribution<> bkg( fMinE, fMaxE );
    std::normal_distribution<> sgn( fMu, fSigma );

    // Generate the energy for s signal events, then check whether they are accepted or not
    for( int s=0; s<fS; s++ )
	{
	    double e = sgn(gen);
	    if( dis(gen) < EvaluateEfficiency( e,
					       fEfficiency,
					       fThreshold,
					       fSigmaEff ) )
		fE.push_back(e);
	}

    // Generate the energy for b signal events, then check whether they are accepted or not
    for( int b=0; b<fB; b++ )
	{
	    double e = bkg(gen);
	    if( dis(gen) < EvaluateEfficiency( e,
					       fEfficiency,
					       fThreshold,
					       fSigmaEff ) )	    
		fE.push_back(e);
	}

    // Compute f_s(E) for all events. Given that we have fixed fMu and fSigma,
    // we can compute this quantities just once to speed up the log-L evaluation.
    for( auto& e: fE )
	fGaussianPDF.push_back( exp( - pow( ( e - fMu ) / fSigma, 2. ) / 2. ) / sqrt(2.*M_PI) / fSigma );
    // Compute f_b(E), which actually does not depend on the energy of the event.
    fBkgPDF = 1. / fDeltaE;

    // Create fit parameters for S, B, efficiency, threshold, and sigma of the efficiency curve.
    // Use flat priors for all of them.
    double errSB = sqrt((double)(fS+fB));
    double nErr = 7.;
    double minS = std::max( 0., fS - nErr * errSB );
    double maxS = fS + nErr * errSB;
    AddParameter( "S", minS, maxS, "S", "[counts]" );
    GetParameters().Back().SetPrior( new BCConstantPrior() );
    GetParameters().Back().SetNbins(100);

    double minB = std::max( 0., fB - nErr * errSB );
    double maxB = fB + nErr * errSB;
    AddParameter( "B", minB, maxB, "B", "[counts]" );
    GetParameters().Back().SetPrior( new BCConstantPrior() );
    GetParameters().Back().SetNbins(100);

    // Add parameters
    double minEff = std::max( 0., 0.95 * fEfficiency );
    double maxEff = std::min( 1., 1.05 * fEfficiency );
    AddParameter( "Efficiency", minEff, maxEff, "#epsilon", "" );
    GetParameters().Back().SetPrior( new BCConstantPrior() );// Flat between 0 and 1 for the efficiency
    GetParameters().Back().SetNbins(100);

    AddParameter( "Threshold_eff", 0.9*fThreshold, 1.1*fThreshold, "E_{thr}", "[keV]" );
    GetParameters().Back().SetPrior( new BCConstantPrior() );
    GetParameters().Back().SetNbins(100);

    AddParameter( "Sigma_eff", 0.75*fSigmaEff, 1.25*fSigmaEff, "#sigma_{#epsilon}", "[keV]" );
    GetParameters().Back().SetPrior( new BCConstantPrior() );
    GetParameters().Back().SetNbins(100);
    
}

// ---------------------------------------------------------
SimultaneousFit::~SimultaneousFit()
{
    // destructor
}

// ---------------------------------------------------------
double SimultaneousFit::LogLikelihood(const std::vector<double>& pars)
{

    double logL = 0.;

    // Retrieve the current value of the fit parameters
    // in the Metropolis-Hastings Markov Chain
    double s        = pars[0];
    double b        = pars[1];
    double eff      = pars[2];
    double thr      = pars[3];
    double sigmaeff = pars[4];

    // Binomial log-L for efficiency data
    double E, k, p;
    for( auto& it: fK )
	{
	    E = it.first;
	    k = it.second;
	    p = EvaluateEfficiency( E, eff, thr, sigmaeff );
	    logL += LogBinomial( k, fNInjected, p );
	}

    // Extended Poisson log-L for physics data.
    // Given that the physics data start at 10 keV where the efficiency curve is already flat,
    // we can directly weight s+b with the asymptotic efficiency value.
    double lambda = eff * ( s + b );
    logL -= lambda;
    for( auto& fs: fGaussianPDF )// Loop over events
    	logL += log( s * eff * fs + b * eff * fBkgPDF );

    return logL;
}

double SimultaneousFit::EvaluateEfficiency( double E,
				       double asympEff,
				       double threshold,
				       double sigma )
{
    return 0.5 * asympEff * ( 1. + erf( ( E-threshold ) / sigma ) );
}

double SimultaneousFit::LogBinomial( double k, double n, double p )
{
    // To save time, we drop the Binomial coefficient, which does not depend on p.
    // This reduces the CPU time by a factor ~20

    double exp_1 = k * log(p);
    double exp_2 = (n-k) * log(1.-p);

    return exp_1 + exp_2;
    /*
    // Full formula, including the Binomial coefficient (sum1-sum2)
    double sum_1 = 0;
    double sum_2 = 0;
    for( double i=(int)(n-k+1); i<=n; i+=1. )
        sum_1 += log(i);
    for( double i=1.; i<=k; i+=1. )
        sum_2 += log(i);
    double exp_1 = k * log(p);
    double exp_2 = (n-k) * log(1.-p);
    return sum_1 - sum_2 + exp_1 + exp_2;
    */
}

double SimultaneousFit::LogFactorial( double n )
{
    if( n <= 1. )
        return 0.;
    
    return log(n) + LogFactorial(n-1.);
} 
