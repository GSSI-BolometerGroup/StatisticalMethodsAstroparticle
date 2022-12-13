#include <random>

#include "BAT/BCConstantPrior.h"

#include "RadioactiveDecayFit.h"


RadioactiveDecayFit::RadioactiveDecayFit(const std::string& name)
    : BCModel(name)
{
    fMethod = FitMethod::kBinomial;

    fN = 1000;
    fHalflife = 138.376;// days
    fDecayRate = log(2.) / fHalflife;
    fDeltaT = 3. * fHalflife;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::exponential_distribution dis(fDecayRate);

    // Create event times
    double t;
    for( int n=0; n<fN; n++ )
	{
	    t = dis(gen);
	    if( t < fDeltaT )
		fT.push_back(t);
	}
    fK = fT.size();
    BCLog::OutSummary( "Number of detected events (fK): " + std::to_string(fK) );

    // Computed P and lambda for Binomial and Poisson distributions, respectively
    fP      = 1. - exp( -fDecayRate * fDeltaT );
    fLambda = fN * fP;

    // Fill histogram with event times. This will not be used for the fit,
    // but just for some final plots.
    fDataHisto = new TH1D( "data", "data", 100, 0., fDeltaT );
    fDataHisto->GetXaxis()->SetTitle( "t [days]" );
    fDataHisto->GetYaxis()->SetTitle( "Events" );
    for( auto& t: fT )
	fDataHisto->Fill(t);

    // Create model parameter N
    double minN = std::max( fK, (size_t)(fN - 7. * sqrt(fN)) );// N cannot be smaller than K
    double maxN = fN + 7. * sqrt(fN);//2. * fN;
    AddParameter( "N", minN, maxN, "N", "" );
    GetParameters().Back().SetPrior( new BCConstantPrior() );
    GetParameters().Back().SetNbins(300);
    
}

RadioactiveDecayFit::~RadioactiveDecayFit()
{
    ;
}

double RadioactiveDecayFit::LogLikelihood(const std::vector<double>& pars)
{
    // We define here an unbinned extended likelihood of type:
    // L(k|N,decayrate) = P(k|N,decayrate) * prod_i[f(t_i)]
    // where:
    // P(k|N,decayrate) is the Binomial or Poisson probability of measuring n events
    //                  given the decay rate and an original number of isotopes N
    // f(t_i)=decayrate*exp(-decayrate*t_i) is the probability for the decay to happen at time t_i
    //
    // Notice that our parameter of interest is the original number of isotopes N,
    // which enters only the Binomial or the Poisson term, and not prod_i[f(t_i)]
    
    double N = pars[0];// Retrieve the parameter
    double logL = 0.;
    // Extended term
    if( fMethod == FitMethod::kBinomial )
	{
	    logL = LogBinomial( fK, N, fP );// N enters directly into the Binomial
	}
    else if( fMethod == FitMethod::kPoisson )
	{
	    double lambda = N * fP;// N enters the binomial through lambda
	    logL = LogPoisson( fK, lambda );
	}
    // Single event terms (which are actually constant in this case, so we could also neglect them)
    logL += fK * log(fDecayRate);
    for( auto& t: fT )
    	logL -= fDecayRate * t;

    
    return logL;
}

void RadioactiveDecayFit::SetFitMethod( FitMethod method )
{
    fMethod = method;
    
    return;
}

double RadioactiveDecayFit::LogBinomial( double k, double n, double p )
{
    // This time we cannot use the shorter formula, because we are interested
    // in the evalutaion of n, which enters the factorials.
    //double exp_1 = k * log(p);
    //double exp_2 = (n-k) * log(1.-p);
    //return exp_1 + exp_2;

    // Full formula, including the Binomial coefficient (sum1-sum2)
    double sum_1 = 0;
    double sum_2 = 0;
    for( double i=(int)(n-k+1); i<=n; i+=1. )
        sum_1 += log(i);
    //for( double i=1.; i<=k; i+=1. )// This term can be skipped because it just depends on k,
    //  sum_2 += log(i);             // which is a constant in any case.
    double exp_1 = k * log(p);
    double exp_2 = (n-k) * log(1.-p);
    return sum_1 - sum_2 + exp_1 + exp_2;


}

// Computes the logarithm of the factorial.
double RadioactiveDecayFit::LogFactorial( double n )
{
    if( n <= 1. )
        return 0.;
    return log(n) + LogFactorial(n-1.);
}

double RadioactiveDecayFit::LogPoisson( double n,
					double lambda )
{
    // Faster formula, neglecting the factorial term which does not depend on lambda.
    return -lambda + n * log(lambda);
    
    // Full formula, including factorial term
    //return -lambda + n * log(lambda) - LogFactorial(n);
}

TH1D* RadioactiveDecayFit::GetData()
{
    return fDataHisto;
}

TF1* RadioactiveDecayFit::GetBestFit()
{
    double binwidth = fDataHisto->GetBinWidth(1);
    double min = fDataHisto->GetXaxis()->GetXmin();
    double max = fDataHisto->GetXaxis()->GetXmax();
    TF1* f = new TF1( "BestFit", "[0]*[1]*[2]*exp(-[2]*x)", min, max );
    f->SetParameter( 0, GetBestFitParameters()[0] );
    f->SetParameter( 1, binwidth );
    f->SetParameter( 2, fDecayRate );

    return f;
}
