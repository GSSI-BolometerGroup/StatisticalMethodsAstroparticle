#include <random>

#include "BAT/BCConstantPrior.h"
#include "BAT/BCGaussianPrior.h"

#include "RadioactiveDecayFit.h"


RadioactiveDecayFit::RadioactiveDecayFit(const std::string& name)
    : BCModel(name)
{
    fMethod = FitMethod::kBinomial;

    fN = 1000;
    fHalflife = 138.376;// days
    fDecayRate = log(2.) / fHalflife;
    fDeltaT = 2. * fHalflife;

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
    double minN = fK;
    double maxN = fN + 10. * sqrt(fN);//2. * fN;
    AddParameter( "N", minN, maxN, "N", "" );
    GetParameters().Back().SetPrior( new BCConstantPrior() );
    GetParameters().Back().SetNbins(300);

    // Create model parameter for the decay rate
    double minR = 0.7 * fDecayRate;// The min and max for the decay rate are just guessed.
    double maxR = 1.3 * fDecayRate;
    AddObservable( "R", minR, maxR, "$Gamma", "d^{-1}" );
    GetObservables().Back().SetNbins(300);


    double minHalflife = log(2.) / maxR;
    double maxHalflife = log(2.) / minR;
    AddObservable( "T1/2", minHalflife, maxHalflife, "T_{1/2}", "d" );
    GetObservables().Back().SetNbins(300);

    AddObservable( "N", minN, maxN, "N", "" );
    GetObservables().Back().SetNbins(300);
}

RadioactiveDecayFit::~RadioactiveDecayFit()
{
    ;
}

double RadioactiveDecayFit::LogLikelihood(const std::vector<double>& pars)
{
    // We define here three a binned likelihood with a product of P(N|k) for each bin.
    // We implement three possible likelihood distributions:
    // 1) Binomial for each bin
    // 2) Poisson
    // 3) Multinomial

    double N = (int)(pars[0]+0.5);// Retrieve the parameter N, parse it to integer
    double P = fK / pars[0];// In this case, we don't parse it to integer to maintain the precision on P
    double R = -log(1.-P) / fDeltaT;// R is the decay rate
    double logL = 0.;
    if( fMethod == FitMethod::kBinomial ||
	fMethod == FitMethod::kPoisson )
	for( int b=1; b<=fDataHisto->GetNbinsX(); b++ )
	    {
		double k  = fDataHisto->GetBinContent(b);
		// First option to compute p_i
		double t1 = fDataHisto->GetXaxis()->GetBinLowEdge(b);
		double t2 = fDataHisto->GetXaxis()->GetBinLowEdge(b+1);
		double p  = exp( -R * t1 ) - exp( -R * t2 );
		// Second option to compute p_i, approximating the integral of f(lambda|t)
		// with f(lambda|t)*dt
		//double p  = R * exp( -R * t1 ) * fDataHisto->GetBinWidth(b);
		// Third optoin to compute p_i. This works if k is large for every bin,
		// but fails miserably if k=0.
		//double p = k/pars[0];
		double lambda = p * N;
		if( fMethod == FitMethod::kBinomial )
		    logL += LogBinomial( k, N, p );
		else if( fMethod == FitMethod::kPoisson )
		    logL += LogPoisson( k, lambda );
	    }
    else if( fMethod == FitMethod::kMultinomial )
	{
	    // First, N!
	    logL += LogFactorial(N);
	    // Second, the term corresponding to undetected decays happening after fDeltaT.
	    // This would be ( (1-P)^(N-fK) )/(N-fK)!
	    logL -= LogFactorial( N - fK );
	    logL += ( N - fK ) * log( 1. - P );

	    // Now the term for each bin
	    for( int b=1; b<=fDataHisto->GetNbinsX(); b++ )
		{
		    double k  = fDataHisto->GetBinContent(b);
		    double t1 = fDataHisto->GetXaxis()->GetBinLowEdge(b);
		    double t2 = fDataHisto->GetXaxis()->GetBinLowEdge(b+1);
		    double p  = exp( -R * t1 ) - exp( -R * t2 );
		    logL += k * log(p);
		    logL -= LogFactorial(k);
		}
	}

    return logL;
}

void RadioactiveDecayFit::CalculateObservables(const std::vector<double>& pars)
{
    double P = fK / pars[0];
    double R = -log(1.-P) / fDeltaT;
    GetObservable(0) = R;
    GetObservable(1) = log(2.) / R;

    return;
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
    for( double i=1.; i<=k; i+=1. )// This term can be skipped because it just depends on k,
    	sum_2 += log(i);           // which is a constant in any case.
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
