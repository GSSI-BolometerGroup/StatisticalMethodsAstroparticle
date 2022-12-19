#include "BAT/BCConstantPrior.h"

#include "LikelihoodRatio.h"

LikelihoodRatio::LikelihoodRatio(const std::string& name)
    : BCModel(name)
{
    // Let's create a Gaussian peak centered at 2039 keV (Qbb value of 76Ge)
    // and a flat backround in the 2030-2050 keV region.
    fMinE   = 2030.;
    fMaxE   = 2050.;
    fDeltaE = fMaxE - fMinE;
    fDE     = fDeltaE / 40.;
    fNBins  = fDeltaE / fDE;
    fK      = std::vector<double>(fNBins);
    for( auto& k: fK )
	k = 0.;
    fS      = 1000;
    fB      = 1000;
    fMu     =  2039.;// keV
    fSigma  = 1.5;// keV;

    // Compute the signal PDF for each bin.
    // This will be used in the LogLikelihood method.
    // We compute it now just once to speed up the fit.
    fGaussianPDF = std::vector<double>(fNBins);
    for( int b=0; b<fNBins; b++ )
	{
	    double E = fMinE + fDE * ( 0.5 + b );
	    fGaussianPDF[b] = exp( - 0.5 * std::pow( (E-fMu) / fSigma, 2. ) ) * fDE / sqrt( 2. * M_PI ) / fSigma;
	}

    // Initialize random number machinery
    fGen = std::mt19937(fRD());
    fBkg = std::uniform_real_distribution<>(fMinE,fMaxE);
    fSgn = std::normal_distribution(fMu,fSigma);

    // Define parameters for number of signal and background counts
    double errSB = sqrt( fS + fB );
    double minS = std::max( 0., fS - 7. * errSB );
    double maxS = fS + 7. * errSB;
    AddParameter( "S", minS, maxS, "S", "[counts]" );
    GetParameters().Back().SetPrior( new BCConstantPrior() );
    GetParameters().Back().SetNbins(300);

    double minB = std::max( 0., fB - 7. * errSB );
    double maxB = fB + 7. * errSB;
    AddParameter( "B", minB, maxB, "B", "[counts]" );
    GetParameters().Back().SetPrior( new BCConstantPrior() );
    GetParameters().Back().SetNbins(300);
	
}

void LikelihoodRatio::ResetData()
{
    for( auto& k: fK )
	k=0;
    
    return;
}

void LikelihoodRatio::SetData()
{
    ResetData();

    // Populate spectrum with signal events
    for( int s=0; s<fS; s++ )
	{
	    double E = fSgn(fGen);
	    int i = ( E - fMinE ) / fDE;
	    fK[i] += 1.;
	}
    // Populate spectrum with backround events
    for( int b=0; b<fB; b++ )
	{
	    double E = fBkg(fGen);
	    int i = ( E - fMinE ) / fDE;
	    fK[i] += 1.;
	}

    return;
}

// ---------------------------------------------------------
LikelihoodRatio::~LikelihoodRatio()
{
    // destructor
}

// ---------------------------------------------------------
double LikelihoodRatio::LogLikelihood(const std::vector<double>& pars)
{

    // Compute -Chi2_lambda, which is equal to 2*log[(likelihood)/(max-likelihood)].
    // We approximate max-likelihood imposing lambda_i = k_i,
    // which means that we impose the expectation value for the bin i
    // to be equal to the number of counts in bin i.
    // This is the method proposed by Baker and Cousins in Nucl. Instrum. Methods 221 (1984) 437-442.
    // It works well if all bins have at least a few counts,
    // but fails if a single bin has zero entries.
    double logL = 0.;
    double S = pars[0];
    double B = pars[1];
    double lambda;
    
    for( int b=0; b<fNBins; b++ )
	{
	    double k = fK[b];
	    lambda = B * fDE / fDeltaE + S * fGaussianPDF[b];
	    logL += k - lambda - k * log( k/lambda );// This will fail if k=0.
	}
    logL *= 2.;

    return logL;
}

TH1D* LikelihoodRatio::GetDataHisto()
{
    fDataHisto = new TH1D( "Data", "Data", fNBins, fMinE, fMaxE );
    for( int b=0; b<fNBins; b++ )
	fDataHisto->SetBinContent( b+1, fK[b] );
    fDataHisto->GetXaxis()->SetTitle("Energy [keV]");
    fDataHisto->GetYaxis()->SetTitle("Counts");
    fDataHisto->GetYaxis()->SetRangeUser( 0., 1.1*fDataHisto->GetMaximum() );
    
    return fDataHisto;
}

TF1* LikelihoodRatio::GetFittingFunction()
{
    TF1* func = new TF1( "BestFit", "[0]*[1]/[2] + [1]*[3]/sqrt(2.*TMath::Pi())/[5]*exp( -0.5 * pow( (x-[4])/[5], 2. ) )", fMinE, fMaxE );
    func->SetParameter( 0, GetBestFitParameters()[1] );
    func->SetParameter( 1, fDE );
    func->SetParameter( 2, fDeltaE );
    func->SetParameter( 3, GetBestFitParameters()[0] );
    func->SetParameter( 4, fMu );
    func->SetParameter( 5, fSigma );

    return func;
}

int LikelihoodRatio::GetNDF()
{
    return fNBins - 2;
}
