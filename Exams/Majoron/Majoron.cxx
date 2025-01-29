#include "TH1D.h"
#include "TRandom3.h"

#include "Majoron.h"

#include "BAT/BCConstantPrior.h"
// #include <BAT/BCMath.h>

Majoron::Majoron(const std::string& name)
    : BCModel(name)
{
    ReadFiles();

    CreateData();
    
    // Define parameters here in the constructor.
    // Also define their priors, if using built-in priors.
    // For example:
    // AddParameter("mu", -2, 1, "#mu", "[GeV]");
    // GetParameters.Back().SetPrior(new BCGaussianPrior(-1, 0.25));

    // Define observables here, too. For example:
    // AddObservable("mu_squared", 1, 4, "#mu^{2}", "[GeV^{2}]");

    double minB = std::max( 0., fB - 7.*sqrt(fB+fS) );
    double maxB = fB + 7.*sqrt(fB+fS);
    AddParameter( "B", minB, maxB, "B", "[counts]" );
    GetParameters().Back().SetPrior(new BCConstantPrior());
    
    double minS = 0.;//std::max( 0., fS - 3.*sqrt(fB+fS) );
    double maxS = 10.*fS;//fS + 3.*sqrt(fB+fS);
    AddParameter( "S", minS, maxS, "S", "[counts]" );
    GetParameters().Back().SetPrior(new BCConstantPrior());

}

void Majoron::SetModel( Model m )
{
    fModel = m;
    if( m == Model::kH0 )
	GetParameter("S").Fix(0.);

    return;
}

void Majoron::CreateData()
{
    fB = 1.e10;//2.96e8;
    fS = 100.;

    TRandom3 rdm(12345);
    
    fData = std::vector<int>(fNBins);
    for( size_t i=0; i<fNBins; i++ )
	{
	    fData[i] = rdm.Poisson( fB * fPDF2nbb[i] ) + rdm.Poisson( fS * fPDFMajoron[i] );
	    //std::cout << i << "\t" << fData[i] << std::endl;
	    //fLogFactorial[fData[i]] = LogFactorial( fData[i] );
	    //std::cout << i << "\t" << fData[i] << "\t" << fLogFactorial[fData[i]] << std::endl;
	}

    FillLogFactorial();
    
    return;
}

void Majoron::FillLogFactorial()
{
    int maxN = 0;
    for( auto& it: fData )
	if( it > maxN ) maxN = it;

    long double logF = 0.;
    fLogFactorial[0] = 0.;
    fLogFactorial[1] = 0.;
    for( int i=2; i<=maxN; i++ )
	{
	    logF += log((long double)i);
	    if( std::find( fData.begin(), fData.end(), i ) != fData.end() )
		fLogFactorial[i] = logF;
	}
    
    return;
}

long double Majoron::LogFactorial( int n )
{
    long double logF = 0;
    for( int i=2; i<=n; i++ )
	logF += log((long double)i);
    return logF;
    /*
    if( n <= 1 )
	return 0.;
    else
	return log((double)n) + LogFactorial(n-1);
    */
}

Majoron::~Majoron()
{
    ;
}

double Majoron::LogLikelihood(const std::vector<double>& pars)
{
    // return the log of the conditional probability p(data|pars).
    // This is where you define your model.
    // BCMath contains many functions you will find helpful.

    double b = pars[0];
    double s = pars[1];
    long double logL = 0;
    long double lambda;
    for( size_t i=0; i<fNBins; i++ )
	{
	    lambda = b * fPDF2nbb[i] + s * fPDFMajoron[i];
	    if( lambda == 0 ) continue;
	    
	    logL += log(lambda) * fData[i] - lambda - fLogFactorial[fData[i]];
	    //std::cout << log(lambda) * fData[i] << "\t"
	    //	      << -lambda << "\t"
	    //	      << -fLogFactorial[fData[i]] << "\t"
	    //	      << log(lambda) * fData[i] - lambda - fLogFactorial[fData[i]] << std::endl;
	}
    //std::cout << "----" << std::endl;
    std::cout << logL << "\t" << exp(logL) << std::endl;
    //std::cout << std::numeric_limits<long double>::min() << std::endl;
    //std::cout << "----" << std::endl;
    return logL;
}

// ---------------------------------------------------------
// double Majoron::LogAPrioriProbability(const std::vector<double>& pars)
// {
//     // return the log of the prior probability p(pars)
//     // If you use built-in priors, leave this function commented out.
// }

// ---------------------------------------------------------
// void Majoron::CalculateObservables(const std::vector<double>& pars)
// {
//     // Calculate and store obvserables. For example:
//     GetObservable(0) = pow(pars[0], 2);
// }

void Majoron::ReadFiles()
{
    // Read 2nbb data
    std::ifstream read2nbb( "/home/gbenato/Downloads/100Mo_ssd_2ds.txt" );
    int ix, iy;
    double ex, ey, p;

    int minIx = 10000;
    int maxIx = 0;
    int minIy = 10000;
    int maxIy = 0;
    while( !read2nbb.eof() )
	{
	    read2nbb >> ix >> iy >> ex >> ey >> p;

	    if( ix < minIx ) minIx = ix;
	    if( ix > maxIx ) maxIx = ix;
	    if( iy < minIy ) minIy = iy;
	    if( iy > maxIy ) maxIy = iy;
	}
    
    read2nbb.clear();
    read2nbb.seekg(0, std::ios::beg);
    int nBinsX = maxIx - minIx + 1;
    TH1D* h1_2nbb = new TH1D( "h1_2nbb", "h1_2nbb", nBinsX, 0.25, 3030.25 );
    while( !read2nbb.eof() )
	{
	    read2nbb >> ix >> iy >> ex >> ey >> p;
	    h1_2nbb->AddBinContent( ix+iy, p);
	}
    double int_2nbb = h1_2nbb->Integral(1,fNBins);
    fBinWidth = h1_2nbb->GetBinWidth(1);
    for( int b=1; b<=nBinsX; b++ )
	fPDF2nbb.push_back( h1_2nbb->GetBinContent(b) / int_2nbb );

    // Read Majoron data
    std::ifstream readMajoron( "/home/gbenato/Downloads/m1_2ds.dat" );

    minIx = 10000;
    maxIx = 0;
    minIy = 10000;
    maxIy = 0;
    while( !readMajoron.eof() )
	{
	    readMajoron >> ix >> iy >> ex >> ey >> p;

	    if( ix < minIx ) minIx = ix;
	    if( ix > maxIx ) maxIx = ix;
	    if( iy < minIy ) minIy = iy;
	    if( iy > maxIy ) maxIy = iy;
	}
    
    readMajoron.clear();
    readMajoron.seekg(0, std::ios::beg);
    nBinsX = maxIx - minIx + 1;
    TH1D* h1_majoron = new TH1D( "h1_majoron", "h1_majoron", nBinsX, 0.25, 3030.25 );
    while( !readMajoron.eof() )
	{
	    readMajoron >> ix >> iy >> ex >> ey >> p;
	    h1_majoron->AddBinContent( ix+iy, p);
	}
    double int_majoron = h1_majoron->Integral(1,fNBins);
    for( int b=1; b<=nBinsX; b++ )
	fPDFMajoron.push_back( h1_majoron->GetBinContent(b) / int_majoron );

    fNBins = fPDFMajoron.size();

    return;
}
