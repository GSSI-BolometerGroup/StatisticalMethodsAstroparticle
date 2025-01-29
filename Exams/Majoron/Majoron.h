#ifndef __BAT__MAJORON__H
#define __BAT__MAJORON__H

#include <BAT/BCModel.h>

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <map>

class Majoron : public BCModel
{
public:
    enum Model{ kH0, kH1 };
    
private:
    void ReadFiles();
    void CreateData();
    long double LogFactorial( int n );
    void FillLogFactorial();
    
    std::vector<long double> fPDF2nbb;
    std::vector<long double> fPDFMajoron;
    std::vector<int> fData;
    std::map<int,long double> fLogFactorial;
    double fBinWidth;
    size_t fNBins;

    double fB;
    double fS;
    Model fModel;
public:
    void SetModel( Model m );
    Majoron(const std::string& name);
    
    ~Majoron();

    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double>& pars);

    // Overload LogAprioriProbability if not using built-in 1D priors
    // double LogAPrioriProbability(const std::vector<double> & pars);

    // Overload CalculateObservables if using observables
    // void CalculateObservables(const std::vector<double> & pars);

};
// ---------------------------------------------------------

#endif
