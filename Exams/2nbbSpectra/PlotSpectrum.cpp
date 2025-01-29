#include <fstream>
#include <iostream>

#include "TApplication.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"

using namespace std;

int main()
{

    ifstream read( "/home/gbenato/Downloads/m1_2ds.dat" );
    //ifstream read( "/home/gbenato/Downloads/100Mo_ssd_2ds.txt" );
  int ix, iy;
  double ex, ey, p;

  int minIx = 10000;
  int maxIx = 0;
  int minIy = 10000;
  int maxIy = 0;
  while( !read.eof() )
    {
      read >> ix >> iy >> ex >> ey >> p;

      if( ix < minIx ) minIx = ix;
      if( ix > maxIx ) maxIx = ix;
      if( iy < minIy ) minIy = iy;
      if( iy > maxIy ) maxIy = iy;
      
    }

  read.clear();
  read.seekg(0, ios::beg);
  int nBinsX = maxIx - minIx + 1;
  int nBinsY = maxIy - minIy + 1;
  TH2D* h2 = new TH2D( "h2", "h2", nBinsX, 0.25, 3030.25, nBinsY, 0.25, 3030.25 );
  TH1D* h1 = new TH1D( "h1", "h1", nBinsX, 0.25, 3030.25 );
  while( !read.eof() )
    {
      read >> ix >> iy >> ex >> ey >> p;

      h2->SetBinContent( ix, iy, p );
      h1->AddBinContent( ix+iy, p);
      //h2->Fill( ex*1000., ey*1000., p );
      //h1->Fill( ex*1000. + ey*1000., p );
    }

  h1->Scale( 1./ h1->Integral(1.,nBinsX ) );
  cout << "Integral: " << h1->Integral(400,nBinsX) << endl;;
  //cout << "Integral: " << h1->Integral(1,nBinsX) << endl;;
  
  TApplication* app = new TApplication( "app", NULL, 0 );
  TCanvas* can = new TCanvas( "can", "can", 1600, 900 );
  can->Divide(2,1);
  can->cd(1);
  h2->Draw("colz");
  can->cd(2);
  h1->Draw();
  app->Run();

  return 0;
}
