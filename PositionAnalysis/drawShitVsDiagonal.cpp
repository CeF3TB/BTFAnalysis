#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TVector2.h"

#include "interface/DrawTools.h"



TGraphErrors* getGraph( const std::string varName, std::vector<std::string> runs, std::vector<float> beamX, std::vector<float> beamY );


int main() {

  DrawTools::setStyle();

  std::vector<std::string> runs02;
  std::vector<float> beamX02;
  std::vector<float> beamY02;

  runs02.push_back("BTF_94_20140430-073300_beam");
  beamX02.push_back(-3.);
  beamY02.push_back(+3.);

  runs02.push_back("BTF_96_20140430-083733_beam");
  beamX02.push_back(-6.);
  beamY02.push_back(+6.);

  runs02.push_back("BTF_98_20140430-092026_beam");
  beamX02.push_back(-9.);
  beamY02.push_back(+9.);

  runs02.push_back("BTF_138_20140430-175224_beam");
  beamX02.push_back(+0.);
  beamY02.push_back(+0.);



  std::vector<std::string> runs13;
  std::vector<float> beamX13;
  std::vector<float> beamY13;

  runs13.push_back("BTF_138_20140430-175224_beam");
  beamX13.push_back(+0.);
  beamY13.push_back(+0.);

  runs13.push_back("BTF_100_20140430-101607_beam");
  beamX13.push_back(-9.);
  beamY13.push_back(-9.);
  
  runs13.push_back("BTF_118_20140430-151237_beam");
  beamX13.push_back(-6.);
  beamY13.push_back(-6.);
  
  runs13.push_back("BTF_136_20140430-171004_beam");
  beamX13.push_back(-3.);
  beamY13.push_back(-3.);
  
  runs13.push_back("BTF_141_20140430-183508_beam");
  beamX13.push_back(+3.);
  beamY13.push_back(+3.);
  
  runs13.push_back("BTF_143_20140430-191455_beam");
  beamX13.push_back(+6.);
  beamY13.push_back(+6.);


  TGraphErrors* gr_ratio02_vs_pos = getGraph( "cef3_corr[0]/cef3_corr[2]", runs02, beamX02, beamY02 );
  TGraphErrors* gr_ratio13_vs_pos = getGraph( "cef3_corr[1]/cef3_corr[3]", runs13, beamX13, beamY13 );
  



  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  float xMin = -20.;
  float xMax = +20.;

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 2. );
  h2_axes->SetXTitle("Distance from Center [mm]");
  h2_axes->SetYTitle("Ratio - Offset");

  h2_axes->Draw();

  gr_ratio02_vs_pos->SetMarkerStyle(20);
  gr_ratio02_vs_pos->SetMarkerSize(1.6);

  gr_ratio13_vs_pos->SetMarkerStyle(24);
  gr_ratio13_vs_pos->SetMarkerSize(2.);
  gr_ratio13_vs_pos->SetMarkerColor(kBlue);



  TF1* f1_d02 = new TF1("poly_d02", "[0] + [1]*x + [2]*x*x", 0., 15. );
  f1_d02->FixParameter(0, 1.);
  f1_d02->SetLineWidth(1);
  gr_ratio02_vs_pos->Fit( f1_d02, "R" );

  TF1* f1_d13 = new TF1("poly_d13", "[0] + [1]*x + [2]*x*x", 0., 15. );
  f1_d13->FixParameter(0, 1.);
  f1_d13->SetLineWidth(1);
  f1_d13->SetLineColor(kBlue);
  gr_ratio13_vs_pos->Fit( f1_d13, "R" );

  gr_ratio02_vs_pos->Draw("p same");
  gr_ratio13_vs_pos->Draw("p same");
  f1_d02->Draw("same");
  f1_d13->Draw("same");

  f1_d13->SetParameter(1, -f1_d13->GetParameter(1));
  f1_d13->SetRange( -15., 0. );
  f1_d13->Draw("same");

  TLine* lineOne = new TLine( xMin, 1., xMax, 1. );
  lineOne->Draw("same");

  c1->SaveAs("ratio_vs_diag.eps");
  c1->SaveAs("ratio_vs_diag.png");

  return 0;

}





TGraphErrors* getGraph( const std::string varName, std::vector<std::string> runs, std::vector<float> beamX, std::vector<float> beamY ) {


  bool is_diag02 = (varName=="cef3_corr[0]/cef3_corr[2]");

  float offset=0.;

  for( unsigned i=0; i<runs.size(); ++i ) {

    if( !(beamX[i]==0. && beamY[i]==0.) ) continue;

    TFile* file = TFile::Open(Form("PosAn_%s.root", runs[i].c_str()));
    TTree* tree = (TTree*)file->Get("tree_passedEvents");

    TH1D* h1_diag = new TH1D("temp", "", 100, 0., 3.);
    tree->Project("temp", varName.c_str() );

    float r = h1_diag->GetMean();

    offset = r-1.;

  }

  std::cout << "offset: " <<  offset << std::endl;;


  TGraphErrors* gr_ratio_vs_pos = new TGraphErrors(0);

  for( unsigned i=0; i<runs.size(); ++i ) {

    TFile* file = TFile::Open(Form("PosAn_%s.root", runs[i].c_str()));
    TTree* tree = (TTree*)file->Get("tree_passedEvents");

    std::string diagName = "diag_" + runs[i];
    TH1D* h1_diag = new TH1D(diagName.c_str(), "", 100, 0., 3.);

    tree->Project(diagName.c_str(), varName.c_str() );

    TVector2 v(beamX[i],beamY[i]);
    TVector2 d = v.Rotate( -3.14159/4. );
    float diag = (is_diag02) ? d.Y() : d.X();

    float r = h1_diag->GetMean() - offset;

    //if( !is_diag02 ) {
    //  r -= 1.;
    //  r*= -1.;
    //  r += 1.;
    //}


    if( !is_diag02 ) {
      diag = -diag;
      r =  1./r;
    }

    if( r<1. ) {
      r -= 1.;
      r*= -1.;
      r += 1.;
    }
      
    gr_ratio_vs_pos->SetPoint( i, diag, r );

  }


  return gr_ratio_vs_pos;

}
