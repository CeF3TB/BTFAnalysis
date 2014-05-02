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






int main() {


  std::vector<std::string> runs;
  std::vector<float> beamX;
  std::vector<float> beamY;

  runs.push_back("BTF_94_20140430-073300_beam");
  beamX.push_back(-3.);
  beamY.push_back(+3.);

  runs.push_back("BTF_96_20140430-083733_beam");
  beamX.push_back(-6.);
  beamY.push_back(+6.);

  runs.push_back("BTF_98_20140430-092026_beam");
  beamX.push_back(-9.);
  beamY.push_back(+9.);

  runs.push_back("BTF_100_20140430-101607_beam");
  beamX.push_back(-9.);
  beamY.push_back(-9.);

  runs.push_back("BTF_118_20140430-151237_beam");
  beamX.push_back(-6.);
  beamY.push_back(-6.);

  runs.push_back("BTF_136_20140430-171004_beam");
  beamX.push_back(-3.);
  beamY.push_back(-3.);

  runs.push_back("BTF_138_20140430-175224_beam");
  beamX.push_back(+0.);
  beamY.push_back(+0.);

  runs.push_back("BTF_141_20140430-183508_beam");
  beamX.push_back(+3.);
  beamY.push_back(+3.);

  //runs.push_back("BTF_143_20140430-191455_beam");
  //beamX.push_back(+6.);
  //beamY.push_back(+6.);


  float offset02=0.;
  float offset13=0.;
  for( unsigned i=0; i<runs.size(); ++i ) {

    if( !(beamX[i]==0. && beamY[i]==0.) ) continue;

    TFile* file = TFile::Open(Form("PosAn_%s.root", runs[i].c_str()));
    TTree* tree = (TTree*)file->Get("tree_passedEvents");

    TH1D* h1_diag02 = new TH1D("temp02", "", 100, 0., 3.);
    tree->Project("temp02", "cef3_corr[0]/cef3_corr[2]");

    TH1D* h1_diag13 = new TH1D("temp13", "", 100, 0., 3.);
    tree->Project("temp13", "cef3_corr[1]/cef3_corr[3]");

    float r02 = h1_diag02->GetMean();
    float r13 = h1_diag13->GetMean();

    offset02 = r02-1.;
    offset13 = r13-1.;

  }

  std::cout << "offset02: " <<  offset02 << std::endl;;
  std::cout << "offset13: " <<  offset13 << std::endl;;


  TGraphErrors* gr_ratio02_vs_pos = new TGraphErrors(0);
  TGraphErrors* gr_ratio13_vs_pos = new TGraphErrors(0);
  TGraphErrors* gr_ratio_vs_pos = new TGraphErrors(0);

  for( unsigned i=0; i<runs.size(); ++i ) {

    TFile* file = TFile::Open(Form("PosAn_%s.root", runs[i].c_str()));
    TTree* tree = (TTree*)file->Get("tree_passedEvents");

    std::string diag02Name = "diag02_" + runs[i];
    TH1D* h1_diag02 = new TH1D(diag02Name.c_str(), "", 100, 0., 3.);

    tree->Project(diag02Name.c_str(), "cef3_corr[0]/cef3_corr[2]");

    std::string diag13Name = "diag13_" + runs[i];
    TH1D* h1_diag13 = new TH1D(diag13Name.c_str(), "", 100, 0., 3.);

    tree->Project(diag13Name.c_str(), "cef3_corr[1]/cef3_corr[3]");

    TVector2 v(beamX[i],beamY[i]);
    TVector2 d = v.Rotate( -3.14159/4. );
    float diag13 = d.X();
    float diag02 = d.Y();

    float r02 = h1_diag02->GetMean() - offset02;
    float r13 = h1_diag13->GetMean() - offset13;

    if( r02>1. ) {
      gr_ratio02_vs_pos->SetPoint( i,  diag02, r02 );
      gr_ratio_vs_pos->SetPoint( i,  diag02, r02 );
    } else {
      gr_ratio02_vs_pos->SetPoint( i, -diag02, 1./r02 );
      gr_ratio_vs_pos->SetPoint( i, -diag02, 1./r02 );
    }

    if( r13>1. ) {
      gr_ratio13_vs_pos->SetPoint( i,  diag13, r13 );
      gr_ratio_vs_pos->SetPoint( i,  diag13, r13 );
    } else {
      gr_ratio13_vs_pos->SetPoint( i, -diag13, 1./r13 );
      gr_ratio_vs_pos->SetPoint( i, -diag13, 1./r13 );
    }

  }


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 15., 10, 0., 2. );
  h2_axes->SetXTitle("Distance from Center [mm]");
  h2_axes->SetYTitle("Channel[0]/Channel[2]");

  h2_axes->Draw();

  gr_ratio02_vs_pos->SetMarkerStyle(20);
  gr_ratio02_vs_pos->SetMarkerSize(1.6);

  gr_ratio13_vs_pos->SetMarkerStyle(20);
  gr_ratio13_vs_pos->SetMarkerSize(1.6);
  gr_ratio13_vs_pos->SetMarkerColor(kRed);

  //gr_ratio_vs_pos->SetMarkerStyle(20);
  //gr_ratio_vs_pos->SetMarkerSize(1.6);
  //gr_ratio_vs_pos->SetMarkerColor(kRed);


  TF1* f1 = new TF1("poly", "[0] + [1]*x + [2]*x*x", -15., 15. );
  f1->FixParameter(0, 1.);
  
  //gr_ratio02_vs_pos->Fit( f1, "R" );
  gr_ratio_vs_pos->Fit( f1, "R" );

  gr_ratio02_vs_pos->Draw("p same");
  gr_ratio13_vs_pos->Draw("p same");
  f1->Draw("same");

  c1->SaveAs("ratio_vs_diag.eps");
  c1->SaveAs("ratio_vs_diag.png");

  return 0;

}
