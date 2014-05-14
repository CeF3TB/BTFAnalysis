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
#include "TMath.h"
#include "TLegend.h"

#include "interface/DrawTools.h"
#include "interface/RunHelper.h"




int main() {


  std::string outputdir = "CaloPositioningPlots/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  DrawTools::setStyle();
 
  std::vector<std::string> runs; 
  std::vector<float> beamEnergy; 


  runs.push_back("BTF_179_20140501-002629_beam");
  runs.push_back("BTF_180_20140501-003051_beam");
  runs.push_back("BTF_181_20140501-003544_beam");
  runs.push_back("BTF_182_20140501-004020_beam");
  runs.push_back("BTF_183_20140501-004521_beam");
  runs.push_back("BTF_184_20140501-005306_beam");
  runs.push_back("BTF_185_20140501-005854_beam");
  runs.push_back("BTF_186_20140501-010447_beam");
  runs.push_back("BTF_187_20140501-011129_beam");
  runs.push_back("BTF_188_20140501-011720_beam");
  runs.push_back("BTF_189_20140501-012157_beam");
  runs.push_back("BTF_190_20140501-012655_beam");
  runs.push_back("BTF_191_20140501-013148_beam");






  TGraphErrors* gr_calibX_vs_X = new TGraphErrors(0);
  TGraphErrors* gr_calibY_vs_Y = new TGraphErrors(0);


  for( unsigned i=0; i<runs.size(); ++i ) {

    TFile* file = TFile::Open(Form("CaloPos_%s.root", runs[i].c_str()));

    TH1D* h1_calibX = (TH1D*)file->Get("cef3CalibX");
    TH1D* h1_calibY = (TH1D*)file->Get("cef3CalibY");

    float xBeam, yBeam;
    RunHelper::getBeamPosition( runs[i], xBeam, yBeam );

    gr_calibX_vs_X->SetPoint( i, xBeam, h1_calibX->GetMean() );
    gr_calibX_vs_X->SetPointError( i, 0., h1_calibX->GetMeanError() );

    gr_calibY_vs_Y->SetPoint( i, yBeam, h1_calibY->GetMean() );
    gr_calibY_vs_Y->SetPointError( i, 0., h1_calibY->GetMeanError() );

  }


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();


  TH2D* h2_axes = new TH2D( "axes", "", 10, -5., 5., 10, -0.5, 0.5 );
  h2_axes->SetXTitle("Distance from Center [mm]");
  h2_axes->SetYTitle("Calibration Constant");
  h2_axes->Draw("");

  gr_calibX_vs_X->SetMarkerStyle(20);
  gr_calibX_vs_X->SetMarkerSize(1.6);
  gr_calibX_vs_X->SetMarkerColor(46);
  gr_calibX_vs_X->Draw("p same");

  gr_calibY_vs_Y->SetMarkerStyle(20);
  gr_calibY_vs_Y->SetMarkerSize(1.6);
  gr_calibY_vs_Y->SetMarkerColor(38);
  gr_calibY_vs_Y->Draw("p same");

  TLegend* legend = new TLegend( 0.2, 0.2, 0.5, 0.3 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( gr_calibX_vs_X, "Cx vs x", "P" );
  legend->AddEntry( gr_calibY_vs_Y, "Cy vs y", "P" );
  legend->Draw("same");

  c1->SaveAs( Form( "%s/calib_vs_position.eps", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/calib_vs_position.png", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/calib_vs_position.pdf", outputdir.c_str() ) );




  return 0;

}



