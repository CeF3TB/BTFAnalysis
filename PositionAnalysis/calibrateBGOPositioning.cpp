#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"

#include "interface/DrawTools.h"







int main( int argc, char* argv[] ) {


  std::string tag = "V02";
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }


  DrawTools::setStyle();


  TFile* file = TFile::Open(Form("PosAnTrees_%s/crossScanFile.root", tag.c_str()));
  TTree* tree = (TTree*)file->Get("posTree");


  bool isSingleEle_scintFront;
  tree->SetBranchAddress( "isSingleEle_scintFront", &isSingleEle_scintFront );

  float xBeam;
  tree->SetBranchAddress( "xBeam", &xBeam );
  float yBeam;
  tree->SetBranchAddress( "yBeam", &yBeam );

  float bgo_corr[8];
  tree->SetBranchAddress( "bgo_corr", bgo_corr );




  TProfile* hp_asymm_x = new TProfile("asymm_x", "", 20, -10., 10.);
  TProfile* hp_asymm_y = new TProfile("asymm_y", "", 20, -10., 10.);


  int nentries=tree->GetEntries();


  for( unsigned i=0; i<nentries; ++i ) {

    tree->GetEntry(i);

    if( !isSingleEle_scintFront ) continue;


    //  0 1 2 
    //  3   4
    //  5 6 7

    float left   = bgo_corr[0]+bgo_corr[3]+bgo_corr[5];
    float right  = bgo_corr[2]+bgo_corr[4]+bgo_corr[7];
    float top    = bgo_corr[0]+bgo_corr[1]+bgo_corr[2];
    float bottom = bgo_corr[5]+bgo_corr[6]+bgo_corr[7];

    float asymm_y = (top-bottom)/(top+bottom);
    float asymm_x = (right-left)/(right+left);

    hp_asymm_x->Fill( xBeam, asymm_x );
    hp_asymm_y->Fill( yBeam, asymm_y );
    
  }

  TF1* line_x = new TF1("line_x", "[0]+[1]*x", -10., 10.);
  line_x->SetLineColor(38);
  hp_asymm_x->SetLineColor(38);
  hp_asymm_x->SetMarkerColor(38);
  hp_asymm_x->SetMarkerSize(1.6);
  hp_asymm_x->SetMarkerStyle(20);
  hp_asymm_x->Fit(line_x, "R");

  TF1* line_y = new TF1("line_y", "[0]+[1]*x", -10., 10.);
  line_y->SetLineColor(46);
  hp_asymm_y->SetLineColor(46);
  hp_asymm_y->SetMarkerColor(46);
  hp_asymm_y->SetMarkerSize(1.6);
  hp_asymm_y->SetMarkerStyle(20);
  hp_asymm_y->Fit(line_y, "R");


  std::cout << std::endl;
  std::cout << "+++  X: " << std::endl;
  std::cout << "p0:   " << line_x->GetParameter(0) << " +/- " << line_x->GetParError(0) << std::endl;
  std::cout << "p1:   " << line_x->GetParameter(1) << " +/- " << line_x->GetParError(1) << std::endl;

  std::cout << std::endl;
  std::cout << "+++  Y: " << std::endl;
  std::cout << "p0:   " << line_y->GetParameter(0) << " +/- " << line_y->GetParError(0) << std::endl;
  std::cout << "p1:   " << line_y->GetParameter(1) << " +/- " << line_y->GetParError(1) << std::endl;

  std::string outdir = "BGOPositioningPlots_" + tag;
  system(Form("mkdir -p %s", outdir.c_str()));


  TCanvas* c1 = new TCanvas("c1", "", 600., 600.); 
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, -10., 10., 10., -0.7, 0.7);
  h2_axes->SetXTitle("Beam Position [mm]");
  h2_axes->SetYTitle("BGO Asymmetry");
  h2_axes->Draw();


  TLegend* legend = new TLegend( 0.7, 0.2, 0.9, 0.4 );
  legend->SetTextSize( 0.035 );
  legend->SetFillColor(0);
  legend->AddEntry( hp_asymm_x, "X", "P" );
  legend->AddEntry( hp_asymm_y, "Y", "P" );
  //legend->Draw("same");

  hp_asymm_x->Draw("P same");
  hp_asymm_y->Draw("P same");

  TPaveText* label_x = new TPaveText( 0.2, 0.7, 0.6, 0.9, "brNDC" );
  label_x->SetTextSize( 0.044 );
  label_x->SetTextColor(38);
  label_x->SetFillColor(0);
  label_x->AddText( Form("A_{x} = %.3fx + %.3f", line_x->GetParameter(1), line_x->GetParameter(0)) );
  label_x->Draw("same");

  TPaveText* label_y = new TPaveText( 0.5, 0.3, 0.9, 0.4, "brNDC" );
  label_y->SetTextSize( 0.044 );
  label_y->SetTextColor(46);
  label_y->SetFillColor(0);
  label_y->AddText( Form("A_{y} = %.3fy + %.3f", line_y->GetParameter(1), line_y->GetParameter(0)) );
  label_y->Draw("same");


  TPaveText* labelTop = DrawTools::getLabelTop();
  labelTop->Draw("same");

  c1->SaveAs(Form("%s/asymm_vs_beam.eps", outdir.c_str()));
  c1->SaveAs(Form("%s/asymm_vs_beam.png", outdir.c_str()));
  c1->SaveAs(Form("%s/asymm_vs_beam.pdf", outdir.c_str()));


  std::string posCalibFileXName = outdir + "/posCalib_x.txt";
  ofstream ofs_x( posCalibFileXName.c_str() );
  ofs_x << line_x->GetParameter(0) << std::endl;
  ofs_x << line_x->GetParameter(1) << std::endl;
  ofs_x.close();

  std::string posCalibFileYName = outdir + "/posCalib_y.txt";
  ofstream ofs_y( posCalibFileYName.c_str() );
  ofs_y << line_y->GetParameter(0) << std::endl;
  ofs_y << line_y->GetParameter(1) << std::endl;
  ofs_y.close();

  std::cout << "-> BGO position calibration constants saved in:" << std::endl;
  std::cout << "    X: " << posCalibFileXName << std::endl;
  std::cout << "    Y: " << posCalibFileYName << std::endl;

  return 0;

}
