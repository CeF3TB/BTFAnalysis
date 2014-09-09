#include <iostream>
#include <fstream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"

#include "interface/DrawTools.h"
#include "interface/PositionTools.h"

#include "fastDQM_CeF3_BTF.h"




struct FitResults {

  float p0_x;
  float p1_x;
  float p0_y;
  float p1_y;

};



FitResults drawAndGetCoeff( const std::string& outputdir, const std::string& name, TProfile* hp_x, TProfile* hp_y );
float sumVector(std::vector<float> v);
std::vector<TH1D*> getPerformanceVector( const std::string& name );



int main( int argc, char* argv[] ) {


  std::string tag = "V02";
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }


  DrawTools::setStyle();


  TFile* file = TFile::Open(Form("PosAnTrees_%s/crossScanFile.root", tag.c_str()));
  TTree* tree = (TTree*)file->Get("posTree");


  unsigned int run;
  tree->SetBranchAddress( "run", &run );
  unsigned int event;
  tree->SetBranchAddress( "event", &event );

  bool isSingleEle_scintFront;
  tree->SetBranchAddress( "isSingleEle_scintFront", &isSingleEle_scintFront );
  bool bgo_corr_ok;
  tree->SetBranchAddress( "bgo_corr_ok", &bgo_corr_ok );

  float xBeam;
  tree->SetBranchAddress( "xBeam", &xBeam );
  float yBeam;
  tree->SetBranchAddress( "yBeam", &yBeam );

  float bgo_corr[8];
  tree->SetBranchAddress( "bgo_corr", bgo_corr );
  float cef3_corr[8];
  tree->SetBranchAddress( "cef3_corr", cef3_corr );




  TProfile* hp_asymm_x = new TProfile("asymm_x", "", 30, -15., 15.);
  TProfile* hp_asymm_y = new TProfile("asymm_y", "", 30, -15., 15.);

  TProfile* hp_asymm_log_x = new TProfile("asymm_log_x", "", 30, -15., 15.);
  TProfile* hp_asymm_log_y = new TProfile("asymm_log_y", "", 30, -15., 15.);

  TProfile* hp_wa_x = new TProfile("wa_x", "", 30, -15., 15.);
  TProfile* hp_wa_y = new TProfile("wa_y", "", 30, -15., 15.);

  TProfile* hp_wa_log_x = new TProfile("wa_log_x", "", 30, -15., 15.);
  TProfile* hp_wa_log_y = new TProfile("wa_log_y", "", 30, -15., 15.);

  TProfile* hp_wa_calo_x = new TProfile("wa_calo_x", "", 30, -15., 15.);
  TProfile* hp_wa_calo_y = new TProfile("wa_calo_y", "", 30, -15., 15.);

  TProfile* hp_wa_calo_log_x = new TProfile("wa_calo_log_x", "", 30, -15., 15.);
  TProfile* hp_wa_calo_log_y = new TProfile("wa_calo_log_y", "", 30, -15., 15.);




  std::vector<float> xbgo, ybgo;
  for( unsigned i=0; i<BGO_CHANNELS; ++i ) {
    float x,y;
    PositionTools::getBGOCoordinates( i, x, y );
    xbgo.push_back( x );
    ybgo.push_back( y );
  }

  int nentries=tree->GetEntries();


  for( unsigned i=0; i<nentries; ++i ) {

    tree->GetEntry(i);

    if( !isSingleEle_scintFront ) continue;
    //if( !bgo_corr_ok ) continue;


    std::vector<float> v_bgo_corr;
    for( unsigned i=0; i<BGO_CHANNELS; ++i ) v_bgo_corr.push_back(bgo_corr[i]);



    //  0 1 2 
    //  3   4
    //  5 6 7


    // FIRST METHOD: LINEAR ASYMM

    float left   = bgo_corr[0]+bgo_corr[3]+bgo_corr[5];
    float right  = bgo_corr[2]+bgo_corr[4]+bgo_corr[7];
    float top    = bgo_corr[0]+bgo_corr[1]+bgo_corr[2];
    float bottom = bgo_corr[5]+bgo_corr[6]+bgo_corr[7];

    float asymm_y = ((top+bottom)>0.) ? (top-bottom)/(top+bottom) : 0.;
    float asymm_x = ((right+left)>0.) ? (right-left)/(right+left) : 0.;

    hp_asymm_x->Fill( xBeam, asymm_x );
    hp_asymm_y->Fill( yBeam, asymm_y );


    // SECOND METHOD: LOG ASYMM

    float logbgo[8];
    logbgo[0] = (bgo_corr[0]>=1.) ? log(bgo_corr[0]) : 0.;
    logbgo[1] = (bgo_corr[1]>=1.) ? log(bgo_corr[1]) : 0.;
    logbgo[2] = (bgo_corr[2]>=1.) ? log(bgo_corr[2]) : 0.;
    logbgo[3] = (bgo_corr[3]>=1.) ? log(bgo_corr[3]) : 0.;
    logbgo[4] = (bgo_corr[4]>=1.) ? log(bgo_corr[4]) : 0.;
    logbgo[5] = (bgo_corr[5]>=1.) ? log(bgo_corr[5]) : 0.;
    logbgo[6] = (bgo_corr[6]>=1.) ? log(bgo_corr[6]) : 0.;
    logbgo[7] = (bgo_corr[7]>=1.) ? log(bgo_corr[7]) : 0.;

    float left_log   = logbgo[0]+logbgo[3]+logbgo[5];
    float right_log  = logbgo[2]+logbgo[4]+logbgo[7];
    float top_log    = logbgo[0]+logbgo[1]+logbgo[2];
    float bottom_log = logbgo[5]+logbgo[6]+logbgo[7];

    float asymm_log_y = ((top_log+bottom_log)>0.) ? (top_log-bottom_log)/(top_log+bottom_log) : 0.;
    float asymm_log_x = ((right_log+left_log)>0.) ? (right_log-left_log)/(right_log+left_log) : 0.;

    hp_asymm_log_x->Fill( xBeam, asymm_log_x );
    hp_asymm_log_y->Fill( yBeam, asymm_log_y );


    // THIRD METHOD: WEIGHTED AVERAGE

    float wa_x = 0.;
    float wa_y = 0.;
    float sumw = 0.;

    for( unsigned i=0; i<v_bgo_corr.size(); ++i ) {
  
      float xbgo, ybgo;
      PositionTools::getBGOCoordinates( i, xbgo, ybgo );
  
      float w = v_bgo_corr[i];
  
      wa_x += w*xbgo;
      wa_y += w*ybgo;
      sumw += w;
  
    }
    
    wa_x /= sumw;
    wa_y /= sumw;

    hp_wa_x->Fill( xBeam, wa_x );
    hp_wa_y->Fill( yBeam, wa_y );


    // FOURTH METHOD: WEIGHTED AVERAGE WITH LOG WEIGHTS

    float wa_log_x = 0.;
    float wa_log_y = 0.;
    float sumw_log = 0.;

    for( unsigned i=0; i<v_bgo_corr.size(); ++i ) {
  
      float xbgo, ybgo;
      PositionTools::getBGOCoordinates( i, xbgo, ybgo );
  
      float w = logbgo[i];
  
      wa_log_x += w*xbgo;
      wa_log_y += w*ybgo;
      sumw_log += w;
  
    }
    
    wa_log_x /= sumw_log;
    wa_log_y /= sumw_log;

    hp_wa_log_x->Fill( xBeam, wa_log_x );
    hp_wa_log_y->Fill( yBeam, wa_log_y );


    // FIFTH METHOD: WEIGHTED AVERAGE ADDING ALSO CEF3

    float etot_cef3 = cef3_corr[0] + cef3_corr[1] + cef3_corr[2] + cef3_corr[3];
    float etot_cef3_calib = 0.79*etot_cef3;
    float sumw_calo = sumw + etot_cef3_calib;

    float wa_calo_x = wa_x * sumw / sumw_calo;
    float wa_calo_y = wa_y * sumw / sumw_calo;

    hp_wa_calo_x->Fill( xBeam, wa_calo_x );
    hp_wa_calo_y->Fill( yBeam, wa_calo_y );




    // SIXTH METHOD: LOG-WEIGHTED AVERAGE ADDING ALSO CEF3


    float wacalolog_x=0.;
    float wacalolog_y=0.;
    float sumwcalolog = 0.;


    for( unsigned i=0; i<v_bgo_corr.size(); ++i ) {

      float xbgo, ybgo;
      PositionTools::getBGOCoordinates( i, xbgo, ybgo );

      float w = (v_bgo_corr[i]>=1.) ? log(v_bgo_corr[i]) : 0.;

      wacalolog_x += w*xbgo;
      wacalolog_y += w*ybgo;
      sumwcalolog += w;

    }

    float etot_cef3_log = (etot_cef3_calib>=1.) ? log(etot_cef3_calib) : 0.;
    
    sumwcalolog += etot_cef3_log;

    wacalolog_x /= sumwcalolog;
    wacalolog_y /= sumwcalolog;


    hp_wa_calo_log_x->Fill( xBeam, wacalolog_x );
    hp_wa_calo_log_y->Fill( yBeam, wacalolog_y );

    
  }

  std::string outputdir = "BGOPositioningPlots_" + tag;
  system(Form("mkdir -p %s", outputdir.c_str()));

  FitResults fr_asymm = drawAndGetCoeff( outputdir, "asymm", hp_asymm_x, hp_asymm_y );
  FitResults fr_asymm_log = drawAndGetCoeff( outputdir, "asymm_log", hp_asymm_log_x, hp_asymm_log_y );
  FitResults fr_wa = drawAndGetCoeff( outputdir, "wa", hp_wa_x, hp_wa_y );
  FitResults fr_wa_log = drawAndGetCoeff( outputdir, "wa_log", hp_wa_log_x, hp_wa_log_y );
  FitResults fr_wa_calo = drawAndGetCoeff( outputdir, "wa_calo", hp_wa_calo_x, hp_wa_calo_y );
  FitResults fr_wa_calo_log = drawAndGetCoeff( outputdir, "wa_calo_log", hp_wa_calo_log_x, hp_wa_calo_log_y );



  return 0;

}





FitResults drawAndGetCoeff( const std::string& outputdir, const std::string& name, TProfile* hp_x, TProfile* hp_y ) {


  TF1* line_x = new TF1(Form("line_%s_x", name.c_str()), "[0]+[1]*x", -10., 10.);
  line_x->SetLineColor(38);
  hp_x->SetLineColor(38);
  hp_x->SetMarkerColor(38);
  hp_x->SetMarkerSize(1.6);
  hp_x->SetMarkerStyle(20);
  hp_x->Fit(line_x, "R");

  TF1* line_y = new TF1(Form("line_%s_y", name.c_str()), "[0]+[1]*x", -10., 10.);
  line_y->SetLineColor(46);
  hp_y->SetLineColor(46);
  hp_y->SetMarkerColor(46);
  hp_y->SetMarkerSize(1.6);
  hp_y->SetMarkerStyle(20);
  hp_y->Fit(line_y, "R");


  //std::cout << std::endl;
  //std::cout << "+++  X: " << std::endl;
  //std::cout << "p0:   " << line_x->GetParameter(0) << " +/- " << line_x->GetParError(0) << std::endl;
  //std::cout << "p1:   " << line_x->GetParameter(1) << " +/- " << line_x->GetParError(1) << std::endl;

  //std::cout << std::endl;
  //std::cout << "+++  Y: " << std::endl;
  //std::cout << "p0:   " << line_y->GetParameter(0) << " +/- " << line_y->GetParError(0) << std::endl;
  //std::cout << "p1:   " << line_y->GetParameter(1) << " +/- " << line_y->GetParError(1) << std::endl;



  TCanvas* c1 = new TCanvas("c2", "", 600., 600.); 
  c1->cd();

  float yMax = hp_x->GetMaximum()*1.2;

  TH2D* h2_axes = new TH2D("axes", "", 10, -15., 15., 10, -yMax, yMax );
  h2_axes->SetXTitle("Beam Position [mm]");
  //h2_axes->SetYTitle("BGO Asymmetry");
  h2_axes->SetYTitle(name.c_str());
  h2_axes->Draw();


  TLegend* legend = new TLegend( 0.7, 0.2, 0.9, 0.4 );
  legend->SetTextSize( 0.035 );
  legend->SetFillColor(0);
  legend->AddEntry( hp_x, "X", "P" );
  legend->AddEntry( hp_y, "Y", "P" );
  //legend->Draw("same");

  hp_x->Draw("P same");
  hp_y->Draw("P same");

  TPaveText* label_x = new TPaveText( 0.2, 0.7, 0.6, 0.9, "brNDC" );
  label_x->SetTextSize( 0.044 );
  label_x->SetTextColor(38);
  label_x->SetFillColor(0);
  label_x->AddText( Form("%.3fx + %.3f", line_x->GetParameter(1), line_x->GetParameter(0)) );
  label_x->Draw("same");

  TPaveText* label_y = new TPaveText( 0.5, 0.3, 0.9, 0.4, "brNDC" );
  label_y->SetTextSize( 0.044 );
  label_y->SetTextColor(46);
  label_y->SetFillColor(0);
  label_y->AddText( Form("%.3fy + %.3f", line_y->GetParameter(1), line_y->GetParameter(0)) );
  label_y->Draw("same");


  TPaveText* labelTop = DrawTools::getLabelTop();
  labelTop->Draw("same");

  c1->SaveAs(Form("%s/%s_vs_beam.eps", outputdir.c_str(), name.c_str()));
  c1->SaveAs(Form("%s/%s_vs_beam.png", outputdir.c_str(), name.c_str()));
  c1->SaveAs(Form("%s/%s_vs_beam.pdf", outputdir.c_str(), name.c_str()));


  std::string posCalibFileXName = outputdir + "/posCalib_" + name + "_x.txt";
  ofstream ofs_x( posCalibFileXName.c_str() );
  ofs_x << line_x->GetParameter(0) << std::endl;
  ofs_x << line_x->GetParameter(1) << std::endl;
  ofs_x.close();

  std::string posCalibFileYName = outputdir + "/posCalib_" + name + "_y.txt";
  ofstream ofs_y( posCalibFileYName.c_str() );
  ofs_y << line_y->GetParameter(0) << std::endl;
  ofs_y << line_y->GetParameter(1) << std::endl;
  ofs_y.close();

  std::cout << "-> BGO position calibration constants saved in:" << std::endl;
  std::cout << "    X: " << posCalibFileXName << std::endl;
  std::cout << "    Y: " << posCalibFileYName << std::endl;

  FitResults fr;
  fr.p0_x = line_x->GetParameter(0);
  fr.p1_x = line_x->GetParameter(1);
  fr.p0_y = line_y->GetParameter(0);
  fr.p1_y = line_y->GetParameter(1);

  delete c1;
  delete h2_axes;

  return fr;

}


float sumVector(std::vector<float> v) {

  float s=0.;

  for( unsigned i=0; i<v.size(); ++i ) s+=v[i];

  return s;

}



std::vector<TH1D*> getPerformanceVector( const std::string& name ) {

  std::vector<TH1D*> v_h1;

  for( unsigned i=-10; i<10; ++i ) {

    TH1D* h1 = new TH1D( Form("%s_%d", name.c_str(), i+10), "", 50, -40., 40. );
    v_h1.push_back(h1);

  }


  return v_h1;

}
