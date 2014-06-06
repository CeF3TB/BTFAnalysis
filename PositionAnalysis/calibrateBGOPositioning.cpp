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

  TProfile* hp_asymm_x_log = new TProfile("asymm_x_log", "", 20, -10., 10.);
  TProfile* hp_asymm_y_log = new TProfile("asymm_y_log", "", 20, -10., 10.);

  TProfile* hp_wa_x = new TProfile("wa_x", "", 20, -10., 10.);
  TProfile* hp_wa_y = new TProfile("wa_y", "", 20, -10., 10.);

  TProfile* hp_wa_x_log = new TProfile("wa_x_log", "", 20, -10., 10.);
  TProfile* hp_wa_y_log = new TProfile("wa_y_log", "", 20, -10., 10.);



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

    float asymm_y = (top-bottom)/(top+bottom);
    float asymm_x = (right-left)/(right+left);

    hp_asymm_x->Fill( xBeam, asymm_x );
    hp_asymm_y->Fill( yBeam, asymm_y );


    // SECOND METHOD: LOG ASYMM

    float left_log   = log(bgo_corr[0])+log(bgo_corr[3])+log(bgo_corr[5]);
    float right_log  = log(bgo_corr[2])+log(bgo_corr[4])+log(bgo_corr[7]);
    float top_log    = log(bgo_corr[0])+log(bgo_corr[1])+log(bgo_corr[2]);
    float bottom_log = log(bgo_corr[5])+log(bgo_corr[6])+log(bgo_corr[7]);

    float asymm_y_log = (top_log-bottom_log)/(top_log+bottom_log);
    float asymm_x_log = (right_log-left_log)/(right_log+left_log);

    hp_asymm_x_log->Fill( xBeam, asymm_x_log );
    hp_asymm_y_log->Fill( yBeam, asymm_y_log );


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

    float wa_x_log = 0.;
    float wa_y_log = 0.;
    float sumw_log = 0.;

    for( unsigned i=0; i<v_bgo_corr.size(); ++i ) {
  
      float xbgo, ybgo;
      PositionTools::getBGOCoordinates( i, xbgo, ybgo );
  
      float w = log(v_bgo_corr[i]);
  
      wa_x_log += w*xbgo;
      wa_y_log += w*ybgo;
      sumw_log += w;
  
    }
    
    wa_x_log /= sumw_log;
    wa_y_log /= sumw_log;

    hp_wa_x_log->Fill( xBeam, wa_x_log );
    hp_wa_y_log->Fill( yBeam, wa_y_log );
    
  }

  std::string outputdir = "BGOPositioningPlots_" + tag;
  system(Form("mkdir -p %s", outputdir.c_str()));

  FitResults fr_asymm = drawAndGetCoeff( outputdir, "asymm", hp_asymm_x, hp_asymm_y );
  FitResults fr_asymm_log = drawAndGetCoeff( outputdir, "asymm_log", hp_asymm_x_log, hp_asymm_y_log );
  FitResults fr_wa = drawAndGetCoeff( outputdir, "wa", hp_wa_x, hp_wa_y );
  FitResults fr_wa_log = drawAndGetCoeff( outputdir, "wa_log", hp_wa_x_log, hp_wa_y_log );




/*
  // second round to test the method performances

  std::vector<TH1D*> vh1_perfx_asymm = getPerformanceVector("asymm_x");
  std::vector<TH1D*> vh1_perfx_asymm_log = getPerformanceVector("asymm_x_log");
  std::vector<TH1D*> vh1_perfx_wa = getPerformanceVector("wa_x");
  std::vector<TH1D*> vh1_perfx_wa_log = getPerformanceVector("wa_x_log");

  std::vector<TH1D*> vh1_perfy_asymm = getPerformanceVector("asymm_y");
  std::vector<TH1D*> vh1_perfy_asymm_log = getPerformanceVector("asymm_y_log");
  std::vector<TH1D*> vh1_perfy_wa = getPerformanceVector("wa_y");
  std::vector<TH1D*> vh1_perfy_wa_log = getPerformanceVector("wa_y_log");
  



  for( unsigned i=0; i<nentries; ++i ) {

    tree->GetEntry(i);

    if( !isSingleEle_scintFront ) continue;


    std::vector<float> v_bgo_corr;
    for( unsigned i=0; i<BGO_CHANNELS; ++i ) v_bgo_corr.push_back(bgo_corr[i]);


    int iBin_x = 

    // FIRST METHOD: LINEAR ASYMM

    float left   = bgo_corr[0]+bgo_corr[3]+bgo_corr[5];
    float right  = bgo_corr[2]+bgo_corr[4]+bgo_corr[7];
    float top    = bgo_corr[0]+bgo_corr[1]+bgo_corr[2];
    float bottom = bgo_corr[5]+bgo_corr[6]+bgo_corr[7];

    float asymm_y = (top-bottom)/(top+bottom);
    float asymm_x = (right-left)/(right+left);

    float xPos_asymm = (asymm_x-fr_asymm.p0x)/fr_asymm.p1x;
    float yPos_asymm = (asymm_y-fr_asymm.p0y)/fr_asymm.p1y;

    


    // SECOND METHOD: LOG ASYMM

    float left_log   = log(bgo_corr[0])+log(bgo_corr[3])+log(bgo_corr[5]);
    float right_log  = log(bgo_corr[2])+log(bgo_corr[4])+log(bgo_corr[7]);
    float top_log    = log(bgo_corr[0])+log(bgo_corr[1])+log(bgo_corr[2]);
    float bottom_log = log(bgo_corr[5])+log(bgo_corr[6])+log(bgo_corr[7]);

    float asymm_y_log = (top_log-bottom_log)/(top_log+bottom_log);
    float asymm_x_log = (right_log-left_log)/(right_log+left_log);

    hp_asymm_x_log->Fill( xBeam, asymm_x_log );
    hp_asymm_y_log->Fill( yBeam, asymm_y_log );


    // THIRD METHOD: WEIGHTED AVERAGE

    std::vector<float> xPosW_bgo;
    std::vector<float> yPosW_bgo;

    xPosW_bgo.push_back(v_bgo_corr[0]*xbgo[0]);
    xPosW_bgo.push_back(v_bgo_corr[1]*xbgo[1]);
    xPosW_bgo.push_back(v_bgo_corr[2]*xbgo[2]);
    xPosW_bgo.push_back(v_bgo_corr[3]*xbgo[3]);
    xPosW_bgo.push_back(v_bgo_corr[4]*xbgo[4]);
    xPosW_bgo.push_back(v_bgo_corr[5]*xbgo[5]);
    xPosW_bgo.push_back(v_bgo_corr[6]*xbgo[6]);
    xPosW_bgo.push_back(v_bgo_corr[7]*xbgo[7]);
    
    yPosW_bgo.push_back(v_bgo_corr[0]*ybgo[0]);
    yPosW_bgo.push_back(v_bgo_corr[1]*ybgo[1]);
    yPosW_bgo.push_back(v_bgo_corr[2]*ybgo[2]);
    yPosW_bgo.push_back(v_bgo_corr[3]*ybgo[3]);
    yPosW_bgo.push_back(v_bgo_corr[4]*ybgo[4]);
    yPosW_bgo.push_back(v_bgo_corr[5]*ybgo[5]);
    yPosW_bgo.push_back(v_bgo_corr[6]*ybgo[6]);
    yPosW_bgo.push_back(v_bgo_corr[7]*ybgo[7]);

    float wa_x = sumVector( xPosW_bgo )/sumVector( v_bgo_corr );
    float wa_y = sumVector( yPosW_bgo )/sumVector( v_bgo_corr );
    
    hp_wa_x->Fill( xBeam, wa_x );
    hp_wa_y->Fill( yBeam, wa_y );


    // FOURTH METHOD: WEIGHTED AVERAGE WITH LOG WEIGHTS

    std::vector<float> xPosW_bgo_log;
    std::vector<float> yPosW_bgo_log;

    xPosW_bgo_log.push_back(log(v_bgo_corr[0])*xbgo[0]);
    xPosW_bgo_log.push_back(log(v_bgo_corr[1])*xbgo[1]);
    xPosW_bgo_log.push_back(log(v_bgo_corr[2])*xbgo[2]);
    xPosW_bgo_log.push_back(log(v_bgo_corr[3])*xbgo[3]);
    xPosW_bgo_log.push_back(log(v_bgo_corr[4])*xbgo[4]);
    xPosW_bgo_log.push_back(log(v_bgo_corr[5])*xbgo[5]);
    xPosW_bgo_log.push_back(log(v_bgo_corr[6])*xbgo[6]);
    xPosW_bgo_log.push_back(log(v_bgo_corr[7])*xbgo[7]);
    
    yPosW_bgo_log.push_back(log(v_bgo_corr[0])*ybgo[0]);
    yPosW_bgo_log.push_back(log(v_bgo_corr[1])*ybgo[1]);
    yPosW_bgo_log.push_back(log(v_bgo_corr[2])*ybgo[2]);
    yPosW_bgo_log.push_back(log(v_bgo_corr[3])*ybgo[3]);
    yPosW_bgo_log.push_back(log(v_bgo_corr[4])*ybgo[4]);
    yPosW_bgo_log.push_back(log(v_bgo_corr[5])*ybgo[5]);
    yPosW_bgo_log.push_back(log(v_bgo_corr[6])*ybgo[6]);
    yPosW_bgo_log.push_back(log(v_bgo_corr[7])*ybgo[7]);

    float wa_x_log = sumVector( xPosW_bgo_log )/sumVector( v_bgo_corr );
    float wa_y_log = sumVector( yPosW_bgo_log )/sumVector( v_bgo_corr );
    


  }  
*/

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



  TCanvas* c1 = new TCanvas("c1", "", 600., 600.); 
  c1->cd();

  float yMax = hp_x->GetMaximum()*1.2;

  TH2D* h2_axes = new TH2D("axes", "", 10, -10., 10., 10., -yMax, yMax );
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
