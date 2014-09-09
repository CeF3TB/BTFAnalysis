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
#include "TVector2.h"

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


  std::string runName;
  if( argc>1 ) {
    std::string run_str(argv[1]);
    runName = run_str;
  }

  std::string tag = "V02";
  if( argc>2 ) {
    std::string tag_str(argv[2]);
    tag = tag_str;
  }


  DrawTools::setStyle();


  TFile* file = TFile::Open(Form("PosAnTrees_%s/%s.root", tag.c_str(), runName.c_str()));
  TTree* tree = (TTree*)file->Get("posTree");


  unsigned int run;
  tree->SetBranchAddress( "run", &run );

  bool isSingleEle_scintFront;
  tree->SetBranchAddress( "isSingleEle_scintFront", &isSingleEle_scintFront );

  float xBeam;
  tree->SetBranchAddress( "xBeam", &xBeam );
  float yBeam;
  tree->SetBranchAddress( "yBeam", &yBeam );

  float xPos_bgo_wa;
  tree->SetBranchAddress( "xPos_bgo_wa", &xPos_bgo_wa );
  float yPos_bgo_wa;
  tree->SetBranchAddress( "yPos_bgo_wa", &yPos_bgo_wa );

  float cef3_corr[4];
  tree->SetBranchAddress( "cef3_corr", cef3_corr );



  float xMax = 15.5;
  int nBins = (int)xMax*2.;

  TProfile* hp_asymm_x_vsBeam = new TProfile("asymm_x_vsBeam", "", nBins, -xMax, xMax);
  TProfile* hp_asymm_y_vsBeam = new TProfile("asymm_y_vsBeam", "", nBins, -xMax, xMax);

  TProfile* hp_asymm_log_x_vsBeam = new TProfile("asymm_log_x_vsBeam", "", nBins, -xMax, xMax);
  TProfile* hp_asymm_log_y_vsBeam = new TProfile("asymm_log_y_vsBeam", "", nBins, -xMax, xMax);

  TProfile* hp_wa_x_vsBeam = new TProfile("wa_x_vsBeam", "", nBins, -xMax, xMax);
  TProfile* hp_wa_y_vsBeam = new TProfile("wa_y_vsBeam", "", nBins, -xMax, xMax);

  TProfile* hp_wa_log_x_vsBeam = new TProfile("wa_log_x_vsBeam", "", nBins, -xMax, xMax);
  TProfile* hp_wa_log_y_vsBeam = new TProfile("wa_log_y_vsBeam", "", nBins, -xMax, xMax);



  TProfile* hp_asymm_x_vsBGO = new TProfile("asymm_x_vsBGO", "", nBins, -xMax, xMax);
  TProfile* hp_asymm_y_vsBGO = new TProfile("asymm_y_vsBGO", "", nBins, -xMax, xMax);

  TProfile* hp_asymm_log_x_vsBGO = new TProfile("asymm_log_x_vsBGO", "", nBins, -xMax, xMax);
  TProfile* hp_asymm_log_y_vsBGO = new TProfile("asymm_log_y_vsBGO", "", nBins, -xMax, xMax);

  TProfile* hp_wa_x_vsBGO = new TProfile("wa_x_vsBGO", "", nBins, -xMax, xMax);
  TProfile* hp_wa_y_vsBGO = new TProfile("wa_y_vsBGO", "", nBins, -xMax, xMax);

  TProfile* hp_wa_log_x_vsBGO = new TProfile("wa_log_x_vsBGO", "", nBins, -xMax, xMax);
  TProfile* hp_wa_log_y_vsBGO = new TProfile("wa_log_y_vsBGO", "", nBins, -xMax, xMax);


  float dMax = 35.;

  TProfile* hp_r13_d_vsBeam = new TProfile("r13_d_vsBeam", "", nBins, 0., dMax);
  TProfile* hp_r20_d_vsBeam = new TProfile("r20_d_vsBeam", "", nBins, 0., dMax);

  TProfile* hp_r13_d_vsBGO = new TProfile("r13_d_vsBGO", "", nBins, 0., dMax);
  TProfile* hp_r20_d_vsBGO = new TProfile("r20_d_vsBGO", "", nBins, 0., dMax);

  TProfile* hp_r1_d_vsBeam = new TProfile("r1_d_vsBeam", "", nBins, 0., dMax);
  TProfile* hp_r2_d_vsBeam = new TProfile("r2_d_vsBeam", "", nBins, 0., dMax);

  TProfile* hp_r1_d_vsBGO = new TProfile("r1_d_vsBGO", "", nBins, 0., dMax);
  TProfile* hp_r2_d_vsBGO = new TProfile("r2_d_vsBGO", "", nBins, 0., dMax);






  int nentries=tree->GetEntries();


  for( unsigned i=0; i<nentries; ++i ) {

    tree->GetEntry(i);

    if( !isSingleEle_scintFront ) continue;
    if( abs(xBeam)>12. || abs(yBeam)>12. ) continue;


    std::vector<float> v_cef3_corr;
    for( unsigned i=0; i<CEF3_CHANNELS; ++i ) v_cef3_corr.push_back(cef3_corr[i]);



    //  0   1  
    //  
    //  3   2


    // FIRST METHOD: LINEAR ASYMM

    float left   = cef3_corr[0]+cef3_corr[3];
    float right  = cef3_corr[1]+cef3_corr[2];
    float top    = cef3_corr[0]+cef3_corr[1];
    float bottom = cef3_corr[2]+cef3_corr[3];

    float asymm_y = (top-bottom)/(top+bottom);
    float asymm_x = (right-left)/(right+left);

    hp_asymm_x_vsBeam->Fill( xBeam, asymm_x );
    hp_asymm_y_vsBeam->Fill( yBeam, asymm_y );

    hp_asymm_x_vsBGO ->Fill( xPos_bgo_wa, asymm_x );
    hp_asymm_y_vsBGO ->Fill( yPos_bgo_wa, asymm_y );



    // SECOND METHOD: LOG ASYMM

    float left_log   = log(cef3_corr[0])+log(cef3_corr[3]);
    float right_log  = log(cef3_corr[1])+log(cef3_corr[2]);
    float top_log    = log(cef3_corr[0])+log(cef3_corr[1]);
    float bottom_log = log(cef3_corr[2])+log(cef3_corr[3]);

    float asymm_log_y = (top_log-bottom_log)/(top_log+bottom_log);
    float asymm_log_x = (right_log-left_log)/(right_log+left_log);

    hp_asymm_log_x_vsBeam->Fill( xBeam, asymm_log_x );
    hp_asymm_log_y_vsBeam->Fill( yBeam, asymm_log_y );

    hp_asymm_log_x_vsBGO ->Fill( xPos_bgo_wa, asymm_log_x );
    hp_asymm_log_y_vsBGO ->Fill( yPos_bgo_wa, asymm_log_y );


    // THIRD METHOD: WEIGHTED AVERAGE

    float wa_x = 0.;
    float wa_y = 0.;
    float sumw = 0.;

    for( unsigned i=0; i<v_cef3_corr.size(); ++i ) {
  
      float x, y;
      PositionTools::getCef3Coordinates( i, x, y );
      float w = v_cef3_corr[i];
  
      wa_x += w*x;
      wa_y += w*y;
      sumw += w;
  
    }
    
    wa_x /= sumw;
    wa_y /= sumw;

    hp_wa_x_vsBeam->Fill( xBeam, wa_x );
    hp_wa_y_vsBeam->Fill( yBeam, wa_y );

    hp_wa_x_vsBGO->Fill( xPos_bgo_wa, wa_x );
    hp_wa_y_vsBGO->Fill( yPos_bgo_wa, wa_y );


    // FOURTH METHOD: WEIGHTED AVERAGE WITH LOG WEIGHTS

    float wa_log_x = 0.;
    float wa_log_y = 0.;
    float sumw_log = 0.;

    for( unsigned i=0; i<v_cef3_corr.size(); ++i ) {
  
      float x, y;
      PositionTools::getCef3Coordinates( i, x, y );
      float w = (v_cef3_corr[i]>=1.) ? log(v_cef3_corr[i]) : 0.;
  
      wa_log_x += w*x;
      wa_log_y += w*y;
      sumw_log += w;
  
    }
    
    wa_log_x /= sumw_log;
    wa_log_y /= sumw_log;

    hp_wa_log_x_vsBeam->Fill( xBeam, wa_log_x );
    hp_wa_log_y_vsBeam->Fill( yBeam, wa_log_y );
    
    hp_wa_log_x_vsBGO->Fill( xPos_bgo_wa, wa_log_x );
    hp_wa_log_y_vsBGO->Fill( yPos_bgo_wa, wa_log_y );


    // now try ratios (study fibres 1 and 2):

    float x1, y1;
    PositionTools::getCef3Coordinates( 1, x1, y1 );
    float x2, y2;
    PositionTools::getCef3Coordinates( 2, x2, y2 );

    float r13 = cef3_corr[1]/cef3_corr[3];
    float r20 = cef3_corr[2]/cef3_corr[0];

    float eTot = cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3];
    float r1 = cef3_corr[1]/eTot;
    float r2 = cef3_corr[2]/eTot;

    TVector2 p1(x1, y1);
    TVector2 p2(x2, y2);
    TVector2 pBeam(xBeam, yBeam);
    TVector2 pBGO(xPos_bgo_wa, yPos_bgo_wa);

    TVector2 v_1Beam = p1-pBeam;
    TVector2 v_2Beam = p2-pBeam;
    TVector2 v_1BGO = p1-pBGO;
    TVector2 v_2BGO = p2-pBGO;

    float d_1Beam = v_1Beam.Mod();
    float d_2Beam = v_2Beam.Mod();
    float d_1BGO  = v_1BGO .Mod();
    float d_2BGO  = v_2BGO .Mod();

    
    hp_r13_d_vsBeam->Fill( d_1Beam, r13 );
    hp_r20_d_vsBeam->Fill( d_2Beam, r20 );

    hp_r13_d_vsBGO->Fill( d_1BGO, r13 );
    hp_r20_d_vsBGO->Fill( d_2BGO, r20 );

    hp_r1_d_vsBeam->Fill( d_1Beam, r1 );
    hp_r2_d_vsBeam->Fill( d_2Beam, r2 );

    hp_r1_d_vsBGO->Fill( d_1BGO, r1 );
    hp_r2_d_vsBGO->Fill( d_2BGO, r2 );


    
  }

  TFile* outfile = TFile::Open( Form("Cef3Positioning_%s_%s.root", runName.c_str(), tag.c_str()), "recreate" );
  outfile->cd();

  hp_asymm_x_vsBeam->Write();
  hp_asymm_y_vsBeam->Write();

  hp_asymm_log_x_vsBeam->Write();
  hp_asymm_log_y_vsBeam->Write();

  hp_wa_x_vsBeam->Write();
  hp_wa_y_vsBeam->Write();

  hp_wa_log_x_vsBeam->Write();
  hp_wa_log_y_vsBeam->Write();



  hp_asymm_x_vsBGO->Write();
  hp_asymm_y_vsBGO->Write();

  hp_asymm_log_x_vsBGO->Write();
  hp_asymm_log_y_vsBGO->Write();

  hp_wa_x_vsBGO->Write();
  hp_wa_y_vsBGO->Write();

  hp_wa_log_x_vsBGO->Write();
  hp_wa_log_y_vsBGO->Write();



  hp_r13_d_vsBeam->Write();
  hp_r20_d_vsBeam->Write();

  hp_r13_d_vsBGO->Write();
  hp_r20_d_vsBGO->Write();

  hp_r1_d_vsBeam->Write();
  hp_r2_d_vsBeam->Write();

  hp_r1_d_vsBGO->Write();
  hp_r2_d_vsBGO->Write();


  outfile->Close();

  

  std::string outputdir = "Cef3PositioningPlots_" + runName + "_" + tag;
  system(Form("mkdir -p %s", outputdir.c_str()));

  FitResults fr_asymm = drawAndGetCoeff( outputdir, "asymm", hp_asymm_x_vsBeam, hp_asymm_y_vsBeam );
  FitResults fr_asymm_log = drawAndGetCoeff( outputdir, "asymm_log", hp_asymm_log_x_vsBeam, hp_asymm_log_y_vsBeam );
  FitResults fr_wa = drawAndGetCoeff( outputdir, "wa", hp_wa_x_vsBeam, hp_wa_y_vsBeam );
  FitResults fr_wa_log = drawAndGetCoeff( outputdir, "wa_log", hp_wa_log_x_vsBeam, hp_wa_log_y_vsBeam );


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

  delete h2_axes;
  delete c1;

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
