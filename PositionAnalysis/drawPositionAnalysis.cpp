#include <iostream>
#include <string>
#include <stdlib.h>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TEllipse.h"
#include "TString.h"

#include "interface/DrawTools.h"
#include "interface/RunHelper.h"

#include "fastDQM_CeF3_BTF.h"




std::string runName;


TGraphErrors* get_xyCenter( TH2D* h2_xyPos );
void drawSinglePositionPlot( const std::string& outputdir, TFile* file, const std::string& runName, const std::string& suffix );
void fitHodoWithBeam( const std::string& outputdir, const std::string& suffix, TH1D* h1, float r, float& pos, float& pos_err );
void drawPositionResolutionXY( const std::string& outputdir, TFile* file, const std::string& runName, const std::string& name1, const std::string& name2="Beam" );



Double_t DoubleGauss(Double_t *x, Double_t *par)
{
  return par[0]*exp(-0.5*(  ((x[0]-par[1])*(x[0]-par[1]))/(par[3]*par[3]) + (x[1]-par[2])*(x[1]-par[2]))/(par[3]*par[3])  );
}



int main( int argc, char* argv[] ) {

  
  runName = "test_10";
  if( argc>1 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
  }

  std::string tag = "V02";
  if( argc>2 ) {
    std::string tag_str(argv[2]);
    tag = tag_str;
  }


  std::string fileName = "PosAnTrees_" + tag + "/PosAn_" + runName + ".root";
  TFile* file = TFile::Open( fileName.c_str() );
  std::cout << "-> Opened file: " << fileName << std::endl;
  
  std::string outputdir = "Plots_" + runName + "_" + tag;
  std::string mkdir_command = "mkdir -p " + outputdir;
  system(mkdir_command.c_str());


  TStyle* style = DrawTools::setStyle();
  style->cd();

  drawSinglePositionPlot( outputdir, file, runName, "" );
  drawSinglePositionPlot( outputdir, file, runName, "_singleEle" );


  drawPositionResolutionXY( outputdir, file, runName, "bgo", "Beam" );

  return 0;

}






void drawSinglePositionPlot( const std::string& outputdir, TFile* file, const std::string& runName, const std::string& suffix ) {

  // manually set beam nominal position for some known runs:
  float beamX = -999.;
  float beamY = -999.;
  float beamRX = 4.;
  float beamRY = 2.;
  RunHelper::getBeamPosition( runName, beamX, beamY );


  bool drawBeam = ((beamX>-999.) && (beamY>-999.));
  //bool beamInsideHodo = ((fabs(beamX)<4.) && (fabs(beamY)<4.));


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  //TH2D* h2_xyPos      = (TH2D*)file->Get(Form("xyPos%s_regrBDTG", suffix.c_str())); 
  TH2D* h2_xyPos      = (TH2D*)file->Get(Form("xyPos_new%s", suffix.c_str())); 
  TH2D* h2_xyPos_hodo = (TH2D*)file->Get(Form("xyPos%s_hodo", suffix.c_str())); 
  //TH2D* h2_xyPos_hodo = (TH2D*)file->Get(Form("xyPos%s_hodo", suffix.c_str())); 
  TH2D* h2_xyPos_bgo  = (TH2D*)file->Get(Form("xyPos%s_bgo", suffix.c_str())); 
  TH2D* h2_xyPos_calo  = (TH2D*)file->Get(Form("xyPos%s_calo", suffix.c_str())); 

  float xySize = 25.;
  float xMax = xySize*3./2.;
  TH2D* h2_axes = new TH2D("axes", "", 10, -xMax, xMax, 10, -xMax, xMax);
  h2_axes->SetXTitle("X Position [mm]");
  h2_axes->SetYTitle("Y Position [mm]");

  h2_axes->Draw();


  TGraphErrors* gr_xyCenter      = get_xyCenter( h2_xyPos );
  TGraphErrors* gr_xyCenter_hodo = get_xyCenter( h2_xyPos_hodo );
  TGraphErrors* gr_xyCenter_bgo  = get_xyCenter( h2_xyPos_bgo  );
  TGraphErrors* gr_xyCenter_calo  = get_xyCenter( h2_xyPos_calo  );

  gr_xyCenter->SetMarkerColor(kRed+2);
  gr_xyCenter->SetLineColor(kRed+2);
  gr_xyCenter->SetMarkerStyle(20);
  gr_xyCenter->SetMarkerSize(1.6);

  gr_xyCenter_hodo->SetMarkerColor(kBlack);
  gr_xyCenter_hodo->SetLineColor(kBlack);
  gr_xyCenter_hodo->SetMarkerStyle(20);
  gr_xyCenter_hodo->SetMarkerSize(1.6);

  gr_xyCenter_bgo->SetMarkerColor(kGreen+3);
  gr_xyCenter_bgo->SetLineColor(kGreen+3);
  gr_xyCenter_bgo->SetMarkerStyle(20);
  gr_xyCenter_bgo->SetMarkerSize(1.6);

  gr_xyCenter_calo->SetMarkerColor(kBlue);
  gr_xyCenter_calo->SetLineColor(kBlue);
  gr_xyCenter_calo->SetMarkerStyle(20);
  gr_xyCenter_calo->SetMarkerSize(1.6);



  // fit hodo points with the expected beam size
  float xPos_hodo_fit, xPos_hodo_fit_err;
  float yPos_hodo_fit, yPos_hodo_fit_err;
  fitHodoWithBeam( outputdir, "X"+suffix, h2_xyPos_hodo->ProjectionX(), beamRX, xPos_hodo_fit, xPos_hodo_fit_err );
  fitHodoWithBeam( outputdir, "Y"+suffix, h2_xyPos_hodo->ProjectionY(), beamRY, yPos_hodo_fit, yPos_hodo_fit_err );


  TGraphErrors* gr_xyCenter_hodo_fit = new TGraphErrors(0);
  gr_xyCenter_hodo_fit->SetPoint(0, xPos_hodo_fit, yPos_hodo_fit);
  gr_xyCenter_hodo_fit->SetPointError(0, xPos_hodo_fit_err, yPos_hodo_fit_err);
  gr_xyCenter_hodo_fit->SetMarkerStyle(20);
  gr_xyCenter_hodo_fit->SetMarkerSize(1.6);





  int lineColor = 17;

  //TLine* line_x1 = new TLine( -xMax, -xySize/2., +xMax, -xySize/2. );
  //line_x1->SetLineColor(lineColor);
  //line_x1->Draw("same");

  //TLine* line_x2 = new TLine( -xMax, +xySize/2., +xMax, +xySize/2. );
  //line_x2->SetLineColor(lineColor);
  //line_x2->Draw("same");

  //TLine* line_y1 = new TLine( -xySize/2., -xMax, -xySize/2., +xMax );
  //line_y1->SetLineColor(lineColor);
  //line_y1->Draw("same");

  //TLine* line_y2 = new TLine( +xySize/2., -xMax, +xySize/2., +xMax );
  //line_y2->SetLineColor(lineColor);
  //line_y2->Draw("same");

  TLegend* legend = new TLegend( 0.75, 0.21, 0.9, 0.39 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( gr_xyCenter_hodo_fit, "Hodo", "P" );
  //legend->AddEntry( gr_xyCenter, "CeF3", "P" );
  legend->AddEntry( gr_xyCenter_bgo, "BGO", "P" );
  //legend->AddEntry( gr_xyCenter_calo, "Calo", "P" );
  legend->Draw("same");

  float bgoFrontSize = 22.;
  std::vector<TBox*> b_bgo;

  for( unsigned i=0; i<BGO_CHANNELS; ++i ) {

    float x,y;
    RunHelper::getBGOCoordinates( i, x, y );
    TBox* b_bgo0 = new TBox( x-bgoFrontSize/2., y-bgoFrontSize/2., x+bgoFrontSize/2., y+bgoFrontSize/2. );
    b_bgo0->SetFillColor(0);
    b_bgo0->SetFillStyle(0);
    b_bgo0->SetLineColor(lineColor);
    b_bgo0->SetLineWidth(1.);
    b_bgo0->Draw("L same");

    b_bgo.push_back(b_bgo0);
 
  }


  TGraph* gr_beamPos = new TGraph(0);
  gr_beamPos->SetMarkerStyle(24);
  gr_beamPos->SetMarkerSize(4);
  TLegend* legend2 = new TLegend( 0.18, 0.225, 0.5, 0.255 );
  legend2->SetFillColor(0);
  legend2->SetFillStyle(0);
  legend2->SetLineColor(0);
  legend2->SetLineWidth(0);
  legend2->SetTextSize(0.035);
  legend2->AddEntry( gr_beamPos, "Beam Position", "P" );
  if( drawBeam ) 
    legend2->Draw("same");



  TPaveText* label_top = DrawTools::getLabelTop();
  TPaveText* label_run = DrawTools::getLabelRun(runName);
  label_top->Draw("same");
  label_run->SetFillStyle(0);
  label_run->SetLineColor(0);
  label_run->SetLineWidth(0);
  label_run->Draw("same");




  h2_xyPos_hodo->SetMarkerColor(14);
  h2_xyPos->SetMarkerColor(46);
  h2_xyPos_bgo->SetMarkerColor(30);
  h2_xyPos_calo->SetMarkerColor(38);

  float hodoSize = 8.;
  TLine* lineHodo_x1 = new TLine( -hodoSize/2., -hodoSize/2., +hodoSize/2., -hodoSize/2. );
  lineHodo_x1->SetLineColor(kBlack);
  lineHodo_x1->SetLineStyle(2);
  lineHodo_x1->Draw("same");

  TLine* lineHodo_x2 = new TLine( -hodoSize/2., +hodoSize/2., +hodoSize/2., +hodoSize/2. );
  lineHodo_x2->SetLineColor(kBlack);
  lineHodo_x2->SetLineStyle(2);
  lineHodo_x2->Draw("same");

  TLine* lineHodo_y1 = new TLine( -hodoSize/2., -hodoSize/2., -hodoSize/2., +hodoSize/2. );
  lineHodo_y1->SetLineColor(kBlack);
  lineHodo_y1->SetLineStyle(2);
  lineHodo_y1->Draw("same");

  TLine* lineHodo_y2 = new TLine( +hodoSize/2., -hodoSize/2., +hodoSize/2., +hodoSize/2. );
  lineHodo_y2->SetLineColor(kBlack);
  lineHodo_y2->SetLineStyle(2);
  lineHodo_y2->Draw("same");


  h2_xyPos_bgo->Draw("same");
  //h2_xyPos->Draw("same");
  h2_xyPos_hodo->Draw("same");
  //h2_xyPos_calo->Draw("same");


  TEllipse* beamPos = new TEllipse( beamX, beamY, beamRX, beamRY );
  beamPos->SetLineColor(kBlack);
  beamPos->SetFillStyle(0);
  beamPos->Draw("same"); 

  //gr_xyCenter_hodo->Draw("p same");  // now using hodo_fit
  gr_xyCenter_hodo_fit->Draw("p same");
  gr_xyCenter_bgo->Draw("p same");
  //gr_xyCenter->Draw("p same"); // don't draw for now
  //gr_xyCenter_calo->Draw("p same");
  //gr_xyPos_fit->Draw("p same");


  c1->SaveAs(Form("%s/xyPos%s.eps", outputdir.c_str(), suffix.c_str()) );
  c1->SaveAs(Form("%s/xyPos%s.pdf", outputdir.c_str(), suffix.c_str()) );
  c1->SaveAs(Form("%s/xyPos%s.png", outputdir.c_str(), suffix.c_str()) );

  c1->Clear();

  xMax = xySize/2.;
  TH2D* h2_axes_zoom = new TH2D("axes_zoom", "", 10, -xMax, xMax, 10, -xMax, xMax);
  h2_axes_zoom->SetXTitle("X Position [mm]");
  h2_axes_zoom->SetYTitle("Y Position [mm]");
  h2_axes_zoom->Draw();

  drawBeam = ((fabs(beamX)<xMax) && (fabs(beamY)<xMax));
  if( drawBeam )
    legend2->Draw("same");
  label_top->Draw("same");
  label_run->Draw("same");
  legend->Draw("same");

  h2_xyPos_bgo->Draw("same");
  h2_xyPos_hodo->Draw("same");
  //h2_xyPos_calo->Draw("same");
  //h2_xyPos->Draw("same");


  lineHodo_x1->Draw("same");
  lineHodo_x2->Draw("same");
  lineHodo_y1->Draw("same");
  lineHodo_y2->Draw("same");

  beamPos->Draw("same"); 

  //gr_xyCenter_hodo->Draw("p same"); // now using hodo_fit
  gr_xyCenter_hodo_fit->Draw("p same");
  gr_xyCenter_bgo->Draw("p same");
  //gr_xyCenter->Draw("p same");
  //gr_xyCenter_calo->Draw("p same");
  //gr_xyPos_fit->Draw("p same");


  c1->SaveAs(Form("%s/xyPos%s_zoom.eps", outputdir.c_str(), suffix.c_str()) );
  c1->SaveAs(Form("%s/xyPos%s_zoom.pdf", outputdir.c_str(), suffix.c_str()) );
  c1->SaveAs(Form("%s/xyPos%s_zoom.png", outputdir.c_str(), suffix.c_str()) );

  delete c1;
  delete h2_axes;
  delete h2_axes_zoom;
  delete legend;

}












TGraphErrors* get_xyCenter( TH2D* h2_xyPos ) {

  float x     = h2_xyPos->ProjectionX()->GetMean();
  float x_err = h2_xyPos->ProjectionX()->GetRMS();

  float y     = h2_xyPos->ProjectionY()->GetMean();
  float y_err = h2_xyPos->ProjectionY()->GetRMS();

  TGraphErrors* gr_point = new TGraphErrors(0);
  gr_point->SetPoint( 0, x, y );
  gr_point->SetPointError( 0, x_err, y_err );

  return gr_point;

}







void fitHodoWithBeam( const std::string& outputdir, const std::string& suffix, TH1D* h1, float r, float& pos, float& pos_err ) {

  h1->Rebin(5);
  float nentries =  h1->GetEntries();

  TF1* f1_gaus = new TF1( "gaus_hodo", "gaus" );
  f1_gaus->SetRange(-4., 4.);
  //f1_gaus->SetParameter(0, nentries);
  f1_gaus->FixParameter(0, nentries);
  f1_gaus->SetParameter(1, 0.);

  //f1_gaus->SetParLimits(0, 0.01*nentries, nentries);
  f1_gaus->SetParLimits(1, -40., 40.);
  f1_gaus->FixParameter(2, r );

  h1->Fit(f1_gaus, "RNL" );
  TCanvas* c1 = new TCanvas("c1_temp", "", 600, 600);
  c1->cd();
  h1->SetLineColor(kRed);
  h1->SetLineWidth(2);
  h1->SetXTitle("Position [mm]");
  h1->SetYTitle("Hits");
  h1->Draw("");
  f1_gaus->SetRange(-40., 40.);
  f1_gaus->Draw("same");
  c1->SaveAs(Form("%s/tmpFit%s.eps", outputdir.c_str(), suffix.c_str()));
  c1->SaveAs(Form("%s/tmpFit%s.pdf", outputdir.c_str(), suffix.c_str()));
  c1->SaveAs(Form("%s/tmpFit%s.png", outputdir.c_str(), suffix.c_str()));
  delete c1;

  pos = f1_gaus->GetParameter(1);
  pos_err = f1_gaus->GetParameter(2);

  delete f1_gaus;

}









void drawPositionResolutionXY( const std::string& outputdir, TFile* file, const std::string& runName, const std::string& name1, const std::string& name2 ) {

  std::string treeName1 = (name1=="Beam") ? name1 : "Pos_" + name1;
  std::string treeName2 = (name2=="Beam") ? name2 : "Pos_" + name2;

  TTree* tree = (TTree*)file->Get("posTree");

  int nBins = 50;
  float xMax = 40.;

  TH1D* h1_resoX = new TH1D("resoX", "", nBins, -xMax, xMax);
  TH1D* h1_resoY = new TH1D("resoY", "", nBins, -xMax, xMax);

  tree->Project( "resoX", Form("x%s-x%s", treeName1.c_str(), treeName2.c_str()), "isSingleEle_scintFront && bgo_corr_ok" );
  tree->Project( "resoY", Form("y%s-y%s", treeName1.c_str(), treeName2.c_str()), "isSingleEle_scintFront && bgo_corr_ok" );

  h1_resoX->SetLineColor(38);
  h1_resoX->SetLineWidth(2);
  h1_resoY->SetLineColor(46);
  h1_resoY->SetLineWidth(2);


  TCanvas* c1 = new TCanvas("cc", "", 600, 600);
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, -xMax, xMax, 10, 0., h1_resoX->GetMaximum()*1.3);
  h2_axes->SetXTitle( Form("%s - %s [mm]", name1.c_str(), name2.c_str()) );  
  h2_axes->SetYTitle( "Entries" );
  h2_axes->Draw();

  h1_resoX->Draw("same");
  h1_resoY->Draw("same");

  TPaveText* labelTop = DrawTools::getLabelTop();
  labelTop->Draw("same");

  TPaveText* label_x = new TPaveText( 0.2, 0.8, 0.5, 0.9, "brNDC" );
  label_x->SetTextSize( 0.038 );
  label_x->SetTextColor(38);
  label_x->SetFillColor(0);
  label_x->AddText( Form("Mean(x): %.2f", h1_resoX->GetMean()) ); 
  label_x->AddText( Form("RMS(x) : %.2f", h1_resoX->GetRMS()) ); 
  label_x->Draw("same");

  TPaveText* label_y = new TPaveText( 0.6, 0.8, 0.9, 0.9, "brNDC" );
  label_y->SetTextSize( 0.038 );
  label_y->SetTextColor(46);
  label_y->SetFillColor(0);
  label_y->AddText( Form("Mean(y): %.2f", h1_resoY->GetMean()) ); 
  label_y->AddText( Form("RMS(y) : %.2f", h1_resoY->GetRMS()) ); 
  label_y->Draw("same");
  
  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/posReso_%s_vs_%s.eps", outputdir.c_str(), name1.c_str(), name2.c_str()) );
  c1->SaveAs( Form("%s/posReso_%s_vs_%s.png", outputdir.c_str(), name1.c_str(), name2.c_str()) );
  c1->SaveAs( Form("%s/posReso_%s_vs_%s.pdf", outputdir.c_str(), name1.c_str(), name2.c_str()) );

  delete c1;
  delete h2_axes;

}
