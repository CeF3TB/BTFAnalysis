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
//void drawSinglePlot( const std::string& outputdir, const std::string& saveName, TFile* file, const std::string& name, const std::string& axisName, int nChannels, float xMin=0, float xMax=4095, int rebin=1, bool plotLog=false );
void drawSinglePlot( const std::string& outputdir, const std::string& saveName, TFile* file, const std::string& varName, const std::string& selection, const std::string& axisName, int nChannels, int nBins, float xMin=0, float xMax=4095, bool plotLog=false );
void fitHodoWithBeam( const std::string& outputdir, const std::string& suffix, TH1D* h1, float r, float& pos, float& pos_err );
//TGraphErrors* getFitPositionCeF3( TFile* file );



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


  std::string fileName = "PosAn_" + runName + ".root";
  TFile* file = TFile::Open( fileName.c_str() );
  std::cout << "-> Opened file: " << fileName << std::endl;
  
  std::string outputdir = "Plots_" + runName;
  std::string mkdir_command = "mkdir -p " + outputdir;
  system(mkdir_command.c_str());


  TStyle* style = DrawTools::setStyle();
  style->cd();

  drawSinglePositionPlot( outputdir, file, runName, "" );
  drawSinglePositionPlot( outputdir, file, runName, "_singleEle" );

  //drawSinglePlot( outputdir, "cef3_spectrum"      , file, "cef3"     , "ADC Counts", 4, 0., 3500., 10, true );
  //drawSinglePlot( outputdir, "cef3_corr_spectrum" , file, "cef3_corr", "ADC Counts", 4, 0., 3500., 10, true );

  if( runName=="BTF_000001_cosmics" ) {

    drawSinglePlot( outputdir, "cef3_spectrum"      , file, "cef3"     , "", "ADC Counts", 4, 50, 0., 200., false );
    drawSinglePlot( outputdir, "cef3_corr_spectrum" , file, "cef3_corr", "", "ADC Counts", 4, 50, 0., 200., false );
  
  } else {

    drawSinglePlot( outputdir, "cef3_spectrum"      , file, "cef3"     , "", "ADC Counts", 4, 200, 0., 3500., false );
    drawSinglePlot( outputdir, "cef3_corr_spectrum" , file, "cef3_corr", "", "ADC Counts", 4, 200, 0., 3500., false );
  
    drawSinglePlot( outputdir, "cef3_spectrum_singleEle"      , file, "cef3"     , "scintFront>500. && scintFront<2000.", "ADC Counts", 4, 200, 0., 3500., false );
    drawSinglePlot( outputdir, "cef3_corr_spectrum_singleEle" , file, "cef3_corr", "scintFront>500. && scintFront<2000.", "ADC Counts", 4, 200, 0., 3500., false );
  
    drawSinglePlot( outputdir, "cef3_spectrum_singleEle_hodo"      , file, "cef3"     , "scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1", "ADC Counts", 4, 200, 0., 3500., false );
    drawSinglePlot( outputdir, "cef3_corr_spectrum_singleEle_hodo" , file, "cef3_corr", "scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1", "ADC Counts", 4, 200, 0., 3500., false );

  }


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


  //TGraphErrors* gr_xyPos_fit = getFitPositionCeF3(file);
  //gr_xyPos_fit->SetMarkerColor(kRed+2);
  //gr_xyPos_fit->SetLineColor(kRed+2);
  //gr_xyPos_fit->SetMarkerStyle(24);
  //gr_xyPos_fit->SetMarkerSize(1.6);

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
  legend->AddEntry( gr_xyCenter_calo, "Calo", "P" );
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
  h2_xyPos->Draw("same");
  h2_xyPos_hodo->Draw("same");
  h2_xyPos_calo->Draw("same");


  TEllipse* beamPos = new TEllipse( beamX, beamY, beamRX, beamRY );
  beamPos->SetLineColor(kBlack);
  beamPos->SetFillStyle(0);
  beamPos->Draw("same"); 

  //gr_xyCenter_hodo->Draw("p same");  // now using hodo_fit
  gr_xyCenter_hodo_fit->Draw("p same");
  gr_xyCenter_bgo->Draw("p same");
  gr_xyCenter->Draw("p same"); // don't draw for now
  gr_xyCenter_calo->Draw("p same");
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
  h2_xyPos_calo->Draw("same");
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
  gr_xyCenter_calo->Draw("p same");
  //gr_xyPos_fit->Draw("p same");


  c1->SaveAs(Form("%s/xyPos%s_zoom.eps", outputdir.c_str(), suffix.c_str()) );
  c1->SaveAs(Form("%s/xyPos%s_zoom.pdf", outputdir.c_str(), suffix.c_str()) );
  c1->SaveAs(Form("%s/xyPos%s_zoom.png", outputdir.c_str(), suffix.c_str()) );

  delete c1;
  delete h2_axes;
  delete h2_axes_zoom;
  delete legend;

}




void drawSinglePlot( const std::string& outputdir, const std::string& saveName, TFile* file, const std::string& varName, const std::string& selection, const std::string& axisName, int nChannels, int nBins, float xMin, float xMax, bool plotLog ) {


  TTree* tree = (TTree*)file->Get("tree_passedEvents");


  std::vector<int> colors;
  colors.push_back( 46 );
  colors.push_back( 38 );
  colors.push_back( 30 );
  colors.push_back( 42 );
  colors.push_back( 29 );
  colors.push_back( kBlack );
  colors.push_back( kGreen );
  colors.push_back( kBlue  );

  float yMax_leg = 0.8;
  float yMin_leg = yMax_leg-0.05*nChannels;
  TLegend* legend = new TLegend( 0.6, yMin_leg, 0.9, yMax_leg );
  legend->SetTextSize(0.035);
  legend->SetFillColor(0);

  float yMax=0.;
  std::vector<TH1D*> histos;

  for(unsigned i=0; i<nChannels; ++i ) {

    std::string histoName( Form("histo_%d", i) );
    TH1D* h1 = new TH1D( histoName.c_str(), "", nBins, xMin, xMax );
    tree->Project( histoName.c_str(), Form("%s[%d]", varName.c_str(), i), selection.c_str() );

    h1->SetLineColor( colors[i] );
    h1->SetLineWidth( 2 );

    TF1* f1 = new TF1( Form("gaus_%d", i), "gaus", 10., 45. );
    f1->SetLineColor( colors[i] );
    f1->SetLineWidth( 2 );
    h1->Fit( f1, "R+" );
  
    histos.push_back(h1);

    legend->AddEntry( h1, Form("Channel %d", i), "L" );


    float thisMax = h1->GetMaximum();
    if( thisMax>yMax ) yMax = thisMax;

  }

  


  TCanvas* c1 = new TCanvas( "c1_new", "", 600, 600 );
  c1->cd();
  if( plotLog ) c1->SetLogy();


  float yScaleFactor = (plotLog) ? 5.  : 1.1;
  float yMin         = (plotLog) ? 1. :  0.;


  TH2D* h2_axes = new TH2D("axes_new", "", 10, xMin, xMax, 10, yMin, yMax*yScaleFactor );
  h2_axes->SetXTitle( axisName.c_str() );
  h2_axes->SetYTitle( "Entries" );

  h2_axes->Draw();

  legend->Draw("Same");

  TPaveText* label_top = DrawTools::getLabelTop();
  TPaveText* label_run = DrawTools::getLabelRun(runName, !plotLog);
  label_top->Draw("same");
  label_run->Draw("same");

  for(unsigned i=0; i<nChannels; ++i ) 
    histos[i]->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/%s.eps", outputdir.c_str(), saveName.c_str()) );
  c1->SaveAs( Form("%s/%s.pdf", outputdir.c_str(), saveName.c_str()) );
  c1->SaveAs( Form("%s/%s.png", outputdir.c_str(), saveName.c_str()) );


  c1->Clear();

  // and now draw also individual channels:
  for( unsigned i=0; i<nChannels; ++i ) {
  
    h2_axes->Draw();
    label_top->Draw("Same");
    label_run->Draw("Same");

    histos[i]->SetLineColor(46);
    histos[i]->Draw("same");

    TPaveText* labelChannel = new TPaveText( 0.65, 0.8, 0.9, 0.85, "brNDC" );
    labelChannel->SetTextSize( 0.035 );
    labelChannel->SetFillColor( 0 );
    labelChannel->AddText( Form("Channel %d", i) );
    labelChannel->Draw("same");
  
    c1->SaveAs( Form("%s/%s_%d.eps", outputdir.c_str(), saveName.c_str(), i) );
    c1->SaveAs( Form("%s/%s_%d.pdf", outputdir.c_str(), saveName.c_str(), i) );
    c1->SaveAs( Form("%s/%s_%d.png", outputdir.c_str(), saveName.c_str(), i) );

    delete histos[i];

  }


  delete c1;
  delete h2_axes;
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








/*
TGraphErrors* getFitPositionCeF3( TFile* file ) {

  TH1D* h1_cef3_corr_0 = (TH1D*)file->Get("cef3_corr_0");
  TH1D* h1_cef3_corr_1 = (TH1D*)file->Get("cef3_corr_1");
  TH1D* h1_cef3_corr_2 = (TH1D*)file->Get("cef3_corr_2");
  TH1D* h1_cef3_corr_3 = (TH1D*)file->Get("cef3_corr_3");
  
  float position = 12. - 0.696;

  int nBins = 80;
  float xMax = 40.;
  TH2D* h2_energyMap = new TH2D("energyMap", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  //   0    1
  //          
  //   3    2
  h2_energyMap->SetBinContent( h2_energyMap->GetXaxis()->FindBin(-position), h2_energyMap->GetYaxis()->FindBin(+position), h1_cef3_corr_0->GetMean() );
  h2_energyMap->SetBinContent( h2_energyMap->GetXaxis()->FindBin(+position), h2_energyMap->GetYaxis()->FindBin(+position), h1_cef3_corr_1->GetMean() );
  h2_energyMap->SetBinContent( h2_energyMap->GetXaxis()->FindBin(+position), h2_energyMap->GetYaxis()->FindBin(-position), h1_cef3_corr_2->GetMean() );
  h2_energyMap->SetBinContent( h2_energyMap->GetXaxis()->FindBin(-position), h2_energyMap->GetYaxis()->FindBin(-position), h1_cef3_corr_3->GetMean() );


  TF2* f2_gaus = new TF2( "gaus2d", DoubleGauss,  -1.1*position, -1.1*position, 1.1*position, 1.1*position, 4 );
  f2_gaus->FixParameter(0, h1_cef3_corr_0->GetEntries() );
  f2_gaus->SetParameter(1, 0.);
  f2_gaus->SetParameter(2, 0.);
  f2_gaus->SetParameter(3, 25.);

  f2_gaus->SetParLimits(1, -30., 30.);
  f2_gaus->SetParLimits(2, -30., 30.);
  f2_gaus->SetParLimits(3, 5., 35.);

  h2_energyMap->Fit( f2_gaus, "RN" );

  TGraphErrors* gr_xyPos_fit = new TGraphErrors(0);
  gr_xyPos_fit->SetPoint(0, f2_gaus->GetParameter(1), f2_gaus->GetParameter(2) );
  gr_xyPos_fit->SetPointError(0, f2_gaus->GetParameter(3), f2_gaus->GetParameter(3) );

  delete h2_energyMap;

  return gr_xyPos_fit;

}
*/
