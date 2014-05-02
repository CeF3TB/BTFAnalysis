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



std::string runName;


TGraphErrors* get_xyCenter( TH2D* h2_xyPos );
void drawSinglePositionPlot( const std::string& outputdir, TFile* file, const std::string& runName, const std::string& suffix );
void drawSinglePlot( const std::string& outputdir, const std::string& saveName, TFile* file, const std::string& name, const std::string& axisName, int nChannels, float xMin=0, float xMax=4095, int rebin=1, bool plotLog=false );
void fitHodoWithBeam( const std::string& outputdir, const std::string& suffix, TH1D* h1, float r, float& pos, float& pos_err );
//TGraphErrors* getFitPositionCeF3( TFile* file );
void getBeamPosition( const std::string& runName, float& beamX, float& beamY );



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

  drawSinglePlot( outputdir, "cef3_spectrum_lin"      , file, "cef3"     , "ADC Counts", 4, 0., 3500., 10, false );
  drawSinglePlot( outputdir, "cef3_corr_spectrum_lin" , file, "cef3_corr", "ADC Counts", 4, 0., 3500., 10, false );


  return 0;

}






void drawSinglePositionPlot( const std::string& outputdir, TFile* file, const std::string& runName, const std::string& suffix ) {

  // manually set beam nominal position for some known runs:
  float beamX = -999.;
  float beamY = -999.;
  float beamRX = 4.;
  float beamRY = 2.;
  getBeamPosition( runName, beamX, beamY );


  bool drawBeam = ((beamX>-999.) && (beamY>-999.));
  //bool beamInsideHodo = ((fabs(beamX)<4.) && (fabs(beamY)<4.));


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TH2D* h2_xyPos      = (TH2D*)file->Get(Form("xyPos_new%s", suffix.c_str())); 
  TH2D* h2_xyPos_hodo = (TH2D*)file->Get(Form("xyPos%s_hodo", suffix.c_str())); 
  //TH2D* h2_xyPos_hodo = (TH2D*)file->Get(Form("xyPos%s_hodo", suffix.c_str())); 
  TH2D* h2_xyPos_bgo  = (TH2D*)file->Get(Form("xyPos%s_bgo", suffix.c_str())); 

  float xySize = 25.;
  float xMax = xySize*3./2.;
  TH2D* h2_axes = new TH2D("axes", "", 10, -xMax, xMax, 10, -xMax, xMax);
  h2_axes->SetXTitle("X Position [mm]");
  h2_axes->SetYTitle("Y Position [mm]");

  h2_axes->Draw();


  TGraphErrors* gr_xyCenter      = get_xyCenter( h2_xyPos );
  TGraphErrors* gr_xyCenter_hodo = get_xyCenter( h2_xyPos_hodo );
  TGraphErrors* gr_xyCenter_bgo  = get_xyCenter( h2_xyPos_bgo  );

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




  TGraph* gr_beamPos = new TGraph(0);
  gr_beamPos->SetMarkerStyle(24);
  gr_beamPos->SetMarkerSize(4);
  TLegend* legend2 = new TLegend( 0.18, 0.21, 0.52, 0.25 );
  legend2->SetFillColor(0);
  legend2->SetTextSize(0.035);
  legend2->AddEntry( gr_beamPos, "Beam Position", "P" );
  if( drawBeam ) 
    legend2->Draw("same");


  TPaveText* label_top = DrawTools::getLabelTop();
  TPaveText* label_run = DrawTools::getLabelRun(runName);
  label_top->Draw("same");
  label_run->Draw("same");


  int lineColor = 17;

  TLine* line_x1 = new TLine( -xMax, -xySize/2., +xMax, -xySize/2. );
  line_x1->SetLineColor(lineColor);
  line_x1->Draw("same");

  TLine* line_x2 = new TLine( -xMax, +xySize/2., +xMax, +xySize/2. );
  line_x2->SetLineColor(lineColor);
  line_x2->Draw("same");

  TLine* line_y1 = new TLine( -xySize/2., -xMax, -xySize/2., +xMax );
  line_y1->SetLineColor(lineColor);
  line_y1->Draw("same");

  TLine* line_y2 = new TLine( +xySize/2., -xMax, +xySize/2., +xMax );
  line_y2->SetLineColor(lineColor);
  line_y2->Draw("same");


  TLegend* legend = new TLegend( 0.75, 0.21, 0.9, 0.39 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( gr_xyCenter, "CeF3", "P" );
  legend->AddEntry( gr_xyCenter_bgo, "BGO", "P" );
  legend->AddEntry( gr_xyCenter_hodo_fit, "Hodo", "P" );
  legend->Draw("same");


  h2_xyPos_hodo->SetMarkerColor(14);
  h2_xyPos->SetMarkerColor(46);
  h2_xyPos_bgo->SetMarkerColor(30);

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


  TEllipse* beamPos = new TEllipse( beamX, beamY, beamRX, beamRY );
  beamPos->SetLineColor(kBlack);
  beamPos->SetFillStyle(0);
  beamPos->Draw("same"); 

  //gr_xyCenter_hodo->Draw("p same");  // now using hodo_fit
  gr_xyCenter_hodo_fit->Draw("p same");
  gr_xyCenter_bgo->Draw("p same");
  gr_xyCenter->Draw("p same");
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
  h2_xyPos->Draw("same");


  lineHodo_x1->Draw("same");
  lineHodo_x2->Draw("same");
  lineHodo_y1->Draw("same");
  lineHodo_y2->Draw("same");

  beamPos->Draw("same"); 

  //gr_xyCenter_hodo->Draw("p same"); // now using hodo_fit
  gr_xyCenter_hodo_fit->Draw("p same");
  gr_xyCenter_bgo->Draw("p same");
  gr_xyCenter->Draw("p same");
  //gr_xyPos_fit->Draw("p same");


  c1->SaveAs(Form("%s/xyPos%s_zoom.eps", outputdir.c_str(), suffix.c_str()) );
  c1->SaveAs(Form("%s/xyPos%s_zoom.pdf", outputdir.c_str(), suffix.c_str()) );
  c1->SaveAs(Form("%s/xyPos%s_zoom.png", outputdir.c_str(), suffix.c_str()) );

  delete c1;
  delete h2_axes;
  delete h2_axes_zoom;
  delete legend;

}




void drawSinglePlot( const std::string& outputdir, const std::string& saveName, TFile* file, const std::string& name, const std::string& axisName, int nChannels, float xMin, float xMax, int rebin, bool plotLog ) {


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

  float yMax;
  std::vector<TH1D*> histos;

  for(unsigned i=0; i<nChannels; ++i ) {

    TH1D* h1 = (TH1D*)file->Get(Form("%s_%d", name.c_str(), i));
    h1->SetLineColor( colors[i] );
    h1->SetLineWidth( 2 );
    h1->Rebin(rebin);
  
    histos.push_back(h1);

    legend->AddEntry( h1, Form("Channel %d", i), "L" );


    float thisMax = h1->GetMaximum();
    if( thisMax>yMax ) yMax = thisMax;

  }

  if( !plotLog ) yMax = 3000.;

  TCanvas* c1 = new TCanvas( "c1_new", "", 600, 600 );
  c1->cd();
  if( plotLog ) c1->SetLogy();


  float yScaleFactor = 1.3;
  //float yScaleFactor = (plotLog) ? 5.  : 1.3;
  float yMin         = (plotLog) ? 10. :  0.;

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





void getBeamPosition( const std::string& runName, float& beamX, float& beamY ) {

  TString runName_tstr(runName.c_str());


/*
  if( runName == "BTF_94_20140430-073300_beam" ) {
    beamX = -3.;
    beamY = +3.;
  } else if( runName == "BTF_96_20140430-083733_beam" ) {
    beamX = -6.;
    beamY = +6.;
  } else if( runName == "BTF_98_20140430-092026_beam" ) {
    beamX = -9.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_144_2014") ) {
    beamX = +12.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_145_2014") ) {
    beamX = +10.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_146_2014") ) {
    beamX = +8.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_147_2014") ) {
    beamX = +6.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_148_2014") ) {
    beamX = +4.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_149_2014") ) {
    beamX = +2.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_150_2014") ) {
    beamX = +0.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_151_2014") ) {
    beamX = -2.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_152_2014") ) {
    beamX = -4.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_153_2014") ) {
    beamX = -6.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_154_2014") ) {
    beamX = -8.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_155_2014") ) {
    beamX = -10.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_156_2014") ) {
    beamX = -12.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_157_2014") ) {
    beamX = +0.;
    beamY = +8.;
  } else if( runName_tstr.BeginsWith("BTF_158_2014") ) {
    beamX = +0.;
    beamY = +6.;
  } else if( runName_tstr.BeginsWith("BTF_159_2014") ) {
    beamX = +0.;
    beamY = +4.;
  } else if( runName_tstr.BeginsWith("BTF_160_2014") ) {
    beamX = +0.;
    beamY = +2.;
  } else if( runName_tstr.BeginsWith("BTF_161_2014") ) {
    beamX = +0.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_162_2014") ) {
    beamX = +0.;
    beamY = -2.;
  } else if( runName_tstr.BeginsWith("BTF_163_2014") ) {
    beamX = +0.;
    beamY = -4.;
  } else if( runName_tstr.BeginsWith("BTF_164_2014") ) {
    beamX = +0.;
    beamY = -6.;
  } else if( runName_tstr.BeginsWith("BTF_165_2014") ) {
    beamX = +0.;
    beamY = -8.;
  } else if( runName_tstr.BeginsWith("BTF_218_2014") ) {
    beamX = +11.3;
    beamY = +11.3;
  } else if( runName_tstr.BeginsWith("BTF_219_2014") ) {
    beamX = +11.3;
    beamY = +11.3;
  } else {
    beamX = -999.;
    beamY = -999.;
  }
*/


  float xTable=-999.;
  float yTable=-999.;
       if( runName_tstr.BeginsWith("BTF_51" ) ) { xTable = 511.0; yTable = 202.5; }
  else if( runName_tstr.BeginsWith("BTF_52" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_53" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_54" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_55" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_56" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_57" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_58" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_59" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_60" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_61" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_62" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_63" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_64" ) ) { xTable = 511.0; yTable = 202.5; }
  else if( runName_tstr.BeginsWith("BTF_65" ) ) { xTable = 511.0; yTable = 205.0; }
  else if( runName_tstr.BeginsWith("BTF_66" ) ) { xTable = 511.0; yTable = 197.5; }
  else if( runName_tstr.BeginsWith("BTF_67" ) ) { xTable = 511.0; yTable = 195.0; }
  else if( runName_tstr.BeginsWith("BTF_68" ) ) { xTable = 511.0; yTable = 199.0; }
  else if( runName_tstr.BeginsWith("BTF_69" ) ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_70" ) ) { xTable = 488.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_71" ) ) { xTable = 488.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_72" ) ) { xTable = 488.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_73" ) ) { xTable = 488.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_74" ) ) { xTable = 488.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_75" ) ) { xTable = 488.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_76" ) ) { xTable = 488.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_77" ) ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_78" ) ) { xTable = 488.0; yTable = 180.0; }
  else if( runName_tstr.BeginsWith("BTF_79" ) ) { xTable = 488.0; yTable = 180.0; }
  else if( runName_tstr.BeginsWith("BTF_80" ) ) { xTable = 510.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_81" ) ) { xTable = 510.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_82" ) ) { xTable = 532.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_83" ) ) { xTable = 532.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_84" ) ) { xTable = 532.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_85" ) ) { xTable = 534.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_86" ) ) { xTable = 534.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_87" ) ) { xTable = 534.0; yTable = 222.0; }
  else if( runName_tstr.BeginsWith("BTF_88" ) ) { xTable = 534.0; yTable = 222.0; }
  else if( runName_tstr.BeginsWith("BTF_89" ) ) { xTable = 512.0; yTable = 224.0; }
  else if( runName_tstr.BeginsWith("BTF_90" ) ) { xTable = 490.0; yTable = 224.0; }
  else if( runName_tstr.BeginsWith("BTF_91" ) ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_92" ) ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_93" ) ) { xTable = 514.0; yTable = 198.0; }
  else if( runName_tstr.BeginsWith("BTF_94" ) ) { xTable = 514.0; yTable = 198.0; }
  else if( runName_tstr.BeginsWith("BTF_95" ) ) { xTable = 517.0; yTable = 195.0; }
  else if( runName_tstr.BeginsWith("BTF_96" ) ) { xTable = 517.0; yTable = 195.0; }
  else if( runName_tstr.BeginsWith("BTF_97" ) ) { xTable = 520.0; yTable = 192.0; }
  else if( runName_tstr.BeginsWith("BTF_98" ) ) { xTable = 520.0; yTable = 192.0; }
  else if( runName_tstr.BeginsWith("BTF_99" ) ) { xTable = 520.0; yTable = 210.0; }
  else if( runName_tstr.BeginsWith("BTF_100") ) { xTable = 520.0; yTable = 210.0; }
  else if( runName_tstr.BeginsWith("BTF_103") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_116") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_117") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_118") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_119") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_120") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_131") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_132") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_133") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_134") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_135") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_136") ) { xTable = 514.0; yTable = 204.0; }
  else if( runName_tstr.BeginsWith("BTF_137") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_138") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_139") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_140") ) { xTable = 508.0; yTable = 198.0; }
  else if( runName_tstr.BeginsWith("BTF_141") ) { xTable = 508.0; yTable = 198.0; }
  else if( runName_tstr.BeginsWith("BTF_142") ) { xTable = 505.0; yTable = 195.0; }
  else if( runName_tstr.BeginsWith("BTF_143") ) { xTable = 505.0; yTable = 195.0; }
  else if( runName_tstr.BeginsWith("BTF_144") ) { xTable = 499.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_145") ) { xTable = 501.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_146") ) { xTable = 503.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_147") ) { xTable = 505.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_148") ) { xTable = 507.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_149") ) { xTable = 509.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_150") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_151") ) { xTable = 513.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_152") ) { xTable = 515.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_153") ) { xTable = 517.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_154") ) { xTable = 519.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_155") ) { xTable = 521.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_156") ) { xTable = 523.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_157") ) { xTable = 511.0; yTable = 193.0; }
  else if( runName_tstr.BeginsWith("BTF_158") ) { xTable = 511.0; yTable = 195.0; }
  else if( runName_tstr.BeginsWith("BTF_159") ) { xTable = 511.0; yTable = 197.0; }
  else if( runName_tstr.BeginsWith("BTF_160") ) { xTable = 511.0; yTable = 199.0; }
  else if( runName_tstr.BeginsWith("BTF_161") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_162") ) { xTable = 511.0; yTable = 203.0; }
  else if( runName_tstr.BeginsWith("BTF_163") ) { xTable = 511.0; yTable = 205.0; }
  else if( runName_tstr.BeginsWith("BTF_164") ) { xTable = 511.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_165") ) { xTable = 511.0; yTable = 209.0; }
  else if( runName_tstr.BeginsWith("BTF_166") ) { xTable = 505.0; yTable = 209.0; }
  else if( runName_tstr.BeginsWith("BTF_167") ) { xTable = 502.0; yTable = 192.0; }
  else if( runName_tstr.BeginsWith("BTF_168") ) { xTable = 502.0; yTable = 192.0; }
  else if( runName_tstr.BeginsWith("BTF_169") ) { xTable = 502.0; yTable = 192.0; }
  else if( runName_tstr.BeginsWith("BTF_170") ) { xTable = 502.0; yTable = 210.0; }
  else if( runName_tstr.BeginsWith("BTF_171") ) { xTable = 502.0; yTable = 210.0; }
  else if( runName_tstr.BeginsWith("BTF_172") ) { xTable = 505.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_173") ) { xTable = 505.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_174") ) { xTable = 505.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_175") ) { xTable = 508.0; yTable = 204.0; }
  else if( runName_tstr.BeginsWith("BTF_176") ) { xTable = 508.0; yTable = 204.0; }
  else if( runName_tstr.BeginsWith("BTF_177") ) { xTable = 508.0; yTable = 204.0; }
  else if( runName_tstr.BeginsWith("BTF_178") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_179") ) { xTable = 511.0; yTable = 198.0; }
  else if( runName_tstr.BeginsWith("BTF_180") ) { xTable = 511.0; yTable = 199.0; }
  else if( runName_tstr.BeginsWith("BTF_181") ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_182") ) { xTable = 511.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_183") ) { xTable = 511.0; yTable = 203.0; }
  else if( runName_tstr.BeginsWith("BTF_184") ) { xTable = 511.0; yTable = 203.0; }
  else if( runName_tstr.BeginsWith("BTF_185") ) { xTable = 511.0; yTable = 204.0; }
  else if( runName_tstr.BeginsWith("BTF_186") ) { xTable = 514.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_187") ) { xTable = 513.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_188") ) { xTable = 512.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_189") ) { xTable = 510.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_190") ) { xTable = 509.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_191") ) { xTable = 508.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_192") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_193") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_194") ) { xTable = 511.0; yTable = 186.0; }
  else if( runName_tstr.BeginsWith("BTF_195") ) { xTable = 511.0; yTable = 189.0; }
  else if( runName_tstr.BeginsWith("BTF_196") ) { xTable = 511.0; yTable = 192.0; }
  else if( runName_tstr.BeginsWith("BTF_197") ) { xTable = 511.0; yTable = 195.0; }
  else if( runName_tstr.BeginsWith("BTF_198") ) { xTable = 511.0; yTable = 198.0; }
  else if( runName_tstr.BeginsWith("BTF_199") ) { xTable = 511.0; yTable = 204.0; }
  else if( runName_tstr.BeginsWith("BTF_200") ) { xTable = 511.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_201") ) { xTable = 511.0; yTable = 210.0; }
  else if( runName_tstr.BeginsWith("BTF_202") ) { xTable = 511.0; yTable = 210.0; }
  else if( runName_tstr.BeginsWith("BTF_203") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_204") ) { xTable = 511.0; yTable = 213.0; }
  else if( runName_tstr.BeginsWith("BTF_205") ) { xTable = 511.0; yTable = 216.0; }
  else if( runName_tstr.BeginsWith("BTF_206") ) { xTable = 526.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_207") ) { xTable = 523.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_208") ) { xTable = 520.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_209") ) { xTable = 517.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_210") ) { xTable = 514.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_211") ) { xTable = 508.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_212") ) { xTable = 508.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_213") ) { xTable = 505.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_214") ) { xTable = 502.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_215") ) { xTable = 499.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_216") ) { xTable = 496.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_217") ) { xTable = 496.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_218") ) { xTable = 499.7; yTable = 189.7; }
  else if( runName_tstr.BeginsWith("BTF_219") ) { xTable = 499.7; yTable = 189.7; }
  else if( runName_tstr.BeginsWith("BTF_220") ) { xTable = 499.2; yTable = 189.2; }
  else if( runName_tstr.BeginsWith("BTF_221") ) { xTable = 498.7; yTable = 188.7; }
  else if( runName_tstr.BeginsWith("BTF_222") ) { xTable = 499.7; yTable = 189.7; }
  else if( runName_tstr.BeginsWith("BTF_223") ) { xTable = 499.7; yTable = 189.7; }
  else if( runName_tstr.BeginsWith("BTF_224") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_225") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_226") ) { xTable = 496.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_227") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_228") ) { xTable = 532.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_229") ) { xTable = 532.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_230") ) { xTable = 510.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_231") ) { xTable = 510.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_232") ) { xTable = 488.0; yTable = 180.0; }
  else if( runName_tstr.BeginsWith("BTF_233") ) { xTable = 488.0; yTable = 180.0; }
  else if( runName_tstr.BeginsWith("BTF_234") ) { xTable = 488.0; yTable = 202.1; }
  else if( runName_tstr.BeginsWith("BTF_235") ) { xTable = 488.0; yTable = 202.1; }
  else if( runName_tstr.BeginsWith("BTF_236") ) { xTable = 534.0; yTable = 200.1; }
  else if( runName_tstr.BeginsWith("BTF_237") ) { xTable = 534.0; yTable = 200.1; }
  else if( runName_tstr.BeginsWith("BTF_238") ) { xTable = 534.0; yTable = 222.1; }
  else if( runName_tstr.BeginsWith("BTF_239") ) { xTable = 534.0; yTable = 222.1; }


  beamX = -xTable + 511.;
  beamY = -yTable + 201.;

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
