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

#include "fastDQM_CeF3_BTF.h"




struct HistoAndName {

  HistoAndName( TH1D* h, const std::string& n ) {
    histo = h;
    name = n;
  }

  TH1D* histo;
  std::string name;

};



struct Histo2AndName {

  Histo2AndName( TH2D* h, const std::string& n ) {
    histo = h;
    name = n;
  }

  TH2D* histo;
  std::string name;

};




TH2D* getSingleHisto( TTree* tree, const std::string& name, const std::string& varName, const std::string& varExpr, const std::string& cuts );
std::pair<TH1D*,TH1D*> getBiasAndResoHistos( TH2D* h2 );
void checkLateralScan( const std::string& outputdir, const std::string& name, TTree* tree, const std::string& cut );
void drawPerformancePlot( const std::string& outputdir, const std::string& name, const std::string& var, std::vector< HistoAndName > hn, const std::string& bias_reso );
void drawProjections( const std::string& outputdir, const std::string& name, std::vector<Histo2AndName> hn, const std::string& var );


int main( int argc, char* argv[] ) {

 // if( argv==1 ) {
 //   std::cout << "USAGE: ./checkBGOPositioning [runSet] [tag=\"V02\"]" << std::endl;
 //   exit(11);
 // }



  std::string tag = "V02";
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }


  TChain* tree = new TChain("posTree");
  tree->Add( Form("PosAnTrees_%s/crossScanFile.root/posTree", tag.c_str()) );
  tree->Add( Form("PosAnTrees_%s/diagonalScanFile.root/posTree", tag.c_str()) );


  DrawTools::setStyle();

  std::string outputdir = "BGOPositioningPlots_" + tag;
  system(Form("mkdir -p %s", outputdir.c_str()));

  checkLateralScan( outputdir, "all"    , tree, "(yBeam==0. && abs(xBeam)>0.2) || (xBeam==0. && abs(yBeam)>0.2) || (xBeam==yBeam && xBeam>0.) || (xBeam==-yBeam && xBeam>0.)" );
  checkLateralScan( outputdir, "horiz"  , tree, "yBeam==0. && abs(xBeam)>0.2" );
  checkLateralScan( outputdir, "vert"   , tree, "xBeam==0. && abs(yBeam)>0.2" );
  checkLateralScan( outputdir, "diag13" , tree, "xBeam==yBeam && xBeam>0." );
  checkLateralScan( outputdir, "diag02" , tree, "xBeam==-yBeam && xBeam>0." );


  return 0;

}



void checkLateralScan( const std::string& outputdir, const std::string& name, TTree* tree, const std::string& cut ) {

  TH2D* h2_asymm_x = getSingleHisto( tree, name, "asymm_x", "(xPos_bgo_asymm-xBeam):xBeam", cut );
  TH2D* h2_asymm_y = getSingleHisto( tree, name, "asymm_y", "(yPos_bgo_asymm-yBeam):yBeam", cut );
  TH2D* h2_asymm_log_x = getSingleHisto( tree, name, "asymm_log_x", "(xPos_bgo_asymmlog-xBeam):xBeam", cut );
  TH2D* h2_asymm_log_y = getSingleHisto( tree, name, "asymm_log_y", "(yPos_bgo_asymmlog-yBeam):yBeam", cut );
  TH2D* h2_wa_x = getSingleHisto( tree, name, "wa_x", "(xPos_bgo_wa-xBeam):xBeam", cut );
  TH2D* h2_wa_y = getSingleHisto( tree, name, "wa_y", "(yPos_bgo_wa-yBeam):yBeam", cut );
  TH2D* h2_wa_log_x = getSingleHisto( tree, name, "wa_log_x", "(xPos_bgo_walog-xBeam):xBeam", cut );
  TH2D* h2_wa_log_y = getSingleHisto( tree, name, "wa_log_y", "(yPos_bgo_walog-yBeam):yBeam", cut );

//TFile* file = TFile::Open("prova.root", "recreate");
//file->cd();
//  h2_asymm_x->Write();
//  h2_asymm_y->Write();
//  h2_asymm_log_x->Write();
//  h2_asymm_log_y->Write();
//  h2_wa_x->Write();
//  h2_wa_y->Write();
//  h2_wa_log_x->Write();
//  h2_wa_log_y->Write();
//file->Close();
//exit(1);

  std::vector< Histo2AndName > hn_x;
  hn_x.push_back( Histo2AndName(h2_asymm_x, "Asymmetry") );
  hn_x.push_back( Histo2AndName(h2_asymm_log_x, "Log-Asymmetry") );
  hn_x.push_back( Histo2AndName(h2_wa_x, "Weighted Average") );
  hn_x.push_back( Histo2AndName(h2_wa_log_x, "Log-Weighted Average") );

  std::vector< Histo2AndName > hn_y;
  hn_y.push_back( Histo2AndName(h2_asymm_y, "Asymmetry") );
  hn_y.push_back( Histo2AndName(h2_asymm_log_y, "Log-Asymmetry") );
  hn_y.push_back( Histo2AndName(h2_wa_y, "Weighted Average") );
  hn_y.push_back( Histo2AndName(h2_wa_log_y, "Log-Weighted Average") );

  drawProjections( outputdir, name, hn_x, "x" );
  drawProjections( outputdir, name, hn_y, "y" );

  std::pair<TH1D*,TH1D*> p_asymm_x = getBiasAndResoHistos(h2_asymm_x);
  std::pair<TH1D*,TH1D*> p_asymm_y = getBiasAndResoHistos(h2_asymm_y);
  std::pair<TH1D*,TH1D*> p_asymm_log_x = getBiasAndResoHistos(h2_asymm_log_x);
  std::pair<TH1D*,TH1D*> p_asymm_log_y = getBiasAndResoHistos(h2_asymm_log_y);
  std::pair<TH1D*,TH1D*> p_wa_x = getBiasAndResoHistos(h2_wa_x);
  std::pair<TH1D*,TH1D*> p_wa_y = getBiasAndResoHistos(h2_wa_y);
  std::pair<TH1D*,TH1D*> p_wa_log_x = getBiasAndResoHistos(h2_wa_log_x);
  std::pair<TH1D*,TH1D*> p_wa_log_y = getBiasAndResoHistos(h2_wa_log_y);

  std::vector< HistoAndName > hn_bias_x;
  hn_bias_x.push_back( HistoAndName(p_asymm_x.first, "Asymmetry") );
  hn_bias_x.push_back( HistoAndName(p_asymm_log_x.first, "Log-Asymmetry") );
  hn_bias_x.push_back( HistoAndName(p_wa_x.first, "Weighted Average") );
  hn_bias_x.push_back( HistoAndName(p_wa_log_x.first, "Log-Weighted Average") );

  std::vector< HistoAndName > hn_bias_y;
  hn_bias_y.push_back( HistoAndName(p_asymm_y.first, "Asymmetry") );
  hn_bias_y.push_back( HistoAndName(p_asymm_log_y.first, "Log-Asymmetry") );
  hn_bias_y.push_back( HistoAndName(p_wa_y.first, "Weighted Ave") );
  hn_bias_y.push_back( HistoAndName(p_wa_log_y.first, "Log-Weighted Ave") );

  std::vector< HistoAndName > hn_reso_x;
  hn_reso_x.push_back( HistoAndName(p_asymm_x.second, "Asymmetry") );
  hn_reso_x.push_back( HistoAndName(p_asymm_log_x.second, "Log-Asymmetry") );
  hn_reso_x.push_back( HistoAndName(p_wa_x.second, "Weighted Average") );
  hn_reso_x.push_back( HistoAndName(p_wa_log_x.second, "Log-Weighted Average") );

  std::vector< HistoAndName > hn_reso_y;
  hn_reso_y.push_back( HistoAndName(p_asymm_y.second, "Asymmetry") );
  hn_reso_y.push_back( HistoAndName(p_asymm_log_y.second, "Log-Asymmetry") );
  hn_reso_y.push_back( HistoAndName(p_wa_y.second, "Weighted Ave") );
  hn_reso_y.push_back( HistoAndName(p_wa_log_y.second, "Log-Weighted Ave") );

  if( name!="vert" ) drawPerformancePlot( outputdir, name, "x", hn_bias_x, "bias" );
  if( name!="horiz") drawPerformancePlot( outputdir, name, "y", hn_bias_y, "bias" );
  if( name!="vert" ) drawPerformancePlot( outputdir, name, "x", hn_reso_x, "reso" );
  if( name!="horiz") drawPerformancePlot( outputdir, name, "y", hn_reso_y, "reso" );

}



std::pair<TH1D*,TH1D*> getBiasAndResoHistos( TH2D* h2 ) {

  int  nBins = h2->GetXaxis()->GetNbins();
  float xMin = h2->GetXaxis()->GetXmin();
  float xMax = h2->GetXaxis()->GetXmax();

  TH1D* h1_bias = new TH1D( Form("bias_%s", h2->GetName()), "", nBins, xMin, xMax );
  TH1D* h1_reso = new TH1D( Form("reso_%s", h2->GetName()), "", nBins, xMin, xMax );

  for( unsigned i=0; i<nBins; ++i ) {

    TH1D* h1_proj = h2->ProjectionY(Form("%s_%d", h2->GetName(), i), i+1, i+1);

    if( h1_proj->GetEntries()>0 ) {
      h1_bias->SetBinContent( i+1, h1_proj->GetMean() );
      h1_bias->SetBinError  ( i+1, h1_proj->GetMeanError() );
      h1_reso->SetBinContent( i+1, h1_proj->GetRMS() );
      h1_reso->SetBinError  ( i+1, h1_proj->GetRMSError() );
    }

  }

  std::pair<TH1D*,TH1D*> returnPair;
  returnPair.first  = h1_bias;
  returnPair.second = h1_reso;

  return returnPair;

}




TH2D* getSingleHisto( TTree* tree, const std::string& name, const std::string& varName, const std::string& varExpr, const std::string& cuts ) {

  std::string histoName = name + "_" + varName;
  TH2D* h2 = new TH2D( histoName.c_str(), "", 81, -40.5, 40.5, 40, -40., 40.);

  std::string fullCuts = cuts + " && isSingleEle_scintFront";
  tree->Project( histoName.c_str(), varExpr.c_str(), fullCuts.c_str() );

  return h2;

}



void drawProjections( const std::string& outputdir, const std::string& name, std::vector<Histo2AndName> hn, const std::string& var ) {

  std::vector<int> colors;
  colors.push_back( 38 );
  colors.push_back( 46 );
  colors.push_back( 29 );
  colors.push_back( 42 );
  colors.push_back( 40 );
  colors.push_back( 41 );
  colors.push_back( 45 );
  colors.push_back( 30 );
  colors.push_back( 16 );  


  int  nBins = hn[0].histo->GetXaxis()->GetNbins();
  float xMin = hn[0].histo->GetXaxis()->GetXmin();
  float xMax = hn[0].histo->GetXaxis()->GetXmax();
  float binSize = (xMax-xMin)/((float)nBins);

  for( unsigned i=0; i<nBins; ++i ) {

    TCanvas* c1 = new TCanvas("cxx", "", 600, 600);
    c1->cd();

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 0.2);
    h2_axes->SetXTitle("Estimate - Beam [mm]");
    h2_axes->SetYTitle("Normalized to Unity");
    h2_axes->Draw();

    float xMin_i = xMin + (float)i*binSize;

    TLegend* legend = new TLegend( 0.2, 0.6, 0.45, 0.9, Form("%.1f < %s < %.1f mm", xMin_i, var.c_str(), xMin_i+binSize) );
    legend->SetTextSize(0.033);
    legend->SetFillColor(0);

    bool saveIt = true;

    for( unsigned j=0; j<hn.size(); j++ ) { 

      TH1D* h1_proj = hn[j].histo->ProjectionY(Form("%s_%d", hn[j].histo->GetName(), i), i+1, i+1);
  
      if( h1_proj->GetEntries()==0 ) 
        saveIt = false;

      if( saveIt ) {
        h1_proj->SetLineColor( colors[j] );
        h1_proj->SetLineWidth(2);
        h1_proj->DrawNormalized("same");
        legend->AddEntry( h1_proj, hn[j].name.c_str(), "L" );
      }
      
    }

    if( saveIt ) {

      legend->Draw("same");

      TPaveText* labelTop = DrawTools::getLabelTop();
      labelTop->Draw("same");

      gPad->RedrawAxis();

      std::string subdir = "proj_" + name;
      system( Form("mkdir -p %s/%s", outputdir.c_str(), subdir.c_str()) );


      c1->SaveAs(Form("%s/%s/proj%s_%s_%d.eps", outputdir.c_str(), subdir.c_str(), var.c_str(), name.c_str(), i) );
      c1->SaveAs(Form("%s/%s/proj%s_%s_%d.png", outputdir.c_str(), subdir.c_str(), var.c_str(), name.c_str(), i) );
      c1->SaveAs(Form("%s/%s/proj%s_%s_%d.pdf", outputdir.c_str(), subdir.c_str(), var.c_str(), name.c_str(), i) );

    }

    delete c1;
    delete h2_axes;

  } // i

}


void drawPerformancePlot( const std::string& outputdir, const std::string& name, const std::string& var, std::vector< HistoAndName > hn, const std::string& bias_reso ) {

  bool isReso;
  if( bias_reso=="bias" ) isReso = false;
  else if( bias_reso=="reso" ) isReso = true;
  else {
    std::cout << "ERROR! bias_reso has to be either 'bias' or 'reso'. Exiting." << std::endl;
    exit(91);
  }


  std::vector<int> colors;
  colors.push_back( 38 );
  colors.push_back( 46 );
  colors.push_back( 29 );
  colors.push_back( 42 );
  colors.push_back( 40 );
  colors.push_back( 41 );
  colors.push_back( 45 );
  colors.push_back( 30 );
  colors.push_back( 16 );  




  TCanvas* c1 = new TCanvas( "cx", "", 600, 600 );
  c1->Clear();
  c1->cd();


  // first: bias
  float xMax = 15.;
  float yMin = (isReso) ? 0. : -1.;
  float yMax = (isReso) ? 15. : 3.;

  TH2D* h2_axes = new TH2D("axes", "", 10, -xMax, xMax, 10, yMin, yMax );
  h2_axes->SetXTitle( Form( "Beam %s Position [mm]", var.c_str() ) );
  if(isReso)
    h2_axes->SetYTitle( "RMS(Estimate - Beam) [mm]" );
  else
    h2_axes->SetYTitle( "Estimate - Beam [mm]" );

  

  TLegend* legend = new TLegend( 0.2, 0.69, 0.4, 0.92 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.033);

  for( unsigned i=0; i<hn.size(); ++i ) {

    hn[i].histo->SetLineColor( colors[i] );
    hn[i].histo->SetMarkerColor( colors[i] );
    hn[i].histo->SetMarkerSize( 1.6 );
    hn[i].histo->SetMarkerStyle( 20+i );

    TF1* f1 = new TF1( Form("f1_%s_%s_%d", name.c_str(), var.c_str(), i), "[0]", -xMax, xMax);
    f1->SetRange(-xMax, xMax);
    f1->SetLineWidth(2);
    f1->SetLineColor(colors[i]);
    hn[i].histo->Fit(f1, "R+");

    legend->AddEntry( hn[i].histo, hn[i].name.c_str(), "P" );
  
  }

  h2_axes->Draw();

  for( unsigned i=0; i<hn.size(); ++i )
    hn[i].histo->Draw("p same");

  legend->Draw("same");

  TLine* lineZero = new TLine(-xMax, 0., xMax, 0.);
  lineZero->Draw("same");

  TPaveText* labelTop = DrawTools::getLabelTop();
  labelTop->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/%s_%s_%s.eps", outputdir.c_str(), name.c_str(), bias_reso.c_str(), var.c_str()) );
  c1->SaveAs( Form("%s/%s_%s_%s.png", outputdir.c_str(), name.c_str(), bias_reso.c_str(), var.c_str()) );
  c1->SaveAs( Form("%s/%s_%s_%s.pdf", outputdir.c_str(), name.c_str(), bias_reso.c_str(), var.c_str()) );

  delete c1;
  delete h2_axes;


}


