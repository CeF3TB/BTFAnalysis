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




TH2D* getSingleHisto( TTree* tree, const std::string& name, const std::string& varName, const std::string& varExpr, const std::string& cuts );
std::pair<TH1D*,TH1D*> getBiasAndResoHistos( TH2D* h2 );
void checkLateralScan( const std::string& outputdir, const std::string& name, TTree* tree, const std::string& cut );
void drawSinglePlot( const std::string& outputdir, const std::string& name, const std::string& var, std::vector< HistoAndName > hn );


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

  checkLateralScan( outputdir, "horiz", tree, "yBeam==0." );


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
  hn_bias_x.push_back( HistoAndName(p_asymm_log_x.first, "Asymmetry (log)") );
  hn_bias_x.push_back( HistoAndName(p_wa_x.first, "Weighted Average") );
  hn_bias_x.push_back( HistoAndName(p_wa_log_x.first, "Weighted Average (log)") );

  std::vector< HistoAndName > hn_bias_y;
  hn_bias_y.push_back( HistoAndName(p_asymm_y.first, "Asymmetry") );
  hn_bias_y.push_back( HistoAndName(p_asymm_log_y.first, "Asymmetry (log)") );
  hn_bias_y.push_back( HistoAndName(p_wa_y.first, "Weighted Ave") );
  hn_bias_y.push_back( HistoAndName(p_wa_log_y.first, "Weighted Ave (log)") );

  std::vector< HistoAndName > hn_reso_x;
  hn_reso_x.push_back( HistoAndName(p_asymm_x.second, "Asymmetry") );
  hn_reso_x.push_back( HistoAndName(p_asymm_log_x.second, "Asymmetry (log)") );
  hn_reso_x.push_back( HistoAndName(p_wa_x.second, "Weighted Average") );
  hn_reso_x.push_back( HistoAndName(p_wa_log_x.second, "Weighted Average (log)") );

  std::vector< HistoAndName > hn_reso_y;
  hn_reso_y.push_back( HistoAndName(p_asymm_y.second, "Asymmetry") );
  hn_reso_y.push_back( HistoAndName(p_asymm_log_y.second, "Asymmetry (log)") );
  hn_reso_y.push_back( HistoAndName(p_wa_y.second, "Weighted Ave") );
  hn_reso_y.push_back( HistoAndName(p_wa_log_y.second, "Weighted Ave (log)") );

  drawSinglePlot( outputdir, name, "x", hn_bias_x );
  drawSinglePlot( outputdir, name, "y", hn_bias_y );

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



void drawSinglePlot( const std::string& outputdir, const std::string& name, const std::string& var, std::vector< HistoAndName > hn ) {

  std::vector<int> colors;
  colors.push_back( 40 );
  colors.push_back( 41 );
  colors.push_back( 42 );
  colors.push_back( 45 );
  colors.push_back( 46 );
  colors.push_back( 30 );
  colors.push_back( 38 );
  colors.push_back( 16 );  

  TCanvas* c1 = new TCanvas( "cx", "", 600, 600 );
  c1->cd();


  // first: bias

  TH2D* h2_axes = new TH2D("axes", "", 10, -40., 40., 10, -5., 7. );
  h2_axes->SetXTitle( Form( "Beam %s Position [mm]", var.c_str() ) );
  h2_axes->SetYTitle( "Estimate - Beam [mm]" );

  h2_axes->Draw();
  

  TLegend* legend = new TLegend( 0.2, 0.65, 0.4, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.033);

  for( unsigned i=0; i<hn.size(); ++i ) {

    hn[i].histo->SetLineColor( colors[i] );
    hn[i].histo->SetMarkerColor( colors[i] );
    hn[i].histo->SetMarkerSize( 1.6 );
    hn[i].histo->SetMarkerStyle( 20+i );

    legend->AddEntry( hn[i].histo, hn[i].name.c_str(), "P" );
  
    hn[i].histo->Draw("P same" );

  }

  legend->Draw("same");

  TPaveText* labelTop = DrawTools::getLabelTop();
  labelTop->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/%s_bias_%s.eps", outputdir.c_str(), name.c_str(), var.c_str()) );
  c1->SaveAs( Form("%s/%s_bias_%s.png", outputdir.c_str(), name.c_str(), var.c_str()) );
  c1->SaveAs( Form("%s/%s_bias_%s.pdf", outputdir.c_str(), name.c_str(), var.c_str()) );

  delete c1;
  delete h2_axes;


}
