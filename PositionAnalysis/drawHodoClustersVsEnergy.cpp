#include <iostream>
#include <string>
#include <stdlib.h>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TString.h"

#include "interface/DrawTools.h"


std::string runName;




int main( int argc, char* argv[] ) {

  
  runName = "BTF_92_20140430-020137_beam";
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


  TTree* tree = (TTree*)file->Get("tree_passedEvents");


  int nHodoClustersX;
  tree->SetBranchAddress( "nHodoClustersX", &nHodoClustersX );
  int nHodoClustersY;
  tree->SetBranchAddress( "nHodoClustersY", &nHodoClustersY );
  float cef3_corr[4];
  tree->SetBranchAddress( "cef3_corr", cef3_corr );


  int nBins = 1500;
  float xMin = 0.;
  float xMax = 15000.;

  TH2D* h2_eVsN_x   = new TH2D("eVsN_x"   , "",   5, -0.5, 4.5, nBins, xMin, xMax );
  TH2D* h2_eVsN_y   = new TH2D("eVsN_y"   , "",   5, -0.5, 4.5, nBins, xMin, xMax );
  TH2D* h2_eVsN_same= new TH2D("eVsN_same", "",   5, -0.5, 4.5, nBins, xMin, xMax );
  TH2D* h2_eVsN_sum = new TH2D("eVsN_sum" , "",   10,-0.5, 9.5, nBins, xMin, xMax );

  int nentries = tree->GetEntries();

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry( iEntry );

    float energy = cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3];
    //float meanHodoClusters = 0.5*(nHodoClustersX+nHodoClustersY) ;

    h2_eVsN_x  ->Fill( nHodoClustersX  , energy );
    h2_eVsN_y  ->Fill( nHodoClustersY  , energy );
    h2_eVsN_sum->Fill( nHodoClustersX+nHodoClustersY, energy );
    if( nHodoClustersX==nHodoClustersY )
      h2_eVsN_same->Fill( nHodoClustersY, energy );

  }

  TFile* outfile = TFile::Open("prova.root", "recreate");
  outfile->cd();

  h2_eVsN_x->Write();
  h2_eVsN_y->Write();
  h2_eVsN_sum->Write();
  h2_eVsN_same->Write();
 
  outfile->Close();


  return 0;

}










// +++++ helper and cosmetic functions:


TStyle* setStyle() {

  // set the TStyle
  TStyle* style_ = new TStyle("DrawBaseStyle", "");
  style_->SetCanvasColor(0);
  style_->SetPadColor(0);
  style_->SetFrameFillColor(0);
  style_->SetStatColor(0);
  style_->SetOptStat(0);
  style_->SetTitleFillColor(0);
  style_->SetCanvasBorderMode(0);
  style_->SetPadBorderMode(0);
  style_->SetFrameBorderMode(0);
  style_->SetPadBottomMargin(0.12);
  style_->SetPadLeftMargin(0.12);
  style_->cd();

  // For the canvas:
  style_->SetCanvasBorderMode(0);
  style_->SetCanvasColor(kWhite);
  style_->SetCanvasDefH(600); //Height of canvas
  style_->SetCanvasDefW(600); //Width of canvas
  style_->SetCanvasDefX(0); //POsition on screen
  style_->SetCanvasDefY(0);

  // For the Pad:
  style_->SetPadBorderMode(0);
  style_->SetPadColor(kWhite);
  style_->SetPadGridX(false);
  style_->SetPadGridY(false);
  style_->SetGridColor(0);
  style_->SetGridStyle(3);
  style_->SetGridWidth(1);

  // For the frame:
  style_->SetFrameBorderMode(0);
  style_->SetFrameBorderSize(1);
  style_->SetFrameFillColor(0);
  style_->SetFrameFillStyle(0);
  style_->SetFrameLineColor(1);
  style_->SetFrameLineStyle(1);
  style_->SetFrameLineWidth(1);


  // Margins:
  style_->SetPadTopMargin(0.05);
  style_->SetPadBottomMargin(0.15);//0.13);
  style_->SetPadLeftMargin(0.15);//0.16);
  style_->SetPadRightMargin(0.05);//0.02);

  // For the Global title:

  style_->SetOptTitle(0);
  style_->SetTitleFont(42);
  style_->SetTitleColor(1);
  style_->SetTitleTextColor(1);
  style_->SetTitleFillColor(10);
  style_->SetTitleFontSize(0.05);

  // For the axis titles:

  style_->SetTitleColor(1, "XYZ");
  style_->SetTitleFont(42, "XYZ");
  style_->SetTitleSize(0.05, "XYZ");
  style_->SetTitleXOffset(1.15);//0.9);
  style_->SetTitleYOffset(1.4); // => 1.15 if exponents

  // For the axis labels:

  style_->SetLabelColor(1, "XYZ");
  style_->SetLabelFont(42, "XYZ");
  style_->SetLabelOffset(0.007, "XYZ");
  style_->SetLabelSize(0.045, "XYZ");

  // For the axis:

  style_->SetAxisColor(1, "XYZ");
  style_->SetStripDecimals(kTRUE);
  style_->SetTickLength(0.03, "XYZ");
  style_->SetNdivisions(510, "XYZ");
  style_->SetPadTickX(1); // To get tick marks on the opposite side of the frame
  style_->SetPadTickY(1);

  style_->cd();

  return style_;

}








