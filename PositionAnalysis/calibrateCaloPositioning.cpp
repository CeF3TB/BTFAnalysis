#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <assert.h>

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TCanvas.h"

#include "fastDQM_CeF3_BTF.h"

#include "interface/PositionTools.h"





float sumVector( std::vector<float> v );



int main( int argc, char* argv[] ) {


  std::string runName = "BTF_246_beam";
  if( argc>1 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
  }

  std::string tag = "V02";
  if( argc>2 ) {
    std::string tag_str(argv[2]);
    tag = tag_str;
  }


  //std::string fileName = "PosAnTrees_" + tag + "/PosAn_" + runName + ".root";
  std::string fileName = "PosAnTrees_" + tag + "/crossScanFile.root";
  TFile* file = TFile::Open(fileName.c_str());
  if( file==0 ) {
    std::cout << "ERROR! Din't find file " << fileName << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  TTree* tree = (TTree*)file->Get("posTree");




  float cef3_corr[CEF3_CHANNELS];
  tree->SetBranchAddress( "cef3_corr", cef3_corr );
  float bgo_corr[BGO_CHANNELS];
  tree->SetBranchAddress( "bgo_corr", bgo_corr );
  float xBeam;
  tree->SetBranchAddress( "xBeam", &xBeam );
  float yBeam;
  tree->SetBranchAddress( "yBeam", &yBeam );
  float xPos_bgo_wa;
  tree->SetBranchAddress( "xPos_bgo_wa", &xPos_bgo_wa );
  float yPos_bgo_wa;
  tree->SetBranchAddress( "yPos_bgo_wa", &yPos_bgo_wa );
  bool isSingleEle_scintFront;
  tree->SetBranchAddress( "isSingleEle_scintFront", &isSingleEle_scintFront );
  int nHodoClustersX;
  tree->SetBranchAddress( "nHodoClustersX", &nHodoClustersX );
  int nHodoClustersY;
  tree->SetBranchAddress( "nHodoClustersY", &nHodoClustersY );

  float pos_hodoClustX[HODOX_CHANNELS];
  float pos_hodoClustY[HODOY_CHANNELS];
  int nFibres_hodoClustX[HODOX_CHANNELS];
  int nFibres_hodoClustY[HODOY_CHANNELS];
  tree->SetBranchAddress( "pos_hodoClustX", pos_hodoClustX );
  tree->SetBranchAddress( "nFibres_hodoClustX", nFibres_hodoClustX );
  tree->SetBranchAddress( "pos_hodoClustY", pos_hodoClustY );
  tree->SetBranchAddress( "nFibres_hodoClustY", nFibres_hodoClustY );




  std::string outfileName = "CaloPos_" + runName + ".root";
  TFile* outfile = TFile::Open( outfileName.c_str(), "recreate" );
  outfile->cd();


  float xMax = 5.5;
  int nBins = (int)2.*xMax;

  TProfile* hp_cef3CalibX = new TProfile("cef3CalibX", "", nBins, -xMax, xMax);
  TProfile* hp_cef3CalibY = new TProfile("cef3CalibY", "", nBins, -xMax, xMax);


  std::vector<float> xbgo, ybgo;
  for( unsigned i=0; i<BGO_CHANNELS; ++i ) {
    float x,y;
    PositionTools::getBGOCoordinates( i, x, y );
    xbgo.push_back( x );
    ybgo.push_back( y );
  }



  int nentries = tree->GetEntries();



  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {


    tree->GetEntry(iEntry);

    if( iEntry % 5000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;


    bool isSingleElectron = (isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1);
    //bool isSingleElectron = (isSingleEle_scintFront);
    if( !isSingleElectron ) continue;
    

    // require that clusters have no more than two fibres each
    // for better position precision
    if( nFibres_hodoClustX[0]>2 || nFibres_hodoClustY[0]>2 ) continue;


    // with one X and one Y clusters in hodo
    // define true beam position as hodo cross center
    float xTarget = pos_hodoClustX[0];  
    float yTarget = pos_hodoClustY[0];  

    //float xTarget = xBeam;  
    //float yTarget = yBeam;  

    std::vector<float> v_bgo_corr;
    for( unsigned i=0; i<BGO_CHANNELS; ++i ) v_bgo_corr.push_back(bgo_corr[i]);

    std::vector<float> xPosW_bgo;
    std::vector<float> yPosW_bgo;

    xPosW_bgo.push_back(bgo_corr[0]*xbgo[0]);
    xPosW_bgo.push_back(bgo_corr[1]*xbgo[1]);
    xPosW_bgo.push_back(bgo_corr[2]*xbgo[2]);
    xPosW_bgo.push_back(bgo_corr[3]*xbgo[3]);
    xPosW_bgo.push_back(bgo_corr[4]*xbgo[4]);
    xPosW_bgo.push_back(bgo_corr[5]*xbgo[5]);
    xPosW_bgo.push_back(bgo_corr[6]*xbgo[6]);
    xPosW_bgo.push_back(bgo_corr[7]*xbgo[7]);
    
    yPosW_bgo.push_back(bgo_corr[0]*ybgo[0]);
    yPosW_bgo.push_back(bgo_corr[1]*ybgo[1]);
    yPosW_bgo.push_back(bgo_corr[2]*ybgo[2]);
    yPosW_bgo.push_back(bgo_corr[3]*ybgo[3]);
    yPosW_bgo.push_back(bgo_corr[4]*ybgo[4]);
    yPosW_bgo.push_back(bgo_corr[5]*ybgo[5]);
    yPosW_bgo.push_back(bgo_corr[6]*ybgo[6]);
    yPosW_bgo.push_back(bgo_corr[7]*ybgo[7]);


    float eTot_bgo_corr = 0.;
    for( unsigned i=0; i<BGO_CHANNELS; ++i )
      eTot_bgo_corr += bgo_corr[i];

    float eTot_cef3_corr = 0.;
    for( unsigned i=0; i<CEF3_CHANNELS; ++i )
      eTot_cef3_corr += cef3_corr[i];


    float xPosBGO_wa = sumVector(xPosW_bgo)/eTot_bgo_corr;
    float yPosBGO_wa = sumVector(yPosW_bgo)/eTot_bgo_corr;


    if( eTot_cef3_corr<=0. ) continue;

    if( xTarget!=0. ) {
      float calibX = ( xPosBGO_wa/xTarget - eTot_bgo_corr )/eTot_cef3_corr;
      hp_cef3CalibX->Fill( xTarget, calibX );
    }

    if( yTarget!=0. ) {
      float calibY = ( yPosBGO_wa/yTarget - eTot_bgo_corr )/eTot_cef3_corr;
      hp_cef3CalibY->Fill( yTarget, calibY );
    }
      

  }

  //std::cout << "X: " << h1_cef3CalibX->GetMean() << " +/- " << h1_cef3CalibX->GetMeanError() << std::endl;
  //std::cout << "Y: " << h1_cef3CalibY->GetMean() << " +/- " << h1_cef3CalibY->GetMeanError() << std::endl;
  

  outfile->cd();
  hp_cef3CalibX->Write();
  hp_cef3CalibY->Write();

  outfile->Close();


  //TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  //c1->cd();
  //
  //h1_cef3CalibX->Draw();
  //h1_cef3CalibY->SetLineColor(kRed);
  //h1_cef3CalibY->Draw("same");

  //c1->SaveAs("prova.eps");


  return 0;

}



float sumVector( std::vector<float> v ) {

  float sum=0.;
  for( unsigned i=0; i<v.size(); ++i ) sum += v[i];

  return sum;

}
