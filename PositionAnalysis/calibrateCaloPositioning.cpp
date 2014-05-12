#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <assert.h>

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

#include "fastDQM_CeF3_BTF.h"
#include "interface/HodoCluster.h"
#include "interface/RunHelper.h"

#include "hodo_efficiency.dat"




float sumVector( std::vector<float> v );



int main( int argc, char* argv[] ) {


  std::string runName = "precalib_BGO_pedestal_noSource";
  if( argc>1 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
  }

  std::string fileName = "PosAn_" + runName + ".root";
  TFile* file = TFile::Open(fileName.c_str());
  if( file==0 ) {
    std::cout << "ERROR! Din't find file " << fileName << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  TTree* tree = (TTree*)file->Get("tree_passedEvents");




  std::vector<float> xbgo, ybgo;
  for( unsigned i=0; i<BGO_CHANNELS; ++i ) {
    float x,y;
    RunHelper::getBGOCoordinates( i, x, y );
    xbgo.push_back( x );
    ybgo.push_back( y );
  }



  float cef3_corr[CEF3_CHANNELS];
  tree->SetBranchAddress( "cef3_corr", cef3_corr );
  float bgo_corr[BGO_CHANNELS];
  tree->SetBranchAddress( "bgo_corr", bgo_corr );
  float scintFront;
  tree->SetBranchAddress( "scintFront", &scintFront );
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





  TH1D* h1_cef3CalibX = new TH1D("cef3CalibX", "", 200, -3., 3.);
  TH1D* h1_cef3CalibY = new TH1D("cef3CalibY", "", 200, -3., 3.);




  int nentries = tree->GetEntries();






  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {


    tree->GetEntry(iEntry);

    if( iEntry % 5000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;


    bool isSingleElectron = (scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1);
    if( !isSingleElectron ) continue;
    

    // require that clusters have no more than two fibres each
    // for better position precision
    if( nFibres_hodoClustX[0]>2 || nFibres_hodoClustY[0]>2 ) continue;


    // with one X and one Y clusters in hodo
    // define true beam position as hodo cross center
    float xTarget = pos_hodoClustX[0];  
    float yTarget = pos_hodoClustY[0];  



   

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

    float calibX = ( sumVector( xPosW_bgo )/xTarget - eTot_bgo_corr ) / eTot_cef3_corr;
    float calibY = ( sumVector( yPosW_bgo )/yTarget - eTot_bgo_corr ) / eTot_cef3_corr;

    h1_cef3CalibX->Fill( calibX );
    h1_cef3CalibY->Fill( calibY );
      

  }


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  
  h1_cef3CalibX->Draw();
  h1_cef3CalibY->SetLineColor(kRed);
  h1_cef3CalibY->Draw("same");

  std::cout << "X: " << h1_cef3CalibX->GetMean() << " +/- " << h1_cef3CalibX->GetMeanError() << std::endl;
  std::cout << "Y: " << h1_cef3CalibY->GetMean() << " +/- " << h1_cef3CalibY->GetMeanError() << std::endl;
  
  c1->SaveAs("prova.eps");


  return 0;

}



float sumVector( std::vector<float> v ) {

  float sum=0.;
  for( unsigned i=0; i<v.size(); ++i ) sum += v[i];

  return sum;

}
