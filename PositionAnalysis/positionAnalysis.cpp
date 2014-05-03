#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <assert.h>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector2.h"

#include "fastDQM_CeF3_BTF.h"
#include "interface/HodoCluster.h"
#include "interface/RunHelper.h"


#include "hodo_efficiency.dat"

std::vector< std::pair<float, float> > getPedestals( const std::string& type, const std::string& fileName );
std::vector<float> subtractPedestals( std::vector<float> raw, std::vector< std::pair<float, float> > pedestals, float nSigma );
float sumVector( std::vector<float> v );
bool checkVector( std::vector<float> v, float theMax=4095. );

float getMeanposHodo( std::vector<HodoCluster*> clusters );
std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo_corr, int nClusterMax );
void getCeF3Position( std::vector<float> cef3, float& xPos, float& yPos );
float getSingleCef3Position( float en, bool takemin=false );
float gethodointercalib(TString axis, int n);








int main( int argc, char* argv[] ) {


  std::string runName = "precalib_BGO_pedestal_noSource";
  if( argc>1 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
  }

  std::string fileName = "data/run_" + runName + ".root";
  TFile* file = TFile::Open(fileName.c_str());
  if( file==0 ) {
    std::cout << "ERROR! Din't find file " << fileName << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  TTree* tree = (TTree*)file->Get("eventRawData");


  std::string pedestalFileName = "data/run_BTF_91_20140430-015540_pedestal.root";
  std::cout << "-> Using pedestal file: " << pedestalFileName << std::endl;

  std::vector<std::pair<float, float> > pedestals = getPedestals( "cef3", pedestalFileName );
  std::cout << std::endl;
  std::cout << "-> Got pedestals of CeF3: " << std::endl;
  for( unsigned i=0; i<CEF3_CHANNELS; ++i )
    std::cout << " CeF3 Channel " << i << ": " << pedestals[i].first << " (+- " << pedestals[i].second << ")" << std::endl;
  std::cout << std::endl;

  std::vector<std::pair<float, float> > pedestals_bgo = getPedestals( "bgo", pedestalFileName );
  std::cout << std::endl;
  std::cout << "-> Got pedestals of BGO: " << std::endl;
  for( unsigned i=0; i<BGO_CHANNELS; ++i )
    std::cout << " BGO Channel " << i << ": " << pedestals_bgo[i].first << " (+- " << pedestals_bgo[i].second << ")" << std::endl;
  std::cout << std::endl;

  std::vector<std::pair<float, float> > pedestals_hodox = getPedestals("hodox", pedestalFileName);
  std::vector<std::pair<float, float> > pedestals_hodoy = getPedestals("hodoy", pedestalFileName);
  std::cout << "-> Got Hodoscope pedestals: " << std::endl;
  std::cout << std::endl;
  for( unsigned i=0; i<HODOX_CHANNELS; ++i )
    std::cout << "Channel " << i << ":  X: " << pedestals_hodox[i].first << " (+- " << pedestals_hodox[i].second << ") Y: " << pedestals_hodoy[i].first << " (+- " << pedestals_hodoy[i].second << ")" << std::endl;

  std::cout << std::endl << std::endl;




  // first: BGO precalibration*calibration:

  std::vector<float> bgo_calibration;

  // good, dont touch
  //bgo_calibration.push_back(1.24251419*1.02018 ); // 0
  //bgo_calibration.push_back(0.78909836*0.91473 ); // 1
  //bgo_calibration.push_back(0.91233889*1.07665 ); // 2
  //bgo_calibration.push_back(0.81254220*1.07697 ); // 3
  //bgo_calibration.push_back(1.19382440*0.956594); // 4
  //bgo_calibration.push_back(1.23486403*0.992715); // 5
  //bgo_calibration.push_back(1.05052378*0.996411); // 6
  //bgo_calibration.push_back(0.99823724*0.987254); // 7

  // test
  bgo_calibration.push_back(1.24251419*1.02018 ); // 0
  bgo_calibration.push_back(0.78909836*0.91473 ); // 1
  bgo_calibration.push_back(0.91233889*1.07665 ); // 2
  bgo_calibration.push_back(0.81254220*1.07697/0.73 ); // 3
  bgo_calibration.push_back(1.19382440*0.956594*0.5); // 4
  bgo_calibration.push_back(1.23486403*0.992715); // 5
  bgo_calibration.push_back(1.05052378*0.996411); // 6
  bgo_calibration.push_back(0.99823724*0.987254); // 7

  //float bgoCalibrationAverage = sumVector(bgo_calibration)/bgo_calibration.size();



  // CeF3 calibration:
  std::vector<float> cef3_calibration;
  cef3_calibration.push_back(848.317);
  cef3_calibration.push_back(995.703);
  cef3_calibration.push_back(891.56 );
  cef3_calibration.push_back(928.443);

  float cef3CalibrationAverage = sumVector(cef3_calibration)/cef3_calibration.size();

  for(unsigned i=0; i<cef3_calibration.size(); ++i )
    cef3_calibration[i] = cef3CalibrationAverage/cef3_calibration[i];





  std::vector<float> xbgo, ybgo;
  for( unsigned i=0; i<BGO_CHANNELS; ++i ) {
    float x,y;
    RunHelper::getBGOCoordinates( i, x, y );
    xbgo.push_back( x );
    ybgo.push_back( y );
  }


  UInt_t evtNumber;
  tree->SetBranchAddress( "evtNumber", &evtNumber );
  UInt_t adcData[40];
  tree->SetBranchAddress( "adcData", adcData );
  UInt_t adcBoard[40];
  tree->SetBranchAddress( "adcBoard", adcBoard );
  UInt_t adcChannel[40];
  tree->SetBranchAddress( "adcChannel", adcChannel );

  float xySize = 25.; // in mm

  int nBins = 500;
  float xMax = xySize*3./2.;
  int nHodoFibersX;
  int nHodoFibersY;
  int nHodoClustersX;
  int nHodoClustersY;


  TH1D* h1_xPos = new TH1D("xPos", "", nBins, -xMax, xMax);
  TH1D* h1_yPos = new TH1D("yPos", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos = new TH2D("xyPos", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_xPos_singleEle = new TH1D("xPos_singleEle", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_singleEle = new TH1D("yPos_singleEle", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_singleEle = new TH2D("xyPos_singleEle", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_xPos_new = new TH1D("xPos_new", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_new = new TH1D("yPos_new", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_new = new TH2D("xyPos_new", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_xPos_new_singleEle = new TH1D("xPos_new_singleEle", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_new_singleEle = new TH1D("yPos_new_singleEle", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_new_singleEle = new TH2D("xyPos_new_singleEle", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_cef3_0   = new TH1D("cef3_0",   "", 5000, 0., 5000.);
  TH1D* h1_cef3_1   = new TH1D("cef3_1",   "", 5000, 0., 5000.);
  TH1D* h1_cef3_2   = new TH1D("cef3_2",   "", 5000, 0., 5000.);
  TH1D* h1_cef3_3   = new TH1D("cef3_3",   "", 5000, 0., 5000.);
  TH1D* h1_cef3_tot = new TH1D("cef3_tot", "", 200, 0., 4.*5000.);
  
  TH1D* h1_cef3_corr_0   = new TH1D("cef3_corr_0",   "", 5000, 0., 5000.);
  TH1D* h1_cef3_corr_1   = new TH1D("cef3_corr_1",   "", 5000, 0., 5000.);
  TH1D* h1_cef3_corr_2   = new TH1D("cef3_corr_2",   "", 5000, 0., 5000.);
  TH1D* h1_cef3_corr_3   = new TH1D("cef3_corr_3",   "", 5000, 0., 5000.);
  TH1D* h1_cef3_corr_tot = new TH1D("cef3_corr_tot", "", 200, 0., 4.*5000.);

  TH1D* h1_bgo_corr_0   = new TH1D("bgo_corr_0",   "", 5000, 0., 5000.);
  TH1D* h1_bgo_corr_1   = new TH1D("bgo_corr_1",   "", 5000, 0., 5000.);
  TH1D* h1_bgo_corr_2   = new TH1D("bgo_corr_2",   "", 5000, 0., 5000.);
  TH1D* h1_bgo_corr_3   = new TH1D("bgo_corr_3",   "", 5000, 0., 5000.);
  TH1D* h1_bgo_corr_4   = new TH1D("bgo_corr_4",   "", 5000, 0., 5000.);
  TH1D* h1_bgo_corr_5   = new TH1D("bgo_corr_5",   "", 5000, 0., 5000.);
  TH1D* h1_bgo_corr_6   = new TH1D("bgo_corr_6",   "", 5000, 0., 5000.);
  TH1D* h1_bgo_corr_7   = new TH1D("bgo_corr_7",   "", 5000, 0., 5000.);
  TH1D* h1_bgo_corr_tot   = new TH1D("bgo_corr_tot",   "", 5000, 0., 8.*5000.);

  TH1D* h1_xPos_bgo = new TH1D("xPos_bgo", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_bgo = new TH1D("yPos_bgo", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_bgo = new TH2D("xyPos_bgo", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_xPos_singleEle_bgo = new TH1D("xPos_singleEle_bgo", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_singleEle_bgo = new TH1D("yPos_singleEle_bgo", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_singleEle_bgo = new TH2D("xyPos_singleEle_bgo", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_xPos_hodo = new TH1D("xPos_hodo", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_hodo = new TH1D("yPos_hodo", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_hodo = new TH2D("xyPos_hodo", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_xPos_hodoClust = new TH1D("xPos_hodoClust", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_hodoClust = new TH1D("yPos_hodoClust", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_hodoClust = new TH2D("xyPos_hodoClust", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_xPos_singleEle_hodo = new TH1D("xPos_singleEle_hodo", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_singleEle_hodo = new TH1D("yPos_singleEle_hodo", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_singleEle_hodo = new TH2D("xyPos_singleEle_hodo", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_xPos_singleEle_hodoClust = new TH1D("xPos_singleEle_hodoClust", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_singleEle_hodoClust = new TH1D("yPos_singleEle_hodoClust", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_singleEle_hodoClust = new TH2D("xyPos_singleEle_hodoClust", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_xPos_calo = new TH1D("xPos_calo", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_calo = new TH1D("yPos_calo", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_calo = new TH2D("xyPos_calo", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_xPos_singleEle_calo = new TH1D("xPos_singleEle_calo", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_singleEle_calo = new TH1D("yPos_singleEle_calo", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_singleEle_calo = new TH2D("xyPos_singleEle_calo", "", nBins, -xMax, xMax, nBins, -xMax, xMax);



  TH2D* h2_correlation_cef3_hodo_xPos = new TH2D("correlation_cef3_hodo_xPos", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);
  TH2D* h2_correlation_cef3_hodo_yPos = new TH2D("correlation_cef3_hodo_yPos", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);

  TH2D* h2_correlation_cef3_bgo_xPos = new TH2D("correlation_cef3_bgo_xPos", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);
  TH2D* h2_correlation_cef3_bgo_yPos = new TH2D("correlation_cef3_bgo_yPos", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);

  TH2D* h2_correlation_hodo_bgo_xPos = new TH2D("correlation_hodo_bgo_xPos", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);
  TH2D* h2_correlation_hodo_bgo_yPos = new TH2D("correlation_hodo_bgo_yPos", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);


  TH2D* h2_correlation_cef3_hodo_xPos_singleEle = new TH2D("correlation_cef3_hodo_xPos_singleEle", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);
  TH2D* h2_correlation_cef3_hodo_yPos_singleEle = new TH2D("correlation_cef3_hodo_yPos_singleEle", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);

  TH2D* h2_correlation_cef3_bgo_xPos_singleEle = new TH2D("correlation_cef3_bgo_xPos_singleEle", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);
  TH2D* h2_correlation_cef3_bgo_yPos_singleEle = new TH2D("correlation_cef3_bgo_yPos_singleEle", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);

  TH2D* h2_correlation_hodo_bgo_xPos_singleEle = new TH2D("correlation_hodo_bgo_xPos_singleEle", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);
  TH2D* h2_correlation_hodo_bgo_yPos_singleEle = new TH2D("correlation_hodo_bgo_yPos_singleEle", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);



  int nentries = tree->GetEntries();

  std::string outfileName = "PosAn_" + runName + ".root";
  TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );

  TTree* outTree = new TTree("tree_passedEvents","tree_passedEvents");
  float cef3_[CEF3_CHANNELS],bgo_[BGO_CHANNELS],hodox_[HODOX_CHANNELS],hodoy_[HODOY_CHANNELS];
  float cef3_corr_[CEF3_CHANNELS],bgo_corr_[BGO_CHANNELS],hodox_corr_[HODOX_CHANNELS],hodoy_corr_[HODOY_CHANNELS];
  float scintFront_;
  int cef3_chan=CEF3_CHANNELS;
  int bgo_chan=BGO_CHANNELS;
  int hodox_chan=HODOX_CHANNELS;
  int hodoy_chan=HODOY_CHANNELS; 

  outTree->Branch( "evtNumber", &evtNumber,"evtNumber/F" );
  outTree->Branch( "adcData", adcData, "adcData/i" );
  outTree->Branch( "adcBoard", adcBoard, "adcBoard/i" );
  outTree->Branch( "adcChannel", adcChannel,"adcChannel/i" );
  outTree->Branch( "nHodoFibersX", &nHodoFibersX, "nHodoFibersX/I" );
  outTree->Branch( "nHodoFibersY", &nHodoFibersY, "nHodoFibersY/I" );
  outTree->Branch( "nHodoClustersX", &nHodoClustersX, "nHodoClustersX/I" );
  outTree->Branch( "nHodoClustersY", &nHodoClustersY, "nHodoClustersY/I" );
  outTree->Branch( "hodox_chan", &hodox_chan, "hodox_chan/I" );
  outTree->Branch( "hodoy_chan", &hodoy_chan, "hodoy_chan/I" );
  outTree->Branch( "cef3_chan", &cef3_chan, "cef3_chan/I" );
  outTree->Branch( "bgo_chan", &bgo_chan, "bgo_chan/I" );
  outTree->Branch( "hodox", hodox_, "hodox_[hodox_chan]/F" );
  outTree->Branch( "hodoy", hodoy_, "hodoy_[hodoy_chan]/F" );
  outTree->Branch( "bgo", bgo_, "bgo_[bgo_chan]/F" );
  outTree->Branch( "cef3", cef3_, "cef3_[cef3_chan]/F" );
  outTree->Branch( "bgo_corr", bgo_corr_, "bgo_corr_[bgo_chan]/F" );
  outTree->Branch( "cef3_corr", cef3_corr_, "cef3_corr_[cef3_chan]/F" );
  outTree->Branch( "hodox_corr", hodox_corr_, "hodox_corr_[hodox_chan]/F" );
  outTree->Branch( "hodoy_corr", hodoy_corr_, "hodoy_corr_[hodoy_chan]/F" );
  outTree->Branch( "scintFront", &scintFront_, "scintFront_/F" );


  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    nHodoFibersX=0;
    nHodoFibersY=0;
    nHodoClustersX=0;
    nHodoClustersY=0;

    tree->GetEntry(iEntry);

    if( iEntry % 5000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;


    // CeF3 fibres
    std::vector<float> cef3;
    for( unsigned i=0; i<CEF3_CHANNELS; ++i ) cef3.push_back(-1.);

    // BGO
    std::vector<float> bgo;
    for( unsigned i=0; i<BGO_CHANNELS; ++i ) bgo.push_back(-1.);

    // hodoX
    std::vector<float> hodox;
    for( unsigned i=0; i<HODOX_CHANNELS; ++i ) hodox.push_back(-1.);

    // hodoY
    std::vector<float> hodoy;
    for( unsigned i=0; i<HODOY_CHANNELS; ++i ) hodoy.push_back(-1.);


    for( unsigned i=0; i<40; ++i ) {

      int board  = adcBoard[i];
      int channel= adcChannel[i];

      if( board==BGO_ADC_BOARD ) {
             if( channel==(BGO_ADC_START_CHANNEL  ) )  bgo[0] = adcData[i];
        else if( channel==(BGO_ADC_START_CHANNEL+1) )  bgo[1] = adcData[i];
        else if( channel==(BGO_ADC_START_CHANNEL+2) )  bgo[2] = adcData[i];
        else if( channel==(BGO_ADC_START_CHANNEL+3) )  bgo[3] = adcData[i];
        else if( channel==(BGO_ADC_START_CHANNEL+4) )  bgo[4] = adcData[i];
        else if( channel==(BGO_ADC_START_CHANNEL+5) )  bgo[5] = adcData[i];
        else if( channel==(BGO_ADC_START_CHANNEL+6) )  bgo[6] = adcData[i];
        else if( channel==(BGO_ADC_START_CHANNEL+7) )  bgo[7] = adcData[i];
      }
      if( board==CEF3_ADC_BOARD ) {
             if( channel==(CEF3_ADC_START_CHANNEL  ) )  cef3[0]  = adcData[i];
        else if( channel==(CEF3_ADC_START_CHANNEL+1) )  cef3[1]  = adcData[i];
        else if( channel==(CEF3_ADC_START_CHANNEL+2) )  cef3[2]  = adcData[i];
        else if( channel==(CEF3_ADC_START_CHANNEL+3) )  cef3[3]  = adcData[i];
      }
      if( board==HODOX_ADC_BOARD ) {
             if( channel==(HODOX_ADC_START_CHANNEL  ) ) hodox[0] = adcData[i];
        else if( channel==(HODOX_ADC_START_CHANNEL+1) ) hodox[1] = adcData[i];
        else if( channel==(HODOX_ADC_START_CHANNEL+2) ) hodox[2] = adcData[i];
        else if( channel==(HODOX_ADC_START_CHANNEL+3) ) hodox[3] = adcData[i];
        else if( channel==(HODOX_ADC_START_CHANNEL+4) ) hodox[4] = adcData[i];
        else if( channel==(HODOX_ADC_START_CHANNEL+5) ) hodox[5] = adcData[i];
        else if( channel==(HODOX_ADC_START_CHANNEL+6) ) hodox[6] = adcData[i];
        else if( channel==(HODOX_ADC_START_CHANNEL+7) ) hodox[7] = adcData[i];
      }
      if( board==HODOY_ADC_BOARD ) { 
             if( channel==(HODOY_ADC_START_CHANNEL  ) ) hodoy[0] = adcData[i];
        else if( channel==(HODOY_ADC_START_CHANNEL+1) ) hodoy[1] = adcData[i];
        else if( channel==(HODOY_ADC_START_CHANNEL+2) ) hodoy[2] = adcData[i];
        else if( channel==(HODOY_ADC_START_CHANNEL+3) ) hodoy[3] = adcData[i];
        else if( channel==(HODOY_ADC_START_CHANNEL+4) ) hodoy[4] = adcData[i];
        else if( channel==(HODOY_ADC_START_CHANNEL+5) ) hodoy[5] = adcData[i];
        else if( channel==(HODOY_ADC_START_CHANNEL+6) ) hodoy[6] = adcData[i];
        else if( channel==(HODOY_ADC_START_CHANNEL+7) ) hodoy[7] = adcData[i];
      }
      if( board==SCINT_FRONT_ADC_BOARD && channel==SCINT_FRONT_ADC_CHANNEL ) scintFront_ = adcData[i];
    }


    float nSigma_hodo = 4.;
    std::vector<float>  bgo_corr = subtractPedestals( bgo , pedestals_bgo, 4. );
    std::vector<float> cef3_corr = subtractPedestals( cef3, pedestals,     4. );
    std::vector<float> hodox_corr = subtractPedestals( hodox, pedestals_hodox, nSigma_hodo );
    std::vector<float> hodoy_corr = subtractPedestals( hodoy, pedestals_hodoy, nSigma_hodo );

    bool cef3_ok = checkVector(cef3);
    bool cef3_corr_ok = checkVector(cef3_corr);



    // FIRST GET POSITION FROM HODOSCOPE:

    int clusterSize=2;
    std::vector<HodoCluster*> hodoxFibres   = getHodoClusters( hodox_corr, 1 ); // fibres are just clusters with size = 1
    std::vector<HodoCluster*> hodoxClusters = getHodoClusters( hodox_corr, clusterSize );

    std::vector<HodoCluster*> hodoyFibres   = getHodoClusters( hodoy_corr, 1 );
    std::vector<HodoCluster*> hodoyClusters = getHodoClusters( hodoy_corr, clusterSize );

    nHodoFibersX = hodoxFibres.size();
    nHodoFibersY = hodoyFibres.size();

    nHodoClustersX = hodoxClusters.size();
    nHodoClustersY = hodoyClusters.size();

    float xPos_hodo = getMeanposHodo(hodoxFibres);
    float yPos_hodo = getMeanposHodo(hodoyFibres);

    float xPos_hodoClust = getMeanposHodo(hodoxClusters);
    float yPos_hodoClust = getMeanposHodo(hodoyClusters);


    //bool isSingleElectron = ((nHodoFibersX==1) && (nHodoFibersY==1));
    bool isSingleElectron = (scintFront_>500. && scintFront_<2000.);

    if( xPos_hodo>-100. )
      h1_xPos_hodo->Fill(xPos_hodo);
    if( yPos_hodo>-100. )
      h1_yPos_hodo->Fill(yPos_hodo);

    if( xPos_hodo>-100. && yPos_hodo>-100. ) 
      h2_xyPos_hodo->Fill(xPos_hodo, yPos_hodo);

    if( xPos_hodoClust>-100. )
      h1_xPos_hodoClust->Fill(xPos_hodoClust);
    if( yPos_hodoClust>-100. )
      h1_yPos_hodoClust->Fill(yPos_hodoClust);

    if( xPos_hodoClust>-100. && yPos_hodoClust>-100. ) 
      h2_xyPos_hodoClust->Fill(xPos_hodoClust, yPos_hodoClust);


    if( isSingleElectron ) {

      if( xPos_hodo>-100. )
        h1_xPos_singleEle_hodo->Fill(xPos_hodo);
      if( yPos_hodo>-100. )
        h1_yPos_singleEle_hodo->Fill(yPos_hodo);

      if( xPos_hodo>-100. && yPos_hodo>-100. ) 
        h2_xyPos_singleEle_hodo->Fill(xPos_hodo, yPos_hodo);

      if( xPos_hodoClust>-100. )
        h1_xPos_singleEle_hodoClust->Fill(xPos_hodoClust);
      if( yPos_hodoClust>-100. )
        h1_yPos_singleEle_hodoClust->Fill(yPos_hodoClust);

      if( xPos_hodoClust>-100. && yPos_hodoClust>-100. ) 
        h2_xyPos_singleEle_hodoClust->Fill(xPos_hodoClust, yPos_hodoClust);

    }



    bool bgo_ok = checkVector(bgo, 4095.);
    bool bgo_corr_ok = checkVector(bgo_corr, 4095.);

    float xPos_bgo;
    float yPos_bgo;

    std::vector<float> xPosW_bgo;
    std::vector<float> yPosW_bgo;

    float eTot_bgo_corr;

    if( bgo_ok && bgo_corr_ok ) {


      for(unsigned i=0; i<bgo_calibration.size(); ++i ) 
        bgo_corr[i] *= bgo_calibration[i]; //correct


      eTot_bgo_corr  = sumVector(bgo_corr);

      h1_bgo_corr_0->Fill( bgo_corr[0] );
      h1_bgo_corr_1->Fill( bgo_corr[1] );
      h1_bgo_corr_2->Fill( bgo_corr[2] );
      h1_bgo_corr_3->Fill( bgo_corr[3] );
      h1_bgo_corr_4->Fill( bgo_corr[4] );
      h1_bgo_corr_5->Fill( bgo_corr[5] );
      h1_bgo_corr_6->Fill( bgo_corr[6] );
      h1_bgo_corr_7->Fill( bgo_corr[7] );
      h1_bgo_corr_tot->Fill( eTot_bgo_corr );



      // then proceed to compute position:

      //   0  1  2
      //   3     4
      //   5  6  7

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
      

      xPos_bgo = sumVector( xPosW_bgo )/eTot_bgo_corr;
      yPos_bgo = sumVector( yPosW_bgo )/eTot_bgo_corr;
      
      h1_xPos_bgo->Fill( xPos_bgo );
      h1_yPos_bgo->Fill( yPos_bgo );
      h2_xyPos_bgo->Fill( xPos_bgo, yPos_bgo );

      h2_correlation_hodo_bgo_xPos->Fill( xPos_hodo, xPos_bgo );
      h2_correlation_hodo_bgo_yPos->Fill( yPos_hodo, yPos_bgo );
      
      if( isSingleElectron ) {

        h1_xPos_singleEle_bgo->Fill( xPos_bgo );
        h1_yPos_singleEle_bgo->Fill( yPos_bgo );
        h2_xyPos_singleEle_bgo->Fill( xPos_bgo, yPos_bgo );

        h2_correlation_hodo_bgo_xPos_singleEle->Fill( xPos_hodo, xPos_bgo );
        h2_correlation_hodo_bgo_yPos_singleEle->Fill( yPos_hodo, yPos_bgo );

      }
      
    }  // if bgo ok



    // THEN USE CeF3 DATA:

    if( cef3_ok ) {

      // intercalibration:

      for(unsigned i=0; i<cef3_calibration.size(); ++i )
        cef3_corr[i] *= cef3_calibration[i]; //correct


      float eTot      = sumVector(cef3);
      float eTot_corr = sumVector(cef3_corr);


      h1_cef3_0->Fill( cef3[0] );
      h1_cef3_1->Fill( cef3[1] );
      h1_cef3_2->Fill( cef3[2] );
      h1_cef3_3->Fill( cef3[3] );
      h1_cef3_tot->Fill( eTot );

      

      if( cef3_corr_ok ) {

        h1_cef3_corr_0->Fill( cef3_corr[0] );
        h1_cef3_corr_1->Fill( cef3_corr[1] );
        h1_cef3_corr_2->Fill( cef3_corr[2] );
        h1_cef3_corr_3->Fill( cef3_corr[3] );
        h1_cef3_corr_tot->Fill( eTot_corr );


        //   0      1
        //          
        //          
        //   3      2


        //float chamfer = 2.1; // in mm
        //float position = xySize/2. - chamfer/4.;
        float position = 12. - 0.696; // using FN's infallible trigonometry

        //std::vector

        std::vector<float> xPosW;
        xPosW.push_back(cef3_corr[0]*(-position));
        xPosW.push_back(cef3_corr[1]*(+position));
        xPosW.push_back(cef3_corr[2]*(+position));
        xPosW.push_back(cef3_corr[3]*(-position));

        std::vector<float> yPosW;
        yPosW.push_back(cef3_corr[0]*(+position));
        yPosW.push_back(cef3_corr[1]*(+position));
        yPosW.push_back(cef3_corr[2]*(-position));
        yPosW.push_back(cef3_corr[3]*(-position));


        float xPos_new, yPos_new;
        getCeF3Position( cef3_corr, xPos_new, yPos_new );

        float xPos = sumVector(xPosW)/eTot_corr;
        float yPos = sumVector(yPosW)/eTot_corr;

        h1_xPos->Fill( xPos );
        h1_yPos->Fill( yPos );

        h2_xyPos->Fill( xPos, yPos );

        h1_xPos_new->Fill( xPos_new );
        h1_yPos_new->Fill( yPos_new );

        h2_xyPos_new->Fill( xPos_new, yPos_new );


        // positioning with all 9 calorimeter channels:
        float xPos_calo = sumVector( xPosW_bgo )/(eTot_bgo_corr + eTot_corr*0.1); // cef3 is in 0,0
        float yPos_calo = sumVector( yPosW_bgo )/(eTot_bgo_corr + eTot_corr*0.1); // so counts only in denominator
        //float xPos_calo = sumVector( xPosW_bgo )/(eTot_bgo_corr + eTot_corr*0.791577); // cef3 is in 0,0
        //float yPos_calo = sumVector( yPosW_bgo )/(eTot_bgo_corr + eTot_corr*0.791577); // so counts only in denominator

  
        if( bgo_ok && bgo_corr_ok ) {

          h1_xPos_calo->Fill( xPos_calo );
          h1_yPos_calo->Fill( yPos_calo );
          h2_xyPos_calo->Fill( xPos_calo, yPos_calo );


          // CORRELATIONS BETWEEN CALO AND HODO:

          h2_correlation_cef3_bgo_xPos->Fill( xPos, xPos_bgo );
          h2_correlation_cef3_bgo_yPos->Fill( yPos, yPos_bgo );

        }

        if( isSingleElectron) {

          h1_xPos_singleEle->Fill( xPos );
          h1_yPos_singleEle->Fill( yPos );
  
          h2_xyPos_singleEle->Fill( xPos, yPos );

          h1_xPos_new_singleEle->Fill( xPos_new );
          h1_yPos_new_singleEle->Fill( yPos_new );
  
          h2_xyPos_new_singleEle->Fill( xPos_new, yPos_new );

          h2_correlation_cef3_hodo_xPos_singleEle->Fill( xPos, xPos_hodo );
          h2_correlation_cef3_hodo_yPos_singleEle->Fill( yPos, yPos_hodo );

  
          if( bgo_ok && bgo_corr_ok ) {

            h1_xPos_singleEle_calo->Fill( xPos_calo );
            h1_yPos_singleEle_calo->Fill( yPos_calo );
            h2_xyPos_singleEle_calo->Fill( xPos_calo, yPos_calo );

            h2_correlation_cef3_bgo_xPos_singleEle->Fill( xPos, xPos_bgo );
            h2_correlation_cef3_bgo_yPos_singleEle->Fill( yPos, yPos_bgo );

          }

        }

        for(int i=0;i<CEF3_CHANNELS;i++){
          cef3_[i]=cef3[i]; 
          cef3_corr_[i]=cef3_corr[i];
        }
        for(int i=0;i<BGO_CHANNELS;i++){
          bgo_[i]=bgo[i];
          bgo_corr_[i]=bgo_corr[i];
        }
        for(int i=0;i<HODOX_CHANNELS;i++){
          hodox_[i]=hodox[i];
          hodox_corr_[i]=hodox_corr[i];
        }
        for(int i=0;i<HODOY_CHANNELS;i++){
          hodoy_[i]=hodoy[i];
          hodoy_corr_[i]=hodoy_corr[i];
        }
        
        outTree->Fill();


      } // if cef3_ok

    }

  }


  std::cout << "-> Events passing overflow cut: " << h1_xPos->GetEntries() << "/" << nentries << " (" << 100.* h1_xPos->GetEntries()/nentries << "%)" << std::endl;

  outfile->cd();

  outTree->Write();

  h1_xPos->Write();
  h1_yPos->Write();
  h2_xyPos->Write();

  h1_xPos_singleEle->Write();
  h1_yPos_singleEle->Write();
  h2_xyPos_singleEle->Write();

  h1_xPos_new->Write();
  h1_yPos_new->Write();
  h2_xyPos_new->Write();

  h1_xPos_new_singleEle->Write();
  h1_yPos_new_singleEle->Write();
  h2_xyPos_new_singleEle->Write();

  h1_cef3_0->Write();
  h1_cef3_1->Write();
  h1_cef3_2->Write();
  h1_cef3_3->Write();
  
  h1_cef3_corr_0->Write();
  h1_cef3_corr_1->Write();
  h1_cef3_corr_2->Write();
  h1_cef3_corr_3->Write();
  
  h1_cef3_tot->Write();
  h1_cef3_corr_tot->Write();
  
  h1_xPos_bgo->Write();
  h1_yPos_bgo->Write();
  h2_xyPos_bgo->Write();
  
  h1_xPos_hodo->Write();
  h1_yPos_hodo->Write();
  h2_xyPos_hodo->Write();

  h1_xPos_hodoClust->Write();
  h1_yPos_hodoClust->Write();
  h2_xyPos_hodoClust->Write();

  h2_correlation_cef3_hodo_xPos->Write();
  h2_correlation_cef3_hodo_yPos->Write();

  h2_correlation_cef3_bgo_xPos->Write();
  h2_correlation_cef3_bgo_yPos->Write();

  h2_correlation_hodo_bgo_xPos->Write();
  h2_correlation_hodo_bgo_yPos->Write();

  h1_xPos_singleEle_bgo->Write();
  h1_yPos_singleEle_bgo->Write();
  h2_xyPos_singleEle_bgo->Write();
  
  h1_xPos_singleEle_hodo->Write();
  h1_yPos_singleEle_hodo->Write();
  h2_xyPos_singleEle_hodo->Write();

  h1_xPos_singleEle_hodoClust->Write();
  h1_yPos_singleEle_hodoClust->Write();
  h2_xyPos_singleEle_hodoClust->Write();

  h1_xPos_calo->Write();
  h1_yPos_calo->Write();
  h2_xyPos_calo->Write();
  
  h1_xPos_singleEle_calo->Write();
  h1_yPos_singleEle_calo->Write();
  h2_xyPos_singleEle_calo->Write();
  

  h2_correlation_cef3_hodo_xPos_singleEle->Write();
  h2_correlation_cef3_hodo_yPos_singleEle->Write();

  h2_correlation_cef3_bgo_xPos_singleEle->Write();
  h2_correlation_cef3_bgo_yPos_singleEle->Write();

  h2_correlation_hodo_bgo_xPos_singleEle->Write();
  h2_correlation_hodo_bgo_yPos_singleEle->Write();

  h1_bgo_corr_0->Write();
  h1_bgo_corr_1->Write();
  h1_bgo_corr_2->Write();
  h1_bgo_corr_3->Write();
  h1_bgo_corr_4->Write();
  h1_bgo_corr_5->Write();
  h1_bgo_corr_6->Write();
  h1_bgo_corr_7->Write();
  h1_bgo_corr_tot->Write();
  
  outfile->Close();
  std::cout << "-> Histograms saved in: " << outfile->GetName() << std::endl;

  return 0;

}






std::vector< std::pair<float, float> > getPedestals( const std::string& type, const std::string& fileName ) {

  int nBoard=-1;
  int nChannels=-1;
  int firstChannel=-1;
  if( type=="cef3" ) {
    nBoard       = CEF3_ADC_BOARD;
    nChannels    = CEF3_CHANNELS;
    firstChannel = CEF3_ADC_START_CHANNEL;
  } else if( type=="bgo" ) {
    nBoard       = BGO_ADC_BOARD;
    nChannels    = BGO_CHANNELS;
    firstChannel = BGO_ADC_START_CHANNEL;
  } else if( type=="hodox" ) {
    nBoard       = HODOX_ADC_BOARD;
    nChannels    = HODOX_CHANNELS;
    firstChannel = HODOX_ADC_START_CHANNEL;
  } else if( type=="hodoy" ) {
    nBoard       = HODOY_ADC_BOARD;
    nChannels    = HODOY_CHANNELS;
    firstChannel = HODOY_ADC_START_CHANNEL;
  } else {
    std::cout << "ERROR! Unkown type '" << type << "'!" << std::endl;
    std::cout << "Don't know what pedestals you're looking for." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(77);
  }

    

  TFile* file = TFile::Open(fileName.c_str());
  TTree* tree = (TTree*)file->Get("eventRawData");

  std::vector< std::pair<float, float> > peds;
  for( unsigned i=0; i<nChannels; ++i ) {
    int iChannel = firstChannel+i;
    TH1D* h1_ped = new TH1D("ped", "", 500, 0., 500.);
    tree->Project( "ped", "adcData", Form("adcBoard==%d && adcChannel==%d", nBoard, iChannel) );
    std::pair<float, float>  thispair;
    thispair.first  = h1_ped->GetMean();
    thispair.second = h1_ped->GetRMS();
    peds.push_back(thispair);
    delete h1_ped;
  }

  return peds;

}


std::vector<float> subtractPedestals( std::vector<float> raw, std::vector< std::pair<float, float> > pedestals, float nSigma ) {

  std::vector<float> corr;

  for(unsigned i=0; i<raw.size(); ++i ) {

    float iCorr = ( raw[i] > pedestals[i].first + nSigma*pedestals[i].second ) ? (raw[i] - pedestals[i].first) : 0.;
    corr.push_back( iCorr );

  }

  return corr;

}


float sumVector( std::vector<float> v ) {

  float sum=0.;
  for( unsigned i=0; i<v.size(); ++i ) sum += v[i];

  return sum;

}


bool checkVector( std::vector<float> v, float theMax ) {

  bool returnBool = true;

  for( unsigned i=0; i<v.size(); ++i ) {
    if( v[i]<=0. ) return false;
    if( v[i]>=theMax ) return false;
  }

  return returnBool;

}



float getMeanposHodo( std::vector<HodoCluster*> clusters ) {

  if( clusters.size()==0 ) return -999.;

  float meanpos = 0.;

  for( unsigned i=0; i<clusters.size(); ++i ) {
    meanpos += clusters[i]->getPosition();
  }

  meanpos /= clusters.size();

  return meanpos;

}




float getMeanposHodo( std::vector<float> hodo_corr, int& nHodoFibers, int& nHodoClusters ) {

  // initialize
  nHodoFibers=0;
  nHodoClusters=0;

  float mean = 0.;
  float eTot = 0.;
  int nClustered  = 0;
  int nClusterMax = 3; // maximal number of clustered fibres

  for( unsigned i=0; i<hodo_corr.size(); ++i ) {

    if( hodo_corr[i] > 0.) {
      mean += -(i-3.5);
      eTot += 1.;
      nHodoFibers++;
      if( nClustered<=nClusterMax ) {
        if( nClustered==0 ) nHodoClusters+=1;
        nClustered+=1;
      } else {
        nClustered = 0;
      }

    } else { // whenever you find a hole: reset
      nClustered=0;
    }

  }

  return mean/eTot;

}

float gethodointercalib(TString axis, int n){

  float res=1;

  if (axis==TString("X")){
    if (n>=0 && n<=7) return hodo_efficiency_vector_X[n];
    else assert(false);
  }
  else if (axis==TString("Y")){
    if (n>=0 && n<=7) return hodo_efficiency_vector_Y[n];
    else assert(false);
  }
  else assert(false);
  
  return 1./res;
  
};


std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo_corr, int nClusterMax ) {

  std::vector<HodoCluster*> clusters;

  HodoCluster* currentCluster = new HodoCluster();

  for( unsigned i=0; i<hodo_corr.size(); ++i ) {

    if( hodo_corr[i] > 0.) { // hit

      if( currentCluster->getSize() < nClusterMax ) {

        currentCluster->addFibre( i );

      } else {

        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster();   // create a new one
        currentCluster->addFibre( i );        // get that fibre!

      }

    } else { // as soon as you find a hole
      
      if( currentCluster->getSize() > 0 ) {
     
        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster();   // create a new one

      }

    }


  } // for fibres


  if( currentCluster->getSize()>0 )
    clusters.push_back( currentCluster ); // store last cluster


  return clusters;

}



void getCeF3Position( std::vector<float> cef3, float& xPos, float& yPos ) {

  xPos=0.;
  yPos=0.;

  float offset02 = 0.;
  //float offset02 = 0.0170062;
  float offset13 = 0.0594743;
  float r02 = cef3[0]/cef3[2] - offset02;
  float r13 = cef3[1]/cef3[3] - offset13;
  float diag02 = (r02>1.) ? getSingleCef3Position( r02, false ) : -getSingleCef3Position( 1./r02, false );
  float diag13 = (r13>1.) ? getSingleCef3Position( r13, false ) : -getSingleCef3Position( 1./r13, false );

  TVector2 v( diag13, diag02 );
  float pi = 3.14159;
  float theta = pi/4.; // 45 degrees 
  TVector2 d = v.Rotate(theta);

  xPos = d.X();
  yPos = d.Y();

}


float getSingleCef3Position( float en, bool takemin ) {

  float c = 1. - en;
  //float b = 1.32340e-02;
  //float a = 1.45209e-03;
  //float b = 1.89382e-02;
  //float a = 1.50157e-03;
  float b = 1.40598e-02;
  float a = 1.82353e-03;

  float theSqrt = b*b - 4.*a*c;

  float x1 = ( theSqrt>0. ) ? (-b + sqrt( theSqrt ))/(2.*a) : 0.;
  float x2 = ( theSqrt>0. ) ? (-b - sqrt( theSqrt ))/(2.*a) : 0.;


  //float returnX = (takemin) ? TMath::Min(x1,x2) : TMath::Max(x1,x2);


  float returnX;

  if( takemin ) {
    if( fabs(x1)<fabs(x2) ) {
      returnX = x1;
    } else {
      returnX = x2;
    }
  } else {
    if( fabs(x1)<fabs(x2) ) {
      returnX = x2;
    } else {
      returnX = x1;
    }
  }

  return returnX;

}

