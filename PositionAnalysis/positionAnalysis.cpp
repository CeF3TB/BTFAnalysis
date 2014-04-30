#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"

#include "fastDQM_CeF3_BTF.h"





std::vector< std::pair<float, float> > getPedestals( const std::string& type, const std::string& fileName );
std::vector<float> subtractPedestals( std::vector<float> raw, std::vector< std::pair<float, float> > pedestals, float nSigma );
float sumVector( std::vector<float> v );
bool checkVector( std::vector<float> v, float theMax=4095. );
float getMeanposHodo( std::vector<float> hodo, std::vector< std::pair<float, float> > pedestals, float nSigma);



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

  TH1D* h1_xPos = new TH1D("xPos", "", nBins, -xMax, xMax);
  TH1D* h1_yPos = new TH1D("yPos", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos = new TH2D("xyPos", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_cef3_0   = new TH1D("cef3_0",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_1   = new TH1D("cef3_1",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_2   = new TH1D("cef3_2",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_3   = new TH1D("cef3_3",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_tot = new TH1D("cef3_tot", "", 200, 0., 4.*4097.);
  
  TH1D* h1_cef3_corr_0   = new TH1D("cef3_corr_0",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_corr_1   = new TH1D("cef3_corr_1",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_corr_2   = new TH1D("cef3_corr_2",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_corr_3   = new TH1D("cef3_corr_3",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_corr_tot = new TH1D("cef3_corr_tot", "", 200, 0., 4.*4097.);

  TH1D* h1_bgo_corr_0   = new TH1D("bgo_corr_0",   "", 4097, 0., 4097.);
  TH1D* h1_bgo_corr_1   = new TH1D("bgo_corr_1",   "", 4097, 0., 4097.);
  TH1D* h1_bgo_corr_2   = new TH1D("bgo_corr_2",   "", 4097, 0., 4097.);
  TH1D* h1_bgo_corr_3   = new TH1D("bgo_corr_3",   "", 4097, 0., 4097.);
  TH1D* h1_bgo_corr_4   = new TH1D("bgo_corr_4",   "", 4097, 0., 4097.);
  TH1D* h1_bgo_corr_5   = new TH1D("bgo_corr_5",   "", 4097, 0., 4097.);
  TH1D* h1_bgo_corr_6   = new TH1D("bgo_corr_6",   "", 4097, 0., 4097.);
  TH1D* h1_bgo_corr_7   = new TH1D("bgo_corr_7",   "", 4097, 0., 4097.);
  TH1D* h1_bgo_corr_tot   = new TH1D("bgo_corr_tot",   "", 4097, 0., 8.*4097.);

  TH1D* h1_xPos_bgo = new TH1D("xPos_bgo", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_bgo = new TH1D("yPos_bgo", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_bgo = new TH2D("xyPos_bgo", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_xPos_hodo = new TH1D("xPos_hodo", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_hodo = new TH1D("yPos_hodo", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_hodo = new TH2D("xyPos_hodo", "", nBins, -xMax, xMax, nBins, -xMax, xMax);



  TH2D* h2_correlation_cef3_hodo_xPos = new TH2D("correlation_cef3_hodo_xPos", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);
  TH2D* h2_correlation_cef3_hodo_yPos = new TH2D("correlation_cef3_hodo_yPos", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);

  TH2D* h2_correlation_cef3_bgo_xPos = new TH2D("correlation_cef3_bgo_xPos", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);
  TH2D* h2_correlation_cef3_bgo_yPos = new TH2D("correlation_cef3_bgo_yPos", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);

  TH2D* h2_correlation_hodo_bgo_xPos = new TH2D("correlation_hodo_bgo_xPos", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);
  TH2D* h2_correlation_hodo_bgo_yPos = new TH2D("correlation_hodo_bgo_yPos", "", 100, -xySize/2., xySize/2.,  100, -xySize/2., xySize/2.);


  int nentries = tree->GetEntries();

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    if( iEntry % 1000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;


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
    }


    std::vector<float>  bgo_corr = subtractPedestals( bgo , pedestals_bgo, 4. );
    std::vector<float> cef3_corr = subtractPedestals( cef3, pedestals,     4. );


    bool cef3_ok = checkVector(cef3);
    bool cef3_corr_ok = checkVector(cef3_corr);



    // FIRST GET POSITION FROM HODOSCOPE:

    bool hodox_ok = checkVector(hodox, 99999.);
    bool hodoy_ok = checkVector(hodoy, 99999.);

    float nSigma_hodo = 4.;
    float xPos_hodo = getMeanposHodo(hodox, pedestals_hodox, nSigma_hodo);
    float yPos_hodo = getMeanposHodo(hodoy, pedestals_hodoy, nSigma_hodo);

    if( hodox_ok )
      h1_xPos_hodo->Fill(xPos_hodo);
    if( hodoy_ok )
      h1_yPos_hodo->Fill(yPos_hodo);

    if( hodox_ok && hodoy_ok ) 
      h2_xyPos_hodo->Fill(xPos_hodo, yPos_hodo);



    bool bgo_ok = checkVector(bgo, 4095.);
    bool bgo_corr_ok = checkVector(bgo_corr, 4095.);

    float xPos_bgo;
    float yPos_bgo;

    if( bgo_ok && bgo_corr_ok ) {


      // first: BGO precalibration:

      std::vector<float> bgo_precalibration;
      bgo_precalibration.push_back(193.7);
      bgo_precalibration.push_back(305.0);
      bgo_precalibration.push_back(263.8);
      bgo_precalibration.push_back(296.2);
      bgo_precalibration.push_back(201.6);
      bgo_precalibration.push_back(194.9);
      bgo_precalibration.push_back(229.1);
      bgo_precalibration.push_back(241.1);

      float bgoCalibrationAverage = sumVector(bgo_precalibration)/bgo_precalibration.size();

      for(unsigned i=0; i<bgo_precalibration.size(); ++i ) {
        bgo_precalibration[i] /= bgoCalibrationAverage;
        bgo_corr[i] *= bgo_precalibration[i]; //correct
      }

      float eTot_bgo_corr  = sumVector(bgo_corr);

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


      float position_bgo = xySize; // in mm
      //float position_bgo = 22.; // in mm
      
      std::vector<float> xPosW_bgo;
      xPosW_bgo.push_back(bgo_corr[0]*(-position_bgo));
      xPosW_bgo.push_back(0.);
      xPosW_bgo.push_back(bgo_corr[2]*(+position_bgo));
      xPosW_bgo.push_back(bgo_corr[3]*(-position_bgo));
      xPosW_bgo.push_back(bgo_corr[4]*(+position_bgo));
      xPosW_bgo.push_back(bgo_corr[5]*(-position_bgo));
      xPosW_bgo.push_back(0.);
      xPosW_bgo.push_back(bgo_corr[7]*(+position_bgo));
      
      std::vector<float> yPosW_bgo;
      yPosW_bgo.push_back(bgo_corr[0]*(+position_bgo));
      yPosW_bgo.push_back(bgo_corr[1]*(+position_bgo));
      yPosW_bgo.push_back(bgo_corr[2]*(+position_bgo));
      yPosW_bgo.push_back(0.);
      yPosW_bgo.push_back(0.);
      yPosW_bgo.push_back(bgo_corr[5]*(-position_bgo));
      yPosW_bgo.push_back(bgo_corr[6]*(-position_bgo));
      yPosW_bgo.push_back(bgo_corr[7]*(-position_bgo));

      xPos_bgo = sumVector( xPosW_bgo )/eTot_bgo_corr;
      yPos_bgo = sumVector( yPosW_bgo )/eTot_bgo_corr;
      
      h1_xPos_bgo->Fill( xPos_bgo );
      h1_yPos_bgo->Fill( yPos_bgo );
      h2_xyPos_bgo->Fill( xPos_bgo, yPos_bgo );

      h2_correlation_hodo_bgo_xPos->Fill( xPos_hodo, xPos_bgo );
      h2_correlation_hodo_bgo_yPos->Fill( yPos_hodo, yPos_bgo );
      
    }  // if bgo ok


    if( cef3_ok ) {

      // THEN USE CeF3 DATA:

      float eTot      = sumVector(cef3);
      float eTot_corr = sumVector(cef3_corr);

      h1_cef3_0->Fill( cef3[0] );
      h1_cef3_1->Fill( cef3[1] );
      h1_cef3_2->Fill( cef3[2] );
      h1_cef3_3->Fill( cef3[3] );
      h1_cef3_tot->Fill( eTot );

      h1_cef3_corr_0->Fill( cef3_corr[0] );
      h1_cef3_corr_1->Fill( cef3_corr[1] );
      h1_cef3_corr_2->Fill( cef3_corr[2] );
      h1_cef3_corr_3->Fill( cef3_corr[3] );
      h1_cef3_corr_tot->Fill( eTot_corr );
      

      if( cef3_corr_ok ) {

        //   0      1
        //          
        //          
        //   3      2


        float chamfer = 2.1; // in mm

        float position = xySize/2. - chamfer/4.;

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


        float xPos = sumVector(xPosW)/eTot_corr;
        float yPos = sumVector(yPosW)/eTot_corr;

        h1_xPos->Fill( xPos );
        h1_yPos->Fill( yPos );

        h2_xyPos->Fill( xPos, yPos );


        // CORRELATIONS BETWEEN CALO AND HODO:
  
        h2_correlation_cef3_hodo_xPos->Fill( xPos, xPos_hodo );
        h2_correlation_cef3_hodo_yPos->Fill( yPos, yPos_hodo );

        if( bgo_ok && bgo_corr_ok ) {
          h2_correlation_cef3_bgo_xPos->Fill( xPos, xPos_bgo );
          h2_correlation_cef3_bgo_yPos->Fill( yPos, yPos_bgo );
        }

      } // if cef3_ok

    }

  }


  std::cout << "-> Events passing overflow cut: " << h1_xPos->GetEntries() << "/" << nentries << " (" << 100.* h1_xPos->GetEntries()/nentries << "%)" << std::endl;

  std::string outfileName = "PosAn_" + runName + ".root";
  TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );
  outfile->cd();

  h1_xPos->Write();
  h1_yPos->Write();
  
  h2_xyPos->Write();

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

  h2_correlation_cef3_hodo_xPos->Write();
  h2_correlation_cef3_hodo_yPos->Write();

  h2_correlation_cef3_bgo_xPos->Write();
  h2_correlation_cef3_bgo_yPos->Write();

  h2_correlation_hodo_bgo_xPos->Write();
  h2_correlation_hodo_bgo_yPos->Write();

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

float getMeanposHodo( std::vector<float> hodo, std::vector<std::pair<float,float> > pedestals, float nSigma ) {

  float mean = 0.;
  float eTot = 0.;
  for( unsigned i=0; i<hodo.size(); ++i ) {
    if( hodo[i] > (pedestals[i].first + nSigma* pedestals[i].second) ) {
      mean += -(i-3.5);
      eTot += 1.;
    }
  }

  return mean/eTot;

}

