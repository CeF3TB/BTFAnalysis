#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector2.h"
#include "TString.h"
#include "TF1.h"


#include "fastDQM_CeF3_BTF.h"
#include "interface/HodoCluster.h"
#include "interface/RunHelper.h"
#include "interface/CalibrationUtility.h"
#include "interface/EnergyCalibration.h"

#include "TApplication.h"



std::vector< std::pair<float, float> > getPedestals( const std::string& type, const std::string& fileName, int runNumber );
std::vector<float> subtractPedestals( std::vector<float> raw, std::vector< std::pair<float, float> > pedestals, float nSigma );
float sumVector( std::vector<float> v );
bool checkVector( std::vector<float> v, float theMax=4095. );

std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo_corr, int nClusterMax );


int main( int argc, char* argv[] ) {

  TApplication* a = new TApplication("a", 0, 0);


  std::string runName = "precalib_BGO_pedestal_noSource";
  std::string tag = "default";
  if( argc>2 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
    std::string tag_str(argv[2]);
    tag = tag_str;
  }else{
    std::cout<<"Usage:"<<std::endl;
    std::cout<<"./makeAnalysisTree BTF_XXX tag"<<std::endl;
    exit(12345);
  }

  TString runName_tstr(runName);
  bool isOnlyRunNumber = !(runName_tstr.BeginsWith("BTF_"));


  TChain* tree = new TChain("eventRawData");
  if( isOnlyRunNumber ) {
    std::cout << "-> We believe you are passing the program only the run number!" << std::endl;
    std::cout << "-> So for instance you are passing '246' for run 'BTF_246_20140501-212512_beam'" << std::endl;
    std::cout << "(if this is not the case this means TROUBLE)" << std::endl;
    std::cout << "-> Will look for runs matching run number: " << runName << std::endl;
    tree->Add(Form("data/run_BTF_%s_2014*_beam.root/eventRawData", runName.c_str()) );
    if( tree->GetEntries()==0 ) {
      std::cout << "WARNING! Didn't find any events matching run: " << runName << std::endl;
      std::cout << "Exiting" << std::endl;
      exit(1913);
    }
  } else {
    std::string fileName = "data/run_" + runName + ".root";
    TFile* file = TFile::Open(fileName.c_str());
    if( file==0 ) {
      std::cout << "ERROR! Din't find file " << fileName << std::endl;
      std::cout << "Exiting." << std::endl;
      exit(11);
    }
    tree = (TChain*)file->Get("eventRawData");
  }



  // get run number with a trick:
  unsigned int runNumber_;
  if( isOnlyRunNumber ) {
    runNumber_ = atoi(runName.c_str());
  } else {
    char runNumber_cstr[3];
    runNumber_cstr[0] = runName.at(4);
    runNumber_cstr[1] = runName.at(5);
    runNumber_cstr[2] = runName.at(6);
    TString runNumber_tstr(runNumber_cstr);
    if(runNumber_tstr.EndsWith("_")) runNumber_tstr.Chop();
    std::string runNumber_str(runNumber_tstr.Data());
    runNumber_ = atoi(runNumber_str.c_str());
  }

  if( isOnlyRunNumber ) {
    // then modify runname in such a way that it's useful for getBeamPosition and outfile:
    runName = "BTF_" + runName + "_beam";
  }




  std::string pedestalFileName = "pedestalFile.root";
  std::cout << "-> Using pedestal file: " << pedestalFileName << std::endl;

  std::vector<std::pair<float, float> > pedestals = getPedestals( "cef3", pedestalFileName, runNumber_ );
  std::cout << std::endl;
  std::cout << "-> Got pedestals of CeF3: " << std::endl;
  for( unsigned i=0; i<CEF3_CHANNELS; ++i )
    std::cout << " CeF3 Channel " << i << ": " << pedestals[i].first << " (+- " << pedestals[i].second << ")" << std::endl;
  std::cout << std::endl;

  std::vector<std::pair<float, float> > pedestals_bgo = getPedestals( "bgo", pedestalFileName, runNumber_ );
  std::cout << std::endl;
  std::cout << "-> Got pedestals of BGO: " << std::endl;
  for( unsigned i=0; i<BGO_CHANNELS; ++i )
    std::cout << " BGO Channel " << i << ": " << pedestals_bgo[i].first << " (+- " << pedestals_bgo[i].second << ")" << std::endl;
  std::cout << std::endl;

  std::vector<std::pair<float, float> > pedestals_hodox = getPedestals("hodox", pedestalFileName, runNumber_);
  std::vector<std::pair<float, float> > pedestals_hodoy = getPedestals("hodoy", pedestalFileName, runNumber_);
  std::cout << "-> Got Hodoscope pedestals: " << std::endl;
  std::cout << std::endl;
  for( unsigned i=0; i<HODOX_CHANNELS; ++i )
    std::cout << "Channel " << i << ":  X: " << pedestals_hodox[i].first << " (+- " << pedestals_hodox[i].second << ") Y: " << pedestals_hodoy[i].first << " (+- " << pedestals_hodoy[i].second << ")" << std::endl;

  std::cout << std::endl << std::endl;



  //set the tag for calibration
  CalibrationUtility calibUtil(tag);
  EnergyCalibration cef3Calib(calibUtil.getCeF3FileName());
  EnergyCalibration bgoCalib(calibUtil.getBGOFileName());



  UInt_t evtNumber;
  tree->SetBranchAddress( "evtNumber", &evtNumber );
  UInt_t adcData[40];
  tree->SetBranchAddress( "adcData", adcData );
  UInt_t adcBoard[40];
  tree->SetBranchAddress( "adcBoard", adcBoard );
  UInt_t adcChannel[40];
  tree->SetBranchAddress( "adcChannel", adcChannel );



  std::string outdir = "analysisTrees_" + tag;
  system( Form("mkdir -p %s", outdir.c_str()) );
  std::string outfileName = outdir + "/Reco_" + runName + ".root";
  TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );


  TTree* outTree = new TTree("recoTree","recoTree");
  float cef3_[CEF3_CHANNELS],cef3_pedSubtracted_[CEF3_CHANNELS],bgo_[BGO_CHANNELS],hodox_[HODOX_CHANNELS],hodoy_[HODOY_CHANNELS];
  float cef3_corr_[CEF3_CHANNELS],bgo_corr_[BGO_CHANNELS],bgo_pedSubtracted_[BGO_CHANNELS],hodox_corr_[HODOX_CHANNELS],hodoy_corr_[HODOY_CHANNELS];
  float scintFront_;
  bool isSingleEle_scintFront_;
  int cef3_chan=CEF3_CHANNELS;
  int bgo_chan=BGO_CHANNELS;
  int hodox_chan=HODOX_CHANNELS;
  int hodoy_chan=HODOY_CHANNELS; 
  float xBeam_, yBeam_;

  int nHodoFibersX;
  int nHodoFibersY;
  int nHodoClustersX;
  int nHodoClustersY;
  float pos_hodoClustX_[HODOX_CHANNELS];
  float pos_hodoClustY_[HODOY_CHANNELS];
  int nFibres_hodoClustX_[HODOX_CHANNELS];
  int nFibres_hodoClustY_[HODOY_CHANNELS];
  bool cef3_ok_;
  bool cef3_corr_ok_;
  bool bgo_ok_;
  bool bgo_corr_ok_;

  outTree->Branch( "runNumber", &runNumber_,"runNumber/i" );
  outTree->Branch( "evtNumber", &evtNumber,"evtNumber/i" );
  outTree->Branch( "adcData", adcData, "adcData/i" );
  outTree->Branch( "adcBoard", adcBoard, "adcBoard/i" );
  outTree->Branch( "adcChannel", adcChannel,"adcChannel/i" );
  outTree->Branch( "nHodoFibersX", &nHodoFibersX, "nHodoFibersX/I" );
  outTree->Branch( "nHodoFibersY", &nHodoFibersY, "nHodoFibersY/I" );
  outTree->Branch( "nHodoClustersX", &nHodoClustersX, "nHodoClustersX/I" );
  outTree->Branch( "pos_hodoClustX", pos_hodoClustX_, "pos_hodoClustX_[nHodoClustersX]/F" );
  outTree->Branch( "nFibres_hodoClustX", nFibres_hodoClustX_, "nFibres_hodoClustX_[nHodoClustersX]/I" );
  outTree->Branch( "nHodoClustersY", &nHodoClustersY, "nHodoClustersY/I" );
  outTree->Branch( "pos_hodoClustY", pos_hodoClustY_, "pos_hodoClustY_[nHodoClustersY]/F" );
  outTree->Branch( "nFibres_hodoClustY", nFibres_hodoClustY_, "nFibres_hodoClustY_[nHodoClustersY]/I" );
  outTree->Branch( "hodox_chan", &hodox_chan, "hodox_chan/I" );
  outTree->Branch( "hodoy_chan", &hodoy_chan, "hodoy_chan/I" );
  outTree->Branch( "cef3_chan", &cef3_chan, "cef3_chan/I" );
  outTree->Branch( "bgo_chan", &bgo_chan, "bgo_chan/I" );
  outTree->Branch( "hodox", hodox_, "hodox_[hodox_chan]/F" );
  outTree->Branch( "hodoy", hodoy_, "hodoy_[hodoy_chan]/F" );
  outTree->Branch( "bgo", bgo_, "bgo_[bgo_chan]/F" );
  outTree->Branch( "cef3", cef3_, "cef3_[cef3_chan]/F" );
  outTree->Branch( "cef3_pedSubtracted", cef3_pedSubtracted_, "cef3_pedSubtracted_[cef3_chan]/F" );
  outTree->Branch( "bgo_corr", bgo_corr_, "bgo_corr_[bgo_chan]/F" );
  outTree->Branch( "bgo_pedSubtracted", bgo_pedSubtracted_, "bgo_pedSubtracted_[bgo_chan]/F" );
  outTree->Branch( "cef3_corr", cef3_corr_, "cef3_corr_[cef3_chan]/F" );
  outTree->Branch( "hodox_corr", hodox_corr_, "hodox_corr_[hodox_chan]/F" );
  outTree->Branch( "hodoy_corr", hodoy_corr_, "hodoy_corr_[hodoy_chan]/F" );
  outTree->Branch( "scintFront", &scintFront_, "scintFront_/F" );
  outTree->Branch( "isSingleEle_scintFront", &isSingleEle_scintFront_, "isSingleEle_scintFront_/O" );
  outTree->Branch( "xBeam", &xBeam_, "xBeam_/F" );
  outTree->Branch( "yBeam", &yBeam_, "yBeam_/F" );
  outTree->Branch( "cef3_ok", &cef3_ok_, "cef3_ok_/O" );
  outTree->Branch( "cef3_corr_ok", &cef3_corr_ok_, "cef3_corr_ok_/O" );
  outTree->Branch( "bgo_ok", &bgo_ok_, "bgo_ok_/O" );
  outTree->Branch( "bgo_corr_ok", &bgo_corr_ok_, "bgo_corr_ok_/O" );


  
  
  
  int nentries = tree->GetEntries();



  RunHelper::getBeamPosition( runName, xBeam_, yBeam_ );


  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    nHodoFibersX=0;
    nHodoFibersY=0;
    nHodoClustersX=0;
    nHodoClustersY=0;
    

    for( unsigned i=0; i<HODOX_CHANNELS; ++i ) {
      nFibres_hodoClustX_[i] = -1.;
      pos_hodoClustX_[i] = -9999.;
    }

    for( unsigned i=0; i<HODOY_CHANNELS; ++i ) {
      nFibres_hodoClustY_[i] = -1.;
      pos_hodoClustY_[i] = -9999.;
    }

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


    std::vector<float>  bgo_corr = subtractPedestals( bgo , pedestals_bgo, 4. );
    std::vector<float> cef3_corr = subtractPedestals( cef3, pedestals,     4. );
    std::vector<float> hodox_corr = subtractPedestals( hodox, pedestals_hodox, 4. );
    std::vector<float> hodoy_corr = subtractPedestals( hodoy, pedestals_hodoy, 4. );

    bgoCalib.applyCalibration(bgo_corr);
    cef3Calib.applyCalibration(cef3_corr);


    std::vector<float> cef3_pedSubtracted = subtractPedestals( cef3, pedestals,    4. );
    std::vector<float>  bgo_pedSubtracted= subtractPedestals( bgo , pedestals_bgo, 4. );

    cef3_ok_ = checkVector(cef3);
    cef3_corr_ok_ = checkVector(cef3_corr);



    // FIRST GET POSITION FROM HODOSCOPE:

    int clusterSize=4;
    std::vector<HodoCluster*> hodoxFibres   = getHodoClusters( hodox_corr, 1 ); // fibres are just clusters with size = 1
    std::vector<HodoCluster*> hodoxClusters = getHodoClusters( hodox_corr, clusterSize );

    std::vector<HodoCluster*> hodoyFibres   = getHodoClusters( hodoy_corr, 1 );
    std::vector<HodoCluster*> hodoyClusters = getHodoClusters( hodoy_corr, clusterSize );

    nHodoFibersX = hodoxFibres.size();
    nHodoFibersY = hodoyFibres.size();

    nHodoClustersX = hodoxClusters.size();
    nHodoClustersY = hodoyClusters.size();

    for( unsigned i=0; i<hodoxClusters.size(); ++i ) {
      nFibres_hodoClustX_[i] = hodoxClusters[i]->getSize();
      pos_hodoClustX_[i] = hodoxClusters[i]->getPosition();
    }

    for( unsigned i=0; i<hodoyClusters.size(); ++i ) {
      nFibres_hodoClustY_[i] = hodoyClusters[i]->getSize();
      pos_hodoClustY_[i] = hodoyClusters[i]->getPosition();
    }


    float scintFrontMin = (runNumber_<=100) ? 110. : 500.;
    float scintFrontMax = (runNumber_<=100) ? 700. : 2000.;
    isSingleEle_scintFront_ = (scintFront_> scintFrontMin && scintFront_< scintFrontMax);




    bgo_ok_ = checkVector(bgo, 4095.);
    bgo_corr_ok_ = checkVector(bgo_corr, 4095.);



    for(int i=0;i<CEF3_CHANNELS;i++){
      cef3_[i]=cef3[i]; 
      cef3_corr_[i]=cef3_corr[i];
      cef3_pedSubtracted_[i]=cef3_pedSubtracted[i];
    }
    for(int i=0;i<BGO_CHANNELS;i++){
      bgo_[i]=bgo[i];
      bgo_corr_[i]=bgo_corr[i];
      bgo_pedSubtracted_[i]=bgo_pedSubtracted[i];
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

  }



  outfile->cd();

  outTree->Write();


  
  outfile->Close();
  std::cout << "-> Analysis Tree saved in: " << outfile->GetName() << std::endl;


  return 0;

}






std::vector< std::pair<float, float> > getPedestals( const std::string& type, const std::string& fileName, int runNumber ) {

  int nChannels=-1;

  if( type=="cef3" ) {
    nChannels    = CEF3_CHANNELS;
  } else if( type=="bgo" ) {
    nChannels    = BGO_CHANNELS;
  } else if( type=="hodox" ) {
    nChannels    = HODOX_CHANNELS;
  } else if( type=="hodoy" ) {
    nChannels    = HODOY_CHANNELS;
  } else {
    std::cout << "ERROR! Unkown type '" << type << "'!" << std::endl;
    std::cout << "Don't know what pedestals you're looking for." << std::endl;
    std::cout << "Exiting." << std::endl;

    exit(77);
  }

    

  TFile* file = TFile::Open(fileName.c_str());

  std::vector< std::pair<float, float> > peds;

  for( unsigned i=0; i<nChannels; ++i ) {

    TH1D* h1_ped = (TH1D*)file->Get(Form("%s_%d", type.c_str(), i));
    
    int iped = runNumber;

    float ped=-1.;
    float pedrms=0.;
    while( ped<0. ) { // get closest run before current one
      ped = h1_ped->GetBinContent(iped);
      
      pedrms = h1_ped->GetBinError(iped);

      iped--;
    }
    std::pair<float, float>  thispair;
    thispair.first  = ped;
    thispair.second = pedrms;
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



