#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector2.h"
#include "TF1.h"

#include "TMVA/Reader.h"

#include "fastDQM_CeF3_BTF.h"
#include "interface/RunHelper.h"

//#include "hodo_efficiency.dat"





float sumVector( std::vector<float> v );
bool checkVector( std::vector<float> v, float theMax=4095. );

float getMeanposHodo( int n, float* pos );
void getCeF3Position( std::vector<float> cef3, float& xPos, float& yPos );
TF1* getCef3Function( const std::string& diagName );
//float gethodointercalib(TString axis, int n);







int main( int argc, char* argv[] ) {


  std::string runName = "precalib_BGO_pedestal_noSource";
  if( argc>1 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
  }

  std::string tag = "V00";
  if( argc>2 ) {
    std::string tag_str(argv[2]);
    tag = tag_str;
  }

  TString runName_tstr(runName);
  bool isOnlyRunNumber = !(runName_tstr.BeginsWith("BTF_"));


  TChain* tree = new TChain("recoTree");
  if( isOnlyRunNumber ) {
    std::cout << "-> We believe you are passing the program only the run number!" << std::endl;
    std::cout << "-> So for instance you are passing '246' for run 'BTF_246_20140501-212512_beam'" << std::endl;
    std::cout << "(if this is not the case this means TROUBLE)" << std::endl;
    std::cout << "-> Will look for runs matching run number: " << runName << std::endl;
    tree->Add(Form("analysisTrees_%s/Reco_BTF_%s_*beam.root/recoTree", tag.c_str(), runName.c_str()) );
    if( tree->GetEntries()==0 ) {
      std::cout << "WARNING! Didn't find any events matching run: " << runName << std::endl;
      std::cout << "Exiting" << std::endl;
      exit(1913);
    }
  } else {
    std::string fileName = "analysisTrees_"+tag+"/Reco_" + runName + ".root";
    TFile* file = TFile::Open(fileName.c_str());
    if( file==0 ) {
      std::cout << "ERROR! Din't find file " << fileName << std::endl;
      std::cout << "Exiting." << std::endl;
      exit(11);
    }
    tree = (TChain*)file->Get("recoTree");
  }





  UInt_t evtNumber;
  tree->SetBranchAddress( "evtNumber", &evtNumber );
  UInt_t adcData[40];
  tree->SetBranchAddress( "adcData", adcData );
  UInt_t adcBoard[40];
  tree->SetBranchAddress( "adcBoard", adcBoard );
  UInt_t adcChannel[40];
  tree->SetBranchAddress( "adcChannel", adcChannel );


  unsigned int runNumber;
  int nHodoFibersX;
  int nHodoFibersY;
  int nHodoClustersX;
  int nHodoClustersY;
  float cef3_corr[CEF3_CHANNELS];
  float bgo_corr[BGO_CHANNELS];
  float scintFront;
  float pos_hodoClustX[HODOX_CHANNELS];
  float pos_hodoClustY[HODOY_CHANNELS];
  int nFibres_hodoClustX[HODOX_CHANNELS];
  int nFibres_hodoClustY[HODOY_CHANNELS];
  float xBeam, yBeam;
  bool isSingleEle_scintFront;
  bool cef3_ok;
  bool cef3_corr_ok;
  bool bgo_ok;
  bool bgo_corr_ok;

  tree->SetBranchAddress( "runNumber", &runNumber );

  tree->SetBranchAddress( "scintFront", &scintFront );
  tree->SetBranchAddress( "cef3_corr", cef3_corr );
  tree->SetBranchAddress( "bgo_corr", bgo_corr );

  tree->SetBranchAddress( "nHodoFibersX", &nHodoFibersX );
  tree->SetBranchAddress( "nHodoFibersY", &nHodoFibersY );
  tree->SetBranchAddress( "nHodoClustersX", &nHodoClustersX );
  tree->SetBranchAddress( "pos_hodoClustX", pos_hodoClustX );
  tree->SetBranchAddress( "nFibres_hodoClustX", nFibres_hodoClustX );
  tree->SetBranchAddress( "nHodoClustersY", &nHodoClustersY );
  tree->SetBranchAddress( "pos_hodoClustY", pos_hodoClustY );
  tree->SetBranchAddress( "nFibres_hodoClustY", nFibres_hodoClustY );
  tree->SetBranchAddress( "scintFront", &scintFront );
  tree->SetBranchAddress( "isSingleEle_scintFront", &isSingleEle_scintFront );
  tree->SetBranchAddress( "xBeam", &xBeam );
  tree->SetBranchAddress( "yBeam", &yBeam );

  tree->SetBranchAddress( "cef3_ok", &cef3_ok );
  tree->SetBranchAddress( "cef3_corr_ok", &cef3_corr_ok );
  tree->SetBranchAddress( "bgo_ok", &bgo_ok );
  tree->SetBranchAddress( "bgo_corr_ok", &bgo_corr_ok );




  int nBins = 500;
  float xMax = 25.*3./2.;


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


  TH1D* h1_xPos_bgo = new TH1D("xPos_bgo", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_bgo = new TH1D("yPos_bgo", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_bgo = new TH2D("xyPos_bgo", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_xPos_singleEle_bgo = new TH1D("xPos_singleEle_bgo", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_singleEle_bgo = new TH1D("yPos_singleEle_bgo", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_singleEle_bgo = new TH2D("xyPos_singleEle_bgo", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_xPos_hodo = new TH1D("xPos_hodo", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_hodo = new TH1D("yPos_hodo", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_hodo = new TH2D("xyPos_hodo", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

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

  TH1D* h1_xPos_calo_vs_hodo = new TH1D("xPos_calo_vs_hodo", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_calo_vs_hodo = new TH1D("yPos_calo_vs_hodo", "", nBins, -xMax, xMax);

  TH1D* h1_xPos_calo_vs_beam = new TH1D("xPos_calo_vs_beam", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_calo_vs_beam = new TH1D("yPos_calo_vs_beam", "", nBins, -xMax, xMax);

  TH1D* h1_xPos_calo_vs_hodo_singleElectron = new TH1D("xPos_calo_vs_hodo_singleElectron", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_calo_vs_hodo_singleElectron = new TH1D("yPos_calo_vs_hodo_singleElectron", "", nBins, -xMax, xMax);

  TH1D* h1_xPos_calo_vs_beam_singleElectron = new TH1D("xPos_calo_vs_beam_singleElectron", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_calo_vs_beam_singleElectron = new TH1D("yPos_calo_vs_beam_singleElectron", "", nBins, -xMax, xMax);

  TH1D* h1_xPos_calo_vs_hodo_singleElectron_HR = new TH1D("xPos_calo_vs_hodo_singleElectron_HR", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_calo_vs_hodo_singleElectron_HR = new TH1D("yPos_calo_vs_hodo_singleElectron_HR", "", nBins, -xMax, xMax);

  TH1D* h1_xPos_calo_vs_beam_singleElectron_HR = new TH1D("xPos_calo_vs_beam_singleElectron_HR", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_calo_vs_beam_singleElectron_HR = new TH1D("yPos_calo_vs_beam_singleElectron_HR", "", nBins, -xMax, xMax);


  TH2D* h2_correlation_cef3_hodo_xPos = new TH2D("correlation_cef3_hodo_xPos", "", 100, -12.5, 12.5,  100, -12.5, 12.5);
  TH2D* h2_correlation_cef3_hodo_yPos = new TH2D("correlation_cef3_hodo_yPos", "", 100, -12.5, 12.5,  100, -12.5, 12.5);

  TH2D* h2_correlation_cef3_bgo_xPos = new TH2D("correlation_cef3_bgo_xPos", "", 100, -12.5, 12.5,  100, -12.5, 12.5);
  TH2D* h2_correlation_cef3_bgo_yPos = new TH2D("correlation_cef3_bgo_yPos", "", 100, -12.5, 12.5,  100, -12.5, 12.5);

  TH2D* h2_correlation_hodo_bgo_xPos = new TH2D("correlation_hodo_bgo_xPos", "", 100, -12.5, 12.5,  100, -12.5, 12.5);
  TH2D* h2_correlation_hodo_bgo_yPos = new TH2D("correlation_hodo_bgo_yPos", "", 100, -12.5, 12.5,  100, -12.5, 12.5);


  TH2D* h2_correlation_cef3_hodo_xPos_singleEle = new TH2D("correlation_cef3_hodo_xPos_singleEle", "", 100, -12.5, 12.5,  100, -12.5, 12.5);
  TH2D* h2_correlation_cef3_hodo_yPos_singleEle = new TH2D("correlation_cef3_hodo_yPos_singleEle", "", 100, -12.5, 12.5,  100, -12.5, 12.5);

  TH2D* h2_correlation_cef3_bgo_xPos_singleEle = new TH2D("correlation_cef3_bgo_xPos_singleEle", "", 100, -12.5, 12.5,  100, -12.5, 12.5);
  TH2D* h2_correlation_cef3_bgo_yPos_singleEle = new TH2D("correlation_cef3_bgo_yPos_singleEle", "", 100, -12.5, 12.5,  100, -12.5, 12.5);

  TH2D* h2_correlation_hodo_bgo_xPos_singleEle = new TH2D("correlation_hodo_bgo_xPos_singleEle", "", 100, -12.5, 12.5,  100, -12.5, 12.5);
  TH2D* h2_correlation_hodo_bgo_yPos_singleEle = new TH2D("correlation_hodo_bgo_yPos_singleEle", "", 100, -12.5, 12.5,  100, -12.5, 12.5);



  int nentries = tree->GetEntries();


  if( isOnlyRunNumber ) {
    // modify runname in such a way that it's useful for getBeamPosition and outfile:
    runName = "BTF_" + runName + "_beam";
  }



  std::string outputdir = "PosAnTrees_"+tag;
  system( Form("mkdir -p %s", outputdir.c_str()) );
  std::string outfileName = outputdir + "/PosAn_" + runName + ".root";
  TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );

  TTree* outTree = new TTree("posTree","posTree");
  float xPos_calo_, yPos_calo_;
  float xPos_bgo_, yPos_bgo_;
  float xPos_new_, yPos_new_;
  float xPos_regr2D_, yPos_regr2D_;
  float r02_, r13_;

  outTree->Branch( "isSingleEle_scintFront", &isSingleEle_scintFront, "isSingleEle_scintFront/O" );
  outTree->Branch( "nHodoClustersX", &nHodoClustersX, "nHodoClustersX/I" );
  outTree->Branch( "nHodoClustersY", &nHodoClustersY, "nHodoClustersY/I" );
  outTree->Branch( "cef3_corr", cef3_corr, "cef3_corr[4]/F" );
  outTree->Branch( "r02", &r02_, "r02_/F" );
  outTree->Branch( "r13", &r13_, "r13_/F" );
  outTree->Branch( "xBeam", &xBeam, "xBeam/F" );
  outTree->Branch( "yBeam", &yBeam, "yBeam/F" );
  outTree->Branch( "xPos_bgo", &xPos_bgo_, "xPos_bgo_/F" );
  outTree->Branch( "yPos_bgo", &yPos_bgo_, "yPos_bgo_/F" );
  outTree->Branch( "xPos_calo", &xPos_calo_, "xPos_calo_/F" );
  outTree->Branch( "yPos_calo", &yPos_calo_, "yPos_calo_/F" );
  outTree->Branch( "xPos_new", &xPos_new_, "xPos_new_/F" );
  outTree->Branch( "yPos_new", &yPos_new_, "yPos_new_/F" );
  outTree->Branch( "xPos_regr2D", &xPos_regr2D_, "xPos_regr2D_/F" );
  outTree->Branch( "yPos_regr2D", &yPos_regr2D_, "yPos_regr2D_/F" );

  float diag02_calo_;
  float diag13_calo_;
  outTree->Branch( "diag02_calo", &diag02_calo_, "diag02_calo_/F" );
  outTree->Branch( "diag13_calo", &diag13_calo_, "diag13_calo_/F" );
  float diag02_new_;
  float diag13_new_;
  outTree->Branch( "diag02_new", &diag02_new_, "diag02_new_/F" );
  outTree->Branch( "diag13_new", &diag13_new_, "diag13_new_/F" );
  float diag02_beam_;
  float diag13_beam_;
  outTree->Branch( "diag02_beam", &diag02_beam_, "diag02_beam_/F" );
  outTree->Branch( "diag13_beam", &diag13_beam_, "diag13_beam_/F" );

  
  std::vector<float> xbgo, ybgo;
  for( unsigned i=0; i<BGO_CHANNELS; ++i ) {
    float x,y;
    RunHelper::getBGOCoordinates( i, x, y );
    xbgo.push_back( x );
    ybgo.push_back( y );
  }


  
  



  float cef3_regr[CEF3_CHANNELS];
  //TMVA::Reader* readerRegrX = new TMVA::Reader( "!Color:!Silent" );
  //readerRegrX->AddVariable("cef3_corr[0]/(cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3])", &cef3_regr[0] );
  //readerRegrX->AddVariable("cef3_corr[1]/(cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3])", &cef3_regr[1] );
  //readerRegrX->AddVariable("cef3_corr[2]/(cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3])", &cef3_regr[2] );
  //readerRegrX->AddVariable("cef3_corr[3]/(cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3])", &cef3_regr[3] );

  //TMVA::Reader* readerRegrY = new TMVA::Reader( "!Color:!Silent" );
  //readerRegrY->AddVariable("cef3_corr[0]/(cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3])", &cef3_regr[0] );
  //readerRegrY->AddVariable("cef3_corr[1]/(cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3])", &cef3_regr[1] );
  //readerRegrY->AddVariable("cef3_corr[2]/(cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3])", &cef3_regr[2] );
  //readerRegrY->AddVariable("cef3_corr[3]/(cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3])", &cef3_regr[3] ); // let try this trick

  //TMVA::Reader* readerRegr2D = new TMVA::Reader( "!Color:!Silent" );
  //readerRegr2D->AddVariable("cef3_corr[0]", &cef3_corr_[0] );
  //readerRegr2D->AddVariable("cef3_corr[1]", &cef3_corr_[1] );
  //readerRegr2D->AddVariable("cef3_corr[2]", &cef3_corr_[2] );
  //readerRegr2D->AddVariable("cef3_corr[3]", &cef3_corr_[3] );
  //readerRegr2D->BookMVA( "MLP", "TMVA/weights/TMVARegression_MLP.weights.xml" );


  std::vector<std::string> methodNames;
  //methodNames.push_back("BDTG");
  ////methodNames.push_back("FDA_MT");
  //methodNames.push_back("LD");
  //methodNames.push_back("MLP");
  ////methodNames.push_back("PDERS");

 
  //
  //std::cout << "-> Booking TMVA Reader" << std::endl;
  //for( unsigned i=0; i<methodNames.size(); ++i ) {
  //  readerRegrX->BookMVA( methodNames[i], Form("TMVA/weights/TMVARegression_%s.weights.xml", methodNames[i].c_str()) ); 
  //  readerRegrY->BookMVA( methodNames[i], Form("TMVA/weights/TMVARegression_%s.weights.xml", methodNames[i].c_str()) ); 
  //}


  std::vector< TH1D* > h1_xPos_regr_vs_calo;
  std::vector< TH1D* > h1_yPos_regr_vs_calo;
  std::vector< TH2D* > h2_xyPos_regr;
  std::vector< TH1D* > h1_xPos_regr_vs_calo_singleEle;
  std::vector< TH1D* > h1_yPos_regr_vs_calo_singleEle;
  std::vector< TH2D* > h2_xyPos_singleEle_regr;
  for( unsigned i=0; i<methodNames.size(); ++i ) {
    TH1D* newHistx = new TH1D( Form("xPos_regr%s_vs_calo", methodNames[i].c_str()), "", nBins, -xMax, xMax);
    h1_xPos_regr_vs_calo.push_back( newHistx );
    TH1D* newHisty = new TH1D( Form("yPos_regr%s_vs_calo", methodNames[i].c_str()), "", nBins, -xMax, xMax);
    h1_yPos_regr_vs_calo.push_back( newHisty );
    TH2D* newHistxy = new TH2D( Form("xyPos_regr%s", methodNames[i].c_str()), "", nBins, -xMax, xMax, nBins, -xMax, xMax);
    h2_xyPos_regr.push_back( newHistxy );
    TH1D* newHistx_singleEle = new TH1D( Form("xPos_regr%s_vs_calo_singleEle", methodNames[i].c_str()), "", nBins, -xMax, xMax);
    h1_xPos_regr_vs_calo_singleEle.push_back( newHistx_singleEle );
    TH1D* newHisty_singleEle = new TH1D( Form("yPos_regr%s_vs_calo_singleEle", methodNames[i].c_str()), "", nBins, -xMax, xMax);
    h1_yPos_regr_vs_calo_singleEle.push_back( newHisty_singleEle );
    TH2D* newHistxy_singleEle = new TH2D( Form("xyPos_singleEle_regr%s", methodNames[i].c_str()), "", nBins, -xMax, xMax, nBins, -xMax, xMax);
    h2_xyPos_singleEle_regr.push_back( newHistxy_singleEle );
  }


  TH2D* h2_xyPos_regr2D = new TH2D("xyPos_regr2D", "", nBins, -xMax, xMax, nBins, -xMax, xMax);
  TH2D* h2_xyPos_singleEle_regr2D = new TH2D("xyPos_singleEle_regr2D", "", nBins, -xMax, xMax, nBins, -xMax, xMax);





  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    xPos_bgo_ = -999.;
    yPos_bgo_ = -999.;

    xPos_calo_ = -999.;
    yPos_calo_ = -999.;


    tree->GetEntry(iEntry);

    if( iEntry % 5000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;

    if( !bgo_corr_ok ) continue;

    r02_ = cef3_corr[0]/cef3_corr[2];
    r13_ = cef3_corr[1]/cef3_corr[3];

    // FIRST GET POSITION FROM HODOSCOPE:


    float xPos_hodo = getMeanposHodo(nHodoClustersX, pos_hodoClustX);
    float yPos_hodo = getMeanposHodo(nHodoClustersY, pos_hodoClustY);



    if( xPos_hodo>-100. )
      h1_xPos_hodo->Fill(xPos_hodo);
    if( yPos_hodo>-100. )
      h1_yPos_hodo->Fill(yPos_hodo);

    if( xPos_hodo>-100. && yPos_hodo>-100. ) 
      h2_xyPos_hodo->Fill(xPos_hodo, yPos_hodo);


    if( isSingleEle_scintFront ) {

      if( xPos_hodo>-100. )
        h1_xPos_singleEle_hodo->Fill(xPos_hodo);
      if( yPos_hodo>-100. )
        h1_yPos_singleEle_hodo->Fill(yPos_hodo);

      if( xPos_hodo>-100. && yPos_hodo>-100. ) 
        h2_xyPos_singleEle_hodo->Fill(xPos_hodo, yPos_hodo);

    }



    std::vector<float> xPosW_bgo;
    std::vector<float> yPosW_bgo;

    std::vector<float> v_bgo_corr;
    for( unsigned i=0; i<BGO_CHANNELS; ++i ) v_bgo_corr.push_back(bgo_corr[i]);

    float eTot_bgo_corr = sumVector(v_bgo_corr);

    if( bgo_ok && bgo_corr_ok ) {

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
      

      xPos_bgo_ = sumVector( xPosW_bgo )/eTot_bgo_corr;
      yPos_bgo_ = sumVector( yPosW_bgo )/eTot_bgo_corr;
      
      h1_xPos_bgo->Fill( xPos_bgo_ );
      h1_yPos_bgo->Fill( yPos_bgo_ );
      h2_xyPos_bgo->Fill( xPos_bgo_, yPos_bgo_ );

      h2_correlation_hodo_bgo_xPos->Fill( xPos_hodo, xPos_bgo_ );
      h2_correlation_hodo_bgo_yPos->Fill( yPos_hodo, yPos_bgo_ );
      
      if( isSingleEle_scintFront ) {

        h1_xPos_singleEle_bgo->Fill( xPos_bgo_ );
        h1_yPos_singleEle_bgo->Fill( yPos_bgo_ );
        h2_xyPos_singleEle_bgo->Fill( xPos_bgo_, yPos_bgo_ );

        h2_correlation_hodo_bgo_xPos_singleEle->Fill( xPos_hodo, xPos_bgo_ );
        h2_correlation_hodo_bgo_yPos_singleEle->Fill( yPos_hodo, yPos_bgo_ );

      }
      
    }  // if bgo ok



    // THEN USE CeF3 DATA:

    if( cef3_ok ) {


      std::vector<float> v_cef3_corr;
      for(unsigned i=0; i<CEF3_CHANNELS; ++i) v_cef3_corr.push_back(cef3_corr[i]);

      float eTot_corr = sumVector(v_cef3_corr);

      

      if( cef3_corr_ok ) {


        //   0      1
        //          
        //          
        //   3      2


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


        getCeF3Position( v_cef3_corr, xPos_new_, yPos_new_ );
        //diag02_new_ = xPos_new_;
        //diag13_new_ = yPos_new_;

        float pi = 3.14159;
        float theta = pi/4.; // 45 degrees 
        TVector2 vnew( xPos_new_, yPos_new_ );
        TVector2 dnew = vnew.Rotate(-theta);
        diag02_new_ = dnew.Y();
        diag13_new_ = dnew.X();


        TVector2 vBeam( xBeam, yBeam );
        TVector2 dBeam = vBeam.Rotate(-theta);
        diag02_beam_ = dBeam.Y();
        diag13_beam_ = dBeam.X();

        float xPos = sumVector(xPosW)/eTot_corr;
        float yPos = sumVector(yPosW)/eTot_corr;

        h1_xPos->Fill( xPos );
        h1_yPos->Fill( yPos );

        h2_xyPos->Fill( xPos, yPos );

        h1_xPos_new->Fill( xPos_new_ );
        h1_yPos_new->Fill( yPos_new_ );

        h2_xyPos_new->Fill( xPos_new_, yPos_new_ );


        // positioning with all 9 calorimeter channels:
        //float xPos_calo = sumVector( xPosW_bgo )/(eTot_bgo_corr + eTot_corr*0.07); // cef3 is in 0,0
        //float yPos_calo = sumVector( yPosW_bgo )/(eTot_bgo_corr + eTot_corr*0.08); // so counts only in denominator
        xPos_calo_ = sumVector( xPosW_bgo )/(eTot_bgo_corr + eTot_corr*0.06); // cef3 is in 0,0
        yPos_calo_ = sumVector( yPosW_bgo )/(eTot_bgo_corr + eTot_corr*0.10); // so counts only in denominator
        //float xPos_calo = sumVector( xPosW_bgo )/(eTot_bgo_corr + eTot_corr*0.791577); // cef3 is in 0,0
        //float yPos_calo = sumVector( yPosW_bgo )/(eTot_bgo_corr + eTot_corr*0.791577); // so counts only in denominator

        TVector2 vcalo( xPos_calo_, yPos_calo_ );
        TVector2 dcalo = vcalo.Rotate(-theta);
        diag02_calo_ = dcalo.Y();
        diag13_calo_ = dcalo.X();

        //xPos_regr2D_ = readerRegr2D->EvaluateRegression( "MLP" )[0];
        //yPos_regr2D_ = readerRegr2D->EvaluateRegression( "MLP" )[1];

        cef3_regr[0] = cef3_corr[0]/(cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]);
        cef3_regr[1] = cef3_corr[1]/(cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]);
        cef3_regr[2] = cef3_corr[2]/(cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]);
        cef3_regr[3] = cef3_corr[3]/(cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]);
  
        if( bgo_ok && bgo_corr_ok ) {

          h1_xPos_calo->Fill( xPos_calo_ );
          h1_yPos_calo->Fill( yPos_calo_ );
          h2_xyPos_calo->Fill( xPos_calo_, yPos_calo_ );

          //for( unsigned i=0; i<methodNames.size(); ++i ) {
          //  Float_t xPos_regr = (readerRegrX->EvaluateRegression( methodNames[i] ))[0];
          //  Float_t yPos_regr = (readerRegrY->EvaluateRegression( methodNames[i] ))[0];
          //  h1_xPos_regr_vs_calo[i]->Fill( xPos_regr-xPos_calo_ );
          //  h1_yPos_regr_vs_calo[i]->Fill( yPos_regr-yPos_calo_ );
          //  h2_xyPos_regr[i]->Fill( xPos_regr, yPos_regr );
          //}
          //h2_xyPos_regr2D->Fill( xPos_regr2D_, yPos_regr2D_ );

          h1_xPos_calo_vs_hodo->Fill( xPos_calo_-xPos_hodo );
          h1_yPos_calo_vs_hodo->Fill( yPos_calo_-yPos_hodo );

          h1_xPos_calo_vs_beam->Fill( xPos_calo_-xBeam );
          h1_yPos_calo_vs_beam->Fill( yPos_calo_-yBeam );

          // CORRELATIONS BETWEEN CALO AND HODO:

          h2_correlation_cef3_bgo_xPos->Fill( xPos, xPos_bgo_ );
          h2_correlation_cef3_bgo_yPos->Fill( yPos, yPos_bgo_ );

        }

        if( isSingleEle_scintFront ) {

          h1_xPos_singleEle->Fill( xPos );
          h1_yPos_singleEle->Fill( yPos );
          h2_xyPos_singleEle->Fill( xPos, yPos );

          h1_xPos_new_singleEle->Fill( xPos_new_ );
          h1_yPos_new_singleEle->Fill( yPos_new_ );
          h2_xyPos_new_singleEle->Fill( xPos_new_, yPos_new_ );

          h1_xPos_calo_vs_hodo_singleElectron->Fill( xPos_calo_-xPos_hodo );
          h1_yPos_calo_vs_hodo_singleElectron->Fill( yPos_calo_-yPos_hodo );

          h1_xPos_calo_vs_beam_singleElectron->Fill( xPos_calo_-xBeam );
          h1_yPos_calo_vs_beam_singleElectron->Fill( yPos_calo_-yBeam );

          if( nHodoClustersX==1 && nHodoClustersY==1 && nFibres_hodoClustX[0]<=2 && nFibres_hodoClustY[0]<=2 ) {

            h1_xPos_calo_vs_hodo_singleElectron_HR->Fill( xPos_calo_-xPos_hodo );
            h1_yPos_calo_vs_hodo_singleElectron_HR->Fill( yPos_calo_-yPos_hodo );
  
            h1_xPos_calo_vs_beam_singleElectron_HR->Fill( xPos_calo_-xBeam );
            h1_yPos_calo_vs_beam_singleElectron_HR->Fill( yPos_calo_-yBeam );

          }


          h2_correlation_cef3_hodo_xPos_singleEle->Fill( xPos, xPos_hodo );
          h2_correlation_cef3_hodo_yPos_singleEle->Fill( yPos, yPos_hodo );

  
          if( bgo_ok && bgo_corr_ok ) {

            h1_xPos_singleEle_calo->Fill( xPos_calo_ );
            h1_yPos_singleEle_calo->Fill( yPos_calo_ );
            h2_xyPos_singleEle_calo->Fill( xPos_calo_, yPos_calo_ );

            //for( unsigned i=0; i<methodNames.size(); ++i ) {
            //  Float_t xPos_regr = (readerRegrX->EvaluateRegression( methodNames[i] ))[0];
            //  Float_t yPos_regr = (readerRegrY->EvaluateRegression( methodNames[i] ))[0];
            //  h1_xPos_regr_vs_calo_singleEle[i]->Fill( xPos_regr-xPos_calo_ );
            //  h1_yPos_regr_vs_calo_singleEle[i]->Fill( yPos_regr-yPos_calo_ );
            //  h2_xyPos_singleEle_regr[i]->Fill( xPos_regr, yPos_regr );
            //}

            //h2_xyPos_singleEle_regr2D->Fill( xPos_regr2D_, yPos_regr2D_ );

            h2_correlation_cef3_bgo_xPos_singleEle->Fill( xPos, xPos_bgo_ );
            h2_correlation_cef3_bgo_yPos_singleEle->Fill( yPos, yPos_bgo_ );

          }

        }

        
        outTree->Fill();


      } // if cef3_ok

    }

  }


  outfile->cd();

  outTree->Write();


  for( unsigned i=0; i<h1_xPos_regr_vs_calo.size(); ++i ) {
    h1_xPos_regr_vs_calo[i]->Write();
    h1_yPos_regr_vs_calo[i]->Write();
    h1_xPos_regr_vs_calo_singleEle[i]->Write();
    h1_yPos_regr_vs_calo_singleEle[i]->Write();
    h2_xyPos_regr[i]->Write();
    h2_xyPos_singleEle_regr[i]->Write();
  }




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

  
  
  h1_xPos_bgo->Write();
  h1_yPos_bgo->Write();
  h2_xyPos_bgo->Write();
  

  h1_xPos_hodo->Write();
  h1_yPos_hodo->Write();
  h2_xyPos_hodo->Write();

  
  h1_xPos_calo_vs_hodo->Write();
  h1_yPos_calo_vs_hodo->Write();

  h1_xPos_calo_vs_beam->Write();
  h1_yPos_calo_vs_beam->Write();
  
  h1_xPos_calo_vs_hodo_singleElectron->Write();
  h1_yPos_calo_vs_hodo_singleElectron->Write();

  h1_xPos_calo_vs_beam_singleElectron->Write();
  h1_yPos_calo_vs_beam_singleElectron->Write();
  
  h1_xPos_calo_vs_hodo_singleElectron_HR->Write();
  h1_yPos_calo_vs_hodo_singleElectron_HR->Write();

  h1_xPos_calo_vs_beam_singleElectron_HR->Write();
  h1_yPos_calo_vs_beam_singleElectron_HR->Write();
  

  std::cout << std::endl;

  

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


  outfile->Close();
  std::cout << "-> Histograms saved in: " << outfile->GetName() << std::endl;


  return 0;

}







float sumVector( std::vector<float> v ) {

  float sum=0.;
  for( unsigned i=0; i<v.size(); ++i ) sum += v[i];

  return sum;

}

float sumVector( int n, float* x ) {

  float sum=0.;
  for( unsigned i=0; i<n; ++i ) sum += x[i];

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



float getMeanposHodo( int n, float* pos ) {

  if( n==0 ) return -999.;

  float meanpos = 0.;

  for( unsigned i=0; i<n; ++i ) {
    meanpos += pos[n];
  }

  meanpos /= (float)n;

  return meanpos;

}




//float gethodointercalib(TString axis, int n){
//
//  float res=1;
//
//  if (axis==TString("X")){
//    if (n>=0 && n<=7) return hodo_efficiency_vector_X[n];
//    else assert(false);
//  }
//  else if (axis==TString("Y")){
//    if (n>=0 && n<=7) return hodo_efficiency_vector_Y[n];
//    else assert(false);
//  }
//  else assert(false);
//  
//  return 1./res;
//  
//};





void getCeF3Position( std::vector<float> cef3, float& xPos, float& yPos ) {

  xPos=0.;
  yPos=0.;

  float r02 = cef3[0]/cef3[2];
  float r13 = cef3[1]/cef3[3];

  TF1* f1_d02 = getCef3Function("diag02");
  TF1* f1_d13 = getCef3Function("diag13");

  float diag02 = f1_d02->GetX(r02, -50., 50.);
  float diag13 = f1_d13->GetX(r13, -50., 50.);

  TVector2 v( diag13, diag02 );
  float pi = 3.14159;
  float theta = pi/4.; // 45 degrees 
  TVector2 d = v.Rotate(+theta);

  xPos = d.X();
  yPos = d.Y();
  //xPos = diag02;
  //yPos = diag13;


  delete f1_d02;
  delete f1_d13;

}


TF1* getCef3Function( const std::string& diagName ) {

  std::string fName = "f1_" + diagName;

  //ifstream ifs(fileName.c_str());
  TF1* f1 = new TF1(fName.c_str(), "pol5", -100., 100.);


  if( diagName=="diag02" ) {
    f1->SetParameter( 0, 1.03747 );
    f1->SetParameter( 1, 0.00895044 );
    f1->SetParameter( 2, 0 );
    f1->SetParameter( 3, 0.000125462 );
    f1->SetParameter( 4, 2.1783e-06 );
    f1->SetParameter( 5, 0 );
  } else if( diagName=="diag13" ) {
    f1->SetParameter( 0, 1.06795 );
    f1->SetParameter( 1, 0.011661 );
    f1->SetParameter( 2, 0.000461872 );
    f1->SetParameter( 3, 0.000142064 );
    f1->SetParameter( 4, 1.90787e-05 );
    f1->SetParameter( 5, 1.0005e-06 );
    //f1->SetParameter( 0, 1.07932 );
    //f1->SetParameter( 1, 0.014342 );
    //f1->SetParameter( 2, 0.000722962 );
    //f1->SetParameter( 3, 0.000192131 );
    //f1->SetParameter( 4, 1.86405e-05 );
    //f1->SetParameter( 5, 8.16208e-07 );
  } else {
    std::cout << "WARNING! Unkown diagName: " << diagName << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11111);
  }

  //int i=0;
  //while( i<4 ) {
  //  float x;
  //  ifs >> x;
  //  f1->SetParameter(i, x);
  //  i++;
  //}

  return f1;

}



