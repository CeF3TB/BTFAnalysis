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
#include "interface/PositionTools.h"


float sumVector( std::vector<float> v );

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

  int nentries = tree->GetEntries();


  if( isOnlyRunNumber ) {
    // modify runname in such a way that it's useful for getBeamPosition and outfile:
    runName = "BTF_" + runName + "_beam";
  }

  std::string outputdir = "SingleElectronSelectionTrees_"+tag;
  system( Form("mkdir -p %s", outputdir.c_str()) );
  std::string outfileName = outputdir + "/SingleEleSelAn_" + runName + ".root";
  TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );

  TTree* outTree = new TTree("singleEleSelTree","singleEleSelTree");

  outTree->Branch( "run", &runNumber, "run/i" );
  outTree->Branch( "scintFront", &scintFront, "scintFront/F" );
  outTree->Branch( "isSingleEle_scintFront", &isSingleEle_scintFront, "isSingleEle_scintFront/O" );
  outTree->Branch( "nHodoClustersX", &nHodoClustersX, "nHodoClustersX/I" );
  outTree->Branch( "nHodoClustersY", &nHodoClustersY, "nHodoClustersY/I" );
  outTree->Branch( "cef3_corr", cef3_corr, "cef3_corr[4]/F" );
  outTree->Branch( "bgo_corr", bgo_corr, "bgo_corr[8]/F" );
  outTree->Branch( "bgo_corr_ok", &bgo_corr_ok, "bgo_corr_ok/O");
  outTree->Branch( "xBeam", &xBeam, "xBeam/F" );
  outTree->Branch( "yBeam", &yBeam, "yBeam/F" );

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {
    tree->GetEntry(iEntry);
    
    if( iEntry % 10000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;
    if( cef3_ok ) {
      
      
      std::vector<float> v_cef3_corr;
      for(unsigned i=0; i<CEF3_CHANNELS; ++i) v_cef3_corr.push_back(cef3_corr[i]);
      
      float eTot_corr = sumVector(v_cef3_corr);

      if( cef3_corr_ok ) {


      }
      outTree->Fill();
    }

  }

  outfile->cd();
  outTree->Write();
  outfile->Close();
  std::cout << "-> Histograms saved in: " << outfile->GetName() << std::endl;


  return 0;

}

float sumVector( std::vector<float> v ) {

  float sum=0.;
  for( unsigned i=0; i<v.size(); ++i ) sum += v[i];

  return sum;

}
