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

//This makes the origninal EEShash*.root trees into Reco_Simulation*.root trees

//Usage:
// ./makeSimulationToAnalysisTree * tag     for normal simulation
// ./makeSimulationToAnalysisTree Ideal_* tag         for "ideal"=5x5 setup

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
    std::cout<<"./makeSimulationToAnalysisTree BeamEnergy(in MeV) tag"<<std::endl;
    exit(12345);
  }

  TString runName_tstr(runName);


  bool isIdeal = false;

  //run name = energy for now, might have to change this for different beams...
    std::string fileName = "OriginalSimulationData/EEShash" + runName + ".root";
    TFile* file = TFile::Open(fileName.c_str());
    if( file==0 ) {
    std::string fileName = "OriginalSimulationData/EEShashIdeal_" + runName + ".root";
    TFile* file = TFile::Open(fileName.c_str());
    isIdeal = true;
    }
    if( file==0 ) {
      std::cout << "ERROR! Din't find file " << fileName << std::endl;
      std::cout << "Exiting." << std::endl;
      exit(11);
    }

  TChain* tree = new TChain("EEShash");
  tree = (TChain*)file->Get("EEShash");





  //std::string outdir = "analysisTrees_" + tag;
  std::string outdir = "OriginalSimulationData";
  system( Form("mkdir -p %s", outdir.c_str()) );
  std::string outfileName;
  if(isIdeal==0){ outfileName = outdir + "/Reco_Simulation" + runName + ".root";}
  else { outfileName = outdir + "/Reco_Simulation_Ideal_" + runName + ".root";}
  TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );
  

  TTree* outTree = new TTree("recoTree","recoTree");

  double Eact_0;
  tree->SetBranchAddress("Eact_0", &Eact_0);
  double Eact_1;
  tree->SetBranchAddress("Eact_1", &Eact_1);
  double Eact_2;
  tree->SetBranchAddress("Eact_2", &Eact_2);
  double Eact_3;
  tree->SetBranchAddress("Eact_3", &Eact_3);
  double Eact_4;
  tree->SetBranchAddress("Eact_4", &Eact_4);
  double Eact_5;
  tree->SetBranchAddress("Eact_5", &Eact_5);
  double Eact_6;
  tree->SetBranchAddress("Eact_6", &Eact_6);
  double Eact_7;
  tree->SetBranchAddress("Eact_7", &Eact_7);
  double Eact_8;
  tree->SetBranchAddress("Eact_8", &Eact_8);
  double Eact_9;
  tree->SetBranchAddress("Eact_9", &Eact_9);

  double Ebgo_0;
  tree->SetBranchAddress( "Ebgo_0", &Ebgo_0);
  double Ebgo_1;
  tree->SetBranchAddress( "Ebgo_1", &Ebgo_1);
  double Ebgo_2;
  tree->SetBranchAddress( "Ebgo_2", &Ebgo_2);
  double Ebgo_3;
  tree->SetBranchAddress("Ebgo_3", &Ebgo_3);
  double Ebgo_4;
  tree->SetBranchAddress("Ebgo_4", &Ebgo_4);
  double Ebgo_5;
  tree->SetBranchAddress("Ebgo_5", &Ebgo_5);
  double Ebgo_6;
  tree->SetBranchAddress("Ebgo_6", &Ebgo_6);
  double Ebgo_7;
  tree->SetBranchAddress("Ebgo_7", &Ebgo_7);

  double Ebgo;
  tree->SetBranchAddress("Ebgo", &Ebgo);

  double Eabs;
  tree->SetBranchAddress("Eabs", &Eabs);
  double Eact;
  tree->SetBranchAddress("Eact", &Eact);

  double Escint;
  tree->SetBranchAddress("Escint", &Escint);

  double Ehodo;
  tree->SetBranchAddress("Ehodo", &Ehodo);



  float cef3_corr_[CEF3_CHANNELS],bgo_corr_[BGO_CHANNELS];
  bool isSingleEle_scintFront_;
  int nHodoClustersY; int nHodoClustersX;
  float E_abs;


  int bgo_chan=BGO_CHANNELS;
  int cef3_chan=CEF3_CHANNELS;

  outTree->Branch( "nHodoClustersX", &nHodoClustersX, "nHodoClustersX/I" );
  outTree->Branch( "nHodoClustersY", &nHodoClustersY, "nHodoClustersY/I" );

  outTree->Branch( "cef3_chan", &cef3_chan, "cef3_chan/I" );
  outTree->Branch( "bgo_chan", &bgo_chan, "bgo_chan/I" );

  outTree->Branch( "bgo_corr", bgo_corr_, "bgo_corr_[bgo_chan]/F" );
  outTree->Branch( "cef3_corr", cef3_corr_, "cef3_corr_[cef3_chan]/F" );

  outTree->Branch( "isSingleEle_scintFront", &isSingleEle_scintFront_, "isSingleEle_scintFront_/O" );

  outTree->Branch( "E_abs", &E_abs, "E_abs/F");

  float LYSF[] = {0.85, 0.94, 0.95, 0.98, 1.00, 1.02, 1.05, 1.05, 1.05, 0.74};
  /*
    std::string fullVarName = "";
    for( unsigned i=0; i<10; ++i ) {
    std::string plusSign = (i==0) ? "" : " + ";
    std::string thisPiece(Form("%s%f*Eact_%d", plusSign.c_str(), LYSF[i], i));
    fullVarName += thisPiece;
    }
  */
  
  
  unsigned nentries = tree->GetEntries();
  
  double hodoEff;
  double scintEff;
  
  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {
    
    tree->GetEntry(iEntry);
    
    if( iEntry % 5000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;
    
    if(Ehodo>0.1){
      nHodoClustersY=1; nHodoClustersX=1;
      hodoEff++;
    } else{   nHodoClustersY=0; nHodoClustersX=0;}
    
    if(Escint>0.25){
      isSingleEle_scintFront_=1;
      scintEff++;
    } else{       isSingleEle_scintFront_=0;}
    
       
    E_abs= Eabs;
    
    std::vector<float> cef3_corr;
    for( int i=0; i<CEF3_CHANNELS; ++i ) cef3_corr.push_back(-1.);
    
    std::vector<float> bgo_corr;
    for( int i=0; i<BGO_CHANNELS; ++i ) bgo_corr.push_back(-1.);
    
    for(int i=0; i<10; ++i){
      cef3_corr[0]= (LYSF[0]*Eact_0+LYSF[1]*Eact_1+LYSF[2]*Eact_2 +LYSF[3]*Eact_3 +LYSF[4]*Eact_4 +LYSF[5]*Eact_5 +LYSF[6]*Eact_6 +LYSF[7]*Eact_7 +LYSF[8]*Eact_8 +LYSF[9]*Eact_9)/4. ; //yes not very elegant, but it didn't want to work otherwise
      // cef3_corr[0] = Eact /4.;
    }
    
    
    bgo_corr[0] = Ebgo/8.;
    bgo_corr[1] = Ebgo/8.;
    bgo_corr[2] = Ebgo/8.;
    bgo_corr[3] = Ebgo/8.;
    bgo_corr[4] = Ebgo/8.;
    bgo_corr[5] = Ebgo/8.;
    bgo_corr[6] = Ebgo/8.;
    bgo_corr[7] = Ebgo/8.;
    
    
    for(int k=0; k<CEF3_CHANNELS; ++k){
      cef3_corr_[k] = cef3_corr[0];
    }
    
    for (int k=0; k<BGO_CHANNELS; ++k){ 
      bgo_corr_[k] = bgo_corr[k];
    }
    
    
    outTree->Fill();
  }
  
  
  outfile->cd();
  outTree->Write();
  outfile->Close();
  
  std::cout << "-> Simulated Analysis Tree saved in: " << outfile->GetName() << std::endl;

  std::cout << "Scintillator efficiency = " << scintEff/nentries << std::endl;
  std::cout << "Hodoscope efficiency = " << hodoEff/nentries << std::endl;
 

  return 0;

}






