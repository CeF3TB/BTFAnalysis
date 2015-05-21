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


//This makes the origninal EEShash*.root trees into Reco_Simulation*.root trees

//Usage:
// ./makeSimulationToAnalysisTree * tag     for normal simulation
// ./makeSimulationToAnalysisTree Ideal_* tag         for "ideal"=5x5 setup if you name the shit that way... you can also make your life easy and just don't do that...

int main( int argc, char* argv[] ) {

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
  double Eact_10;
  tree->SetBranchAddress("Eact_10", &Eact_10);

  double Eact_11;
  tree->SetBranchAddress("Eact_11", &Eact_11);
  double Eact_12;
  tree->SetBranchAddress("Eact_12", &Eact_12);
  double Eact_13;
  tree->SetBranchAddress("Eact_13", &Eact_13);
  double Eact_14;
  tree->SetBranchAddress("Eact_14", &Eact_14);


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

  double Escint1;
  tree->SetBranchAddress("Escint1", &Escint1);

  double Ehodo11;
  tree->SetBranchAddress("Ehodo11", &Ehodo11);

  double Ehodo12;
  tree->SetBranchAddress("Ehodo12", &Ehodo12);


  double Escint;
  tree->SetBranchAddress("Escint", &Escint);

  double Ehodo;
  tree->SetBranchAddress("Ehodo", &Ehodo);

  double xPosition;
  tree->SetBranchAddress("xPosition",&xPosition);
  double yPosition;
  tree->SetBranchAddress("yPosition",&yPosition);


  float cef3_corr_[CEF3_CHANNELS],bgo_corr_[BGO_CHANNELS];
  bool isSingleEle_scintFront_=1;
  int nHodoClustersY=1; int nHodoClustersX=1;
  int nHodoClusters1Y=1; int nHodoClusters1X=1;
  float E_abs;

  float xPos;
  float yPos;

  int bgo_chan=BGO_CHANNELS;
  int cef3_chan=CEF3_CHANNELS;

  outTree->Branch( "nHodoClustersX", &nHodoClustersX, "nHodoClustersX/I" );
  outTree->Branch( "nHodoClustersY", &nHodoClustersY, "nHodoClustersY/I" );
  outTree->Branch( "nHodoClusters1X", &nHodoClusters1X, "nHodoClusters1X/I" );
  outTree->Branch( "nHodoClusters1Y", &nHodoClusters1Y, "nHodoClusters1Y/I" );

  outTree->Branch( "cef3_chan", &cef3_chan, "cef3_chan/I" );
  outTree->Branch( "bgo_chan", &bgo_chan, "bgo_chan/I" );

  outTree->Branch( "bgo_corr", bgo_corr_, "bgo_corr_[bgo_chan]/F" );
  outTree->Branch( "cef3_corr", cef3_corr_, "cef3_corr_[cef3_chan]/F" );

  outTree->Branch( "isSingleEle_scintFront", &isSingleEle_scintFront_, "isSingleEle_scintFront_/O" );

  outTree->Branch( "E_abs", &E_abs, "E_abs/F");

  outTree->Branch( "xPos", &xPos, "xPos/F");
  outTree->Branch( "yPos", &yPos, "yPos/F");


  float LYSF[] = {0.85, 0.94, 0.95, 0.98, 1.00, 1.02, 1.05, 1.05, 1.05, 0.74};

  // float LYSF[] = {0.85, 0.94, 0.95, 0.98, 1.00, 0.7, 1.05, 1.05, 1.05, 0.74};

  // float LYSF[] =  {0.87, 0.96, 0.97, 1.00, 1.02, 1.04, 1.07, 1.07, 1.08, 0.75};
  //as a test with the best being the best one
  //  float LYSF[] =  {0.81, 0.89, 0.90, 0.93, 0.95, 0.97, 1.00, 1.00, 1.00, 0.70};


  //Test if 6.crystal fails
  //  float LYSFH[] = {0.87, 0.96, 0.97, 1.00, 1.02, 0.50, 1.07, 1.07, 1.08, 1.07, 1.08, 1.09, 1.09, 1.04, 0.75};
  ///



   float LYSFH[] = {0.87, 0.96, 0.97, 1.00, 1.02, 1.04, 1.07, 1.07, 1.08, 1.07, 1.08, 1.09, 1.09, 1.04, 0.75};

  /*
    std::string fullVarName = "";
    for( unsigned i=0; i<15; ++i ) {
    std::string plusSign = (i==0) ? "" : " + ";
    std::string thisPiece(Form("%s%f*Eact_%d", plusSign.c_str(), LYSFH[i], i));
    fullVarName += thisPiece;
    }
  
  */
  
  unsigned nentries = tree->GetEntries();
  
  double hodoEff;
  double scintEff;
  
  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {
    
    tree->GetEntry(iEntry);
    
    if( iEntry % 5000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;
 
 
        
    xPos = xPosition;
    yPos = yPosition;
    
    if(Ehodo11>0.1){
       nHodoClustersX=1;
      hodoEff++;
    } else{   nHodoClusters1X=0;}

    if(Ehodo12>0.1){
      nHodoClustersY=1; 
      hodoEff++;
    } else{   nHodoClusters1Y=0; }
   
    
    if(Escint1>0.25){
      isSingleEle_scintFront_=1;
      scintEff++;
    } else{       isSingleEle_scintFront_=0;}
    
   
 
   


    if(Ehodo>0.1){
       nHodoClustersX=1; nHodoClustersY=1;
      hodoEff++;
    } else{   nHodoClustersX=0;  nHodoClustersY=0;}
  
    /*
    if(Escint>0.25 && Escint< 1.5){
      isSingleEle_scintFront_=1;
      scintEff++;
    } else{       isSingleEle_scintFront_=0;}     
    */ 
    
    E_abs= Eabs;


    
    std::vector<float> cef3_corr;
    for( int i=0; i<CEF3_CHANNELS; ++i ) cef3_corr.push_back(-1.);
    
    std::vector<float> bgo_corr;
    for( int i=0; i<BGO_CHANNELS; ++i ) bgo_corr.push_back(-1.);
 
    /*
     
     for(int i=0; i<10; ++i){
     cef3_corr[0]= (LYSF[0]*Eact_0+LYSF[1]*Eact_1+LYSF[2]*Eact_2 +LYSF[3]*Eact_3 +LYSF[4]*Eact_4 +LYSF[5]*Eact_5 +LYSF[6]*Eact_6 +LYSF[7]*Eact_7 +LYSF[8]*Eact_8 +LYSF[9]*Eact_9)/4. ; //yes not very elegant, but it didn't want to work otherwise
     // cef3_corr[0] = Eact /4.;
     }
    */
      
  
    
    for(int i=0; i<15; ++i){
      cef3_corr[0]= (LYSFH[0]*Eact_0+LYSFH[1]*Eact_1+LYSFH[2]*Eact_2 +LYSFH[3]*Eact_3 +LYSFH[4]*Eact_4 +LYSFH[5]*Eact_5 +LYSFH[6]*Eact_6 +LYSFH[7]*Eact_7 +LYSFH[8]*Eact_8 +LYSFH[9]*Eact_9 +LYSFH[10]*Eact_10  +LYSFH[11]*Eact_11 +LYSFH[12]*Eact_12 +LYSFH[13]*Eact_13  +LYSFH[14]*Eact_14 )/4. ; //yes not very elegant, but it didn't want to work otherwise

      //  cef3_corr[0] = Eact /4.;
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






