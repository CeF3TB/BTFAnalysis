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
#include "TCanvas.h"
#include "TString.h"

#include "TMVA/Reader.h"

#include "fastDQM_CeF3_BTF.h"
#include "interface/RunHelper.h"
#include "interface/PositionTools.h"
#include "interface/DrawTools.h"



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

  std::string cuts="";
  if( argc>3 ) {
    std::string cuts_str(argv[3]);
    cuts = cuts_str;
  }
  std::cout<<"cuts="<<cuts<<std::endl;

  TString cuts_tstr(cuts);

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



   UInt_t          runNumber;
   UInt_t          evtNumber;
   UInt_t          adcData;
   UInt_t          adcBoard;
   UInt_t          adcChannel;
   Int_t           nHodoFibersX;
   Int_t           nHodoFibersY;
   Int_t           nHodoClustersX;
   Float_t         pos_hodoClustX[CEF3_CHANNELS];   //[nHodoClustersX]
   Int_t           nFibres_hodoClustX[CEF3_CHANNELS];   //[nHodoClustersX]
   Int_t           nHodoClustersY;
   Float_t         pos_hodoClustY[CEF3_CHANNELS];   //[nHodoClustersY]
   Int_t           nFibres_hodoClustY[CEF3_CHANNELS];   //[nHodoClustersY]
   Int_t           hodox_chan;
   Int_t           hodoy_chan;
   Int_t           cef3_chan;
   Int_t           bgo_chan;
   Float_t         hodox[BGO_CHANNELS];   //[hodox_chan]
   Float_t         hodoy[BGO_CHANNELS];   //[hodoy_chan]
   Float_t         bgo[BGO_CHANNELS];   //[bgo_chan]
   Float_t         cef3[CEF3_CHANNELS];   //[cef3_chan]
   Float_t         cef3_pedSubtracted[CEF3_CHANNELS];   //[cef3_chan]
   Float_t         bgo_corr[BGO_CHANNELS];   //[bgo_chan]
   Float_t         bgo_pedSubtracted[BGO_CHANNELS];   //[bgo_chan]
   Float_t         cef3_corr[CEF3_CHANNELS];   //[cef3_chan]
   Float_t         hodox_corr[BGO_CHANNELS];   //[hodox_chan]
   Float_t         hodoy_corr[BGO_CHANNELS];   //[hodoy_chan]
   Float_t         scintFront;
   Bool_t          isSingleEle_scintFront;
   Float_t         xBeam;
   Float_t         yBeam;
   Bool_t          cef3_ok;
   Bool_t          cef3_corr_ok;
   Bool_t          bgo_ok;
   Bool_t          bgo_corr_ok;

   tree->SetBranchAddress("runNumber", &runNumber);
   tree->SetBranchAddress("evtNumber", &evtNumber);
   tree->SetBranchAddress("adcData", &adcData);
   tree->SetBranchAddress("adcBoard", &adcBoard);
   tree->SetBranchAddress("adcChannel", &adcChannel);
   tree->SetBranchAddress("nHodoFibersX", &nHodoFibersX);
   tree->SetBranchAddress("nHodoFibersY", &nHodoFibersY);
   tree->SetBranchAddress("nHodoClustersX", &nHodoClustersX);
   tree->SetBranchAddress("pos_hodoClustX", pos_hodoClustX);
   tree->SetBranchAddress("nFibres_hodoClustX", nFibres_hodoClustX);
   tree->SetBranchAddress("nHodoClustersY", &nHodoClustersY);
   tree->SetBranchAddress("pos_hodoClustY", pos_hodoClustY);
   tree->SetBranchAddress("nFibres_hodoClustY", nFibres_hodoClustY);
   tree->SetBranchAddress("hodox_chan", &hodox_chan);
   tree->SetBranchAddress("hodoy_chan", &hodoy_chan);
   tree->SetBranchAddress("cef3_chan", &cef3_chan);
   tree->SetBranchAddress("bgo_chan", &bgo_chan);
   tree->SetBranchAddress("hodox", hodox);
   tree->SetBranchAddress("hodoy", hodoy);
   tree->SetBranchAddress("bgo", bgo);
   tree->SetBranchAddress("cef3", cef3);
   tree->SetBranchAddress("cef3_pedSubtracted", cef3_pedSubtracted);
   tree->SetBranchAddress("bgo_corr", bgo_corr);
   tree->SetBranchAddress("bgo_pedSubtracted", bgo_pedSubtracted);
   tree->SetBranchAddress("cef3_corr", cef3_corr);
   tree->SetBranchAddress("hodox_corr", hodox_corr);
   tree->SetBranchAddress("hodoy_corr", hodoy_corr);
   tree->SetBranchAddress("scintFront", &scintFront);
   tree->SetBranchAddress("isSingleEle_scintFront", &isSingleEle_scintFront);
   tree->SetBranchAddress("xBeam", &xBeam);
   tree->SetBranchAddress("yBeam", &yBeam);
   tree->SetBranchAddress("cef3_ok", &cef3_ok);
   tree->SetBranchAddress("cef3_corr_ok", &cef3_corr_ok);
   tree->SetBranchAddress("bgo_ok", &bgo_ok);
   tree->SetBranchAddress("bgo_corr_ok", &bgo_corr_ok);


   int nentries = tree->GetEntries();

  if( isOnlyRunNumber ) {
    // modify runname in such a way that it's useful for getBeamPosition and outfile:
    runName = "BTF_" + runName + "_beam";
  }

  std::string outputdir = "2DAnTrees_"+tag;

  system( Form("mkdir -p %s", outputdir.c_str()) );
  std::string outfileName = outputdir + "/2DAn_" + runName + ".root";
  TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );

  TStyle* style = DrawTools::setStyle();
  //  tree->Print();

  TH2F* cef3VsScintFront= new TH2F("cef3VsScintFront","cef3VsScintFront",40,0.,4000.,40,0.,12000.);
  TH2F* cef3TopVsBottom= new TH2F("cefTopVsBottom","cefTopVsBottom",40,0.,8000.,40,0.,8000.);
  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {
    tree->GetEntry(iEntry);
    if( iEntry % 10000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;


    if(cuts=="hodo"){
      if(!(nHodoClustersX>=1 && nHodoClustersY>=1))continue;
    }else if (cuts=="signalFibers"){
      if(!(cef3_corr[0]>0 && cef3_corr[1]>0 && cef3_corr[2]>0 && cef3_corr[3]>0))continue;
    }else if (cuts=="bgo"){
      if(!((bgo_corr[0]+bgo_corr[1]+bgo_corr[2]+bgo_corr[3]+bgo_corr[4]+bgo_corr[5]+bgo_corr[6]+bgo_corr[7])<1000))continue;
    }

    cef3VsScintFront->Fill(scintFront,cef3[0]+cef3[1]+cef3[2]+cef3[3]);
    cef3TopVsBottom->Fill(cef3_corr[2]+cef3_corr[3],cef3_corr[0]+cef3_corr[1]);
  }

  style->SetPadTopMargin(0.05);
  style->SetPadBottomMargin(0.13);//0.13);                                                                                                                   
  style->SetPadLeftMargin(0.17);//0.16);                                                                                                                     
  style->SetPadRightMargin(0.13);//0.02); 

  TCanvas* c1= new TCanvas("c1"," ", 600,600);
  c1->cd();
  cef3VsScintFront->SetYTitle("#Sigma CeF_{3} Fibers");
  cef3VsScintFront->SetTitleOffset(1.7,"Y");
  cef3VsScintFront->SetXTitle("scintFront");
  cef3VsScintFront->Draw("colz");  
  TPaveText* labelTop = DrawTools::getLabelTop();
  labelTop->Draw("same");
  c1->SetLogz();
  
  c1->SaveAs("cef3VsScintFront_"+cuts_tstr+".png");
  c1->SaveAs("cef3VsScintFront_"+cuts_tstr+".eps");

  c1->Clear();
  c1->cd();
  cef3TopVsBottom->SetYTitle("CeF_{3} Fibers Top");
  cef3TopVsBottom->SetTitleOffset(1.7,"Y");
  cef3TopVsBottom->SetXTitle("CeF_{3} Fibers Bottom");
  cef3TopVsBottom->Draw("colz");  
  labelTop->Draw("same");
  c1->SetLogz();
  
  c1->SaveAs("cef3TopVsBottom_"+cuts_tstr+".png");
  c1->SaveAs("cef3TopVsBottom_"+cuts_tstr+".eps");

}
