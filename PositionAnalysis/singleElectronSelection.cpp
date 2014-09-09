#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <stdlib.h>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLegend.h"

#include "TMath.h"
#include "TTree.h"

#include "TMVA/Reader.h"

#include "fastDQM_CeF3_BTF.h"
#include "interface/RunHelper.h"
#include "interface/PositionTools.h"

struct FitResults{

   float ped_mu;
   float ped_mu_err;
   float ped_sigma;
   float ped_sigma_err;

   float mu;
   float mu_err;
   float offset;
   float offset_err;
   float Q1;
   float Q1_err;
   float sigma;
   float sigma_err;

};


float sumVector( std::vector<float> v );
Double_t cef3Function(Double_t *x, Double_t *par);
FitResults fitSingleHisto( TH1D* histo, double pedMin, double pedMax, double xMin, double xMax );

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

  TH1D* h1_cef3_corr_tot = new TH1D("cef3_corr_tot", "", 1500, 0., 15000);

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {
    tree->GetEntry(iEntry);
    
    if( iEntry % 10000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;
    if( cef3_ok ) {
      
      
      std::vector<float> v_cef3_corr;
      for(unsigned i=0; i<CEF3_CHANNELS; ++i) v_cef3_corr.push_back(cef3_corr[i]);
      
      float eTot_corr = sumVector(v_cef3_corr);

      if( cef3_corr_ok ) {
	h1_cef3_corr_tot->Fill(eTot_corr);
      }
      outTree->Fill();
    }

  }

  FitResults fr_0 = fitSingleHisto( h1_cef3_corr_tot, 0.,0., 1000., 11500. );

  outfile->cd();
  h1_cef3_corr_tot->Write();
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

Double_t cef3Function(Double_t *x, Double_t *par){

   float N = par[0];
   float mu = par[1];
   float Q1 = par[2];
   float sigma = par[3];
   float offset = par[4];
   float sigmaoffset = par[5];
   float alpha = par[6];
   float w = par[7];

   float xx = x[0];
   double value = 0.;
   int NORDERS=6;
   for( unsigned i=1; i<NORDERS; ++i ) {
     double sigma_n = sqrt( (double)(i)*sigma*sigma + sigmaoffset*sigmaoffset);
     value = value + N*(TMath::Poisson( i, mu ) * TMath::Gaus( xx, (double)i*Q1 + offset, sigma_n) );
   }

   return value;


}

FitResults fitSingleHisto( TH1D* histo, double pedMin, double pedMax, double xMin, double xMax ) {

  float integral = histo->Integral();
  TF1* f1_ped = new TF1( "ped", "gaus", pedMin, pedMax );
  f1_ped->SetParameter(0, integral);
  f1_ped->SetParameter(1, 110.);
  f1_ped->SetParameter(2, 10.);

  f1_ped->SetLineColor(kRed+2);
  if(pedMin>10){
  histo->Fit( f1_ped, "RQN" );

  int nSteps = 2;
  for( unsigned iStep=0; iStep<nSteps; iStep++ ) {

    float ped_mean = f1_ped->GetParameter(1);
    float ped_sigma = f1_ped->GetParameter(2);

    float nSigma = 2.;
    float newMin = ped_mean-nSigma*ped_sigma;
    float newMax = ped_mean+nSigma*ped_sigma;

    f1_ped->SetRange( newMin, newMax );
 
    std::string option = (iStep<(nSteps-1)) ? "RQN" : "RQ+";
    histo->Fit( f1_ped, option.c_str() );

  }
  }

  TF1* f1 = new TF1( "func", cef3Function, xMin, xMax, 8 );
  f1->SetParameter( 0, integral ); //normalization
  f1->SetParameter( 1, 1. ); //poiss mu
  f1->SetParameter( 2, 3000. ); //gauss step
  f1->SetParameter( 3, 500. ); //gauss sigma
  f1->SetParameter( 4, 1000 ); //offset
  f1->SetParameter( 5, 30. ); //sigmaoffset
  f1->SetParameter( 6, 0.3 ); //alpha
  f1->SetParameter( 7, 4 ); //w

  f1->FixParameter( 5, 0. ); //sigmaoffset
  f1->FixParameter( 6, 0. ); //alpha
  f1->FixParameter( 7, 0. ); //w

//  f1->SetParLimits( 1, 0.5, 2.5 ); //poiss mu
//  f1->SetParLimits( 2, 10., 40. ); //gauss step
//  f1->SetParLimits( 3, 3., 12. ); //gauss sigma
  //f1->SetParLimits( 4, 90., 110.); //offset
  //f1->SetParLimits( 5, 0., 8. ); //gauss sigma
  f1->SetLineColor(kRed);
  f1->SetLineWidth(3);

  if(pedMin<10){
    f1->FixParameter( 4, 0 );
  }

  histo->Fit( f1, "R+" );

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();

  c1->SetLogy();
  TLine *lineLow1 = new TLine(f1->GetParameter(2)-2*f1->GetParameter(3),0,f1->GetParameter(2)-2*f1->GetParameter(3),1.2*histo->GetMaximum());
  TLine *lineLow2 = new TLine(2*f1->GetParameter(2)-2*sqrt(2)*f1->GetParameter(3),0,2*f1->GetParameter(2)-2*sqrt(2)*f1->GetParameter(3),1.2*histo->GetMaximum());
  TLine *lineUp2 = new TLine(2*f1->GetParameter(2)+2*sqrt(2)*f1->GetParameter(3),0,2*f1->GetParameter(2)+2*sqrt(2)*f1->GetParameter(3),1.2*histo->GetMaximum());
  lineLow1->SetLineWidth(2);
  lineLow2->SetLineWidth(2);
  lineUp2->SetLineWidth(2);
  lineLow1->SetLineColor(kRed);
  lineLow2->SetLineColor(kBlue);
  lineUp2->SetLineColor(kBlue);

  std::cout<<f1->GetParameter(2)-2*f1->GetParameter(3)<<" "<<2*f1->GetParameter(2)-2*sqrt(2)*f1->GetParameter(3)<<std::endl;
  histo->SetXTitle( "ADC Counts" );
  histo->Draw();
  TString histoName(histo->GetName());
  c1->SaveAs( histoName + ".eps" );
  c1->SaveAs( histoName + ".png" );
  lineLow1->Draw("same");
  lineLow2->Draw("same");
  lineUp2->Draw("same");





  c1->SaveAs( histoName + "_sigma.eps" );
  c1->SaveAs( histoName + "_sigma.png" );

  FitResults fr;
  fr.ped_mu = f1_ped->GetParameter(1);
  fr.ped_mu_err = f1_ped->GetParError(1);
  fr.ped_sigma = f1_ped->GetParameter(2);
  fr.ped_sigma_err = f1_ped->GetParError(2);

  fr.mu = f1->GetParameter(1);
  fr.mu_err = f1->GetParError(1);
  fr.offset = f1->GetParameter(4);
  fr.offset_err = f1->GetParError(4);
  fr.Q1 = f1->GetParameter(2);
  fr.Q1_err = f1->GetParError(2);
  fr.sigma = f1->GetParameter(3);
  fr.sigma_err = f1->GetParError(3);

  delete c1;
  delete f1;
  delete f1_ped;

  return fr;



}
