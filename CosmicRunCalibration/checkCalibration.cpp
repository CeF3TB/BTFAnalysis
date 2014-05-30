#include <TStyle.h>
#include <TCanvas.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "../PositionAnalysis/fastDQM_CeF3_BTF.h"

void setStyle(){
  // set the TStyle
  TStyle* style_ = new TStyle("DrawBaseStyle", "");
  style_->SetCanvasColor(0);
  style_->SetPadColor(0);
  style_->SetFrameFillColor(0);
  style_->SetStatColor(0);
  style_->SetOptStat(0);
  style_->SetTitleFillColor(0);
  style_->SetCanvasBorderMode(0);
  style_->SetPadBorderMode(0);
  style_->SetFrameBorderMode(0);
  style_->SetPadBottomMargin(0.12);
  style_->SetPadLeftMargin(0.12);
  style_->cd();

  // For the canvas:
  style_->SetCanvasBorderMode(0);
  style_->SetCanvasColor(kWhite);
  style_->SetCanvasDefH(600); //Height of canvas
  style_->SetCanvasDefW(600); //Width of canvas
  style_->SetCanvasDefX(0); //POsition on screen
  style_->SetCanvasDefY(0);

  // For the Pad:
  style_->SetPadBorderMode(0);
  style_->SetPadColor(kWhite);
  style_->SetPadGridX(false);
  style_->SetPadGridY(false);
  style_->SetGridColor(0);
  style_->SetGridStyle(3);
  style_->SetGridWidth(1);

  // For the frame:
  style_->SetFrameBorderMode(0);
  style_->SetFrameBorderSize(1);
  style_->SetFrameFillColor(0);
  style_->SetFrameFillStyle(0);
  style_->SetFrameLineColor(1);
  style_->SetFrameLineStyle(1);
  style_->SetFrameLineWidth(1);


  // Margins:
  style_->SetPadTopMargin(0.05);
  style_->SetPadBottomMargin(0.15);//0.13);
  style_->SetPadLeftMargin(0.15);//0.16);
  style_->SetPadRightMargin(0.05);//0.02);

  // For the Global title:

  style_->SetOptTitle(0);
  style_->SetTitleFont(42);
  style_->SetTitleColor(1);
  style_->SetTitleTextColor(1);
  style_->SetTitleFillColor(10);
  style_->SetTitleFontSize(0.05);

  // For the axis titles:

  style_->SetTitleColor(1, "XYZ");
  style_->SetTitleFont(42, "XYZ");
  style_->SetTitleSize(0.05, "XYZ");
  style_->SetTitleXOffset(1.15);//0.9);
  style_->SetTitleYOffset(1.4); // => 1.15 if exponents

  // For the axis labels:

  style_->SetLabelColor(1, "XYZ");
  style_->SetLabelFont(42, "XYZ");
  style_->SetLabelOffset(0.007, "XYZ");
  style_->SetLabelSize(0.045, "XYZ");

  // For the axis:

  style_->SetAxisColor(1, "XYZ");
  style_->SetStripDecimals(kTRUE);
  style_->SetTickLength(0.03, "XYZ");
  style_->SetNdivisions(510, "XYZ");
  style_->SetPadTickX(1); // To get tick marks on the opposite side of the frame
  style_->SetPadTickY(1);

  style_->cd();
}




int main( int argc, char* argv[] ) {

  setStyle();
  
  std::string runName = "test_10";
  if( argc>1 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
  }
  
  
  std::string fileName = "../PositionAnalysis/PosAn_" + runName + ".root";
  TFile* file = TFile::Open( fileName.c_str() );
  std::cout << "-> Opened file: " << fileName << std::endl;

  TTree* tree = (TTree*)file->Get("tree_passedEvents");

  float hodox_corr[HODOX_CHANNELS], hodoy_corr[HODOY_CHANNELS],cef3[CEF3_CHANNELS],cef3_corr[CEF3_CHANNELS],cef3_pedSubtracted[CEF3_CHANNELS],scintFront;
  int nHodoClustersX,nHodoClustersY;

  tree->SetBranchAddress( "nHodoClustersX", &nHodoClustersX );
  tree->SetBranchAddress( "nHodoClustersY", &nHodoClustersY );
  tree->SetBranchAddress( "hodox_corr", hodox_corr );
  tree->SetBranchAddress( "hodoy_corr", hodoy_corr );
  tree->SetBranchAddress( "cef3", cef3 );
  tree->SetBranchAddress( "cef3_corr", cef3_corr );
  tree->SetBranchAddress( "cef3_pedSubtracted", cef3_pedSubtracted );
  tree->SetBranchAddress( "scintFront", &scintFront );

  int nentries = tree->GetEntries();

  std::string outfileName = "PhotoElectronAn_" + runName + ".root";
  TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );

  TH1D* h1_singleEleEnergyDistribution = new TH1D("singleEleEnergyDistribution", "",1000 , 0 , 10000);
  TH1D* h1_singleEleEnergyDistributionCalib = new TH1D("singleEleEnergyDistribution", "",1000 , 0 , 10000);
  TH1D* h1_singleEleEnergyDistribution_ch[4];
  TH1D* h1_singleEleEnergyDistributionCalib_ch[4];

//  Q1 calibration
//  std::vector<float> correctionFactors;
//  correctionFactors.push_back(1.11228);
//  correctionFactors.push_back(0.855333);
//  correctionFactors.push_back(0.9802);
//  correctionFactors.push_back(1.08781);

//  mu*Q1 calibration
//  std::vector<float> correctionFactors;
//  correctionFactors.push_back(1.19514);
//  correctionFactors.push_back(0.860177);
//  correctionFactors.push_back(0.952363);
//  correctionFactors.push_back(0.977169);

  //mumean*Q1 calibration
  std::vector<float> correctionFactors;
  correctionFactors.push_back(1.09716);
  correctionFactors.push_back(0.843707);
  correctionFactors.push_back(0.966413);
  correctionFactors.push_back(1.07241);



  for (int i=0;i<4;++i){
    TString nameHisto="singleEnergyDistribution";
    TString nameHistoCalib="singleEnergyDistributionCalib";
    nameHisto+=i;
    nameHistoCalib+=i;
    h1_singleEleEnergyDistribution_ch[i] = new TH1D(nameHisto, "",100 , 0 , 4000.);
    h1_singleEleEnergyDistributionCalib_ch[i] = new TH1D(nameHistoCalib, "",100 , 0 , 4000);
  }



  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {
    
    tree->GetEntry(iEntry);

    if( iEntry % 1000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;

    bool isSingleEle=scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1;
    float cef3Total=0.;    
    float cef3TotalCalib=0.;   
 
    for (int i=0;i<CEF3_CHANNELS;i++){

      float energy=cef3_pedSubtracted[i];

      cef3Total+=energy;
      cef3TotalCalib+=energy*correctionFactors[i];

      if(isSingleEle){
	h1_singleEleEnergyDistribution_ch[i]->Fill(energy);
	h1_singleEleEnergyDistributionCalib_ch[i]->Fill(energy*correctionFactors[i]);
      }
    }

    if(isSingleEle) {
      h1_singleEleEnergyDistribution->Fill(cef3Total);
      h1_singleEleEnergyDistributionCalib->Fill(cef3TotalCalib);
    }
  }

  TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();



  TH2D* h2_axes = new TH2D("axes", "", 10, 100., 2000., 10, 0., 1.15*h1_singleEleEnergyDistributionCalib_ch[2]->GetMaximum() );
  h2_axes->SetXTitle( "ADC Counts" );
  h2_axes->Draw();

  h1_singleEleEnergyDistributionCalib_ch[0]->SetLineColor(kBlack);
  h1_singleEleEnergyDistributionCalib_ch[1]->SetLineColor(kRed);
  h1_singleEleEnergyDistributionCalib_ch[2]->SetLineColor(kBlue);
  h1_singleEleEnergyDistributionCalib_ch[3]->SetLineColor(kMagenta);
     					  
  h1_singleEleEnergyDistributionCalib_ch[0]->SetLineWidth(2);
  h1_singleEleEnergyDistributionCalib_ch[1]->SetLineWidth(2);
  h1_singleEleEnergyDistributionCalib_ch[2]->SetLineWidth(2);
  h1_singleEleEnergyDistributionCalib_ch[3]->SetLineWidth(2);

  TLegend* legend = new TLegend( 0.7, 0.7, 0.90, 0.9 );
  legend->SetLineColor(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( h1_singleEleEnergyDistributionCalib_ch[0], "channel 0", "L" );
  legend->AddEntry( h1_singleEleEnergyDistributionCalib_ch[1], "channel 1", "L" );
  legend->AddEntry( h1_singleEleEnergyDistributionCalib_ch[2], "channel 2", "L" );
  legend->AddEntry( h1_singleEleEnergyDistributionCalib_ch[3], "channel 3", "L" );
  legend->Draw("same");

  h1_singleEleEnergyDistributionCalib_ch[2]->Draw("same");
  h1_singleEleEnergyDistributionCalib_ch[0]->Draw("same");
  h1_singleEleEnergyDistributionCalib_ch[1]->Draw("same");
  h1_singleEleEnergyDistributionCalib_ch[3]->Draw("same");

  c1->SaveAs("calibratedChannels.png");
  c1->SaveAs("calibratedChannels.eps");

  c1->Clear();

  h2_axes->Draw();

  h1_singleEleEnergyDistribution_ch[0]->SetLineColor(kBlack);
  h1_singleEleEnergyDistribution_ch[1]->SetLineColor(kRed);
  h1_singleEleEnergyDistribution_ch[2]->SetLineColor(kBlue);
  h1_singleEleEnergyDistribution_ch[3]->SetLineColor(kMagenta);
     					  
  h1_singleEleEnergyDistribution_ch[0]->SetLineWidth(2);
  h1_singleEleEnergyDistribution_ch[1]->SetLineWidth(2);
  h1_singleEleEnergyDistribution_ch[2]->SetLineWidth(2);
  h1_singleEleEnergyDistribution_ch[3]->SetLineWidth(2);

  h1_singleEleEnergyDistribution_ch[0]->Draw("same");
  h1_singleEleEnergyDistribution_ch[1]->Draw("same");
  h1_singleEleEnergyDistribution_ch[2]->Draw("same");
  h1_singleEleEnergyDistribution_ch[3]->Draw("same");

  legend->Draw("same");

  c1->SaveAs("uncalibratedChannels.png");
  c1->SaveAs("uncalibratedChannels.eps");


  //plot mu,rms of uncalib distribution
  TH1D* h1_mean = new TH1D("mean", "", 4, -0.5, 3.5);
  TH1D* h1_rms = new TH1D("rms", "", 4, -0.5, 3.5);

  TF1* f1_gaus = new TF1( "f1_gaus", "gaus",0,2000 );

  f1_gaus->SetRange(h1_singleEleEnergyDistribution_ch[0]->GetMean()-3*h1_singleEleEnergyDistribution_ch[0]->GetRMS(),h1_singleEleEnergyDistribution_ch[0]->GetMean()+3*h1_singleEleEnergyDistribution_ch[0]->GetRMS() );
  h1_singleEleEnergyDistribution_ch[0]->Fit(f1_gaus,"R");
  h1_mean->SetBinContent( 1, f1_gaus->GetParameter(1) );
  h1_mean->SetBinError( 1, f1_gaus->GetParError(1) );
  h1_rms->SetBinContent( 1, f1_gaus->GetParameter(2) );
  h1_rms->SetBinError( 1, f1_gaus->GetParError(2) );


  f1_gaus->SetRange(h1_singleEleEnergyDistribution_ch[1]->GetMean()-3*h1_singleEleEnergyDistribution_ch[1]->GetRMS(),h1_singleEleEnergyDistribution_ch[1]->GetMean()+3*h1_singleEleEnergyDistribution_ch[1]->GetRMS() );
  h1_singleEleEnergyDistribution_ch[1]->Fit(f1_gaus,"R");
  h1_mean->SetBinContent( 2, f1_gaus->GetParameter(1) );
  h1_mean->SetBinError( 2, f1_gaus->GetParError(1) );
  h1_rms->SetBinContent( 2, f1_gaus->GetParameter(2) );
  h1_rms->SetBinError( 2, f1_gaus->GetParError(2) );


  f1_gaus->SetRange(h1_singleEleEnergyDistribution_ch[2]->GetMean()-3*h1_singleEleEnergyDistribution_ch[2]->GetRMS(),h1_singleEleEnergyDistribution_ch[2]->GetMean()+3*h1_singleEleEnergyDistribution_ch[2]->GetRMS() );
  h1_singleEleEnergyDistribution_ch[2]->Fit(f1_gaus,"R");
  h1_mean->SetBinContent( 3, f1_gaus->GetParameter(1) );
  h1_mean->SetBinError( 3, f1_gaus->GetParError(2) );
  h1_rms->SetBinContent( 3, f1_gaus->GetParameter(2) );
  h1_rms->SetBinError( 3, f1_gaus->GetParError(2) );

  f1_gaus->SetRange(h1_singleEleEnergyDistribution_ch[3]->GetMean()-3*h1_singleEleEnergyDistribution_ch[3]->GetRMS(),h1_singleEleEnergyDistribution_ch[3]->GetMean()+3*h1_singleEleEnergyDistribution_ch[3]->GetRMS() );
  h1_singleEleEnergyDistribution_ch[3]->Fit(f1_gaus,"R");
  h1_mean->SetBinContent( 4, f1_gaus->GetParameter(1) );
  h1_mean->SetBinError( 4, f1_gaus->GetParError(1) );
  h1_rms->SetBinContent( 4, f1_gaus->GetParameter(2) );
  h1_rms->SetBinError( 4, f1_gaus->GetParError(2) );


  c1->Clear();
  h1_mean->SetXTitle( "Channel Number");
  h1_mean->SetMarkerStyle(20);
  h1_mean->SetMarkerSize(1.6);

  h1_mean->GetYaxis()->SetRangeUser(720.,890.);
  h1_mean->Draw();

  c1->SaveAs("muUncalibratedChannels.png");
  c1->SaveAs("muUncalibratedChannels.eps");

  c1->Clear();
  h1_rms->SetXTitle( "Channel Number");
  h1_rms->SetMarkerStyle(20);
  h1_rms->SetMarkerSize(1.6);


  h1_rms->Draw();

  c1->SaveAs("rmsUncalibratedChannels.png");
  c1->SaveAs("rmsUncalibratedChannels.eps");


  //plot mu,rms of uncalib distribution
  TH1D* h1_meanCalib = new TH1D("meanCalib", "", 4, -0.5, 3.5);
  TH1D* h1_rmsCalib = new TH1D("rmsCalib", "", 4, -0.5, 3.5);


  f1_gaus->SetRange(h1_singleEleEnergyDistributionCalib_ch[0]->GetMean()-3*h1_singleEleEnergyDistributionCalib_ch[0]->GetRMS(),h1_singleEleEnergyDistributionCalib_ch[0]->GetMean()+3*h1_singleEleEnergyDistributionCalib_ch[0]->GetRMS() );
  h1_singleEleEnergyDistributionCalib_ch[0]->Fit(f1_gaus,"R");
  h1_meanCalib->SetBinContent( 1, f1_gaus->GetParameter(1) );
  h1_meanCalib->SetBinError( 1, f1_gaus->GetParError(1) );
  h1_rmsCalib->SetBinContent( 1, f1_gaus->GetParameter(2) );
  h1_rmsCalib->SetBinError( 1, f1_gaus->GetParError(2) );


  f1_gaus->SetRange(h1_singleEleEnergyDistributionCalib_ch[1]->GetMean()-3*h1_singleEleEnergyDistributionCalib_ch[1]->GetRMS(),h1_singleEleEnergyDistributionCalib_ch[1]->GetMean()+3*h1_singleEleEnergyDistributionCalib_ch[1]->GetRMS() );
  h1_singleEleEnergyDistributionCalib_ch[1]->Fit(f1_gaus,"R");
  h1_meanCalib->SetBinContent( 2, f1_gaus->GetParameter(1) );
  h1_meanCalib->SetBinError( 2, f1_gaus->GetParError(1) );
  h1_rmsCalib->SetBinContent( 2, f1_gaus->GetParameter(2) );
  h1_rmsCalib->SetBinError( 2, f1_gaus->GetParError(2) );


  f1_gaus->SetRange(h1_singleEleEnergyDistributionCalib_ch[2]->GetMean()-3*h1_singleEleEnergyDistributionCalib_ch[2]->GetRMS(),h1_singleEleEnergyDistributionCalib_ch[2]->GetMean()+3*h1_singleEleEnergyDistributionCalib_ch[2]->GetRMS() );
  h1_singleEleEnergyDistributionCalib_ch[2]->Fit(f1_gaus,"R");
  h1_meanCalib->SetBinContent( 3, f1_gaus->GetParameter(1) );
  h1_meanCalib->SetBinError( 3, f1_gaus->GetParError(2) );
  h1_rmsCalib->SetBinContent( 3, f1_gaus->GetParameter(2) );
  h1_rmsCalib->SetBinError( 3, f1_gaus->GetParError(2) );

  f1_gaus->SetRange(h1_singleEleEnergyDistributionCalib_ch[3]->GetMean()-3*h1_singleEleEnergyDistributionCalib_ch[3]->GetRMS(),h1_singleEleEnergyDistributionCalib_ch[3]->GetMean()+3*h1_singleEleEnergyDistributionCalib_ch[3]->GetRMS() );
  h1_singleEleEnergyDistributionCalib_ch[3]->Fit(f1_gaus,"R");
  h1_meanCalib->SetBinContent( 4, f1_gaus->GetParameter(1) );
  h1_meanCalib->SetBinError( 4, f1_gaus->GetParError(1) );
  h1_rmsCalib->SetBinContent( 4, f1_gaus->GetParameter(2) );
  h1_rmsCalib->SetBinError( 4, f1_gaus->GetParError(2) );


  c1->Clear();
  h1_meanCalib->GetYaxis()->SetRangeUser(720.,890.);
  h1_meanCalib->SetXTitle( "Channel Number");
  h1_meanCalib->SetMarkerStyle(20);
  h1_meanCalib->SetMarkerSize(1.6);


  h1_meanCalib->Draw();

  c1->SaveAs("muCalibratedChannels.png");
  c1->SaveAs("muCalibratedChannels.eps");

  c1->Clear();
  h1_rmsCalib->SetXTitle( "Channel Number");
  h1_rmsCalib->SetMarkerStyle(20);
  h1_rmsCalib->SetMarkerSize(1.6);


  h1_rmsCalib->Draw();

  c1->SaveAs("rmsCalibratedChannels.png");
  c1->SaveAs("rmsCalibratedChannels.eps");

  float mean=0.,rms=0.;
  float meanCalib=0.,rmscalib=0.;
  for (int i=1;i<=4;i++){
    mean+=h1_mean->GetBinContent(i);
    meanCalib+=h1_meanCalib->GetBinContent(i);
  }
  mean=mean/4.;
  meanCalib=meanCalib/4.;

  for (int i=1;i<=4;i++){
    rms+= (h1_mean->GetBinContent(i)-mean)*(h1_mean->GetBinContent(i)-mean);
    rmscalib+= (h1_meanCalib->GetBinContent(i)-meanCalib)*(h1_meanCalib->GetBinContent(i)-meanCalib);
  }
  rms=sqrt(rms/4);
  rmscalib=sqrt(rmscalib/4);
  std::cout<<"rms: "<<rms<<" rms_calib: "<<rmscalib<<std::endl;

  outfile->Write();
  outfile->Close();

  return 0;

    
}

