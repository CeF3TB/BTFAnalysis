#include <iostream>
#include <string>
#include <stdlib.h>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLegend.h"

#include "TMath.h"
#include "/Users/micheli/work/BTFAnalysis_20140515/PositionAnalysis/fastDQM_CeF3_BTF.h"
#include "TTree.h"

struct FitResults {

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


int NORDERS;


Double_t PMTFunction(Double_t *x, Double_t *par)
{

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

   for( unsigned i=1; i<NORDERS; ++i ) {

     //double Qn = offset + (double)(i)*Q1;
     double sigma_n = sqrt( (double)(i)*sigma*sigma + sigmaoffset*sigmaoffset);

     //double poisson = TMath::Poisson( i, mu );
     //double gauss   = TMath::Gaus( xx, Qn, sigma_n );

     //double xxp = xx     - Qn - alpha*sigma_n*sigma_n;
     //double Q0p = offset - Qn - alpha*sigma_n*sigma_n;
     //double bg      = 0.5*alpha * TMath::Exp(-alpha*xxp)* (
     //                     TMath::Erf( abs(Q0p)/(sigma_n*sqrt(2) ) ) +  xxp/abs(xxp) * TMath::Erf( abs(xxp)/(sigma_n*sqrt(2)) ) );
     //value = value + N*( poisson * ( (1.-w)*gauss + w*bg ) );
     value = value + N*(TMath::Poisson( i, mu ) * TMath::Gaus( xx, (double)i*Q1 + offset, sigma_n) );
     //value = value + N*(TMath::Poisson( i, mu ) * TMath::Gaus( xx, (double)i*Q1 + offset, sqrt((double)i)*sigma ));
   }

   return value;

}

std::vector<float> intercalibrateFibers(FitResults fr_0,FitResults fr_1,FitResults fr_2,FitResults fr_3){

  std::vector<float> correctionFactors;

  float mean=0.;

  mean+=(fr_0.offset+fr_0.Q1)/(fr_0.Q1_err*fr_0.Q1_err+fr_0.offset_err*fr_0.offset_err);
  mean+=(fr_1.offset+fr_1.Q1)/(fr_1.Q1_err*fr_1.Q1_err+fr_1.offset_err*fr_1.offset_err);
  mean+=(fr_2.offset+fr_2.Q1)/(fr_2.Q1_err*fr_2.Q1_err+fr_2.offset_err*fr_2.offset_err);
  mean+=(fr_3.offset+fr_3.Q1)/(fr_3.Q1_err*fr_3.Q1_err+fr_3.offset_err*fr_3.offset_err);

  mean=mean/(1/(fr_0.Q1_err*fr_0.Q1_err+fr_0.offset_err*fr_0.offset_err)+1/(fr_1.Q1_err*fr_1.Q1_err+fr_1.offset_err*fr_1.offset_err)+1/(fr_2.Q1_err*fr_2.Q1_err+fr_2.offset_err*fr_2.offset_err)+1/(fr_3.Q1_err*fr_3.Q1_err+fr_3.offset_err*fr_3.offset_err));

  std::cout<<"---mean is "<<mean<<std::endl;

  correctionFactors.push_back(mean/(fr_0.Q1+fr_0.offset));
  correctionFactors.push_back(mean/(fr_1.Q1+fr_1.offset));
  correctionFactors.push_back(mean/(fr_2.Q1+fr_2.offset));
  correctionFactors.push_back(mean/(fr_3.Q1+fr_3.offset));

  return correctionFactors;

}







FitResults fitSingleHisto( TH1D* histo, double pedMin, double pedMax, double xMin, double xMax ) {

  float integral = histo->Integral();
  TF1* f1_ped = new TF1( "ped", "gaus", pedMin, pedMax );
  f1_ped->SetParameter(0, integral);
  f1_ped->SetParameter(1, 110.);
  f1_ped->SetParameter(2, 10.);

  f1_ped->SetLineColor(kRed+2);

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

  TF1* f1 = new TF1( "func", PMTFunction, xMin, xMax, 8 );
  f1->SetParameter( 0, integral ); //normalization
  f1->SetParameter( 1, 1. ); //poiss mu
  f1->SetParameter( 2, 25. ); //gauss step
  f1->SetParameter( 3, 10. ); //gauss sigma
  f1->SetParameter( 4, 100 ); //offset
  f1->SetParameter( 5, 3. ); //sigmaoffset
  f1->SetParameter( 6, 0.03 ); //alpha
  f1->SetParameter( 7, 0.4 ); //w

  f1->FixParameter( 5, 0. ); //sigmaoffset
  f1->FixParameter( 6, 0. ); //alpha
  f1->FixParameter( 7, 0. ); //w

  f1->SetParLimits( 1, 0.5, 2.5 ); //poiss mu
  f1->SetParLimits( 2, 10., 40. ); //gauss step
  f1->SetParLimits( 3, 3., 12. ); //gauss sigma
  //f1->SetParLimits( 4, 90., 110.); //offset
  //f1->SetParLimits( 5, 0., 8. ); //gauss sigma
  f1->SetLineColor(kRed+2);

  histo->Fit( f1, "R+" );

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();

  c1->SetLogy();

  TH2D* h2_axes = new TH2D("axes", "", 10, 100., 350., 10, 9., 7.*histo->GetMaximum() );
  h2_axes->SetXTitle( "ADC Counts" );
  h2_axes->Draw();

  histo->Draw("same");


  TString histoName(histo->GetName());



  c1->SaveAs( histoName + ".eps" );
  c1->SaveAs( histoName + ".png" );

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
  delete h2_axes;

  return fr;



}


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



void doSummaryPlots(FitResults fr_0,FitResults fr_1,FitResults fr_2,FitResults fr_3,FitResults fr_corr_0,FitResults fr_corr_1,FitResults fr_corr_2,FitResults fr_corr_3){

  TH1D* h1_Q1 = new TH1D("Q1", "", 4, -0.5, 3.5);
  h1_Q1->SetXTitle( "Channel Number");
  h1_Q1->SetBinContent( 1, fr_0.Q1 );
  h1_Q1->SetBinContent( 2, fr_1.Q1 );
  h1_Q1->SetBinContent( 3, fr_2.Q1 );
  h1_Q1->SetBinContent( 4, fr_3.Q1 );

  h1_Q1->SetBinError( 1, fr_0.Q1_err );
  h1_Q1->SetBinError( 2, fr_1.Q1_err );
  h1_Q1->SetBinError( 3, fr_2.Q1_err );
  h1_Q1->SetBinError( 4, fr_3.Q1_err );



  h1_Q1->SetMarkerStyle(20);
  h1_Q1->SetMarkerSize(1.6);




  TH1D* h1_offset = new TH1D("offset", "", 4, -0.5, 3.5);
  h1_offset->SetXTitle( "Channel Number");
  h1_offset->SetBinContent( 1, fr_0.offset );
  h1_offset->SetBinContent( 2, fr_1.offset );
  h1_offset->SetBinContent( 3, fr_2.offset );
  h1_offset->SetBinContent( 4, fr_3.offset );

  h1_offset->SetBinError( 1, fr_0.offset_err );
  h1_offset->SetBinError( 2, fr_1.offset_err );
  h1_offset->SetBinError( 3, fr_2.offset_err );
  h1_offset->SetBinError( 4, fr_3.offset_err );

  h1_offset->SetMarkerStyle(21);
  h1_offset->SetMarkerSize(1.6);
  h1_offset->SetMarkerColor(kRed);

  TH1D* h1_ped_mu = new TH1D("ped_mu", "", 4, -0.5, 3.5);
  h1_ped_mu->SetXTitle( "Channel Number");
  h1_ped_mu->SetBinContent( 1, fr_0.ped_mu );
  h1_ped_mu->SetBinContent( 2, fr_1.ped_mu );
  h1_ped_mu->SetBinContent( 3, fr_2.ped_mu );
  h1_ped_mu->SetBinContent( 4, fr_3.ped_mu );

  h1_ped_mu->SetBinError( 1, fr_0.ped_mu_err );
  h1_ped_mu->SetBinError( 2, fr_1.ped_mu_err );
  h1_ped_mu->SetBinError( 3, fr_2.ped_mu_err );
  h1_ped_mu->SetBinError( 4, fr_3.ped_mu_err );

  h1_ped_mu->SetMarkerStyle(20);
  h1_ped_mu->SetMarkerSize(1.6);
  h1_ped_mu->SetMarkerColor(38);





  TCanvas* c2 = new TCanvas( "c2", "", 600, 600 );

  c2->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, -0.5, 3.5, 10, 0., 200.);
  h2_axes->SetXTitle( "Channel Number" );
  h2_axes->SetYTitle( "ADC Counts" );
  h2_axes->Draw();

  TLegend* legend = new TLegend( 0.4, 0.7, 0.68, 0.9 );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( h1_Q1, "Q1", "P" );
  legend->AddEntry( h1_offset, "Offset", "P" );
  legend->AddEntry( h1_ped_mu, "Pedestal", "P" );
  legend->Draw("same");

  TPaveText* label_top = new TPaveText(0.4,0.953,0.975,0.975, "brNDC");
  label_top->SetFillColor(kWhite);
  label_top->SetTextSize(0.038);
  label_top->SetTextAlign(31); // align right
  label_top->SetTextFont(62);
  label_top->AddText("Cosmic Run");
  label_top->Draw("same");

  gStyle->SetErrorX(0);

  h1_Q1->Draw("P same");
  h1_offset->Draw("P same");
  h1_ped_mu->Draw("P same");

  c2->SaveAs("summaryPlot.eps");
  c2->SaveAs("summaryPlot.png");


  c2->Clear();



  TH1D* h1_mu = new TH1D("mu", "", 4, -0.5, 3.5);
  h1_mu->SetXTitle( "Channel Number");
  h1_mu->SetBinContent( 1, fr_0.mu );
  h1_mu->SetBinContent( 2, fr_1.mu );
  h1_mu->SetBinContent( 3, fr_2.mu );
  h1_mu->SetBinContent( 4, fr_3.mu );

  h1_mu->SetBinError( 1, fr_0.mu_err );
  h1_mu->SetBinError( 2, fr_1.mu_err );
  h1_mu->SetBinError( 3, fr_2.mu_err );
  h1_mu->SetBinError( 4, fr_3.mu_err );

  h1_mu->SetMarkerStyle(20);
  h1_mu->SetMarkerSize(1.6);
  h1_mu->SetMarkerColor(46);

  TH2D* h2_axes_2 = new TH2D("axes_2", "", 10, -0.5, 3.5, 10, 0., 3.);
  h2_axes_2->SetXTitle( "Channel Number" );
  h2_axes_2->SetYTitle( "Poisson #mu" );
  h2_axes_2->Draw();


  label_top->Draw("same");
  h1_mu->Draw("same");

  c2->SaveAs("poissonMu.eps");
  c2->SaveAs("poissonMu.png");

  c2->Clear();

  //first peak position
  TH1D* h1_1stPeak = new TH1D("1stPeak", "", 4, -0.5, 3.5);
  h1_1stPeak->SetXTitle( "Channel Number");
  h1_1stPeak->SetBinContent( 1, fr_0.Q1 + fr_0.ped_mu );
  h1_1stPeak->SetBinContent( 2, fr_1.Q1 + fr_1.ped_mu);
  h1_1stPeak->SetBinContent( 3, fr_2.Q1 + fr_2.ped_mu);
  h1_1stPeak->SetBinContent( 4, fr_3.Q1 + fr_3.ped_mu);
  h1_1stPeak->SetBinError( 1, sqrt(fr_0.Q1_err*fr_0.Q1_err + fr_0.ped_mu_err*fr_0.ped_mu_err));
  h1_1stPeak->SetBinError( 2, sqrt(fr_1.Q1_err*fr_1.Q1_err + fr_1.ped_mu_err*fr_1.ped_mu_err));
  h1_1stPeak->SetBinError( 3, sqrt(fr_2.Q1_err*fr_2.Q1_err + fr_2.ped_mu_err*fr_2.ped_mu_err));
  h1_1stPeak->SetBinError( 4, sqrt(fr_3.Q1_err*fr_3.Q1_err + fr_3.ped_mu_err*fr_3.ped_mu_err));
  h1_1stPeak->SetMarkerStyle(20);
  h1_1stPeak->SetMarkerSize(1.6);
  h1_1stPeak->SetMarkerColor(kBlue);

  TH2D* h2_axes_3 = new TH2D("axes_3", "", 10, -0.5, 3.5, 10, 130., 160.);
  h2_axes_3->SetXTitle( "Channel Number" );
  h2_axes_3->SetYTitle( "ADC Counts" );


  TPaveText* label_top_2 = new TPaveText(0.4,0.953,0.975,0.975, "brNDC");
  label_top_2->SetFillColor(kWhite);
  label_top_2->SetTextSize(0.038);
  label_top_2->SetTextAlign(31); // align right
  label_top_2->SetTextFont(62);
  label_top_2->AddText("1st peak position");


  //first peak position
  TH1D* h1_1stPeak_corr = new TH1D("1stPeak_corr", "", 4, -0.5, 3.5);
  h1_1stPeak_corr->SetXTitle( "Channel Number");

  std::cout<<"-----1st peak "<<fr_corr_0.Q1+fr_corr_0.ped_mu<<" "<<fr_corr_1.Q1+fr_corr_1.ped_mu<<" "<<fr_corr_2.Q1+fr_corr_2.ped_mu<<" "<<fr_corr_3.Q1+fr_corr_3.ped_mu<<std::endl;
  h1_1stPeak_corr->SetBinContent( 1, fr_corr_0.Q1 + fr_corr_0.ped_mu );
  h1_1stPeak_corr->SetBinContent( 2, fr_corr_1.Q1 + fr_corr_1.ped_mu);
  h1_1stPeak_corr->SetBinContent( 3, fr_corr_2.Q1 + fr_corr_2.ped_mu);
  h1_1stPeak_corr->SetBinContent( 4, fr_corr_3.Q1 + fr_corr_3.ped_mu);
  h1_1stPeak_corr->SetBinError( 1, sqrt(fr_corr_0.Q1_err*fr_corr_0.Q1_err + fr_corr_0.ped_mu_err*fr_corr_0.ped_mu_err));
  h1_1stPeak_corr->SetBinError( 2, sqrt(fr_corr_1.Q1_err*fr_corr_1.Q1_err + fr_corr_1.ped_mu_err*fr_corr_1.ped_mu_err));
  h1_1stPeak_corr->SetBinError( 3, sqrt(fr_corr_2.Q1_err*fr_corr_2.Q1_err + fr_corr_2.ped_mu_err*fr_corr_2.ped_mu_err));
  h1_1stPeak_corr->SetBinError( 4, sqrt(fr_corr_3.Q1_err*fr_corr_3.Q1_err + fr_corr_3.ped_mu_err*fr_corr_3.ped_mu_err));
  h1_1stPeak_corr->SetMarkerStyle(20);
  h1_1stPeak_corr->SetMarkerSize(1.6);
  h1_1stPeak_corr->SetMarkerColor(kRed);

  h2_axes_3->Draw();
  h1_1stPeak->Draw("same");
  h1_1stPeak_corr->Draw("same");
  label_top_2->Draw("same");

  c2->SaveAs("peak1Position.eps");
  c2->SaveAs("peak1Position.png");





}



int main( int argc, char* argv[] ) {

  setStyle();

  std::string runName = "precalib_BGO_pedestal_noSource";
  if( argc>1 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
  }

  NORDERS = 6;
  if( argc>2 ) {
    NORDERS = atoi(argv[2]);
    std::cout << "-> NORDER is set to: " << NORDERS << std::endl;
  }


  std::string fileName = "run_" + runName + ".root";
  TFile* file = TFile::Open(fileName.c_str());
  
  
  if( argc>1 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
  }
  if( file==0 ) {
    std::cout << "ERROR! Din't find file " << fileName << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  
  TTree* tree = (TTree*)file->Get("eventRawData");
  UInt_t evtNumber;
  tree->SetBranchAddress( "evtNumber", &evtNumber );
  UInt_t adcData[40];
  tree->SetBranchAddress( "adcData", adcData );
  UInt_t adcBoard[40];
  tree->SetBranchAddress( "adcBoard", adcBoard );
  UInt_t adcChannel[40];
  tree->SetBranchAddress( "adcChannel", adcChannel );
  
  int nentries = tree->GetEntries();

  std::string outfileName = "calibAn_" + runName + ".root";
  TFile* outFile = TFile::Open( outfileName.c_str(), "RECREATE" );

  TH1D* h1_cef3_0   = new TH1D("cef3_0",   "", 400, 0., 400.);
  TH1D* h1_cef3_1   = new TH1D("cef3_1",   "", 400, 0., 400.);
  TH1D* h1_cef3_2   = new TH1D("cef3_2",   "", 400, 0., 400.);
  TH1D* h1_cef3_3   = new TH1D("cef3_3",   "", 400, 0., 400.);
  TH1D* h1_cef3_tot = new TH1D("cef3_tot", "", 400, 0., 4.*400.);

  TH1D* h1_cef3_corr_0   = new TH1D("cef3_corr_0",   "", 400, 0., 400.);
  TH1D* h1_cef3_corr_1   = new TH1D("cef3_corr_1",   "", 385, 0., 400.4);
  TH1D* h1_cef3_corr_2   = new TH1D("cef3_corr_2",   "", 400, 0., 404.);
  TH1D* h1_cef3_corr_3   = new TH1D("cef3_corr_3",   "", 400, 0., 400.);
  TH1D* h1_cef3_corr_tot = new TH1D("cef3_corr_tot", "", 400, 0., 4.*400.);

  
  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {
    
    tree->GetEntry(iEntry);
    
    if( iEntry % 5000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;
    
    for( unsigned i=0; i<40; ++i ) {
      
      int board  = adcBoard[i];
      int channel= adcChannel[i];
      
      
      float cef3=0;
      
      
      if( board==CEF3_ADC_BOARD ) {
	if( channel==(CEF3_ADC_START_CHANNEL  ) ){
	  h1_cef3_0->Fill(adcData[i]);
	  cef3+=adcData[i];
	}
        else if( channel==(CEF3_ADC_START_CHANNEL+1) ){
	  h1_cef3_1->Fill(adcData[i]);
	  cef3+=adcData[i];
	}
        else if( channel==(CEF3_ADC_START_CHANNEL+2) ) {
	  h1_cef3_2->Fill(adcData[i]);
	  cef3+=adcData[i];
	}
        else if( channel==(CEF3_ADC_START_CHANNEL+3) ) {
	  h1_cef3_3->Fill(adcData[i]);
	  cef3+=adcData[i];
	}
	h1_cef3_tot->Fill(cef3);
      }
      
    }
  }




  h1_cef3_0->SetLineWidth(2);
  h1_cef3_1->SetLineWidth(2);
  h1_cef3_2->SetLineWidth(2);
  h1_cef3_3->SetLineWidth(2);

  h1_cef3_0->SetLineColor(kBlack);
  h1_cef3_1->SetLineColor(kRed);
  h1_cef3_2->SetLineColor(kBlue);
  h1_cef3_3->SetLineColor(kMagenta);

  h1_cef3_corr_0->SetLineWidth(2);
  h1_cef3_corr_1->SetLineWidth(2);
  h1_cef3_corr_2->SetLineWidth(2);
  h1_cef3_corr_3->SetLineWidth(2);

  h1_cef3_corr_0->SetLineColor(kBlack);
  h1_cef3_corr_1->SetLineColor(kRed);
  h1_cef3_corr_2->SetLineColor(kBlue);
  h1_cef3_corr_3->SetLineColor(kMagenta);

  //plot uncorr energy (before intercalibration)
  TCanvas* cuncorr = new TCanvas( "cuncorr", "", 600, 600 );



  cuncorr->cd();
  cuncorr->SetLogy();

  TH1D* histo_axes = new TH1D("cef3_0",   "", 400, 80., 350.);


  histo_axes->GetXaxis()->SetRangeUser(80.,350.);
  histo_axes->GetYaxis()->SetRangeUser(10.,h1_cef3_0->GetMaximum()+h1_cef3_0->GetMaximum()*0.10);
  histo_axes->SetXTitle( "ADC Counts" );
  histo_axes->Draw();


  h1_cef3_0->Draw("same");
  h1_cef3_1->Draw("same");
  h1_cef3_2->Draw("same");
  h1_cef3_3->Draw("same");

  cuncorr->SaveAs("uncorrEnergyAllChannels.png");
  cuncorr->SaveAs("uncorrEnergyAllChannels.eps");

  TCanvas* cuncorr_2 = new TCanvas( "cuncorr_2", "", 600, 600 );

  cuncorr_2->cd();
  cuncorr_2->SetLogy();

  histo_axes->GetXaxis()->SetRangeUser(120.,200.);
  histo_axes->GetYaxis()->SetRangeUser(10.,h1_cef3_0->GetMaximum()+h1_cef3_0->GetMaximum()*0.10);
  histo_axes->SetXTitle( "ADC Counts" );
  histo_axes->Draw();


  h1_cef3_0->Draw("same");
  h1_cef3_1->Draw("same");
  h1_cef3_2->Draw("same");
  h1_cef3_3->Draw("same");

  cuncorr_2->SaveAs("uncorrEnergyAllChannels_zoom.png");
  cuncorr_2->SaveAs("uncorrEnergyAllChannels_zoom.eps");



  //intercalibration of single fibers with photoelectrons
  FitResults fr_0 = fitSingleHisto( h1_cef3_0, 110., 135., 138., 185. );
  FitResults fr_1 = fitSingleHisto( h1_cef3_1, 100., 125., 125., 190. );
  FitResults fr_2 = fitSingleHisto( h1_cef3_2, 100., 125., 128., 190. );
  FitResults fr_3 = fitSingleHisto( h1_cef3_3, 100., 125., 137., 200. );

  std::vector<float> correctionFactors = intercalibrateFibers(fr_0,fr_1,fr_2,fr_3);

  for(int i=0;i<4;i++){
    std::cout<<correctionFactors[i]<<std::endl;
  }

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {
    
    tree->GetEntry(iEntry);
    
    if( iEntry % 5000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;
    
    for( unsigned i=0; i<40; ++i ) {
      
      int board  = adcBoard[i];
      int channel= adcChannel[i];
      
      
      float cef3_corr=0;
      
      
      if( board==CEF3_ADC_BOARD ) {
	if( channel==(CEF3_ADC_START_CHANNEL  ) ){
	  h1_cef3_corr_0->Fill(adcData[i]*correctionFactors[0]);
	  cef3_corr+=adcData[i]*correctionFactors[0];
	}
        else if( channel==(CEF3_ADC_START_CHANNEL+1) ){
	  h1_cef3_corr_1->Fill(adcData[i]*correctionFactors[1]);
	  cef3_corr+=adcData[i]*correctionFactors[1];
	}
        else if( channel==(CEF3_ADC_START_CHANNEL+2) ) {
	  h1_cef3_corr_2->Fill(adcData[i]*correctionFactors[2]);
	  cef3_corr+=adcData[i]*correctionFactors[2];
	}
        else if( channel==(CEF3_ADC_START_CHANNEL+3) ) {
	  h1_cef3_corr_3->Fill(adcData[i]*correctionFactors[3]);
	  cef3_corr+=adcData[i]*correctionFactors[3];
	}
	h1_cef3_tot->Fill(cef3_corr);
      }
      
    }
  }




  //plot corr energy (after intercalibration)
  TCanvas* ccorr = new TCanvas( "ccorr", "", 600, 600 );

  ccorr->cd();
  ccorr->SetLogy();

  TH2D* h2_axes_corr = new TH2D("axes_corr", "", 10, 80., 350., 10, 10., 1.1*h1_cef3_corr_0->GetMaximum() );
  h2_axes_corr->SetXTitle( "ADC Counts" );
  h2_axes_corr->Draw();


  h1_cef3_corr_0->Draw("same");
  h1_cef3_corr_1->Draw("same");
  h1_cef3_corr_2->Draw("same");
  h1_cef3_corr_3->Draw("same");

  ccorr->SaveAs("corrEnergyAllChannels.png");
  ccorr->SaveAs("corrEnergyAllChannels.eps");

  TCanvas* ccorr_2 = new TCanvas( "ccorr_2", "", 600, 600 );

  ccorr_2->cd();
  ccorr_2->SetLogy();

  h2_axes_corr->GetXaxis()->SetRangeUser(120.,200.);

  h2_axes_corr->Draw();
  h1_cef3_corr_0->Draw("same");
  h1_cef3_corr_1->Draw("same");
  h1_cef3_corr_2->Draw("same");
  h1_cef3_corr_3->Draw("same");

  ccorr_2->SaveAs("corrEnergyAllChannels_zoom.png");
  ccorr_2->SaveAs("corrEnergyAllChannels_zoom.eps");


  //fit corrected histos
  FitResults fr_corr_0 = fitSingleHisto( h1_cef3_corr_0, 0.966779*120., 0.966779*130., 0.966779*138., 0.966779*185. );
  FitResults fr_corr_1 = fitSingleHisto( h1_cef3_corr_1, 1.04148*100., 1.04148*125., 1.04148*125., 1.04148*189. );
  FitResults fr_corr_2 = fitSingleHisto( h1_cef3_corr_2, 1.01311*100., 1.01311*125., 1.01311*128., 1.01311*189. );
  FitResults fr_corr_3 = fitSingleHisto( h1_cef3_corr_3, 0.984112*100., 0.984112*125., 0.984112*137., 0.984112*200. );



  h1_cef3_0->Write(); 
  h1_cef3_1->Write(); 
  h1_cef3_2->Write(); 
  h1_cef3_3->Write(); 
  h1_cef3_tot->Write(); 
                 
  h1_cef3_corr_0->Write();    
  h1_cef3_corr_1->Write();    
  h1_cef3_corr_2->Write();    
  h1_cef3_corr_3->Write();    
  h1_cef3_corr_tot->Write();  

  outFile->Write(); 


  doSummaryPlots(fr_0,fr_1,fr_2,fr_3,fr_corr_0,fr_corr_1,fr_corr_2,fr_corr_3);



  return 0;

}




