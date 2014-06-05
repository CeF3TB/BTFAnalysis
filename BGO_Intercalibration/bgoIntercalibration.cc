#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TLatex.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"

#include "../PositionAnalysis/fastDQM_CeF3_BTF.h"
#include "../PositionAnalysis/interface/DrawTools.h"

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooAbsReal.h"
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "RooPlot.h"
using namespace RooFit ;

#define BGO_CHANNELS 8
enum {ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8};
Double_t bgo_index    [BGO_CHANNELS] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5};
Double_t bgo_index_err[BGO_CHANNELS] = {0.,0.,0.,0.,0.,0.,0.,0.};
TString BGOchannel[BGO_CHANNELS] = {"0", "1", "2", "3", "4", "5", "6", "7" };
TString BGOnumber [BGO_CHANNELS] = {"1", "2", "3", "4", "5", "6", "7", "8" };

std::string BGOrun[BGO_CHANNELS] = {
  "BTF_530_20140505-014707_beam", // CHANNEL 0
  "BTF_538_20140505-053339_beam", // CHANNEL 1 
  "BTF_537_20140505-050526_beam", // CHANNEL 2 
  "BTF_531_20140505-020502_beam", // CHANNEL 3 
  "BTF_536_20140505-043018_beam", // CHANNEL 4 
  "BTF_533_20140505-024712_beam", // CHANNEL 5 
  "BTF_534_20140505-031520_beam", // CHANNEL 6  
  "BTF_535_20140505-035012_beam", // CHANNEL 7 
};

TH1D* fitSingleChannelBGO( const std::string& outputdir, const std::string& name, const std::string& runName, int iChannel, double &mean_output, double &meanError_output, float calibConst=1.);
TH1D* fitSingleChannel   ( const std::string& outputdir, const std::string& name, const std::string& runName, int iChannel, const std::string& varName, const std::string& plotName, 
			   int nBins, float xMin, float xMax, double &mean_output, double &meanError_output );

void drawHistos( const std::string& outputdir, std::vector<TH1D*> histos, const std::string& name, float yMax, float xMax = 4000. );
float sumVector( std::vector<float> v );

std::pair<float,float> getMeanRMS(double mean[BGO_CHANNELS]);

// -- Crystal Ball Function
Double_t CrystalBall(Double_t *x, Double_t *par){
  
  Double_t t = (x[0]-par[2])/par[3];
  if (par[0] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[0]);
  if (t >= -absAlpha) { return par[4]*exp(-0.5*t*t); }
  else {
    Double_t a = TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[1]/absAlpha - absAlpha;
    return par[4]*(a/TMath::Power(b - t, par[1]));
  }

}

int main(){

  // -- create folder to store the txt files and plots
  std::string outputdir = "";
  outputdir = "BGOInterCalibration_CB_28May_bis";

  std::string mkdir_command = "mkdir -p " + outputdir;
  system(mkdir_command.c_str());
  
  // -- define histograms
  std::vector<TH1D*> rawHistos;
  std::vector<TH1D*> calibHistos;
  std::vector<float> calibConstants;
  float yMax = 0.;
  
  // -- Raw distributions
  std::string meansRawFileName;
  meansRawFileName = outputdir + "/meansRaw_CB.txt";
  ofstream ofsmr(meansRawFileName.c_str());

  double rawMean    [BGO_CHANNELS];
  double rawMeanErr [BGO_CHANNELS];
  double channelRawMean    = 0.0;
  double channelRawMeanErr = 0.0;
  for( unsigned i=0; i<BGO_CHANNELS; ++i ) {
    TH1D* h1_raw = fitSingleChannelBGO( outputdir, "raw", BGOrun[i], i, channelRawMean, channelRawMeanErr, 1. );
    rawHistos.push_back(h1_raw);
    
    double mean = channelRawMean;
    double err  = channelRawMeanErr;

    rawMean   [i] = mean;
    rawMeanErr[i] = err;
    ofsmr << i << "\t " << mean << "\t +- " << err << std::endl;

    calibConstants.push_back(mean);
    
    float thisMax = h1_raw->GetMaximum()/h1_raw->Integral();
    if( thisMax>yMax )
      yMax = thisMax;
  }
  
  ofsmr.close();
  
  float calibAve = sumVector(calibConstants)/calibConstants.size();
  
  // -- create txt with constants for intercalibration
  std::string constantsFileName;
  constantsFileName = outputdir + "/constants_V0.txt";
  ofstream ofs(constantsFileName.c_str());

  // -- Calibrated distributions
  std::string meansCalibFileName;
  meansCalibFileName = outputdir + "/meansCalib_CB.txt";
  ofstream ofsmc(meansCalibFileName.c_str());

  double calibratedMean    [BGO_CHANNELS];
  double calibratedMeanErr [BGO_CHANNELS];
  double channelCalibMean    = 0.0;
  double channelCalibMeanErr = 0.0;
  for( unsigned i=0; i<BGO_CHANNELS; ++i ) {
    float thisCalib = calibAve/calibConstants[i];
    calibHistos.push_back(fitSingleChannelBGO( outputdir, "calib", BGOrun[i], i, channelCalibMean, channelCalibMeanErr, thisCalib ));
   
    // ofs << i << "\t " << thisCalib << std::endl; // Old format, with channel
    ofs << thisCalib << std::endl; // New VXX format, without channel

    double mean = channelCalibMean;
    double err  = channelCalibMeanErr;

    calibratedMean   [i] = mean;
    calibratedMeanErr[i] = err;
    ofsmc << i << "\t " << mean << "\t +- " << err << std::endl;

  }

  ofs.close();
  ofsmc.close();

  drawHistos( outputdir, rawHistos,   "rawSpectra"  , yMax );
  drawHistos( outputdir, calibHistos, "calibSpectra", yMax );

  std::cout << std::endl;
  std::cout << "-> Calibration constants saved in: " << constantsFileName << std::endl;
  std::cout << "Calibration average for BGO: " << calibAve << std::endl;
  std::cout << std::endl;

  // -- estimate the RMS for giving the intercalibration estimate
  std::pair<float,float>        rawValues = getMeanRMS(        rawMean );
  std::pair<float,float> calibratedValues = getMeanRMS( calibratedMean );

  double rawIntercalibration = 100*rawValues.second/rawValues.first;
  double calIntercalibration = 100*calibratedValues.second/calibratedValues.first;

  // -- print the intercalibration values
  std::cout << " All: RAW        (mean +- rms) = ( " << rawValues.first        << " +- " << rawValues.second        << " ) " << std::endl; 
  std::cout << " ---> Intercalibration: " << rawIntercalibration << " % " << std::endl;
  std::cout << "" << std::endl;
  std::cout << " All: CALIBRATED (mean +- rms) = ( " << calibratedValues.first << " +- " << calibratedValues.second << " ) " << std::endl; 
  std::cout << " ---> Intercalibration: " << calIntercalibration << " % " << std::endl;
  
  TCanvas *cg = new TCanvas("mean graph RAW","mean graph RAW");
  TGraphErrors *grRaw = new TGraphErrors(BGO_CHANNELS, bgo_index, rawMean, bgo_index_err, rawMeanErr);
  TGraphErrors *grCal = new TGraphErrors(BGO_CHANNELS, bgo_index, calibratedMean, bgo_index_err, calibratedMeanErr);
  cg->cd();

  TH2D* h2_axes = new TH2D("axes", "", 8, 0., 8., 10, 2300., 3000. );
  h2_axes->SetXTitle( "BGO channel" );
  h2_axes->SetYTitle( "mean ADC counts" );
  h2_axes->SetTitleOffset(1.3,"Y");
  h2_axes->Draw("");
  h2_axes->SetStats(0);

  grRaw->SetMarkerStyle(20);
  grRaw->SetMarkerSize(1.0);
  grRaw->SetMarkerColor(kBlack);
  grRaw->Draw("P");

  grCal->SetMarkerStyle(21);
  grCal->SetMarkerSize(1.0);
  grCal->SetMarkerColor(kRed);
  grCal->Draw("P");
  
  TLegend* legend = new TLegend( 0.55, 0.70, 0.85, 0.85 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->SetTextFont(40);
  legend->SetHeader("mean value from fit");
  legend->AddEntry( grRaw, "Raw", "p");
  legend->AddEntry( grCal, "Calibrated", "p");
  legend->Draw("same");

  cg->SaveAs( Form("%s/graphRaw.png", outputdir.c_str() ) );



}

TH1D* fitSingleChannelBGO( const std::string& outputdir, const std::string& name, const std::string& runName, int iChannel, double &mean_output, double &meanError_output, float calibConst ) {

  return fitSingleChannel( outputdir, name, runName, iChannel, Form("bgo_pedSubtracted[%d]*%f", iChannel, calibConst), Form("BGO Channel %d", iChannel), 50, 1000., 3500., mean_output, meanError_output );

}

TH1D* fitSingleChannel( const std::string& outputdir, const std::string& name, const std::string& runName,  int iChannel, const std::string& varName, const std::string& plotName, 
			int nBins, float xMin, float xMax, double &mean_output, double &meanError_output ) {

  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << "##################################################################### " << std::endl;
  std::cout << name << " running file " << runName << std::endl;
  std::cout << "##################################################################### " << std::endl;
  std::cout << " " << std::endl;
  // -- open the roofiles with the BGOs with pedestal substracted
  TFile* file = TFile::Open(Form("analysisTrees/Reco_BTF_%s_beam.root", runName.c_str()));
  TTree* tree = (TTree*)file->Get("recoTree");

  // -- define TH1D histos to be used for the fit
  std::string histoName = runName;
  TH1D* hBGO = new TH1D(histoName.c_str(), "", nBins, xMin, xMax );
  hBGO->Sumw2();

  // -- select single electron events: front scintillator entries in range [500, 2000]
  tree->Project( histoName.c_str(), varName.c_str(), "scintFront>500. && scintFront<2000.");

  // -- fit first a gaussian to get the inizialitation values for mean and sigma of the peak
  TF1* faux = new TF1( Form("gaussianAux_%s", runName.c_str()), "gaus(0)" );
  faux->SetRange(1800., 3500.);
  faux->SetParameter(1, 2500.);
  faux->SetParameter(2, 200.);
  hBGO->Fit( faux, "RQN");
  
  // -- initialize crystal ball parameters, using the gaussian fit and the histo normalization
  double mean  = faux->GetParameter(1);
  double sigma = faux->GetParameter(2);
  double norm  = hBGO->Integral();
  
  // -- data
  RooRealVar     x(   "x", "ADC counts", 1000., 3500.);
  RooDataHist data("data",    "dataset", x, hBGO);
  
  // -- parameters for Crystal Ball: mean, sigma, normalization, n and alpha
  RooRealVar cbmean (  "cbmean",  "cbmean" ,  mean, 2000., 3000.); // inizialite from gauss fit
  RooRealVar cbsigma( "cbsigma", "cbsigma" , sigma,  100.,  300.); // inizialite from gauss fit
  RooRealVar cbsig  (   "cbsig", "cbsignal",  norm,    0., 5000.); // inizialite from histogram
  RooRealVar n      (       "n",         "",  10.0,    0.,  200.);
  RooRealVar alpha  (   "alpha",         "",   1.0,    0.,   10.);

  // -- pdf: Crystal Ball
  RooCBShape cball  (  "cball", "crystal ball", x, cbmean, cbsigma, alpha, n);

  // -- fit
  cball.fitTo(data);
  
  // -- plot
  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  
  RooPlot* xframe = x.frame();
  //data.plotOn(xframe);
  //cball.plotOn(xframe,LineColor(kBlue));    

  data.plotOn (xframe, Name("data"));
  cball.plotOn(xframe, LineColor(kBlue), Name("cball"));    

  xframe->SetName("");
  xframe->SetTitle("");
  xframe->SetTitleOffset(1.4,"Y");
  xframe->Draw();
  // double chi2 = xframe->chiSquare();
  double chi2 = xframe->chiSquare("cball", "data", 4);
  double ndof = xframe->GetNbinsX();
  
  // -- draw BGO channel
  TPaveText* labelChan = new TPaveText( 0.2, 0.7, 0.3, 0.8, "brNDC");
  labelChan->SetFillColor(0);
  labelChan->SetTextSize(0.038);
  labelChan->SetTextFont(42);
  labelChan->SetTextAlign(10);
  TString header = "";
  if (name.find("raw")  !=std::string::npos) header = "Raw";
  if (name.find("calib")!=std::string::npos) header = "Calibrated";
  labelChan->AddText(header);
  labelChan->AddText(plotName.c_str() );
  labelChan->Draw("same");
  
  // -- draw run number
  TPaveText* labelChan2 = new TPaveText( 0.1, 0.9, 0.7, 1.0, "brNDC");
  labelChan2->SetFillColor(0);
  labelChan2->SetTextSize(0.038);
  labelChan2->SetTextFont(42);
  labelChan2->SetTextAlign(10);
  labelChan2->AddText( runName.c_str() );
  // labelChan2->Draw("same");
  
  // -- draw fitted CB parameters and chi2
  TPaveText* labelChan3 = new TPaveText(0.2, 0.35, 0.3, 0.65, "brNDC");
  labelChan3->SetFillColor(0);
  labelChan3->SetTextSize(0.032);
  labelChan3->SetTextFont(42);
  labelChan3->SetTextAlign(10);
  labelChan3->AddText( "Crystal Ball" );
  labelChan3->AddText( "" );
  labelChan3->AddText( TString::Format(   "\\mu = %.1f \\pm %.1f",  cbmean.getVal(),  cbmean.getError()) );
  labelChan3->AddText( TString::Format("\\sigma = %.0f \\pm %.0f", cbsigma.getVal(), cbsigma.getError()) );
  labelChan3->AddText( TString::Format("\\alpha = %.2f \\pm %.2f",   alpha.getVal(),   alpha.getError()) );
  labelChan3->AddText( "" );
  labelChan3->AddText( TString::Format("\\chi^{2} = %.1f", chi2) );
  labelChan3->Draw("same");
  
  // -- fill the mean and its error
  mean_output      = cbmean.getVal();
  meanError_output = cbmean.getError();

  // -- some prints on screen, to check the fit
  std::cout << "" << std::endl;
  std::cout << " Mean value for distribution is: (" << mean_output << " +- " << meanError_output << " ), with Chi2 for the fit: " << chi2 << std::endl;
  std::cout << "------------------------------------------------------------------" << std::endl;
  std::cout << plotName.c_str() << " " << runName.c_str() << std::endl;
  std::cout << "------------------------------------------------------------------" << std::endl; std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::endl;
  std::cout << " ... [AFTER ] Fitted parameters for CB " << std::endl;
  std::cout << "  alfa  = ( " << alpha.getVal()   << "   +- " << alpha.getError()   << " ) " << std::endl;
  std::cout << "  n     = ( " << n.getVal()       << "   +- " << n.getError()       << " ) " << std::endl; std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::endl;
  std::cout << "  mean  = ( " << mean_output      << " +- "   << meanError_output   << " ) " << std::endl;
  std::cout << "  sigma = ( " << cbsigma.getVal() << "  +- "  << cbsigma.getError() << " ) " << std::endl; std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::endl;
  std::cout << "  N     = ( " << cbsig.getVal()   << " +- "   << cbsig.getError()   << " ) " << std::endl;

  // -- save plots
  c1->SaveAs( Form("%s/fit_%s_%s.png", outputdir.c_str(), name.c_str(), runName.c_str()) );
  
  delete c1;
  
  // -- return the BGO adc counts distribution
  return hBGO;

}

float sumVector( std::vector<float> v ) {

  float sum=0.;
  for( unsigned i=0; i<v.size(); ++i ) sum += v[i];

  return sum;

}

void drawHistos( const std::string& outputdir, std::vector<TH1D*> histos, const std::string& name, float yMax, float xMax ) {

  std::vector<int> colors;
  
  colors.push_back( 40 );
  colors.push_back( 41 );
  colors.push_back( 42 );
  colors.push_back( 45 );
  colors.push_back( 46 );
  colors.push_back( 30 );
  colors.push_back( 38 );
  colors.push_back( 16  );

  TCanvas* c2 = new TCanvas("c2", "", 600, 600);
  c2->cd();

  TLegend* legend = new TLegend( 0.15, 0.45, 0.35, 0.85 );
  TString header = "";
  if (name.find("raw")  !=std::string::npos) header = "Raw";
  if (name.find("calib")!=std::string::npos) header = "Calibrated";
  legend->SetHeader(header);
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->SetTextFont(40);
  legend->SetBorderSize(0);

  TLatex* luminosity = new TLatex(0.6, 0.92, "No hodoscopes / Electron Beam");
  luminosity->SetNDC();
  luminosity->SetTextAlign(32);
  luminosity->SetTextFont(42);
  luminosity->SetTextSize(0.035);

  
  double ymax_ = 0.0;
  for( unsigned i=0; i<histos.size(); ++i ) {
    // double thisMax = histos[i]->GetMaximum()/histos[i]->Integral();
    double thisMax = histos[i]->GetMaximum()*1; // -- multiply by rebin size
    if( thisMax>ymax_ )
      ymax_ = thisMax;
  }

  TH2D* h2_axes = new TH2D("axes", "", 10, 1000., 3500., 10, 0., 1.1*ymax_ );
  h2_axes->SetXTitle( "BGO ADC counts" );
  h2_axes->SetYTitle( "events" );
  h2_axes->SetTitleOffset(1.3,"Y");
  h2_axes->Draw("");
  h2_axes->SetStats(0);

  for( unsigned i=0; i<histos.size(); ++i ) {
    // histos[i]->Rebin(2);
    histos[i]->SetLineColor( colors[i] );
    histos[i]->SetLineWidth( 2 );
    histos[i]->GetXaxis()->SetRangeUser(1000., 3500.);
    // histos[i]->DrawNormalized( "histo same" );
    histos[i]->Draw( "histo same" );
    histos[i]->SetStats(0);
    legend->AddEntry( histos[i], Form("Ch %d", i), "l" );
  }

  legend->Draw("same");
  luminosity->Draw("same");

  c2->SaveAs( Form("%s/%s.png", outputdir.c_str(), name.c_str()) );

  delete c2;
  delete h2_axes;

}

std::pair<float,float> getMeanRMS(double mean[BGO_CHANNELS]){
  
  std::pair<float, float>  thispair;

  double mymean = 0.0;
  for(int i=0;i<BGO_CHANNELS;i++){
    mymean += mean[i];
  }

  mymean = mymean / BGO_CHANNELS;
  
  double myrms = 0.0;
  for(int i=0;i<BGO_CHANNELS;i++){
    myrms += (mymean-mean[i])*(mymean-mean[i]);
  }

  myrms = myrms/BGO_CHANNELS;
  myrms = TMath::Sqrt(myrms);
  
  thispair.first  = mymean;
  thispair.second = myrms;
  
  return thispair;

}
