#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TVector2.h"
#include "TMath.h"

#include "interface/DrawTools.h"


void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name );
TF1* fitSingleElectronPeak( const std::string& outputdir, const std::string& runName, TTree* tree );


int main() {


  std::string outputdir = "EnergyScanPlots/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  DrawTools::setStyle();
 
  std::vector<std::string> runs; 
  std::vector<float> beamEnergy; 


  //runs.push_back("BTF_314_20140503-024715_beam");
  //beamEnergy.push_back(98.3);

  //runs.push_back("BTF_308_20140503-002534_beam");
  //beamEnergy.push_back(147.4);

  //runs.push_back("BTF_293_20140502-180258_beam");
  //beamEnergy.push_back(196.5);

  runs.push_back("BTF_286_287");
  beamEnergy.push_back(294.8-2.);

  runs.push_back("BTF_259_20140502-012847_beam");
  beamEnergy.push_back(491.4-2.);

  //runs.push_back("BTF_246_20140501-212512_beam");
  //beamEnergy.push_back(491.4);




  TGraphErrors* gr_resp_vs_energy = new TGraphErrors(0);
  TGraphErrors* gr_reso_vs_energy = new TGraphErrors(0);
  TGraphErrors* gr_Sres_vs_energy = new TGraphErrors(0);


  for( unsigned i=0; i<runs.size(); ++i ) {

    TFile* file = TFile::Open(Form("PosAn_%s.root", runs[i].c_str()));
    TTree* tree = (TTree*)file->Get("tree_passedEvents");

    TF1* thisFunc = fitSingleElectronPeak( outputdir, runs[i], tree );

    float energy = beamEnergy[i];
    float energyErr = (energy>200.) ? 0.01*energy : 0.05*energy;
    //float energyErr = (energy>200.) ? 0.005*energy : 0.05*energy;

    float mean = thisFunc->GetParameter(1);
    float meanErr = thisFunc->GetParError(1);

    float rms = thisFunc->GetParameter(2);
    float rmsErr = thisFunc->GetParError(2);

    float reso = 100.* rms/mean; //in percent
    float resoErr = 100.* sqrt( rmsErr*rmsErr/(mean*mean) + rms*rms*meanErr*meanErr/(mean*mean*mean*mean) );

    gr_resp_vs_energy->SetPoint( i, energy, mean );
    gr_resp_vs_energy->SetPointError( i, energyErr, meanErr );

    gr_reso_vs_energy->SetPoint( i, energy, reso );
    gr_reso_vs_energy->SetPointError( i, energyErr, resoErr );

    float energyGeV = energy/1000.;
    float energyGeVErr = energyErr/1000.;
    float Sres = reso*sqrt(energyGeV);
    //float Sres_error = sqrt( rmsErr*rmsErr*energyGeV );
    //float Sres_error = sqrt( rms*rms*0.5*0.5*energyGeVErr*energyGeVErr/energyGeV );
    float Sres_error = sqrt( resoErr*resoErr*energyGeV + reso*reso*0.5*0.5*energyGeVErr*energyGeVErr/energyGeV );
    gr_Sres_vs_energy->SetPoint( i, energy, Sres );
    gr_Sres_vs_energy->SetPointError( i, energyErr, Sres_error );

  }


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  float xMax = 550.;

  TH2D* h2_axes = new TH2D( "axes", "", 10, 0., xMax, 10, 0., 4000. );
  h2_axes->SetXTitle("Electron Beam Energy [MeV]");
  h2_axes->SetYTitle("CeF3 Response [ADC Counts]");
  h2_axes->Draw("");

  gr_resp_vs_energy->SetMarkerStyle(20);
  gr_resp_vs_energy->SetMarkerSize(1.6);
  gr_resp_vs_energy->SetMarkerColor(46);
  gr_resp_vs_energy->Draw("p same");

  TF1* f1_line = new TF1("line", "[0] + [1]*x", 0., xMax );
  //f1_line->FixParameter(0, 0.);
  gr_resp_vs_energy->Fit(f1_line, "RN");
  f1_line->SetLineWidth(1.);
  f1_line->Draw("L same");

  std::cout << "Chi2: " << f1_line->GetChisquare() << std::endl;
  std::cout << "Prob: " << TMath::Prob(f1_line->GetChisquare(), 2 ) << std::endl;

  //TF1* f1_line2 = new TF1("line2", "[0] + [1]*x", 0., xMax );
  //gr_resp_vs_energy->Fit(f1_line2, "RN");
  //f1_line2->SetLineWidth(1.);
  //f1_line2->SetLineColor(46);
  //f1_line2->SetRange(0., xMax);
  //f1_line2->Draw("L same");

  c1->SaveAs( Form( "%s/resp_vs_energy.eps", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/resp_vs_energy.png", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/resp_vs_energy.pdf", outputdir.c_str() ) );


  c1->Clear();

  TH2D* h2_axes2 = new TH2D( "axes", "", 10, 0., xMax, 10, 0., 50. );
  h2_axes2->SetXTitle("Electron Beam Energy [MeV]");
  h2_axes2->SetYTitle("CeF3 Resolution [%]");
  h2_axes2->Draw("");

  gr_reso_vs_energy->SetMarkerStyle(20);
  gr_reso_vs_energy->SetMarkerSize(1.6);
  gr_reso_vs_energy->SetMarkerColor(46);
  gr_reso_vs_energy->Draw("p same");

  gr_Sres_vs_energy->SetMarkerStyle(21);
  gr_Sres_vs_energy->SetMarkerSize(1.6);
  gr_Sres_vs_energy->SetMarkerColor(38);
  gr_Sres_vs_energy->Draw("p same");

  TF1* f1_const = new TF1("const", "[0]", 0., xMax );
  gr_Sres_vs_energy->Fit(f1_const, "RN");
  f1_const->SetLineWidth(1.);
  f1_const->Draw("L same");


  c1->SaveAs( Form( "%s/reso_vs_energy.eps", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/reso_vs_energy.png", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/reso_vs_energy.pdf", outputdir.c_str() ) );


  return 0;

}



TF1* fitSingleElectronPeak( const std::string& outputdir, const std::string& runName, TTree* tree ) {

  std::string histoName(Form("h1_%s", runName.c_str()));
  TH1D* h1 = new TH1D(histoName.c_str(), "", 200, 0., 6000.);
  tree->Project( histoName.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1)");

  TF1* f1 = new TF1( Form("gaus_%s", runName.c_str()), "gaus", 400., 6000.);
  f1->SetParameter(0, h1->Integral() );
  f1->SetParameter(1, h1->GetMean() );
  f1->SetParameter(2, h1->GetRMS() );

  f1->SetParError(1, h1->GetMeanError() );
  f1->SetParError(2, h1->GetRMSError() );

  doSingleFit( h1, f1, outputdir, runName );

  return f1;

}




void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name ) {

  h1->Fit( f1, "RQN" );

  int niter = 4.;
  float nSigma = 1.5;

  for( unsigned iter=0; iter<niter; iter++ ) {

    float mean  = f1->GetParameter(1);
    float sigma = f1->GetParameter(2);
    float fitMin = mean - nSigma*sigma;
    float fitMax = mean + nSigma*sigma;
    f1->SetRange( fitMin, fitMax );
    if( iter==(niter-1) )
      h1->Fit( f1, "RQN" );
    else
      h1->Fit( f1, "RQ+" );
  }


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  h1->Draw();

  c1->SaveAs( Form("%s/fit_%s.eps", outputdir.c_str(), name.c_str()) );
  c1->SaveAs( Form("%s/fit_%s.png", outputdir.c_str(), name.c_str()) );

  delete c1;

}
