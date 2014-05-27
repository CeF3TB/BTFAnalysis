#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h> 

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"

#include "interface/DrawTools.h"



void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name );
TF1* fitSingleElectronPeak( const std::string& outputdir, int i, TTree* tree );
TF1* checkTotalResolution( const std::string& outputdir, TTree* tree );
float sumVector( std::vector<float> v );
bool savePlots=false;

int main( int argc, char* argv[] ) {

  std::string inputDir = "./analysisTrees";
  std::string runName = "BTF_92_20140430-020137_beam_uncalib";
  std::string tag = "default";

  if( argc == 3 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
    std::string tag_str(argv[2]);
    tag = tag_str;
  }else if (argc==4){
    std::string inputDir_str(argv[2]);
    inputDir =inputDir_str;
  } else{
    std::cout<<"Usage:"<<std::endl;
    std::cout<<"./calibrateCef3 BTF_XXX tag [inputDir]"<<std::endl;
    exit(12345);
  }


  DrawTools::setStyle();


  TFile* file = TFile::Open(Form("%s/Reco_%s_%s.root", inputDir.c_str(),runName.c_str(),tag.c_str()));
  std::cout<<"opening file:"<<file->GetName();


  TTree* tree = (TTree*)file->Get("recoTree");

  std::string outputdir = "CeF3Calibration/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  std::string ofsName = outputdir + "/constants.txt";
  ofstream ofs(ofsName.c_str());

  int nentries = tree->GetEntries();

  std::vector<float> cef3_calibration;

  for( unsigned i=0; i<4; ++i ) {
    TF1* f1 = fitSingleElectronPeak( outputdir, i, tree );
    float mean  = f1->GetParameter(1);
    float sigma = f1->GetParameter(2);
    std::cout << std::endl;
    std::cout << "Channel " << i << std::endl;
    std::cout << "  Mean       : " << mean << std::endl;
    std::cout << "  Sigma      : " << sigma << std::endl;
    std::cout << "  Resolution : " << sigma/mean << std::endl;

    cef3_calibration.push_back(mean);
  }


  float cef3CalibrationAverage = sumVector(cef3_calibration)/cef3_calibration.size();

  for(unsigned i=0; i<cef3_calibration.size(); ++i ){
    cef3_calibration[i] = cef3CalibrationAverage/cef3_calibration[i];
    ofs << i << "\t" << cef3_calibration[i] << std::endl;
  }


  ofs.close();



  std::cout << "-> Saved constants in: " << ofsName << std::endl;


  checkTotalResolution( outputdir, tree );

  return 0;

}



TF1* fitSingleElectronPeak( const std::string& outputdir, int i, TTree* tree ) {

  std::string histoName(Form("h1_%d", i));
  TH1D* h1 = new TH1D(histoName.c_str(), "", 100, 0., 3000.);
  //tree->Project( histoName.c_str(), Form("cef3_corr[%d]", i), "");
  tree->Project( histoName.c_str(), Form("cef3_corr[%d]", i), "(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1)");

  TF1* f1 = new TF1( Form("gaus_%d", i), "gaus", 400., 1200.);
  f1->SetParameter(0, 3000.);
  f1->SetParameter(1, 800.);
  f1->SetParameter(2, 150.);

  doSingleFit( h1, f1, outputdir, Form("%d", i) );

  return f1;

}



TF1* checkTotalResolution( const std::string& outputdir, TTree* tree ) {

  std::string histoName("h1_tot");
  TH1D* h1 = new TH1D(histoName.c_str(), "", 500, 0., 12000.);
  //tree->Project( histoName.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "");
  tree->Project( histoName.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1)");

  TF1* f1 = new TF1("gaus_tot", "gaus", 1600., 4800.);
  f1->SetParameter(0, 3000.);
  f1->SetParameter(1, 3000.);
  f1->SetParameter(2, 600.);

  doSingleFit( h1, f1, outputdir, "tot" );

  std::cout << std::endl;
  std::cout << std::endl;
  float mean  = f1->GetParameter(1);
  float sigma = f1->GetParameter(2);
  float reso = sigma/mean;
  std::cout << "Total " << std::endl;
  std::cout << "  Mean       : " << mean << std::endl;
  std::cout << "  Sigma      : " << sigma << std::endl;
  std::cout << "  Resolution : " << sigma/mean << std::endl;
  std::cout << "Corresponds to a stochastic term of: " << 100.*reso*sqrt(0.5) << " %" << std::endl;

  return f1;

}


void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name ) {

  h1->Fit( f1, "RQN" );

  int niter = 4.;
  float nSigma = 2.;

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

  if(savePlots){
    c1->SaveAs( Form("%s/fit_%s.eps", outputdir.c_str(), name.c_str()) );
    c1->SaveAs( Form("%s/fit_%s.png", outputdir.c_str(), name.c_str()) );
  }

  delete c1;

}

float sumVector( std::vector<float> v ) {

  float sum=0.;
  for( unsigned i=0; i<v.size(); ++i ) sum += v[i];

  return sum;

}
