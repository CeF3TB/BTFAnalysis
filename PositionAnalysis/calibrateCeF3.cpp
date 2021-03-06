#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h> 

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h" 
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLine.h"

#include "interface/DrawTools.h"

#include "TApplication.h"


void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name );
TF1* fitSingleElectronPeak( const std::string& outputdir, int i, TTree* tree );
TF1* fitSingleElectronPeakCentral( const std::string& outputdir, int i, TTree* tree );
TF1* checkTotalResolution( const std::string& outputdir, TTree* tree );
void checkIntercalibration(std::vector<float> constant, std::vector<float> const_uncert, std::vector<float> const_central, std::vector<float>  const_central_uncert,const std::string& outputdir, const std::string& runName, const std::string& tag);
float sumVector( std::vector<float> v );
bool savePlots=true;
bool checkIntercal = false;

int main( int argc, char* argv[] ) {

  TApplication* a = new TApplication("a", 0, 0);
  TStyle* style = DrawTools::setStyle();
  style->cd();


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

  if(tag!="default") inputDir = inputDir + "_"+tag;


  DrawTools::setStyle();


  TFile* file = TFile::Open(Form("%s/Reco_%s.root", inputDir.c_str(),runName.c_str()));
  std::cout<<"opening file:"<<file->GetName();


  TTree* tree = (TTree*)file->Get("recoTree");

  std::string outputdir = "CeF3Calibration/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  std::string ofsName = outputdir +Form( "/constants_%s_%s.txt", runName.c_str(), tag.c_str() );
  ofstream ofs(ofsName.c_str());
  std::string ofsNameU = outputdir + Form( "/constants_uncert_%s_%s.txt", runName.c_str(), tag.c_str() );
  ofstream ofsU(ofsNameU.c_str());

  //for central constants
  std::string ofsName2= outputdir + Form( "/constants_central_%s_%s.txt", runName.c_str(), tag.c_str() );
  ofstream ofs2(ofsName2.c_str());
  std::string ofsName2U = outputdir + Form("/constants_central_uncert_%s_%s.txt",runName.c_str(), tag.c_str() ) ;
  ofstream ofs2U(ofsName2U.c_str());

  std::vector<float> cef3_calibration_central;
  std::vector<float> cef3_calib_central_uncert;

  int nentries = tree->GetEntries();

  std::vector<float> cef3_calibration;
  std::vector<float> cef3_calib_uncert;



  for( unsigned i=0; i<4; ++i ) {
    TF1* f1 = fitSingleElectronPeak( outputdir, i, tree );
    float mean  = f1->GetParameter(1);
    float sigma = f1->GetParameter(2);
    float mean_err = f1->GetParError(1);
    std::cout << std::endl;
    std::cout << "Channel " << i << std::endl;
    std::cout << "  Mean       : " << mean << std::endl;
    std::cout << "  Sigma      : " << sigma << std::endl;
    std::cout << "  Resolution : " << sigma/mean << std::endl;

    cef3_calibration.push_back(mean);
    cef3_calib_uncert.push_back(mean_err);

    TF1* f2= fitSingleElectronPeakCentral( outputdir, i, tree);
    float central_mean = f2->GetParameter(1);
    float central_mean_uncert = f2->GetParError(1);
    
    cef3_calibration_central.push_back(central_mean);
    cef3_calib_central_uncert.push_back(central_mean_uncert);
  }


  float cef3CalibrationAverage = sumVector(cef3_calibration)/cef3_calibration.size();

  for(unsigned i=0; i<cef3_calibration.size(); ++i ){
    cef3_calib_uncert[i] = abs(cef3CalibrationAverage-cef3_calibration[i]/4.)/(cef3_calibration[i]*cef3_calibration[i])*cef3_calib_uncert[i];
    ofsU << cef3_calib_uncert[i] << std::endl;

    cef3_calibration[i] = cef3CalibrationAverage/cef3_calibration[i];
    ofs << cef3_calibration[i] << std::endl;
  }



///////////////////////For central one
  float cef3CalibrationAverage_central = sumVector(cef3_calibration_central)/cef3_calibration_central.size();


  for(unsigned i=0; i<cef3_calibration_central.size(); ++i ){
    cef3_calib_central_uncert[i]= abs(cef3CalibrationAverage_central-cef3_calibration_central[i]/4)/(cef3_calibration_central[i]*cef3_calibration_central[i])*cef3_calib_central_uncert[i];
    ofs2U << cef3_calib_central_uncert[i] << std::endl;

    cef3_calibration_central[i] = cef3CalibrationAverage_central/cef3_calibration_central[i];
    ofs2 << cef3_calibration_central[i] << std::endl;
  }

  ofs.close();
  ofs2.close();
  ofsU.close();
  ofs2U.close();

  if(checkIntercal == true){
  checkIntercalibration(cef3_calibration, cef3_calib_uncert, cef3_calibration_central, cef3_calib_central_uncert, outputdir, runName, tag);
  }

  std::cout << "-> Saved constants in: " << ofsName << std::endl;


  checkTotalResolution( outputdir, tree );

  return 0;

}



TF1* fitSingleElectronPeak( const std::string& outputdir, int i, TTree* tree ) {

  std::string histoName(Form("h1_%d", i));
  TH1D* h1 = new TH1D(histoName.c_str(), "", 100, 0., 3000.);
  //tree->Project( histoName.c_str(), Form("cef3_corr[%d]", i), "");
  tree->Project( histoName.c_str(), Form("cef3_corr[%d]", i), "(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1)");

  //tree->Project( histoName.c_str(), Form("cef3_corr[%d]", i), "(isSingleEle_scintFront==1 && nHodoClustersX==1 && nHodoClustersY==1 )");

  TF1* f1 = new TF1( Form("gaus_%d", i), "gaus", 0., 1200.);
  //  TF1* f1 = new TF1( Form("gaus_%d", i), "gaus", 400., 1200.);
  f1->SetParameter(0, 3000.);
  f1->SetParameter(1, 800.);
  f1->SetParameter(2, 150.);

  doSingleFit( h1, f1, outputdir, Form("%d", i) );

  return f1;

}

//same as fitSingleElectronPeak but with |pos_hodoClust|<2 and nFibres<3
TF1* fitSingleElectronPeakCentral( const std::string& outputdir, int i, TTree* tree ) {

  std::string histoName(Form("h1_%d", i));
  TH1D* h1 = new TH1D(histoName.c_str(), "", 100, 0., 3000.);
  //tree->Project( histoName.c_str(), Form("cef3_corr[%d]", i), "");
  tree->Project( histoName.c_str(), Form("cef3_corr[%d]", i), "(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1 && -2.<pos_hodoClustX<2. && -2.<pos_hodoClustY<2. && (nFibres_hodoClustX< 3) && (nFibres_hodoClustY< 3))");
  //tree->Project( histoName.c_str(), Form("cef3_corr[%d]", i), "(isSingleEle_scintFront==1 && nHodoClustersX==1 && nHodoClustersY==1 && -2.<pos_hodoClustX<2. && -2.<pos_hodoClustY<2. && (nFibres_hodoClustX< 3) && (nFibres_hodoClustY< 3))");

  TF1* f1 = new TF1( Form("gaus_%d", i), "gaus", 0., 1200.);
  //TF1* f1 = new TF1( Form("gaus_%d", i), "gaus", 400., 1200.);
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

  //  TF1* f1 = new TF1("gaus_tot", "gaus", 1600., 4800.);
  TF1* f1 = new TF1("gaus_tot", "gaus", 100., 4800.);
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
  float nSigma =1.5 ;

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





void checkIntercalibration(std::vector<float> constant, std::vector<float> const_uncert, std::vector<float> const_central, std::vector<float>  const_central_uncert,const std::string& outputdir, const std::string& runName, const std::string& tag){

  TCanvas* canny = new TCanvas("canny", "",200,200);

  int n=4;
  float vv[4] = { -0.1,0.9,1.9,2.9 };
  std::vector<float> x(&vv[0], &vv[0]+4);
  std::vector<float> xerr (4,0);

  TGraphErrors* graf = new TGraphErrors(n,&(x[0]),&(constant[0]),&(xerr[0]),&(const_uncert[0]));
  graf->SetMarkerStyle(9);
  graf->SetMarkerSize(0.5);
  graf->SetMarkerColor(4);

  float v2[4] = { 0.1,1.1,2.1 ,3.1 };
std::vector<float> xx(&v2[0], &v2[0]+4);

  TGraphErrors* graf2 = new TGraphErrors(n,&(xx[0]),&(const_central[0]),&(xerr[0]),&(const_central_uncert[0]));
  graf2->SetMarkerStyle(9);
  graf2->SetMarkerSize(0.5);
  graf2->SetMarkerColor(kGreen+1);

  TMultiGraph *multi= new TMultiGraph();
  multi->Add(graf);
  multi->Add(graf2);
  multi->SetTitle(" ;Channel Nr.; Correction Factor");
  multi->Draw("AP");



multi->GetYaxis()->SetRangeUser(0.975,1.025);
//multi->GetYaxis()->SetRangeUser(0.9,1.1);
multi->Draw("AP");
canny->Update();


  TLine* lin = new TLine(-0.25,1.,3.25,1.);
  lin->SetLineColor(kRed);
  lin->Draw();

  TLegend* leg = new TLegend(0.6, 0.75, 0.9, 0.9); 
  leg->AddEntry(graf,"single Electron","P");
  leg->AddEntry(graf2,"+central + nFibres #\leq 2 ","P");
  leg->SetFillColor(0);
  leg->Draw("same");

  if(savePlots==1){  canny->SaveAs( Form( "%s/corrPlot_%s_%s.pdf", outputdir.c_str(), runName.c_str(), tag.c_str()  ));
}

  delete canny;

}
