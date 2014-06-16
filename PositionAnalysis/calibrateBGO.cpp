#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include <vector>
#include <string>
#include <fstream>

#include "interface/DrawTools.h"



//TH1D* fitSingleChannel( const std::string& outputdir, const std::string& name, const std::string& runName, int iChannel, float calibConst=1. );
TH1D* fitSingleChannelBGO( const std::string& outputdir, const std::string& name, const std::string& tag, const std::string& runName, int iChannel, float calibConst=1. );
TH1D* fitCeF3( const std::string& outputdir, const std::string& name, const std::string& tag, const std::string& runName );
TH1D* fitSingleChannel( const std::string& outputdir, const std::string& name, const std::string& tag, const std::string& runName,  const std::string& varName, const std::string& plotName, int nBins, float xMin, float xMax );
void drawHistos( const std::string& outputdir, std::vector<TH1D*> histos, const std::string& name, float yMax, float xMax = 4000. );
float sumVector( std::vector<float> v );



int main() {


  std::string tag = "V02";


  DrawTools::setStyle();

  std::vector<std::string> runs;
  runs.push_back( "BTF_229_20140501-155939_beam" ); // 0
  runs.push_back( "BTF_231_20140501-163742_beam" ); // 1
  runs.push_back( "BTF_233_20140501-171741_beam" ); // 2
  runs.push_back( "BTF_237_20140501-183950_beam" ); // 3
  runs.push_back( "BTF_235_20140501-175948_beam" ); // 4
  runs.push_back( "BTF_239_20140501-191908_beam" ); // 5
  runs.push_back( "BTF_241_20140501-200053_beam" ); // 6
  runs.push_back( "BTF_243_20140501-203838_beam" ); // 7


  std::string outputdir = "BGOCalibration_" + tag;
  std::string mkdir_command = "mkdir -p " + outputdir;
  system(mkdir_command.c_str());


  std::vector<TH1D*> rawHistos;
  std::vector<TH1D*> calibHistos;
  std::vector<float> calibConstants;
  float yMax = 0.;

  for( unsigned i=0; i<runs.size(); ++i ) {
    TH1D* h1_raw = fitSingleChannelBGO( outputdir, "raw", tag, runs[i], i, 1. );
    rawHistos.push_back(h1_raw);
    TF1* thisFunc = (TF1*)(h1_raw->GetListOfFunctions()->FindObject(Form("gaussian_%s", runs[i].c_str())));
    calibConstants.push_back(thisFunc->GetParameter(1));
    float thisMax = h1_raw->GetMaximum()/h1_raw->Integral();
    if( thisMax>yMax )
      yMax = thisMax;
  }


  float calibAve = sumVector(calibConstants)/calibConstants.size();


  std::string constantsFileName = outputdir + "/constants.txt";
  ofstream ofs(constantsFileName.c_str());

  for( unsigned i=0; i<runs.size(); ++i ) {
    float thisCalib = calibAve/calibConstants[i];
    calibHistos.push_back(fitSingleChannelBGO( outputdir, "calib", tag, runs[i], i, thisCalib ));
    ofs << i << "\t" << thisCalib << std::endl;
  }

  ofs.close();

  drawHistos( outputdir, rawHistos,   "rawSpectra"  , yMax );
  drawHistos( outputdir, calibHistos, "calibSpectra", yMax );


  std::cout << std::endl;
  std::cout << "-> Calibration constants saved in: " << constantsFileName << std::endl;

  std::cout << "Calibration average for BGO: " << calibAve << std::endl;


  std::string run_cef3 = "BTF_246_beam";
  TH1D* h1_cef3 = fitCeF3( outputdir, "cef3_raw", tag, run_cef3 );
  TF1* f1_cef3 = (TF1*)(h1_cef3->GetListOfFunctions()->FindObject(Form("gaussian_%s", run_cef3.c_str())));

  std::cout << "BGO/CeF3 relative calibration: " << calibAve/f1_cef3->GetParameter(1) << std::endl;
  
  return 0;

}




void drawHistos( const std::string& outputdir, std::vector<TH1D*> histos, const std::string& name, float yMax, float xMax ) {

  std::vector<int> colors;
  colors.push_back( kRed );
  colors.push_back( 38 );
  colors.push_back( 30 );
  colors.push_back( 42 );
  colors.push_back( 29 );
  colors.push_back( kBlack );
  colors.push_back( kGreen );
  colors.push_back( kBlue  );

  TCanvas* c2 = new TCanvas("c2", "", 600, 600);
  c2->cd();

  TLegend* legend = new TLegend( 0.2, 0.55, 0.5, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);

  TH2D* h2_axes = new TH2D("axes", "", 10, 0., xMax, 10, 0., 1.1*yMax );
  h2_axes->SetXTitle( "BGO ADC Channel" );
  h2_axes->SetYTitle( "Normalized to Unity" );
  h2_axes->Draw("");

  for( unsigned i=0; i<histos.size(); ++i ) {

    histos[i]->SetLineColor( colors[i] );
    histos[i]->SetLineWidth( 2 );
    histos[i]->DrawNormalized( "histo same" );
    legend->AddEntry( histos[i], Form("Channel %d", i), "L" );

  }


  TPaveText* labelTop = DrawTools::getLabelTop(); 
  labelTop->Draw("same");

  legend->Draw("same");

  c2->SaveAs( Form("%s/%s.eps", outputdir.c_str(), name.c_str()) );
  c2->SaveAs( Form("%s/%s.png", outputdir.c_str(), name.c_str()) );
  c2->SaveAs( Form("%s/%s.pdf", outputdir.c_str(), name.c_str()) );

  delete c2;
  delete h2_axes;

}







TH1D* fitSingleChannelBGO( const std::string& outputdir, const std::string& name, const std::string& tag, const std::string& runName, int iChannel, float calibConst ) {

  return fitSingleChannel( outputdir, name, tag, runName, Form("bgo_corr[%d]*%f", iChannel, calibConst), Form("Channel %d", iChannel), 50, 0., 4000. );

}

TH1D* fitCeF3( const std::string& outputdir, const std::string& name, const std::string& tag, const std::string& runName ) {

  return fitSingleChannel( outputdir, name, tag, runName, "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "CeF3", 200, 0., 16000. );

}



TH1D* fitSingleChannel( const std::string& outputdir, const std::string& name, const std::string& tag, const std::string& runName,  const std::string& varName, const std::string& plotName, int nBins, float xMin, float xMax ) {


  TFile* file = TFile::Open(Form("PosAnTrees_%s/PosAn_%s.root", tag.c_str(), runName.c_str()));
  std::cout << "-> Opened file: " << file->GetName() << std::endl;
  TTree* tree = (TTree*)file->Get("posTree");

  std::string histoName = runName;
  TH1D* h1_bgo = new TH1D(histoName.c_str(), "", nBins, xMin, xMax );

  //tree->Project( histoName.c_str(), Form("bgo_corr[%d]*%f", iChannel, calibConst), "scintFront>500. && scintFront<2000.");
  tree->Project( histoName.c_str(), varName.c_str(), "isSingleEle_scintFront" );

  TF1* f1 = new TF1( Form("gaussian_%s", runName.c_str()), "gaus(0)" );
  f1->SetRange(1800., 3200.);
  f1->SetParameter(1, 2500.);
  f1->SetParameter(2, 200.);

  h1_bgo->Fit( f1, "RQN" );

  int niter = 6;
  float nSigmaPlus  = 1.5;
  float nSigmaMinus = 0.7;

  for( unsigned iter=0; iter<niter; iter++ ) {

    float mean  = f1->GetParameter(1);
    float sigma = f1->GetParameter(2);
    float fitMin = mean - nSigmaMinus*sigma;
    float fitMax = mean + nSigmaPlus *sigma;
    f1->SetRange( fitMin, fitMax );
    if( iter==(niter-1) )
      h1_bgo->Fit( f1, "RQ+" );
    else
      h1_bgo->Fit( f1, "RQN" );
  }


  TCanvas* c1 = new TCanvas("c2", "", 600, 600);
  c1->cd();

  h1_bgo->Draw();
  f1->SetLineColor(kRed);
  f1->Draw("same");

  TPaveText* label = DrawTools::getLabelTop();
  label->Draw("same");

  float xMin_label = (plotName=="CeF3") ? 0.6 : 0.2;
  float xMax_label = (plotName=="CeF3") ? 0.9 : 0.5;
  TPaveText* labelChan = new TPaveText( xMin_label, 0.55, xMax_label, 0.9, "brNDC");
  labelChan->SetFillColor(0);
  labelChan->SetTextSize(0.038);
  labelChan->AddText( plotName.c_str() );
  labelChan->Draw("same");

  c1->SaveAs( Form("%s/fit_%s_%s.eps", outputdir.c_str(), name.c_str(), runName.c_str()) );
  c1->SaveAs( Form("%s/fit_%s_%s.png", outputdir.c_str(), name.c_str(), runName.c_str()) );
  c1->SaveAs( Form("%s/fit_%s_%s.pdf", outputdir.c_str(), name.c_str(), runName.c_str()) );

  delete c1;


  return h1_bgo;

}



float sumVector( std::vector<float> v ) {

  float sum=0.;
  for( unsigned i=0; i<v.size(); ++i ) sum += v[i];

  return sum;

}


