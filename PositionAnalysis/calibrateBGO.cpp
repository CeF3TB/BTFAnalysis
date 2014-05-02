#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include <vector>
#include <string>



TH1D* fitSingleChannel( const std::string& outputdir, const std::string& runName, int iChannel );



int main() {


  std::vector<std::string> runs;
  runs.push_back( "BTF_229_20140501-155939_beam" ); // 0
  runs.push_back( "BTF_231_20140501-163742_beam" ); // 1
  runs.push_back( "BTF_233_20140501-171741_beam" ); // 2
  runs.push_back( "BTF_237_20140501-183950_beam" ); // 3
  runs.push_back( "BTF_235_20140501-175948_beam" ); // 4
  runs.push_back( "BTF_239_20140501-191908_beam" ); // 5
  runs.push_back( "BTF_241_20140501-200053_beam" ); // 6
  runs.push_back( "BTF_243_20140501-203838_beam" ); // 7


  std::string outputdir = "BGOCalibration/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system(mkdir_command.c_str());

  std::vector<TH1D*> calib;

  for( unsigned i=0; i<runs.size(); ++i )
     calib.push_back(fitSingleChannel( outputdir, runs[i], i ) );




  return 0;

}




TH1D* fitSingleChannel( const std::string& outputdir, const std::string& runName, int iChannel ) {


  TFile* file = TFile::Open(Form("PosAn_%s.root", runName.c_str()));
  TTree* tree = (TTree*)file->Get("tree_passedEvents");

  std::string histoName = "bgo_"+runName;
  TH1D* h1_bgo = new TH1D(histoName.c_str(), "", 100, 0., 4000.);

  tree->Project( histoName.c_str(), Form("bgo_corr[%d]", iChannel), "scintFront>500. && scintFront<2000.");

  TF1* f1 = new TF1( "gaussian", "gaus(0)" );
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


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  h1_bgo->Draw();
  f1->SetLineColor(kRed);
  f1->Draw("same");

  c1->SaveAs( Form("%s/fit_%s.eps", outputdir.c_str(), runName.c_str()) );
  c1->SaveAs( Form("%s/fit_%s.png", outputdir.c_str(), runName.c_str()) );

  delete c1;


  return h1_bgo;

}
