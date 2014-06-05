#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TVector2.h"
#include "TProfile.h"
#include "TLegend.h"

#include "interface/DrawTools.h"
#include "interface/RunHelper.h"








int main() {


  DrawTools::setStyle();

  
  TFile* file = TFile::Open("TMVA/trainingFile.root");
  TTree* tree = (TTree*)file->Get("tree_passedEvents");



  unsigned int runNumber;
  tree->SetBranchAddress( "runNumber", &runNumber );
  bool bgo_corr_ok;
  tree->SetBranchAddress( "bgo_corr_ok", &bgo_corr_ok );
  float scintFront;
  tree->SetBranchAddress( "scintFront", &scintFront );
  int nHodoClustersX;
  tree->SetBranchAddress( "nHodoClustersX", &nHodoClustersX );
  int nHodoClustersY;
  tree->SetBranchAddress( "nHodoClustersY", &nHodoClustersY );
  int nFibres_hodoClustX[8];
  tree->SetBranchAddress( "nFibres_hodoClustX", nFibres_hodoClustX );
  float pos_hodoClustX[8];
  tree->SetBranchAddress( "pos_hodoClustX", pos_hodoClustX );
  int nFibres_hodoClustY[8];
  tree->SetBranchAddress( "nFibres_hodoClustY", nFibres_hodoClustY );
  float pos_hodoClustY[8];
  tree->SetBranchAddress( "pos_hodoClustY", pos_hodoClustY );
  float xBeam;
  tree->SetBranchAddress( "xBeam", &xBeam );
  float yBeam;
  tree->SetBranchAddress( "yBeam", &yBeam );
  float xCalo;
  tree->SetBranchAddress( "xPos_calo", &xCalo );
  float yCalo;
  tree->SetBranchAddress( "yPos_calo", &yCalo );
  float cef3[4];
  tree->SetBranchAddress( "cef3_corr", cef3 );


  TFile* outfile = TFile::Open( "CeF3PositionCalibration.root", "recreate" );
  outfile->cd();

  TProfile* hp_diag02_vs_dBeam = new TProfile( "diag02_vs_dBeam", "", 50, -20., 20.);
  TProfile* hp_diag13_vs_dBeam = new TProfile( "diag13_vs_dBeam", "", 50, -20., 20.);

  TProfile* hp_diag02_vs_dCalo = new TProfile( "diag02_vs_dCalo", "", 50, -20., 20.);
  TProfile* hp_diag13_vs_dCalo = new TProfile( "diag13_vs_dCalo", "", 50, -20., 20.);


  int nentries = tree->GetEntries();

  int cachedRun = 0;
  TGraph* gr_run_xyMap = new TGraph(0);


  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    if( iEntry % 100000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;

    tree->GetEntry(iEntry); 

    if( runNumber!=cachedRun ) {
      cachedRun=runNumber;
      gr_run_xyMap->SetPoint( gr_run_xyMap->GetN(), xBeam, yBeam );
    }

    bool isScintFrontSingleEle = ( runNumber<100 ) ? (scintFront>120. && scintFront<700.) : (scintFront>500. && scintFront<2000.);

    if( !bgo_corr_ok ) continue;
    if( !isScintFrontSingleEle ) continue;

    float r02 = cef3[0]/cef3[2];
    float r13 = cef3[1]/cef3[3];

    TVector2 v_Beam(xBeam, yBeam);
    TVector2 d_Beam = v_Beam.Rotate( -3.14159/4. );
    float diag13_Beam = d_Beam.X();
    float diag02_Beam = d_Beam.Y();

    TVector2 v_Calo(xCalo, yCalo);
    TVector2 d_Calo = v_Calo.Rotate( -3.14159/4. );
    float diag13_Calo = d_Calo.X();
    float diag02_Calo = d_Calo.Y();

    hp_diag02_vs_dBeam->Fill(diag02_Beam, r02 );
    hp_diag13_vs_dBeam->Fill(diag13_Beam, r13 );
    hp_diag02_vs_dCalo->Fill(diag02_Calo, r02 );
    hp_diag13_vs_dCalo->Fill(diag13_Calo, r13 );


  } // for entries


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  hp_diag02_vs_dBeam->SetLineColor(kBlack);
  hp_diag13_vs_dBeam->SetLineColor(kRed);
  hp_diag02_vs_dCalo->SetLineColor(kBlue);
  hp_diag13_vs_dCalo->SetLineColor(kGreen);


  TF1* f1_d02 = new TF1("f1_d02", "pol5", -17., 17.);
  f1_d02->SetLineColor( kGray);
  f1_d02->FixParameter( 2, 0. );
  f1_d02->FixParameter( 5, 0. );
  //f1_d02->SetParLimits( 2, 0., 0.1 );
  //f1_d02->SetParameter( 5, 0.05 );
  //f1_d02->SetParLimits( 5, 0., 0.1 );
  hp_diag02_vs_dCalo->Fit( f1_d02, "RB" );

  TF1* f1_d13 = new TF1("f1_d13", "pol5", -17., 17.);
  f1_d13->SetLineColor(kRed);
  hp_diag13_vs_dCalo->Fit( f1_d13, "R+" );




  TH2D* h2_axes = new TH2D("axes", "", 10, -20., 20., 10, 0., 3.);
  h2_axes->Draw();

  TLine* lineone = new TLine( -20., 1., 20., 1. );

  TLegend* legend = new TLegend( 0.2, 0.55, 0.5, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( f1_d02, "02", "L" );
  legend->AddEntry( f1_d13, "13", "L" );
  legend->Draw("same");

  hp_diag02_vs_dCalo->Draw("same");
  lineone->Draw("same");
  hp_diag13_vs_dCalo->Draw("same");
  hp_diag02_vs_dBeam->Draw("same");
  hp_diag13_vs_dBeam->Draw("same");

  c1->SaveAs("prova.eps"); 


  outfile->cd();

  gr_run_xyMap->Write();

  hp_diag02_vs_dBeam->Write();
  hp_diag13_vs_dBeam->Write();

  hp_diag02_vs_dCalo->Write();
  hp_diag13_vs_dCalo->Write(); 

  outfile->Close();

  ofstream ofs02("diag02_parameters.txt");
  ofs02 << f1_d02->GetParameter(0) << std::endl;
  ofs02 << f1_d02->GetParameter(1) << std::endl;
  ofs02 << f1_d02->GetParameter(2) << std::endl;
  ofs02 << f1_d02->GetParameter(3) << std::endl;
  ofs02 << f1_d02->GetParameter(4) << std::endl;
  ofs02 << f1_d02->GetParameter(5) << std::endl;
  ofs02.close();
  std::cout << "Diag02 fit parameters saved in: diag02_parameters.txt" << std::endl;


  ofstream ofs13("diag13_parameters.txt");
  ofs13 << f1_d13->GetParameter(0) << std::endl;
  ofs13 << f1_d13->GetParameter(1) << std::endl;
  ofs13 << f1_d13->GetParameter(2) << std::endl;
  ofs13 << f1_d13->GetParameter(3) << std::endl;
  ofs13 << f1_d13->GetParameter(4) << std::endl;
  ofs13 << f1_d13->GetParameter(5) << std::endl;
  ofs13.close();
  std::cout << "Diag13 fit parameters saved in: diag13_parameters.txt" << std::endl;

  return 0;

}

