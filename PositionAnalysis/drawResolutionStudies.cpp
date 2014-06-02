#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"


#include "interface/DrawTools.h"
#include "interface/FitTools.h"



struct ResoStruct {

  float resp;
  float resp_error;

  float reso;
  float reso_error;

  float Sres;
  float Sres_error;

};


struct DiagonalScanStruct {

 DiagonalScanStruct( float d, TTree* t_data, TTree* t_mc ) {
   diag = d;
   tree_mc = t_mc;
   tree_data = t_data;
 }

 float diag;
 TTree* tree_data;
 TTree* tree_mc;

};


ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, float LYSF[], const std::string& name );
std::string getVarName( float LYSF[] );
ResoStruct getRespAndReso( TF1* f1, float energyErrorPercent );
float getRatioError( float num, float denom, float numErr, float denomErr );
ResoStruct getRespResoFromHisto( TH1D* h1 );




int main() {

  DrawTools::setStyle();

  TFile* file_mc = TFile::Open("EEShash_491MeV_10000ev_smear.root");
  //TFile* file_mc = TFile::Open("EEShash_491MeV_10000ev.root");

  TFile* file_data = TFile::Open("analysisTrees_V00/Reco_BTF_259_20140502-012847_beam.root");

  TTree* tree_data = (TTree*)file_data->Get("recoTree");
  TTree* tree_mc = (TTree*)file_mc->Get("EEShash");

  std::string outputdir = "ResolutionStudiesPlots";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );
  
  TF1* f1_data = FitTools::fitSingleElectronPeak( outputdir, "data", tree_data );

  ResoStruct rs_data = getRespAndReso(f1_data, 0.01);



  float LYSF_ideal[] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
  float LYSF_real [] = {0.85, 0.94, 0.95, 0.98, 1.00, 1.02, 1.05, 1.05, 1.05, 0.74};
  float LYSF_hole [] = {0.85, 0.94, 0.95, 0.98, 1.00, 0.00, 1.05, 1.05, 1.05, 0.74};

  ResoStruct rs_ideal = getResponseResolutionMC( outputdir, tree_mc, LYSF_ideal, "ideal" );
  ResoStruct rs_real  = getResponseResolutionMC( outputdir, tree_mc, LYSF_real, "real" );
  ResoStruct rs_hole  = getResponseResolutionMC( outputdir, tree_mc, LYSF_hole, "hole" );
  
  std::cout << std::endl;
  std::cout << "* DATA: " << std::endl;
  std::cout << "reso: " << rs_data.reso << " +/- " << rs_data.reso_error << std::endl;
  std::cout << "Sres: " << rs_data.Sres << " +/- " << rs_data.Sres_error << std::endl;
  
  std::cout << "MC ideal: " << rs_ideal.Sres << " +/- " << rs_ideal.Sres_error << std::endl;
  std::cout << "MC real: " << rs_real.Sres << " +/- " << rs_real.Sres_error << std::endl;
  std::cout << "MC hole: " << rs_hole.Sres << " +/- " << rs_hole.Sres_error << std::endl;


  TFile* file_mc_3x3y = TFile::Open("EEShash_491MeV_10000ev_smear_3x3y.root");
  TFile* file_mc_6x6y = TFile::Open("EEShash_491MeV_10000ev_smear_6x6y.root");
  TFile* file_mc_9x9y = TFile::Open("EEShash_491MeV_10000ev_smear_9x9y.root");
  TFile* file_mc_11p3x11p3y = TFile::Open("EEShash_491MeV_10000ev_smear_11p3x11p3y.root");

  TTree* tree_mc_3x3y = (TTree*)file_mc_3x3y->Get("EEShash");
  TTree* tree_mc_6x6y = (TTree*)file_mc_6x6y->Get("EEShash");
  TTree* tree_mc_9x9y = (TTree*)file_mc_9x9y->Get("EEShash");
  TTree* tree_mc_11p3x11p3y = (TTree*)file_mc_11p3x11p3y->Get("EEShash");

  TFile* file_data_3x3y       = TFile::Open("PosAnTrees_V00/PosAn_BTF_141_beam.root");
  TFile* file_data_6x6y       = TFile::Open("PosAnTrees_V00/PosAn_BTF_143_beam.root");
  TFile* file_data_9x9y       = TFile::Open("PosAnTrees_V00/PosAn_BTF_167_beam.root");
  TFile* file_data_11p3x11p3y = TFile::Open("PosAnTrees_V00/PosAn_BTF_219_beam.root");

  TTree* tree_data_3x3y = (TTree*)file_data_3x3y->Get("posTree");
  TTree* tree_data_6x6y = (TTree*)file_data_6x6y->Get("posTree");
  TTree* tree_data_9x9y = (TTree*)file_data_9x9y->Get("posTree");
  TTree* tree_data_11p3x11p3y = (TTree*)file_data_11p3x11p3y->Get("posTree");

  //TFile* file_data_3x3y       = TFile::Open("AnalysisTrees_V00/Reco_BTF_141_20140430-183508_beam.root");
  //TFile* file_data_6x6y       = TFile::Open("AnalysisTrees_V00/Reco_BTF_143_20140430-191455_beam.root");
  //TFile* file_data_9x9y       = TFile::Open("AnalysisTrees_V00/Reco_BTF_167_20140430-210839_beam.root");
  //TFile* file_data_11p3x11p3y = TFile::Open("AnalysisTrees_V00/Reco_BTF_219_20140501-092151_beam.root");

  //TTree* tree_data_3x3y = (TTree*)file_data_3x3y->Get("recoTree");
  //TTree* tree_data_6x6y = (TTree*)file_data_6x6y->Get("recoTree");
  //TTree* tree_data_9x9y = (TTree*)file_data_9x9y->Get("recoTree");
  //TTree* tree_data_11p3x11p3y = (TTree*)file_data_11p3x11p3y->Get("recoTree");

  std::vector<DiagonalScanStruct> dss;
  dss.push_back( DiagonalScanStruct(0., tree_data, tree_mc) );
  dss.push_back( DiagonalScanStruct(3., tree_data_3x3y, tree_mc_3x3y) );
  dss.push_back( DiagonalScanStruct(6., tree_data_6x6y, tree_mc_6x6y) );
  dss.push_back( DiagonalScanStruct(9., tree_data_9x9y, tree_mc_9x9y) );
  dss.push_back( DiagonalScanStruct(11.3, tree_data_11p3x11p3y, tree_mc_11p3x11p3y) );


  std::string fullVarName_mc = getVarName(LYSF_hole);

  ResoStruct rs_ref_mc;
  ResoStruct rs_ref_data;
  TGraphErrors* gr_RespVsDiag_data = new TGraphErrors(0);
  TGraphErrors* gr_ResoVsDiag_data = new TGraphErrors(0);
  TGraphErrors* gr_RespVsDiag_mc   = new TGraphErrors(0);
  TGraphErrors* gr_ResoVsDiag_mc   = new TGraphErrors(0);

  for( unsigned i=0; i<dss.size(); ++i ) {

    std::string histoName_data(Form("data_%.0f", dss[i].diag));
    TH1D* h1_data = new TH1D( histoName_data.c_str(), "", 200, 0., 5000. );
    dss[i].tree_data->Project( histoName_data.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "isSingleEle_scintFront" );

    std::string histoName_mc(Form("mc_%.0f", dss[i].diag));
    TH1D* h1_mc = new TH1D( histoName_mc.c_str(), "", 100, 0., 500. );
    dss[i].tree_mc->Project( histoName_mc.c_str(), fullVarName_mc.c_str() );

    if( i==0 ) {

      rs_ref_mc   = getRespResoFromHisto( h1_mc );
      rs_ref_data = getRespResoFromHisto( h1_data );

    }

    ResoStruct rs_mc = getRespResoFromHisto( h1_mc );
    ResoStruct rs_data = getRespResoFromHisto( h1_data );

    gr_RespVsDiag_data->SetPoint( i, dss[i].diag, rs_data.resp/rs_ref_data.resp );
    gr_ResoVsDiag_data->SetPoint( i, dss[i].diag, rs_data.reso/rs_ref_data.reso );
    gr_RespVsDiag_mc  ->SetPoint( i, dss[i].diag, rs_mc.resp/rs_ref_mc.resp );
    gr_ResoVsDiag_mc  ->SetPoint( i, dss[i].diag, rs_mc.reso/rs_ref_mc.reso );
    
    gr_RespVsDiag_data->SetPointError( i, 0., getRatioError( rs_data.resp, rs_ref_data.resp, rs_data.resp_error, rs_ref_data.resp_error) );
    gr_ResoVsDiag_data->SetPointError( i, 0., getRatioError( rs_data.reso, rs_ref_data.reso, rs_data.reso_error, rs_ref_data.reso_error) );
    gr_RespVsDiag_mc  ->SetPointError( i, 0., getRatioError( rs_mc  .resp, rs_ref_mc  .resp, rs_mc  .resp_error, rs_ref_mc  .resp_error) );
    gr_ResoVsDiag_mc  ->SetPointError( i, 0., getRatioError( rs_mc  .reso, rs_ref_mc  .reso, rs_mc  .reso_error, rs_ref_mc  .reso_error) );
    
  }


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  float xMin = -1.;
  float xMax = 13.;
  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.1);
  h2_axes->SetXTitle( "Diagonal Distance From Center [mm]");
  h2_axes->SetYTitle( "Response Ratio" );

  h2_axes->Draw();

  TLine* line_one = new TLine( xMin, 1., xMax, 1. );
  line_one->Draw("same");

  gr_RespVsDiag_data->SetMarkerSize(1.6);
  gr_RespVsDiag_mc->SetMarkerSize(1.6);

  gr_RespVsDiag_data->SetMarkerStyle(20);
  gr_RespVsDiag_mc->SetMarkerStyle(24);

  gr_RespVsDiag_data->SetMarkerColor(46);
  gr_RespVsDiag_mc->SetMarkerColor(kBlack);

  gr_RespVsDiag_data->Draw("p same");
  gr_RespVsDiag_mc->Draw("p same");

  TLegend* legend = new TLegend( 0.23, 0.27, 0.5, 0.47 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( gr_RespVsDiag_data, "Data", "P" );
  legend->AddEntry( gr_RespVsDiag_mc, "Geant4", "P" );
  legend->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/resp_vs_diag.eps", outputdir.c_str()) );
  c1->SaveAs( Form("%s/resp_vs_diag.png", outputdir.c_str()) );
  c1->SaveAs( Form("%s/resp_vs_diag.pdf", outputdir.c_str()) );

  return 0;

}




ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, float LYSF[], const std::string& name ) {

  std::string fullVarName = getVarName(LYSF);

  TH1D* h1 = new TH1D( name.c_str(), "", 500, 0., 500. );

  tree->Project( name.c_str(), fullVarName.c_str() );

  TF1* f1 = new TF1( Form("f1_%s", name.c_str()), "gaus", 100., 500.);

  FitTools::doSingleFit( h1, f1, outputdir, name );


  ResoStruct rs = getRespAndReso( f1, 0. );

  return rs;

}


std::string getVarName( float LYSF[] ) {

  std::string fullVarName = "";
  for( unsigned i=0; i<10; ++i ) {
    std::string plusSign = (i==0) ? "" : " + ";
    std::string thisPiece(Form("%s%f*Eact_%d", plusSign.c_str(), LYSF[i], i));
    fullVarName += thisPiece;
  }

  return fullVarName;

}


ResoStruct getRespAndReso( TF1* f1, float energyErrorPercent ) {

  float energy = 491.;
  float energyErr = energyErrorPercent*energy;

  float mean = f1->GetParameter(1);
  float meanErr = f1->GetParError(1);

  float rms = f1->GetParameter(2);
  float rmsErr = f1->GetParError(2);

  float reso = 100.* rms/mean; //in percent
  float resoErr = 100.* getRatioError( mean, rms, meanErr, rmsErr );

  float energyGeV = energy/1000.;
  float energyGeVErr = energyErr/1000.;
  float Sres = reso*sqrt(energyGeV);
  float Sres_error = sqrt( resoErr*resoErr*energyGeV + reso*reso*0.5*0.5*energyGeVErr*energyGeVErr/energyGeV );

  ResoStruct rs;
  rs.resp = mean;
  rs.resp_error = meanErr;
  rs.reso = reso;
  rs.reso_error = resoErr;
  rs.Sres = Sres;
  rs.Sres_error = Sres_error;

  return rs;

}


float getRatioError( float num, float denom, float numErr, float denomErr ) {

  return sqrt( denomErr*denomErr/(num*num) + denom*denom*numErr*numErr/(num*num*num*num) );

}

ResoStruct getRespResoFromHisto( TH1D* h1 ) {

  float mean = h1->GetMean();
  float meanErr = h1->GetMeanError();

  float rms = h1->GetRMS();
  float rmsErr = h1->GetRMS();

  float reso = 100.* rms/mean; //in percent
  float resoErr = 100.* getRatioError(mean, rms, meanErr, rmsErr);

  ResoStruct rs;
  rs.resp = mean;
  rs.resp_error = meanErr;
  rs.reso = reso;
  rs.reso_error = resoErr;

  return rs;

}
