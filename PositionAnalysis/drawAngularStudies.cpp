#include <stdio.h>
#include <cmath>
#include <math.h>

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


struct LateralScanStruct {

 LateralScanStruct( float d, TTree* t_data, TTree* t_mc ) {
   offset = d;
   tree_mc = t_mc;
   tree_data = t_data;
 }

 float offset;
 TTree* tree_data;
 TTree* tree_mc;

};


ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, float LYSF[], const std::string& name );
std::string getVarName( float LYSF[] );
ResoStruct getRespAndReso( TF1* f1, float energyErrorPercent );
float getRatioError( float num, float denom, float numErr, float denomErr );
ResoStruct getRespResoFromHisto( TH1D* h1 );
ResoStruct addPhotoStatistics( ResoStruct rs );




int main( int argc, char* argv[]) {


  std::string tag = "V00";
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }

  DrawTools::setStyle();


  std::string outputdir = "ResolutionStudiesPlots_"+tag;
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );


  float LYSF_real [] = {0.85, 0.94, 0.95, 0.98, 1.00, 1.02, 1.05, 1.05, 1.05, 0.74};

  std::vector< float > angles;
  angles.push_back( 0. );
  angles.push_back( 0.5);
  angles.push_back( 1. );
  angles.push_back( 1.5);
  angles.push_back( 2. );
  angles.push_back( 2.5);
  angles.push_back( 3. );
  angles.push_back( 4. );
  angles.push_back( 5. );
  angles.push_back( 6. );
  angles.push_back( 7. );
  angles.push_back( 8. );
  angles.push_back( 10. );


  TGraphErrors* gr_reso = new TGraphErrors(0);
  TGraphErrors* gr_reso_ps = new TGraphErrors(0);

  for( unsigned i=0; i< angles.size(); ++i ) {

    double intPart;
    double fracPart = modf( angles[i], &intPart );

    char angleText[100];
    if( fracPart == 0. ) 
      sprintf( angleText, "%.0fdeg", floor(angles[i]) );
    else
      sprintf( angleText, "%.0fp%.0fdeg", floor(angles[i]), fracPart*10. );


    TFile* file = TFile::Open( Form("EEShash_%s.root", angleText) );
    if( file==0 ) {
      std::cout << "WARNING: didn't find file for " << angleText << std::endl;
      continue;
    }
    TTree* tree = (TTree*)file->Get("EEShash");
    if( tree==0 ) {
      std::cout << "WARNING: didn't find tree in file: " << file->GetName() << std::endl;
      continue;
    }
     
    ResoStruct rs = getResponseResolutionMC( outputdir, tree, LYSF_real, angleText );
    ResoStruct rs_ps = addPhotoStatistics( rs );

    gr_reso->SetPoint(i, angles[i], rs.reso);
    gr_reso->SetPointError(i, 0., rs.reso_error);

    gr_reso_ps->SetPoint(i, angles[i], rs_ps.reso);
    gr_reso_ps->SetPointError(i, 0., rs_ps.reso_error);

  }



  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, -0.5, 10.5, 10, 10., 30.);
  h2_axes->SetXTitle("Beam Impact Angle [deg]");
  h2_axes->SetYTitle("Single-Channel Energy Resolution [%]");
  h2_axes->Draw("same");

  TLegend* legend = new TLegend( 0.2, 0.75, 0.45, 0.9 );
  legend->SetTextSize(0.035);
  legend->SetFillColor(0);
  legend->AddEntry( gr_reso, "Pure Sampling", "P" );
  legend->AddEntry( gr_reso_ps, "+ Photostatistics", "P" );
  legend->Draw("same");

  gr_reso->SetMarkerStyle(24);
  gr_reso->SetMarkerColor(kBlack);
  gr_reso->SetMarkerSize(2.);

  gr_reso_ps->SetMarkerStyle(20);
  gr_reso_ps->SetMarkerColor(kBlack);
  gr_reso_ps->SetMarkerSize(2.);

  gr_reso->Draw("p same");
  gr_reso_ps->Draw("p same");

  TPaveText* labelTop = DrawTools::getLabelTop("MC Simulation, 491 MeV Electrons");
  labelTop->Draw("same");

  c1->SaveAs("reso_vs_angle.eps");
  c1->SaveAs("reso_vs_angle.png");


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
  float resoErr = 100.* getRatioError( rms, mean, meanErr, rmsErr );

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

  return sqrt( numErr*numErr/(denom*denom) + denomErr*denomErr*num*num/(denom*denom*denom*denom) );

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





ResoStruct addPhotoStatistics( ResoStruct rs ) {

  // MC response is already in MeV, 0.49 is the number of p.e. per MeV
  // RM = 90% of energy, so assume 90% of energy (0.9*491) is deposited in central channel
  float nADC = rs.resp/185.*3200.;
  float nPhotoElectrons = nADC/35.;

  float poissonError = 100./sqrt( nPhotoElectrons ); // in percent

  rs.reso = sqrt( rs.reso*rs.reso + poissonError*poissonError );
  rs.Sres = rs.reso*sqrt(0.491); // E = 491 MeV

  return rs;

} 
