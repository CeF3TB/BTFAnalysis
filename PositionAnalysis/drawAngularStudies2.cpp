#include <stdio.h>
#include <cmath>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "TApplication.h"

#include "interface/DrawTools.h"
#include "interface/FitTools.h"



struct ResoStruct {

  float resp;
  float resp_error;

  float reso;
  float reso_error;

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


ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, const std::string& name );
ResoStruct getResponseResolutionMCdata( const std::string& outputdir, TTree* tree, const std::string& name );
std::string getVarName( float LYSF[] );
ResoStruct getRespAndReso( TF1* f1, float energyErrorPercent );
float getRatioError( float num, float denom, float numErr, float denomErr );
ResoStruct getRespResoFromHisto( TH1D* h1 );
ResoStruct addPhotoStatistics( ResoStruct rs );




int main( int argc, char* argv[]) {

  TApplication* a = new TApplication("a", 0, 0);

  std::string tag = "V00";
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }

  DrawTools::setStyle();
  gStyle->SetOptFit(1);

  std::string outputdir = "AngularResolutionStudiesPlots";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );


//Data for the tilte BOX (beam centred at hodoscope):
  std::vector<std::string> runName;
  std::vector<float> runAngle;

  // runName.push_back("BTF_259_20140502-012847_beam");
  // runAngle.push_back(0.0);
  // runName.push_back("BTF_321_20140503-053105_beam"); //dont have it atm
  //  runAngle.push_back(0.0);
  runName.push_back("BTF_334_20140503-103227_beam");
  runAngle.push_back(1.);
  runName.push_back("BTF_335_20140503-105047_beam");
  runAngle.push_back(1.0);
  runName.push_back("BTF_341_20140503-115519_beam");
  runAngle.push_back(3.0);
  runName.push_back("BTF_344_20140503-124045_beam");
  runAngle.push_back(5.0);




//Data for the tilte BEAM (beam centred at ?):
  std::vector<std::string> runName2;
  std::vector<float> runAngle2;

  // runName.push_back("BTF_259_20140502-012847_beam");
  // runAngle.push_back(0.0);
  runName2.push_back("BTF_321_20140503-053105_beam"); 
  runAngle2.push_back(0.0);
  runName2.push_back("BTF_323_20140503-070141_beam");
  runAngle2.push_back(0.45);
  runName2.push_back("BTF_324_20140503-074214_beam");
  runAngle2.push_back(0.9);
  runName2.push_back("BTF_329_20140503-085651_beam");
  runAngle2.push_back(1.12);






  std::vector< float > angles;
  std::vector<std::string> anglesText; 
  angles.push_back(0.0);
  angles.push_back( 0.5);
  angles.push_back( 1.0 ); 
  angles.push_back( 1.5);
  angles.push_back( 2.0 );
  angles.push_back( 2.5);
  angles.push_back( 3. );
  angles.push_back( 3.5 );
  angles.push_back( 4. );
  angles.push_back( 4.5 );
  angles.push_back( 5. );
  angles.push_back( 5.5 );
  angles.push_back( 6. );
  angles.push_back( 7. );
  angles.push_back( 8. );
  angles.push_back( 9. );
  angles.push_back( 10.0 );

  anglesText.push_back("0.0");
  anglesText.push_back("0.5");
  anglesText.push_back("1.0"); 
  anglesText.push_back("1.5");
  anglesText.push_back("2.0"); 
  anglesText.push_back("2.5");
  anglesText.push_back("3.0");
  anglesText.push_back("3.5");
  anglesText.push_back("4.0");
  anglesText.push_back("4.5");
  anglesText.push_back("5.0");
  anglesText.push_back("5.5");
  anglesText.push_back( "6.0" );
  anglesText.push_back( "7.0" );
  anglesText.push_back( "8.0" );
  anglesText.push_back( "9.0" );
  anglesText.push_back("10.0"); 

/*
  angles.push_back( 2.5);
  angles.push_back( 3. );
  angles.push_back( 4. );
  angles.push_back( 5. );
  angles.push_back( 6. );
  angles.push_back( 7. );
  angles.push_back( 8. );
  angles.push_back( 10. );
*/

  TGraphErrors* gr_reso = new TGraphErrors(0);
  TGraphErrors* gr_reso_ps = new TGraphErrors(0);

  TGraphErrors* gr_reso_data = new TGraphErrors(0);
  TGraphErrors* gr_reso_data2 = new TGraphErrors(0);



  for( unsigned i=0; i< runName.size(); ++i ) {

    TFile* file = TFile::Open( Form("analysisTrees_V03/Reco_%s.root", runName[i].c_str()  ) );
    if( file==0 ) {
      std::cout << "WARNING: didn't find file for " << runName[i] << std::endl;
      continue;
    }
    TTree* tree = (TTree*)file->Get("recoTree");
    if( tree==0 ) {
      std::cout << "WARNING: didn't find tree in file: " << file->GetName() << std::endl;
      continue;
    }
     
    std::string name = Form("%0f",runAngle[i]);
    ResoStruct rs = getResponseResolutionMCdata( outputdir, tree, name.c_str() );
    ResoStruct rs_ps = addPhotoStatistics( rs );

    gr_reso_data->SetPoint(i, runAngle[i], rs.reso);
    gr_reso_data->SetPointError(i, 0., rs.reso_error);
  }


  for( unsigned i=0; i< runName2.size(); ++i ) {

    TFile* file = TFile::Open( Form("analysisTrees_V03/Reco_%s.root", runName2[i].c_str()  ) );
    if( file==0 ) {
      std::cout << "WARNING: didn't find file for " << runName2[i] << std::endl;
      continue;
    }
    TTree* tree = (TTree*)file->Get("recoTree");
    if( tree==0 ) {
      std::cout << "WARNING: didn't find tree in file: " << file->GetName() << std::endl;
      continue;
    }
     
    std::string name = Form("%2f",runAngle2[i]);
    ResoStruct rs = getResponseResolutionMCdata( outputdir, tree, name.c_str() );
    ResoStruct rs_ps = addPhotoStatistics( rs );

    gr_reso_data2->SetPoint(i, runAngle2[i], rs.reso);
    gr_reso_data2->SetPointError(i, 0., rs.reso_error);
  }



  for( unsigned i=0; i< angles.size(); ++i ) {

    double intPart;
    double fracPart = modf( angles[i], &intPart );

    char angleText[100];
    //    if( fracPart == 0. ) 
    //  sprintf( angleText, "%.0f", floor(angles[i]) );
    //else
      sprintf( angleText, "%.0f.%.0f", floor(angles[i]), fracPart*10. );


TFile* file = TFile::Open( Form("OriginalSimulationData/angularStudiesBTF/Reco_Simulation%s.root", angleText ) );
    if( file==0 ) {
      std::cout << "WARNING: didn't find file for " << angleText << std::endl;
      continue;
    }
    TTree* tree = (TTree*)file->Get("recoTree");
    if( tree==0 ) {
      std::cout << "WARNING: didn't find tree in file: " << file->GetName() << std::endl;
      continue;
    }
     
    ResoStruct rs = getResponseResolutionMC( outputdir, tree, angleText );
    ResoStruct rs_ps = addPhotoStatistics( rs );

    gr_reso->SetPoint(i, angles[i], rs.reso);
    gr_reso->SetPointError(i, 0., rs.reso_error);

    gr_reso_ps->SetPoint(i, angles[i], rs_ps.reso);
    gr_reso_ps->SetPointError(i, 0., rs_ps.reso_error);

  }



  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, -0.5, 10.5, 10, 10., 40.);
  h2_axes->SetXTitle("Beam Impact Angle [deg]");
  h2_axes->SetYTitle("CeF_{3} Energy Resolution [%]");
  h2_axes->Draw("same");

  TLegend* legend = new TLegend( 0.2, 0.75, 0.45, 0.9 );
  legend->SetTextSize(0.035);
  legend->SetFillColor(0);
  legend->AddEntry( gr_reso_data, "Data Tilted Box", "P" );
  legend->AddEntry( gr_reso_data2, "Data Tilted Beam", "P" );
  legend->AddEntry( gr_reso_ps, "MC Centred on Hodo", "P" );
  legend->Draw("same");

//  gr_reso->SetMarkerStyle(24);
//  gr_reso->SetMarkerColor(kBlack);
//  gr_reso->SetMarkerSize(2.);

  gr_reso_data->SetMarkerStyle(20);
  gr_reso_data->SetMarkerColor(46);
  gr_reso_data->SetMarkerSize(1.6);
  gr_reso_data->Draw("p same");

  gr_reso_data2->SetMarkerStyle(20);
  gr_reso_data2->SetMarkerColor(38);
  gr_reso_data2->SetMarkerSize(1.6);
  gr_reso_data2->Draw("p same");

  gr_reso_ps->SetMarkerStyle(24);
  gr_reso_ps->SetMarkerColor(kBlack);
  gr_reso_ps->SetMarkerSize(1.4);
  gr_reso_ps->Draw("p same");




  //  gr_reso->Draw("p same");


  TPaveText* labelTop = DrawTools::getLabelTop("491 MeV Electron Beam");
  labelTop->Draw("same");

// c1->SaveAs("reso_vs_angle.eps");
// c1->SaveAs("reso_vs_angle.png");
c1->SaveAs(Form("%s/reso_vs_angle.pdf",outputdir.c_str()) );

  return 0;

}




ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree,  const std::string& name ) {



  TH1D* h1 = new TH1D( name.c_str(), "", 500, 0., 500. );

tree->Project( name.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]");

  TF1* f1 = new TF1( Form("f1_%s", name.c_str()), "gaus", 100., 500.);

  FitTools::doSingleFit( h1, f1, outputdir, name );


  ResoStruct rs = getRespAndReso( f1, 0. );

  return rs;

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

  ResoStruct rs;
  rs.resp = mean;
  rs.resp_error = meanErr;
  rs.reso = reso;
  rs.reso_error = resoErr;


  return rs;

}


float getRatioError( float num, float denom, float numErr, float denomErr ) {

  return sqrt( numErr*numErr/(denom*denom) + denomErr*denomErr*num*num/(denom*denom*denom*denom) );

}




ResoStruct addPhotoStatistics( ResoStruct rs ) {

  // MC response is already in MeV, 0.49 is the number of p.e. per MeV
  // RM = 90% of energy, so assume 90% of energy (0.9*491) is deposited in central channel
  //  float nADC = rs.resp/185.*3200.;
  float  nADC = rs.resp/169.445*3224.36;
  float nPhotoElectrons = nADC/27.3;
  //  float nPhotoElectrons = nADC/35.;

  float poissonError = 100./sqrt( nPhotoElectrons ); // in percent

  float resoUnsmeared = rs.reso;
  
  float resoSmeared =  sqrt( rs.reso*rs.reso + poissonError*poissonError );
  
  rs.reso = resoSmeared;

  rs.reso_error = rs.reso_error * resoSmeared / resoUnsmeared; 

  return rs;

} 




ResoStruct getResponseResolutionMCdata( const std::string& outputdir, TTree* tree,  const std::string& name ) {



  TH1D* h1 = new TH1D( name.c_str(), "", 500, 0., 6000. );

  tree->Project( name.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]");

  TF1* f1 = new TF1( Form("f1_%s", name.c_str()), "gaus", 100., 6000.);

  FitTools::doSingleFit( h1, f1, outputdir, name );


  ResoStruct rs = getRespAndReso( f1, 0. );

  return rs;

}
