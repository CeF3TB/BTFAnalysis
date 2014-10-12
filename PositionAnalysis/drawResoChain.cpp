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
#include "TLegend.h"

#include "interface/DrawTools.h"
#include "interface/FitTools.h"

#include "TApplication.h"

//the bulky thingy that lets you draw the resolution at 491 MeV
//for different setups of the simulation with increasing "reality"

struct ResoStruct {

  float resp;
  float resp_error;

  float reso;
  float reso_error;

};

ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, const std::string& name, float energyS );
ResoStruct getRespAndReso( TF1* f1, float energyS );
float getRatioError( float num, float denom, float numErr, float denomErr );
ResoStruct addPhotoStatistics( ResoStruct rs );


int main( int argc, char* argv[] ) {

  TApplication* a = new TApplication("a",0,0);

  std::string outputdir = "ResoChain/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  DrawTools::setStyle();

  std::vector<std::string> simulationdir; 
  float  beamEnergy = 491.; 

  //  runs.push_back("BTF_259_20140502-012847_beam");
  //  beamEnergy.push_back(491.4);

  simulationdir.push_back("Ideal20Layers1m");
  simulationdir.push_back("RealWOAirSmear");
  simulationdir.push_back("RealAirSmear");
  //simulationdir.push_back("Ideal20Layers1mAir");
  //simulationdir.push_back("Ideal20Layers1mAirScint");
  simulationdir.push_back("RealScint");
  simulationdir.push_back("BTF_259_20140502-012847_beam");



  TGraphErrors* gr_resp = new TGraphErrors(0);
  TGraphErrors* gr_reso = new TGraphErrors(0);
  TGraphErrors* gr_data = new TGraphErrors(0);


  TFile* fileD = TFile::Open(Form("analysisTrees_V03/Reco_BTF_259_20140502-012847_beam.root"));
  TTree* treeD= (TTree*)fileD->Get("recoTree");
  TH1D* h1 = new TH1D( "data", "", 500, 50., 6000. );
  treeD->Project( "data","cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 )" );
  TF1* f1 = new TF1( Form("f1_%s", "data"), "gaus", 100., 6000.);
  FitTools::doSingleFit( h1, f1, outputdir, "data" , 4 ,1.5 );

  ResoStruct rs5 = getRespAndReso( f1, beamEnergy );
  ResoStruct rs_ps5 = addPhotoStatistics( rs5 );





  TFile* fileS0 = TFile::Open(Form("OriginalSimulationData/%s/Reco_SimulationIdeal_491.root", simulationdir[0].c_str()));
  TTree* treeS0 = (TTree*)fileS0->Get("recoTree");
  
  ResoStruct rs0 = getResponseResolutionMC( outputdir, treeS0, simulationdir[0], beamEnergy);
  
  gr_resp->SetPoint( 0+1, 0+1, rs0.resp );
  gr_resp->SetPointError( 0+1,0, rs0.resp_error);
  
  gr_reso->SetPoint( 0+1, 0+1, rs0.reso );
  gr_reso->SetPointError( 0+1,0,  rs0.reso_error );
  



  for( unsigned i=0; i<simulationdir.size()-1; ++i ) {

    TFile* fileS = TFile::Open(Form("OriginalSimulationData/%s/Reco_SimulationIdeal_491.root", simulationdir[i].c_str()));
    TTree* treeS = (TTree*)fileS->Get("recoTree");

    ResoStruct rs = getResponseResolutionMC( outputdir, treeS, simulationdir[i], beamEnergy );

    ResoStruct rs_ps = addPhotoStatistics( rs );

    gr_resp->SetPoint(i+1+1,  i+1+1,  rs.resp  );
    gr_resp->SetPointError( i+1+1,0, rs.resp_error);

    gr_reso->SetPoint( i+1+1, i+1+1, rs_ps.reso );
    gr_reso->SetPointError( i+1+1,0,  rs_ps.reso_error);

    gr_data->SetPoint( i, i, rs_ps5.reso );
    gr_data->SetPointError(i,0,  rs_ps5.reso_error );

  }


  TFile* fileS = TFile::Open("OriginalSimulationData/RealWHodo/Reco_Simulation491.root");
  TTree* treeS = (TTree*)fileS->Get("recoTree");
  
  ResoStruct rs = getResponseResolutionMC( outputdir, treeS, "RealWHodo", beamEnergy );
  ResoStruct rs_ps = addPhotoStatistics( rs );
  
  gr_resp->SetPoint(6,  6,  rs.resp  );
  gr_resp->SetPointError(6,0, rs.resp_error);
  
  gr_reso->SetPoint( 6, 6, rs_ps.reso );
  gr_reso->SetPointError( 6 ,0,  rs_ps.reso_error);

  TFile* fileE = TFile::Open("OriginalSimulationData/RealEnergy/Reco_Simulation491.root");
  TTree* treeE = (TTree*)fileE->Get("recoTree");
  
  ResoStruct rsE = getResponseResolutionMC( outputdir, treeE, "RealEnergy", beamEnergy );
  ResoStruct rs_psE = addPhotoStatistics( rsE );

  gr_resp->SetPoint(7,  7,  rsE.resp  );
  gr_resp->SetPointError(7,0, rsE.resp_error);
  gr_reso->SetPoint( 7, 7, rs_psE.reso );
  gr_reso->SetPointError( 7 ,0,  rs_psE.reso_error);
    

  gr_data->SetPoint( 4, 4, rs_ps5.reso );
  gr_data->SetPointError(4,0,  rs_ps5.reso_error );
  gr_data->SetPoint( 5, 5, rs_ps5.reso );
  gr_data->SetPointError(5,0,  rs_ps5.reso_error );    
  gr_data->SetPoint( 6, 6, rs_ps5.reso );
  gr_data->SetPointError(6,0,  rs_ps5.reso_error ); 
  gr_data->SetPoint( 7, 7, rs_ps5.reso );
  gr_data->SetPointError(7,0,  rs_ps5.reso_error );
  gr_data->SetPoint( 8, 8, rs_ps5.reso );
  gr_data->SetPointError(8,0,  rs_ps5.reso_error );
  
  
  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();



  TH2D* h2_axes = new TH2D( "axes", "", 10, -0.49, 5.49, 10, 0., 250. );
  h2_axes->SetXTitle("Increasing Reality");
  h2_axes->SetYTitle("CeF3 Response [MeV]");
  h2_axes->Draw("");

  gr_resp->SetMarkerStyle(20);
  gr_resp->SetMarkerSize(1.6);
  gr_resp->SetMarkerColor(46);
  gr_resp->Draw("p same");

  /*  TF1* f1_line = new TF1("line", "[0] + [1]*x", 0., xMax );
    gr_resp->Fit(f1_line, "RN");
  f1_line->SetLineColor(46);
  f1_line->SetLineWidth(1.);
  f1_line->Draw("L same");
  */

  TLegend* leg0 = new TLegend(0.55, 0.2, 0.85, 0.4);
  leg0->SetTextSize(0.038);
  leg0->AddEntry(gr_resp,"simulation","p");
  /*
  leg0->AddEntry(f1_line,Form("offset=%.0f\n #pm %.0f\n",f1_line->GetParameter(0), f1_line->GetParError(0) ),"L");
  leg0->AddEntry( (TObject*)0, Form("chi2/ndf=%.2f\n / %d",f1_line->GetChisquare(), f1_line->GetNDF() ), "");
  leg0->AddEntry(gr_resp_vs_energy_simul,"simulation","p");
  leg0->AddEntry(f1_lines,Form("offset=%.0f\n #pm %.0f\n",f1_lines->GetParameter(0), f1_lines->GetParError(0) ),"L");
  leg0->AddEntry( (TObject*)0, Form("chi2/ndf=%.2f\n / %d",f1_lines->GetChisquare(), f1_lines->GetNDF() ), "");
  */
  leg0->SetFillColor(0);
  leg0->Draw("same");

  c1->SaveAs( Form( "%s/resp_vs_reality.pdf", outputdir.c_str() ) );

  c1->Clear();



  ///////////////////////RESOLUTION///////////////////////////////
  TH2D* h2_axes2 = new TH2D( "axes", "", 5, 0.49,7.49, 50, 0., 20. );
  h2_axes2->SetXTitle("Increasing Reality [Arbitrary Units]");
  h2_axes2->SetYTitle("CeF3 Resolution [%]");
  h2_axes2->Draw("");

  std::string names[5] = {"1mx1m","+P.S.","24mmX24mm","+air&spread","+scint"};
  h2_axes->GetXaxis()->SetNdivisions(5);

  gr_reso->SetFillColor(46);
  gr_reso->SetFillStyle(3005);

  gr_reso->SetMarkerStyle(20);
  gr_reso->SetMarkerSize(1.6);
  gr_reso->SetMarkerColor(46);
  gr_reso->Draw("p 3 same");

  gr_data->SetFillColor(38);
  // gr_data->SetFillStyle(3005);

  gr_data->SetMarkerStyle(20);
  gr_data->SetMarkerSize(1.6);
  gr_data->SetMarkerColor(38);
  gr_data->Draw("3 same");

  TLegend* leg4 = new TLegend(0.55, 0.465, 0.9, 0.2);
  leg4->SetTextSize(0.038);
  leg4->AddEntry(gr_reso,"simulations","p");
  leg4->AddEntry(gr_data,"data","p");
  leg4->SetFillColor(0);
  leg4->Draw("same");

  c1->SaveAs( Form( "%s/reso_vs_reality.pdf", outputdir.c_str() ) );



  return 0;
}






ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, const std::string& name, float energyS ) {

  TH1D* h1 = new TH1D( name.c_str(), "", 500, 0., 500. );

  tree->Project( name.c_str(),"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(isSingleEle_scintFront && nHodoClustersX==1 && nHodoClustersY==1 )" );

  TF1* f1 = new TF1( Form("f1_%s", name.c_str()), "gaus", 100., 500.);

  FitTools::doSingleFit( h1, f1, outputdir, name,4,1.5 );

  ResoStruct rs = getRespAndReso( f1, energyS );

  return rs;

}



ResoStruct getRespAndReso( TF1* f1, float energyS) {

  // float energy = energyS;
  float energyErrS = 0;

  float meanS = f1->GetParameter(1);
  float meanErrS = f1->GetParError(1);

  float rmsS = f1->GetParameter(2);
  float rmsErrS = f1->GetParError(2);

   float resoS = 100.* rmsS/meanS; //in percent
  //float resoS = 100.* rmsS/meanS; //in percent
  float resoErrS = 100.* getRatioError( rmsS, meanS, meanErrS, rmsErrS );

  ResoStruct rs;
  rs.resp = meanS;
  rs.resp_error = meanErrS;
  rs.reso = resoS;
  rs.reso_error = resoErrS;

  return rs;

}


float getRatioError( float num, float denom, float numErr, float denomErr ) {

  return sqrt( numErr*numErr/(denom*denom) + denomErr*denomErr*num*num/(denom*denom*denom*denom) );

}


ResoStruct addPhotoStatistics( ResoStruct rs ) {

  // MC response is already in MeV, 0.49 is the number of p.e. per MeV
  // RM = 90% of energy, so assume 90% of energy (0.9*491) is deposited in central channel
  float nADC = rs.resp/173.1*3236.14;
  //   float nADC = rs.resp/169.784*3236.14;
  float nPhotoElectrons = nADC/27.3;
  //float nADC = rs.resp/185.*3200.;
  //float nPhotoElectrons = nADC/35.;

  float poissonError = 100./sqrt( nPhotoElectrons ); // in percent
 
  rs.reso = sqrt( rs.reso*rs.reso + poissonError*poissonError );
 
  return rs;

} 
